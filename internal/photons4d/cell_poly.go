package photons4d

import "math"


// Polytope										Vertices		Edges		Faces (2D)			Cells (3D)
// 5-cell 	(simplex)					5						10			10 triangles		5 tetrahedra
// 8-cell 	(tesseract)				16					32			24 squares			8 cubes
// 16-cell 	(cross-polytope)	8						24			32 triangles		16 tetrahedra
// 24-cell										24					96			96 triangles		24 octahedra
// 120-cell										600					1200		720 pentagons		120 dodecahedra
// 600-cell										120					720			1200 triangles	600 tetrahedra

// Convex 4D polytope expressed as intersection of planes N[i]·x <= D[i] (N need not be unit).
// Built by taking unit-radius local normals 'nl[i]' and an offset 'Radius' (scalar), then
// transforming by Scale (anisotropic), Rotation, and Center translation.
type cellPoly struct {
	Center Point4
	Scale  Vector4
	R      Mat4
	RT     Mat4

	Color   RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	// planes in world space
	N []Vector4 // not necessarily unit
	U []Vector4 // corresponding unit normals (outward)
	D []Real

	// aabb
	AABBMin Point4
	AABBMax Point4

	// material caches (same layout as others)
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func (cp *cellPoly) buildPlanes(localNormals []Vector4, radius Real) {
	// world N = R * S^{-T} * n
	sx, sy, sz, sw := cp.Scale.X, cp.Scale.Y, cp.Scale.Z, cp.Scale.W
	inv := func(x Real) Real { return 1 / x }
	for _, n := range localNormals {
		// S^{-T}*n
		m := Vector4{n.X * inv(sx), n.Y * inv(sy), n.Z * inv(sz), n.W * inv(sw)}
		// R * (...)
		m = cp.R.MulVec(m)
		u := m
		l := u.Len()
		if l > 0 {
			u = u.Mul(1 / l)
		}
		d := radius + m.X*cp.Center.X + m.Y*cp.Center.Y + m.Z*cp.Center.Z + m.W*cp.Center.W
		cp.N = append(cp.N, m)
		cp.U = append(cp.U, u)
		cp.D = append(cp.D, d)
	}
	// AABB via support function along axes
	cp.computeAABB()
}

func (cp *cellPoly) computeAABB() {
	minV := [4]Real{+1e300, +1e300, +1e300, +1e300}
	maxV := [4]Real{-1e300, -1e300, -1e300, -1e300}

	for axis := 0; axis < 4; axis++ {
		// max along +axis
		tmax := 1e300
		for i := range cp.N {
			a := cp.N[i]
			ai := []Real{a.X, a.Y, a.Z, a.W}[axis]
			if ai > 1e-18 {
				if v := cp.D[i] / ai; v < tmax {
					tmax = v
				}
			}
		}
		// max along -axis => min coord
		tminp := 1e300
		for i := range cp.N {
			a := cp.N[i]
			ai := []Real{a.X, a.Y, a.Z, a.W}[axis]
			if ai < -1e-18 {
				if v := cp.D[i] / (-ai); v < tminp {
					tminp = v
				}
			}
		}

		// Fallback if something degenerated numerically on this axis.
		if tmax == 1e300 {
			tmax = []Real{cp.Center.X, cp.Center.Y, cp.Center.Z, cp.Center.W}[axis]
		}
		if tminp == 1e300 {
			tminp = []Real{cp.Center.X, cp.Center.Y, cp.Center.Z, cp.Center.W}[axis]
			if tminp < 0 {
				tminp = -tminp
			}
			// ensure nonzero thickness
			if tminp == 0 {
				tminp = 1
			}
		}

		maxV[axis] = tmax
		minV[axis] = -tminp
	}

	// Ensure the box is valid and contains the center (conservative fixup).
	cx := []Real{cp.Center.X, cp.Center.Y, cp.Center.Z, cp.Center.W}
	for axis := 0; axis < 4; axis++ {
		if minV[axis] > maxV[axis] {
			minV[axis], maxV[axis] = maxV[axis], minV[axis]
		}
		if cx[axis] < minV[axis] {
			minV[axis] = cx[axis]
		}
		if cx[axis] > maxV[axis] {
			maxV[axis] = cx[axis]
		}
	}

	cp.AABBMin = Point4{minV[0], minV[1], minV[2], minV[3]}
	cp.AABBMax = Point4{maxV[0], maxV[1], maxV[2], maxV[3]}
}

func (cp *cellPoly) materialFrom(color, refl, refr, ior RGB) {
	cp.Color, cp.Reflect, cp.Refract, cp.IOR = color, refl, refr, ior
	cp.refl = [3]Real{refl.R, refl.G, refl.B}
	cp.refr = [3]Real{refr.R, refr.G, refr.B}
	cp.colorArr = [3]Real{color.R, color.G, color.B}
	cp.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - cp.refl[i] - cp.refr[i]
		if p < 0 {
			p = 0
		}
		cp.pAbs[i] = p
		cp.iorInv[i] = 1 / cp.iorArr[i]
		n := cp.iorArr[i]
		r0 := (n - 1) / (n + 1)
		cp.f0[i] = r0 * r0
	}
}

// intersection by half-space clipping (like cell5)
func intersectRayCellPoly(O Point4, D Vector4, cp *cellPoly) (objectHit, bool) {
	tEnter := -1e300
	tExit := 1e300
	enterIdx := -1
	exitIdx := -1

	for i := range cp.N {
		n := cp.N[i]
		d := cp.D[i]
		nO := n.X*O.X + n.Y*O.Y + n.Z*O.Z + n.W*O.W
		nD := n.X*D.X + n.Y*D.Y + n.Z*D.Z + n.W*D.W
		rhs := d - nO

		const eps = 1e-12
		if math.Abs(nD) < eps {
			if rhs < -1e-12 {
				return objectHit{}, false
			}
			continue
		}
		t := rhs / nD
		if nD > 0 {
			if t < tExit {
				tExit = t
				exitIdx = i
			}
		} else {
			if t > tEnter {
				tEnter = t
				enterIdx = i
			}
		}
	}

	if tExit < 1e-12 && tEnter < 1e-12 {
		return objectHit{}, false
	}
	if tEnter > tExit {
		return objectHit{}, false
	}

	inv := tEnter <= 1e-12
	var t Real
	var idx int
	if !inv {
		t, idx = tEnter, enterIdx
	} else {
		t, idx = tExit, exitIdx
	}
	if t <= 1e-12 || idx < 0 {
		return objectHit{}, false
	}
	Nw := cp.U[idx]
	return objectHit{t: t, Nw: Nw, poly: cp, inv: inv}, true
}
