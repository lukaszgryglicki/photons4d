package photons4d

import (
	"fmt"
	"math"
)

// Cell5 represents a (possibly anisotropically scaled) 4D regular cell5 (5-cell).
// Build pipeline (local space):
//  1. start from a centered regular 5-cell with canonical coordinates,
//  2. uniformly scale to match requested side length,
//  3. apply per-axis scale (Scale.X/Y/Z/W),
//  4. rotate by R about origin,
//  5. translate to Center.
//
// We compute world-space facet planes (outward normals) from transformed vertices,
// so normals are correct even under anisotropic scaling.
type Cell5 struct {
	Center Point4
	Scale  Vector4 // per-axis scaling
	R      Mat4    // local->world rotation
	RT     Mat4    // world->local (R^T)

	Diffuse RGB
	Color   RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	// cached geometry
	Verts   [5]Point4 // world-space vertices
	AABBMin Point4
	AABBMax Point4

	// facet planes: N·x <= d is "inside" (N is outward)
	N [5]Vector4
	D [5]Real

	// material caches
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

// ---- math helpers -----------------------------------------------------------

// 4D normal to a 3D facet spanned by (u,v,w): generalized cross product.
// n is orthogonal to u,v,w. Signs per standard cofactor expansion.
func cross3_4D(u, v, w Vector4) Vector4 {
	det3 := func(a1, a2, a3, b1, b2, b3, c1, c2, c3 Real) Real {
		return a1*(b2*c3-b3*c2) - a2*(b1*c3-b3*c1) + a3*(b1*c2-b2*c1)
	}
	nx := det3(u.Y, u.Z, u.W, v.Y, v.Z, v.W, w.Y, w.Z, w.W)
	ny := -det3(u.X, u.Z, u.W, v.X, v.Z, v.W, w.X, w.Z, w.W)
	nz := det3(u.X, u.Y, u.W, v.X, v.Y, v.W, w.X, w.Y, w.W)
	nw := -det3(u.X, u.Y, u.Z, v.X, v.Y, v.Z, w.X, w.Y, w.Z)
	return Vector4{nx, ny, nz, nw}
}

func min4(a, b, c, d Real) Real { return math.Min(math.Min(a, b), math.Min(c, d)) }
func max4(a, b, c, d Real) Real { return math.Max(math.Max(a, b), math.Max(c, d)) }

// canonical regular 5-cell vertices in R^4 centered at origin.
// Constructed from 5D {e_i - centroid} projected onto the 4D
// subspace orthogonal to (1,1,1,1,1) via an orthonormal basis.
func canonicalCell5() [5]Vector4 {
	// Orthonormal rows (b1..b4) ⟂ (1,1,1,1,1).
	// Each bi is length 1 and rows are mutually orthogonal.
	B := [4][5]Real{
		{1 / math.Sqrt2, -1 / math.Sqrt2, 0, 0, 0},                                                                                     // (1,-1,0,0,0)/√2
		{1 / Real(math.Sqrt(6)), 1 / Real(math.Sqrt(6)), -2 / Real(math.Sqrt(6)), 0, 0},                                                // (1,1,-2,0,0)/√6
		{1 / Real(math.Sqrt(12)), 1 / Real(math.Sqrt(12)), 1 / Real(math.Sqrt(12)), -3 / Real(math.Sqrt(12)), 0},                       // (1,1,1,-3,0)/√12
		{1 / Real(math.Sqrt(20)), 1 / Real(math.Sqrt(20)), 1 / Real(math.Sqrt(20)), 1 / Real(math.Sqrt(20)), -4 / Real(math.Sqrt(20))}, // (1,1,1,1,-4)/√20
	}
	// e_i - centroid
	w := [5][5]Real{}
	for i := 0; i < 5; i++ {
		for k := 0; k < 5; k++ {
			w[i][k] = -0.2
		}
		w[i][i] = 0.8
	}
	var V [5]Vector4
	for i := 0; i < 5; i++ {
		V[i] = Vector4{
			B[0][0]*w[i][0] + B[0][1]*w[i][1] + B[0][2]*w[i][2] + B[0][3]*w[i][3] + B[0][4]*w[i][4],
			B[1][0]*w[i][0] + B[1][1]*w[i][1] + B[1][2]*w[i][2] + B[1][3]*w[i][3] + B[1][4]*w[i][4],
			B[2][0]*w[i][0] + B[2][1]*w[i][1] + B[2][2]*w[i][2] + B[2][3]*w[i][3] + B[2][4]*w[i][4],
			B[3][0]*w[i][0] + B[3][1]*w[i][1] + B[3][2]*w[i][2] + B[3][3]*w[i][3] + B[3][4]*w[i][4],
		}
	}
	return V
}

// ---- construction -----------------------------------------------------------

func NewCell5(
	center Point4,
	scale Vector4,
	angles Rot4, // radians
	color, diffuse, reflectivity, refractivity, ior RGB,
) (*Cell5, error) {
	// validate scale
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("cell5 per-axis scale must be > 0 on all axes, got %+v", scale)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	type ch struct {
		n       string
		r, t, d Real
	}
	for _, c := range []ch{
		{"R", reflectivity.R, refractivity.R, diffuse.R},
		{"G", reflectivity.G, refractivity.G, diffuse.G},
		{"B", reflectivity.B, refractivity.B, diffuse.B},
	} {
		if !in01(c.r) || !in01(c.t) || !in01(c.d) {
			return nil, fmt.Errorf("reflect/refract/diffuse must be in [0,1]; channel %s got reflect=%.6g refract=%.6g diffuse=%.6g", c.n, c.r, c.t, c.d)
		}
		if c.r+c.t+c.d > 1+1e-12 {
			return nil, fmt.Errorf("per-channel reflect+refract+diffuse must be ≤1; channel %s got %.6g", c.n, c.r+c.t+c.d)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	R := rotFromAngles(angles)
	sx, sy, sz, sw := scale.X, scale.Y, scale.Z, scale.W

	// canonical vertices (centered at origin)
	V := canonicalCell5()
	// base edge length
	base := V[0].Sub(V[1]).Len()
	u := 1.0 / base

	// transform to world
	var Wv [5]Point4
	for i := 0; i < 5; i++ {
		p := V[i].Mul(u)                                                                   // uniform side length
		p = Vector4{p.X * sx, p.Y * sy, p.Z * sz, p.W * sw}                                // per-axis scale
		pw := R.MulVec(p)                                                                  // rotate
		Wv[i] = Point4{center.X + pw.X, center.Y + pw.Y, center.Z + pw.Z, center.W + pw.W} // translate
	}

	// AABB
	minX, maxX := Wv[0].X, Wv[0].X
	minY, maxY := Wv[0].Y, Wv[0].Y
	minZ, maxZ := Wv[0].Z, Wv[0].Z
	minW, maxW := Wv[0].W, Wv[0].W
	for i := 1; i < 5; i++ {
		minX = math.Min(minX, Wv[i].X)
		maxX = math.Max(maxX, Wv[i].X)
		minY = math.Min(minY, Wv[i].Y)
		maxY = math.Max(maxY, Wv[i].Y)
		minZ = math.Min(minZ, Wv[i].Z)
		maxZ = math.Max(maxZ, Wv[i].Z)
		minW = math.Min(minW, Wv[i].W)
		maxW = math.Max(maxW, Wv[i].W)
	}

	sx4 := &Cell5{
		Center:  center,
		Scale:   scale,
		R:       R,
		RT:      R.Transpose(),
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,

		Verts:   Wv,
		AABBMin: Point4{minX, minY, minZ, minW},
		AABBMax: Point4{maxX, maxY, maxZ, maxW},
	}

	// material caches
	sx4.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	sx4.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	sx4.colorArr = [3]Real{color.R, color.G, color.B}
	sx4.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	sx4.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - sx4.refl[i] - sx4.refr[i] - sx4.diff[i]
		if p < 0 {
			p = 0
		}
		sx4.pAbs[i] = p
		sx4.iorInv[i] = 1 / sx4.iorArr[i]
		n := sx4.iorArr[i]
		r0 := (n - 1) / (n + 1)
		sx4.f0[i] = r0 * r0
	}

	// facet planes: for each i, plane through the other four vertices.
	for i := 0; i < 5; i++ {
		// collect the 4 vertices not equal to i
		idx := make([]int, 0, 4)
		for j := 0; j < 5; j++ {
			if j != i {
				idx = append(idx, j)
			}
		}
		// pick p0 and spanning vectors (in world space)
		p0 := sx4.Verts[idx[0]]
		u1 := Vector4{
			sx4.Verts[idx[1]].X - p0.X,
			sx4.Verts[idx[1]].Y - p0.Y,
			sx4.Verts[idx[1]].Z - p0.Z,
			sx4.Verts[idx[1]].W - p0.W,
		}
		u2 := Vector4{
			sx4.Verts[idx[2]].X - p0.X,
			sx4.Verts[idx[2]].Y - p0.Y,
			sx4.Verts[idx[2]].Z - p0.Z,
			sx4.Verts[idx[2]].W - p0.W,
		}
		u3 := Vector4{
			sx4.Verts[idx[3]].X - p0.X,
			sx4.Verts[idx[3]].Y - p0.Y,
			sx4.Verts[idx[3]].Z - p0.Z,
			sx4.Verts[idx[3]].W - p0.W,
		}
		n := cross3_4D(u1, u2, u3).Norm()
		d := n.X*p0.X + n.Y*p0.Y + n.Z*p0.Z + n.W*p0.W

		// With "inside: n·x ≤ d" and n being the outward normal,
		// the excluded vertex must lie on the INSIDE side (n·v_excl ≤ d).
		vout := sx4.Verts[i]
		if n.X*vout.X+n.Y*vout.Y+n.Z*vout.Z+n.W*vout.W > d {
			n = n.Mul(-1)
			d = -d
		}
		sx4.N[i], sx4.D[i] = n, d
	}

	// DebugLog("Created 5-cell: center=%+v, scale=%+v, AABB=[%+v .. %+v]", center, scale, sx4.AABBMin, sx4.AABBMax)
	DebugLog("Created 5-cell: %+v", sx4)
	return sx4, nil
}

// intersect ray with a convex 4D polytope defined by planes N·x <= d (N outward).
// Returns first positive intersection if starting outside, or the exit if starting inside.
func intersectRayCell5(O Point4, D Vector4, s *Cell5) (hit objectHit, ok bool) {
	// Slab-like half-space clipping across planes
	tEnter := -1e300
	tExit := 1e300
	enterIdx := -1
	exitIdx := -1

	for i := 0; i < 5; i++ {
		n := s.N[i]
		d := s.D[i]
		nO := n.X*O.X + n.Y*O.Y + n.Z*O.Z + n.W*O.W
		nD := n.X*D.X + n.Y*D.Y + n.Z*D.Z + n.W*D.W
		rhs := d - nO

		const eps = 1e-12
		if nD > -eps && nD < eps {
			// parallel to plane
			if rhs < -1e-12 {
				return objectHit{}, false // outside and never entering
			}
			// inside for this plane; no bounds update
			continue
		}

		t := rhs / nD
		if nD > 0 {
			// t <= t_max
			if t < tExit {
				tExit = t
				exitIdx = i
			}
		} else {
			// nD < 0 ⇒ t >= t_min
			if t > tEnter {
				tEnter = t
				enterIdx = i
			}
		}
	}

	if tExit < 1e-12 && tEnter < 1e-12 {
		// cell5 is entirely behind or we're touching at origin; treat as no hit
		return objectHit{}, false
	}
	if tEnter > tExit {
		return objectHit{}, false
	}

	// Decide whether starting inside
	inv := tEnter <= 1e-12
	var t Real
	var idx int
	if !inv {
		t = tEnter
		idx = enterIdx
	} else {
		t = tExit
		idx = exitIdx
	}
	if t <= 1e-12 || idx < 0 {
		return objectHit{}, false
	}
	Nw := s.N[idx] // outward normal of the hit facet

	return objectHit{t: t, Nw: Nw, mat: s, inv: inv}, true
}

func (c *Cell5) PAbsCh(i int) Real   { return c.pAbs[i] }
func (c *Cell5) DiffCh(i int) Real   { return c.diff[i] }
func (c *Cell5) ColorCh(i int) Real  { return c.colorArr[i] }
func (c *Cell5) F0Ch(i int) Real     { return c.f0[i] }
func (c *Cell5) ReflCh(i int) Real   { return c.refl[i] }
func (c *Cell5) RefrCh(i int) Real   { return c.refr[i] }
func (c *Cell5) IORCh(i int) Real    { return c.iorArr[i] }
func (c *Cell5) IORInvCh(i int) Real { return c.iorInv[i] }
