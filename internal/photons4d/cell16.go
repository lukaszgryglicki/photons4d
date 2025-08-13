package photons4d

import (
	"fmt"
	"math"
)

// Cell16: regular 4D cross-polytope (16 facets, 8 vertices).
// Build: canonical ±e_i -> uniform edge scale -> per-axis scale -> rotate -> translate.
// Facets are given by sign-vectors s∈{±1}^4 with planes s·x = 1 in canonical space.
// We transform planes via A^{-T} with A = R * diag(scale), and set d by max over world verts.
type Cell16 struct {
	Center  Point4
	Scale   Vector4
	R       Mat4
	RT      Mat4
	Color   RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	Verts   [8]Point4 // world-space vertices
	AABBMin Point4
	AABBMax Point4

	N [16]Vector4 // outward facet normals (world)
	D [16]Real    // N·x <= D is inside

	// material caches
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func canonical16Verts() [8]Vector4 {
	return [8]Vector4{
		{+1, 0, 0, 0},
		{-1, 0, 0, 0},
		{0, +1, 0, 0},
		{0, -1, 0, 0},
		{0, 0, +1, 0},
		{0, 0, -1, 0},
		{0, 0, 0, +1},
		{0, 0, 0, -1},
	}
}

// 16 sign-vectors for facets: all ±1 choices.
func signVectors16() [16]Vector4 {
	out := [16]Vector4{}
	idx := 0
	for sx := -1; sx <= 1; sx += 2 {
		for sy := -1; sy <= 1; sy += 2 {
			for sz := -1; sz <= 1; sz += 2 {
				for sw := -1; sw <= 1; sw += 2 {
					out[idx] = Vector4{Real(sx), Real(sy), Real(sz), Real(sw)}
					idx++
				}
			}
		}
	}
	return out
}

func NewCell16(
	center Point4,
	scale Vector4,
	angles Rot4,
	color, reflectivity, refractivity, ior RGB,
) (*Cell16, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("16-cell per-axis scale must be > 0 on all axes, got %+v", scale)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	for _, p := range []struct {
		n    string
		r, t Real
	}{
		{"R", reflectivity.R, refractivity.R},
		{"G", reflectivity.G, refractivity.G},
		{"B", reflectivity.B, refractivity.B},
	} {
		if !in01(p.r) || !in01(p.t) {
			return nil, fmt.Errorf("reflect/refract must be in [0,1]; channel %s got reflect=%.6g refract=%.6g", p.n, p.r, p.t)
		}
		if p.r+p.t > 1+1e-12 {
			return nil, fmt.Errorf("per-channel reflect+refract must be ≤1; channel %s got %.6g", p.n, p.r+p.t)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	R := rotFromAngles(angles)
	sx, sy, sz, sw := scale.X, scale.Y, scale.Z, scale.W

	// canonical vertices
	V := canonical16Verts()

	// base edge = min nonzero pair distance in canonical (sqrt(2) for ±e_i)
	base := 1e300
	for i := 0; i < 8; i++ {
		for j := i + 1; j < 8; j++ {
			d := V[i].Sub(V[j]).Len()
			if d > 1e-12 && d < base {
				base = d
			}
		}
	}
	u := 1.0 / base

	// transform vertices to world
	var Wv [8]Point4
	minX, maxX := +1e300, -1e300
	minY, maxY := +1e300, -1e300
	minZ, maxZ := +1e300, -1e300
	minW, maxW := +1e300, -1e300
	for i := 0; i < 8; i++ {
		p := V[i].Mul(u)
		p = Vector4{p.X * sx, p.Y * sy, p.Z * sz, p.W * sw}
		pw := R.MulVec(p)
		w := Point4{center.X + pw.X, center.Y + pw.Y, center.Z + pw.Z, center.W + pw.W}
		Wv[i] = w
		if w.X < minX {
			minX = w.X
		}
		if w.X > maxX {
			maxX = w.X
		}
		if w.Y < minY {
			minY = w.Y
		}
		if w.Y > maxY {
			maxY = w.Y
		}
		if w.Z < minZ {
			minZ = w.Z
		}
		if w.Z > maxZ {
			maxZ = w.Z
		}
		if w.W < minW {
			minW = w.W
		}
		if w.W > maxW {
			maxW = w.W
		}
	}

	c := &Cell16{
		Center:  center,
		Scale:   scale,
		R:       R,
		RT:      R.Transpose(),
		Color:   color,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
		Verts:   Wv,
		AABBMin: Point4{minX, minY, minZ, minW},
		AABBMax: Point4{maxX, maxY, maxZ, maxW},
	}

	// materials caches
	c.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	c.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	c.colorArr = [3]Real{color.R, color.G, color.B}
	c.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - c.refl[i] - c.refr[i]
		if p < 0 {
			p = 0
		}
		c.pAbs[i] = p
		c.iorInv[i] = 1 / c.iorArr[i]
		n := c.iorArr[i]
		r0 := (n - 1) / (n + 1)
		c.f0[i] = r0 * r0
	}

	// planes: 16 sign-vectors s; transform normals by A^{-T} (A=R*S), set d as max(N·v_world)
	Sinv := Mat4{M: [4][4]Real{
		{1 / sx, 0, 0, 0},
		{0, 1 / sy, 0, 0},
		{0, 0, 1 / sz, 0},
		{0, 0, 0, 1 / sw},
	}}
	AinvT := R.Mul(Sinv) // A^{-T} = R * S^{-1} (since R is orthonormal)
	signs := signVectors16()
	for i := 0; i < 16; i++ {
		nw := AinvT.MulVec(signs[i]) // not unit yet; scale cancels with d
		// compute d as support over vertices
		d := -1e300
		for j := 0; j < 8; j++ {
			val := nw.X*c.Verts[j].X + nw.Y*c.Verts[j].Y + nw.Z*c.Verts[j].Z + nw.W*c.Verts[j].W
			if val > d {
				d = val
			}
		}
		// normalize N (makes clipping numerically nicer)
		L := math.Sqrt(nw.Dot(nw))
		if L > 0 {
			nw = nw.Mul(1 / L)
			d /= L
		}
		c.N[i], c.D[i] = nw, d
	}

	// DebugLog("Created 16-cell: center=%+v, scale=%+v, AABB=[%+v .. %+v]", center, scale, c.AABBMin, c.AABBMax)
	DebugLog("Created 16-cell: %+v", c)
	return c, nil
}

func intersectRayCell16(O Point4, D Vector4, c *Cell16) (hit objectHit, ok bool) {
	// identical half-space clip as cell5
	tEnter, tExit := -1e300, 1e300
	enterIdx, exitIdx := -1, -1
	for i := 0; i < 16; i++ {
		n := c.N[i]
		d := c.D[i]
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
	return objectHit{t: t, Nw: c.N[idx], c16: c, inv: inv}, true
}
