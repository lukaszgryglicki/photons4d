package photons4d

import (
	"fmt"
	"math"
)

// TwentyFourCell: regular 24-cell (self-dual). We'll use the D4 root set:
// all permutations of (±1, ±1, 0, 0). That’s 24 canonical vertices.
// Planes are also 24 in number (self-dual). For each canonical normal seed n_l (same set),
// world-space normal n_w = (R*S^{-1}) n_l (scale cancels in d via support).
type TwentyFourCell struct {
	Center  Point4
	Side    Real
	Scale   Vector4
	R       Mat4
	RT      Mat4
	Color   RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	Verts   [24]Point4
	AABBMin Point4
	AABBMax Point4

	N [24]Vector4
	D [24]Real

	// material caches
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func canonical24Verts() [24]Vector4 {
	out := [24]Vector4{}
	idx := 0
	// choose positions of the two ±1's
	pos := [][2]int{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}
	for _, p := range pos {
		i, j := p[0], p[1]
		for si := -1; si <= 1; si += 2 {
			for sj := -1; sj <= 1; sj += 2 {
				v := [4]Real{0, 0, 0, 0}
				v[i] = Real(si)
				v[j] = Real(sj)
				out[idx] = Vector4{v[0], v[1], v[2], v[3]}
				idx++
			}
		}
	}
	return out
}

func NewTwentyFourCell(
	center Point4,
	side Real,
	scale Vector4,
	angles Rot4,
	color, reflectivity, refractivity, ior RGB,
) (*TwentyFourCell, error) {
	if !(side > 0) {
		return nil, fmt.Errorf("24-cell side must be > 0, got %.6g", side)
	}
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("24-cell per-axis scale must be > 0 on all axes, got %+v", scale)
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

	V := canonical24Verts()

	// base edge = min nonzero distance between canonical vertices
	base := 1e300
	for i := 0; i < 24; i++ {
		for j := i + 1; j < 24; j++ {
			d := V[i].Sub(V[j]).Len()
			if d > 1e-12 && d < base {
				base = d
			}
		}
	}
	u := side / base

	var Wv [24]Point4
	minX, maxX := +1e300, -1e300
	minY, maxY := +1e300, -1e300
	minZ, maxZ := +1e300, -1e300
	minW, maxW := +1e300, -1e300
	for i := 0; i < 24; i++ {
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

	c := &TwentyFourCell{
		Center:  center,
		Side:    side,
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

	// material caches
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

	// planes: seeds = same canonical set (self-dual), transform normals by A^{-T}
	Sinv := Mat4{M: [4][4]Real{
		{1 / sx, 0, 0, 0},
		{0, 1 / sy, 0, 0},
		{0, 0, 1 / sz, 0},
		{0, 0, 0, 1 / sw},
	}}
	AinvT := R.Mul(Sinv)
	seeds := canonical24Verts()
	for i := 0; i < 24; i++ {
		nw := AinvT.MulVec(seeds[i]) // not unit yet
		// support value over verts
		d := -1e300
		for j := 0; j < 24; j++ {
			val := nw.X*c.Verts[j].X + nw.Y*c.Verts[j].Y + nw.Z*c.Verts[j].Z + nw.W*c.Verts[j].W
			if val > d {
				d = val
			}
		}
		L := math.Sqrt(nw.Dot(nw))
		if L > 0 {
			nw = nw.Mul(1 / L)
			d /= L
		}
		c.N[i], c.D[i] = nw, d
	}

	DebugLog("Created 24-cell: center=%+v, side=%.6g, scale=%+v, AABB=[%+v .. %+v]", center, side, scale, c.AABBMin, c.AABBMax)
	return c, nil
}

func intersectRayTwentyFour(O Point4, D Vector4, c *TwentyFourCell) (hit objectHit, ok bool) {
	tEnter, tExit := -1e300, 1e300
	enterIdx, exitIdx := -1, -1
	for i := 0; i < 24; i++ {
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
	return objectHit{t: t, Nw: c.N[idx], c24: c, inv: inv}, true
}
