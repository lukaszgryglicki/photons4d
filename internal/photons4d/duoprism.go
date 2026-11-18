package photons4d

import (
	"fmt"
	"math"
)

// Duoprism: the Cartesian product of a regular P-gon in the local XY plane
// and a regular Q-gon in the local ZW plane, both with circumradius 1 before
// Scale is applied (so Scale.X/Y stretch the first polygon and Scale.Z/W the
// second one). It has P·Q vertices and P + Q cells: P prisms over the Q-gon
// and Q prisms over the P-gon. The 4-4 duoprism is the 8-cell, and for large
// P and Q the shape tends to the duocylinder of the same Scale.
type Duoprism struct {
	cellPoly
	P, Q int
}

// Largest polygon accepted for a duoprism (beyond it use a duocylinder).
const maxDuoprismSides = 360

func NewDuoprism(center Point4, scale Vector4, p, q int, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Duoprism, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("duoprism scale must be > 0 on all axes, got %+v", scale)
	}
	if p < 3 || q < 3 || p > maxDuoprismSides || q > maxDuoprismSides {
		return nil, fmt.Errorf("duoprism needs polygons with 3..%d sides, got p=%d q=%d", maxDuoprismSides, p, q)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	for _, c := range []struct {
		n       string
		r, t, d Real
	}{
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

	// The side k of a polygon with n sides lies between the vertices at the
	// angles 2πk/n and 2π(k+1)/n; its normal is scaled so that the side is
	// the line n·x = 1.
	normals := make([]Vector4, 0, p+q)
	for k := 0; k < p; k++ {
		a := (2*Real(k) + 1) * math.Pi / Real(p)
		l := 1 / math.Cos(math.Pi/Real(p))
		normals = append(normals, Vector4{l * math.Cos(a), l * math.Sin(a), 0, 0})
	}
	for k := 0; k < q; k++ {
		a := (2*Real(k) + 1) * math.Pi / Real(q)
		l := 1 / math.Cos(math.Pi/Real(q))
		normals = append(normals, Vector4{0, 0, l * math.Cos(a), l * math.Sin(a)})
	}
	verts := make([]Vector4, 0, p*q)
	for i := 0; i < p; i++ {
		a := 2 * math.Pi * Real(i) / Real(p)
		for j := 0; j < q; j++ {
			b := 2 * math.Pi * Real(j) / Real(q)
			verts = append(verts, Vector4{math.Cos(a), math.Sin(a), math.Cos(b), math.Sin(b)})
		}
	}

	R := rotFromAngles(angles)
	cp := cellPoly{
		Center: center,
		Scale:  scale,
		R:      R,
		RT:     R.Transpose(),
	}
	cp.materialFrom(color, diffuse, reflectivity, refractivity, ior)
	cp.buildPlanes(normals, 1.0)
	cp.computeAABBFromLocalVerts(verts)

	DebugLog("Created %d-%d duoprism: %+v", p, q, cp)
	return &Duoprism{cellPoly: cp, P: p, Q: q}, nil
}
