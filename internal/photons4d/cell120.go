package photons4d

import "fmt"

// 120-cell as intersection of 120 planes whose normals are the 120 vertices of the dual (600-cell).
type Cell120 struct{ cellPoly }

func NewCell120(center Point4, radius Real, scale Vector4, angles Rot4, color, reflectivity, refractivity, ior RGB) (*Cell120, error) {
	if !(radius > 0) {
		return nil, fmt.Errorf("radius must be > 0, got %.6g", radius)
	}
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("scale must be > 0 on all axes, got %+v", scale)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	for _, c := range []struct {
		n    string
		r, t Real
	}{
		{"R", reflectivity.R, refractivity.R},
		{"G", reflectivity.G, refractivity.G},
		{"B", reflectivity.B, refractivity.B},
	} {
		if !in01(c.r) || !in01(c.t) || c.r+c.t > 1+1e-12 {
			return nil, fmt.Errorf("reflect/refract invalid (channel %s)", c.n)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	R := rotFromAngles(angles)
	cp := cellPoly{
		Center: center,
		Scale:  scale,
		R:      R,
		RT:     R.Transpose(),
	}
	cp.materialFrom(color, reflectivity, refractivity, ior)
	cp.buildPlanes(verts600Unit(), radius)

	DebugLog("Created 120-cell: center=%+v, radius=%.6g, scale=%+v, AABB=[%+v..%+v]", center, radius, scale, cp.AABBMin, cp.AABBMax)
	return &Cell120{cp}, nil
}
