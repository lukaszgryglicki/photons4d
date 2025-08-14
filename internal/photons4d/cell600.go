package photons4d

import "fmt"

// 600-cell as intersection of 600 planes whose normals are the 600 vertices of the dual (120-cell).
type Cell600 struct{ cellPoly }

func NewCell600(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Cell600, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("scale must be > 0 on all axes, got %+v", scale)
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
		if !in01(c.r) || !in01(c.t) || !in01(c.d) || c.r+c.t+c.d > 1+1e-12 {
			return nil, fmt.Errorf("diffuse/reflect/refract invalid (channel %s)", c.n)
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
	cp.materialFrom(color, diffuse, reflectivity, refractivity, ior)
	cp.buildPlanes(verts120Unit(), 1.0)

	// DebugLog("Created 600-cell: center=%+v, scale=%+v, AABB=[%+v..%+v]", center, scale, cp.AABBMin, cp.AABBMax)
	DebugLog("Created 600-cell: %+v", cp)
	return &Cell600{cp}, nil
}
