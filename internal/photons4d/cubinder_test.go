package photons4d

import (
	"math"
	"testing"
)

func TestNewCubinderAndIntersect(t *testing.T) {
	obj, err := NewCubinder(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.3, 0.25, 0.15},
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewCubinder: %v", err)
	}

	cases := []struct {
		name string
		O    Point4
		D    Vector4
		t    Real
		inv  bool
		N    Vector4
	}{
		{"W face, entering", Point4{0, 0, 0, 0}, Vector4{0, 0, 0, 1}, 0.45, false, Vector4{0, 0, 0, -1}},
		{"W face, leaving", Point4{0, 0, 0, 0.6}, Vector4{0, 0, 0, 1}, 0.15, true, Vector4{0, 0, 0, 1}},
		{"Z face, entering", Point4{0, 0, -1, 0.6}, Vector4{0, 0, 1, 0}, 0.75, false, Vector4{0, 0, -1, 0}},
		{"Z face, leaving", Point4{0.1, 0.1, 0, 0.6}, Vector4{0, 0, -1, 0}, 0.25, true, Vector4{0, 0, -1, 0}},
		{"round side along X, entering", Point4{-1, 0, 0, 0.6}, Vector4{1, 0, 0, 0}, 0.8, false, Vector4{-1, 0, 0, 0}},
		{"round side along Y, leaving", Point4{0, 0, 0.1, 0.65}, Vector4{0, 1, 0, 0}, 0.3, true, Vector4{0, 1, 0, 0}},
	}
	for _, c := range cases {
		h, ok := intersectRayCubinder(c.O, c.D, obj)
		if !ok {
			t.Fatalf("%s: expected a hit", c.name)
		}
		if h.inv != c.inv || math.Abs(float64(h.t-c.t)) > 1e-9 || h.Nw.Sub(c.N).Len() > 1e-9 {
			t.Fatalf("%s: got t=%.12g inv=%v N=%+v, want t=%.12g inv=%v N=%+v", c.name, h.t, h.inv, h.Nw, c.t, c.inv, c.N)
		}
	}

	// beside the round side, inside both slabs
	if _, ok := intersectRayCubinder(Point4{0.21, 0, 0, 0}, Vector4{0, 0, 0, 1}, obj); ok {
		t.Fatalf("ray beside the round side must miss")
	}
	// inside the ellipse, beside the Z slab
	if _, ok := intersectRayCubinder(Point4{0, 0, 0.26, 0}, Vector4{0, 0, 0, 1}, obj); ok {
		t.Fatalf("ray beside the Z faces must miss")
	}

	if _, err := NewCubinder(Point4{}, Vector4{1, 1, 0, 1}, Rot4{}, RGB{1, 1, 1}, RGB{}, RGB{}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2}); err == nil {
		t.Fatalf("zero scale accepted")
	}
}

func TestSceneAABBWithCubinder(t *testing.T) {
	obj, err := NewCubinder(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		Rot4{XY: 0.2, XZ: -0.15, XW: 0.4, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8},
		RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}

	// The box must hold the boundary and touch it on every side.
	lo := [4]Real{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
	hi := [4]Real{math.Inf(-1), math.Inf(-1), math.Inf(-1), math.Inf(-1)}
	rx, ry, hz, hw := obj.Scale.X, obj.Scale.Y, obj.Scale.Z, obj.Scale.W
	for i := 0; i < 3600; i++ {
		a := 2 * math.Pi * Real(i) / 3600
		for _, sz := range []Real{-hz, hz} {
			for _, sw := range []Real{-hw, hw} {
				p := obj.R.MulVec(Vector4{rx * math.Cos(a), ry * math.Sin(a), sz, sw})
				w := [4]Real{obj.Center.X + p.X, obj.Center.Y + p.Y, obj.Center.Z + p.Z, obj.Center.W + p.W}
				for k := 0; k < 4; k++ {
					lo[k] = math.Min(lo[k], w[k])
					hi[k] = math.Max(hi[k], w[k])
				}
			}
		}
	}
	min := [4]Real{obj.AABBMin.X, obj.AABBMin.Y, obj.AABBMin.Z, obj.AABBMin.W}
	max := [4]Real{obj.AABBMax.X, obj.AABBMax.Y, obj.AABBMax.Z, obj.AABBMax.W}
	for k := 0; k < 4; k++ {
		if lo[k] < min[k]-1e-12 || hi[k] > max[k]+1e-12 {
			t.Fatalf("axis %d: boundary [%g, %g] outside AABB [%g, %g]", k, lo[k], hi[k], min[k], max[k])
		}
		if lo[k]-min[k] > 1e-6 || max[k]-hi[k] > 1e-6 {
			t.Fatalf("axis %d: AABB [%g, %g] is not tight around [%g, %g]", k, min[k], max[k], lo[k], hi[k])
		}
	}
}
