package photons4d

import (
	"math"
	"testing"
)

func TestNewDitorusAndIntersect(t *testing.T) {
	obj, err := NewDitorus(
		Point4{0, 0, 0, 0.6},
		Vector4{1, 1, 1, 1},
		0.3, 0.12, 0.04,
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewDitorus: %v", err)
	}

	// On the X axis of its own hyperplane the solid is |x| in [0.14, 0.22] and [0.38, 0.46].
	cases := []struct {
		name string
		O    Point4
		D    Vector4
		t    Real
		inv  bool
		N    Vector4
	}{
		{"outer wall, entering", Point4{-1, 0, 0, 0.6}, Vector4{1, 0, 0, 0}, 0.54, false, Vector4{-1, 0, 0, 0}},
		{"outer wall, leaving", Point4{-0.42, 0, 0, 0.6}, Vector4{1, 0, 0, 0}, 0.04, true, Vector4{1, 0, 0, 0}},
		{"inner wall, entering", Point4{-0.3, 0, 0, 0.6}, Vector4{1, 0, 0, 0}, 0.08, false, Vector4{-1, 0, 0, 0}},
		{"inner wall, leaving", Point4{-0.18, 0, 0, 0.6}, Vector4{1, 0, 0, 0}, 0.04, true, Vector4{1, 0, 0, 0}},
		{"across the hole, entering", Point4{0, 0, 0, 0.6}, Vector4{0, 1, 0, 0}, 0.14, false, Vector4{0, -1, 0, 0}},
		{"along W, entering", Point4{0.42, 0, 0, 0}, Vector4{0, 0, 0, 1}, 0.56, false, Vector4{0, 0, 0, -1}},
		{"along W, leaving", Point4{0.42, 0, 0, 0.6}, Vector4{0, 0, 0, 1}, 0.04, true, Vector4{0, 0, 0, 1}},
		{"along Z over the ring, entering", Point4{0.3, 0, -1, 0.6}, Vector4{0, 0, 1, 0}, 0.84, false, Vector4{0, 0, -1, 0}},
	}
	for _, c := range cases {
		h, ok := intersectRayDitorus(c.O, c.D, obj)
		if !ok {
			t.Fatalf("%s: expected a hit", c.name)
		}
		if h.inv != c.inv || math.Abs(float64(h.t-c.t)) > 1e-9 || h.Nw.Sub(c.N).Len() > 1e-9 {
			t.Fatalf("%s: got t=%.12g inv=%v N=%+v, want t=%.12g inv=%v N=%+v", c.name, h.t, h.inv, h.Nw, c.t, c.inv, c.N)
		}
	}

	// through the middle of the hole along W, and between the two walls
	for _, O := range []Point4{{0, 0, 0, 0}, {0.3, 0, 0, 0}} {
		if _, ok := intersectRayDitorus(O, Vector4{0, 0, 0, 1}, obj); ok {
			t.Fatalf("ray from %+v along W must miss", O)
		}
	}

	for _, bad := range [][3]Real{{0.3, 0.12, 0.12}, {0.3, 0.2, 0.1}, {0.3, 0.12, 0}} {
		if _, err := NewDitorus(Point4{}, Vector4{1, 1, 1, 1}, bad[0], bad[1], bad[2], Rot4{}, RGB{1, 1, 1}, RGB{}, RGB{}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2}); err == nil {
			t.Fatalf("radii %v accepted", bad)
		}
	}
}

func TestSceneAABBWithDitorus(t *testing.T) {
	obj, err := NewDitorus(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{1.1, 0.9, 1.05, 0.8},
		0.3, 0.11, 0.045,
		Rot4{XY: 0.2, XZ: -0.15, XW: 0.4, YW: -0.5, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8},
		RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}

	// The box must hold the surface and touch it on every side.
	lo := [4]Real{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
	hi := [4]Real{math.Inf(-1), math.Inf(-1), math.Inf(-1), math.Inf(-1)}
	const n = 90
	for i := 0; i < n; i++ {
		a := 2 * math.Pi * Real(i) / n
		for j := 0; j < n; j++ {
			b := 2 * math.Pi * Real(j) / n
			for k := 0; k < n; k++ {
				g := 2 * math.Pi * Real(k) / n
				tube := obj.MiddleRadius + obj.MinorRadius*math.Cos(g)
				ring := obj.MajorRadius + tube*math.Cos(b)
				p := obj.R.MulVec(Vector4{
					ring * math.Cos(a) * obj.Scale.X,
					ring * math.Sin(a) * obj.Scale.Y,
					tube * math.Sin(b) * obj.Scale.Z,
					obj.MinorRadius * math.Sin(g) * obj.Scale.W,
				})
				w := [4]Real{obj.Center.X + p.X, obj.Center.Y + p.Y, obj.Center.Z + p.Z, obj.Center.W + p.W}
				for c := 0; c < 4; c++ {
					lo[c] = math.Min(lo[c], w[c])
					hi[c] = math.Max(hi[c], w[c])
				}
			}
		}
	}
	min := [4]Real{obj.AABBMin.X, obj.AABBMin.Y, obj.AABBMin.Z, obj.AABBMin.W}
	max := [4]Real{obj.AABBMax.X, obj.AABBMax.Y, obj.AABBMax.Z, obj.AABBMax.W}
	for c := 0; c < 4; c++ {
		if lo[c] < min[c]-1e-12 || hi[c] > max[c]+1e-12 {
			t.Fatalf("axis %d: surface [%g, %g] outside AABB [%g, %g]", c, lo[c], hi[c], min[c], max[c])
		}
		if lo[c]-min[c] > 2e-3 || max[c]-hi[c] > 2e-3 {
			t.Fatalf("axis %d: AABB [%g, %g] is not tight around [%g, %g]", c, min[c], max[c], lo[c], hi[c])
		}
	}
}
