package photons4d

import (
	"math"
	"testing"
)

func TestNewSpherinderAndIntersect(t *testing.T) {
	s, err := NewSpherinder(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.2, 0.2, 0.15},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewSpherinder: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRaySpherinder(O, D, s)
	if !ok {
		t.Fatalf("expected spherinder hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.45)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.W+1)) > 1e-9 {
		t.Fatalf("expected -W cap normal, got %+v", h.Nw)
	}

	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRaySpherinder(O2, D2, s)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.15)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W cap normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithSpherinder(t *testing.T) {
	s, err := NewSpherinder(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(s.AABBMin.X <= s.AABBMax.X && s.AABBMin.Y <= s.AABBMax.Y && s.AABBMin.Z <= s.AABBMax.Z && s.AABBMin.W <= s.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", s.AABBMin, s.AABBMax)
	}

	samples := []Point4{
		{s.Center.X + s.R.M[0][0]*s.Scale.X, s.Center.Y + s.R.M[1][0]*s.Scale.X, s.Center.Z + s.R.M[2][0]*s.Scale.X, s.Center.W + s.R.M[3][0]*s.Scale.X},
		{s.Center.X - s.R.M[0][0]*s.Scale.X, s.Center.Y - s.R.M[1][0]*s.Scale.X, s.Center.Z - s.R.M[2][0]*s.Scale.X, s.Center.W - s.R.M[3][0]*s.Scale.X},
		{s.Center.X + s.R.M[0][1]*s.Scale.Y, s.Center.Y + s.R.M[1][1]*s.Scale.Y, s.Center.Z + s.R.M[2][1]*s.Scale.Y, s.Center.W + s.R.M[3][1]*s.Scale.Y},
		{s.Center.X - s.R.M[0][1]*s.Scale.Y, s.Center.Y - s.R.M[1][1]*s.Scale.Y, s.Center.Z - s.R.M[2][1]*s.Scale.Y, s.Center.W - s.R.M[3][1]*s.Scale.Y},
		{s.Center.X + s.R.M[0][2]*s.Scale.Z, s.Center.Y + s.R.M[1][2]*s.Scale.Z, s.Center.Z + s.R.M[2][2]*s.Scale.Z, s.Center.W + s.R.M[3][2]*s.Scale.Z},
		{s.Center.X - s.R.M[0][2]*s.Scale.Z, s.Center.Y - s.R.M[1][2]*s.Scale.Z, s.Center.Z - s.R.M[2][2]*s.Scale.Z, s.Center.W - s.R.M[3][2]*s.Scale.Z},
		{s.Center.X + s.R.M[0][3]*s.Scale.W, s.Center.Y + s.R.M[1][3]*s.Scale.W, s.Center.Z + s.R.M[2][3]*s.Scale.W, s.Center.W + s.R.M[3][3]*s.Scale.W},
		{s.Center.X - s.R.M[0][3]*s.Scale.W, s.Center.Y - s.R.M[1][3]*s.Scale.W, s.Center.Z - s.R.M[2][3]*s.Scale.W, s.Center.W - s.R.M[3][3]*s.Scale.W},
	}
	for i, p := range samples {
		if p.X < s.AABBMin.X-1e-12 || p.X > s.AABBMax.X+1e-12 ||
			p.Y < s.AABBMin.Y-1e-12 || p.Y > s.AABBMax.Y+1e-12 ||
			p.Z < s.AABBMin.Z-1e-12 || p.Z > s.AABBMax.Z+1e-12 ||
			p.W < s.AABBMin.W-1e-12 || p.W > s.AABBMax.W+1e-12 {
			t.Fatalf("support sample %d outside AABB: %+v not in [%+v .. %+v]", i, p, s.AABBMin, s.AABBMax)
		}
	}
}
