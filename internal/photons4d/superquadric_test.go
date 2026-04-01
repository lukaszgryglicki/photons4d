package photons4d

import (
	"math"
	"testing"
)

func TestNewSuperquadricAndIntersect(t *testing.T) {
	obj, err := NewSuperquadric(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.3, 0.25, 0.15},
		4,
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewSuperquadric: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRaySuperquadric(O, D, obj)
	if !ok {
		t.Fatalf("expected superquadric hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.45)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.W+1)) > 1e-9 {
		t.Fatalf("expected -W normal, got %+v", h.Nw)
	}

	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRaySuperquadric(O2, D2, obj)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.15)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithSuperquadric(t *testing.T) {
	power := Real(4)
	obj, err := NewSuperquadric(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		power,
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(obj.AABBMin.X <= obj.AABBMax.X && obj.AABBMin.Y <= obj.AABBMax.Y &&
		obj.AABBMin.Z <= obj.AABBMax.Z && obj.AABBMin.W <= obj.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", obj.AABBMin, obj.AABBMax)
	}

	c := math.Pow(0.25, 1/float64(power))
	localSamples := []Vector4{
		{+obj.Scale.X, 0, 0, 0},
		{-obj.Scale.X, 0, 0, 0},
		{0, +obj.Scale.Y, 0, 0},
		{0, -obj.Scale.Y, 0, 0},
		{0, 0, +obj.Scale.Z, 0},
		{0, 0, -obj.Scale.Z, 0},
		{0, 0, 0, +obj.Scale.W},
		{0, 0, 0, -obj.Scale.W},
		{obj.Scale.X * c, obj.Scale.Y * c, obj.Scale.Z * c, obj.Scale.W * c},
		{-obj.Scale.X * c, obj.Scale.Y * c, -obj.Scale.Z * c, obj.Scale.W * c},
	}
	for i, v := range localSamples {
		p := obj.R.MulVec(v)
		w := Point4{obj.Center.X + p.X, obj.Center.Y + p.Y, obj.Center.Z + p.Z, obj.Center.W + p.W}
		if w.X < obj.AABBMin.X-1e-12 || w.X > obj.AABBMax.X+1e-12 ||
			w.Y < obj.AABBMin.Y-1e-12 || w.Y > obj.AABBMax.Y+1e-12 ||
			w.Z < obj.AABBMin.Z-1e-12 || w.Z > obj.AABBMax.Z+1e-12 ||
			w.W < obj.AABBMin.W-1e-12 || w.W > obj.AABBMax.W+1e-12 {
			t.Fatalf("sample %d outside AABB: %+v not in [%+v .. %+v]", i, w, obj.AABBMin, obj.AABBMax)
		}
	}
}
