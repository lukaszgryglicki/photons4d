package photons4d

import (
	"math"
	"testing"
)

func TestNewDuotorusAndIntersect(t *testing.T) {
	obj, err := NewDuotorus(
		Point4{0, 0, 0, 0.6},
		Vector4{1, 1, 1, 1},
		0.3,
		0.2,
		0.08,
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewDuotorus: %v", err)
	}

	// Through the tube cross-section centered at local (x=Rxy, z=Rzw), along +X.
	O := Point4{0, 0, 0.2, 0.6}
	D := Vector4{1, 0, 0, 0}
	h, ok := intersectRayDuotorus(O, D, obj)
	if !ok {
		t.Fatalf("expected duotorus hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.22)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.X+1)) > 1e-9 {
		t.Fatalf("expected -X normal, got %+v", h.Nw)
	}

	O2 := Point4{0.3, 0, 0.2, 0.6}
	D2 := Vector4{1, 0, 0, 0}
	h2, ok2 := intersectRayDuotorus(O2, D2, obj)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.08)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.X-1)) > 1e-9 {
		t.Fatalf("expected +X normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithDuotorus(t *testing.T) {
	obj, err := NewDuotorus(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{1.1, 0.9, 1.0, 0.95},
		0.28,
		0.22,
		0.07,
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8},
		RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(obj.AABBMin.X <= obj.AABBMax.X && obj.AABBMin.Y <= obj.AABBMax.Y &&
		obj.AABBMin.Z <= obj.AABBMax.Z && obj.AABBMin.W <= obj.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", obj.AABBMin, obj.AABBMax)
	}

	Rxy := obj.MajorRadiusXY
	Rzw := obj.MajorRadiusZW
	r := obj.MinorRadius
	s := obj.Scale

	localSamples := []Vector4{
		{s.X * (Rxy + r), 0, s.Z * Rzw, 0},
		{s.X * (Rxy - r), 0, s.Z * Rzw, 0},
		{-s.X * (Rxy + r), 0, s.Z * Rzw, 0},
		{0, s.Y * (Rxy + r), s.Z * Rzw, 0},
		{0, s.Y * (Rxy - r), s.Z * Rzw, 0},
		{s.X * Rxy, 0, s.Z * (Rzw + r), 0},
		{s.X * Rxy, 0, s.Z * (Rzw - r), 0},
		{s.X * Rxy, 0, 0, s.W * (Rzw + r)},
		{s.X * Rxy, 0, 0, s.W * (Rzw - r)},
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
