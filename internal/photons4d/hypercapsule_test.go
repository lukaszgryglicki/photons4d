package photons4d

import (
	"math"
	"testing"
)

func TestNewHyperCapsuleAndIntersect(t *testing.T) {
	c, err := NewHyperCapsule(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.2, 0.2, 0.1},
		0.15,
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewHyperCapsule: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayHyperCapsule(O, D, c)
	if !ok {
		t.Fatalf("expected hypercapsule hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.35)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.W+1)) > 1e-9 {
		t.Fatalf("expected -W normal at bottom tip, got %+v", h.Nw)
	}

	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayHyperCapsule(O2, D2, c)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.25)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W normal at top tip, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithHyperCapsule(t *testing.T) {
	c, err := NewHyperCapsule(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.12},
		0.16,
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(c.AABBMin.X <= c.AABBMax.X && c.AABBMin.Y <= c.AABBMax.Y &&
		c.AABBMin.Z <= c.AABBMax.Z && c.AABBMin.W <= c.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", c.AABBMin, c.AABBMax)
	}

	localSamples := []Vector4{
		{0, 0, 0, -(c.HalfLength + c.Scale.W)},
		{0, 0, 0, +(c.HalfLength + c.Scale.W)},
		{+c.Scale.X, 0, 0, 0},
		{-c.Scale.X, 0, 0, 0},
		{0, +c.Scale.Y, 0, 0},
		{0, -c.Scale.Y, 0, 0},
		{0, 0, +c.Scale.Z, 0},
		{0, 0, -c.Scale.Z, 0},
	}
	for i, v := range localSamples {
		p := c.R.MulVec(v)
		w := Point4{c.Center.X + p.X, c.Center.Y + p.Y, c.Center.Z + p.Z, c.Center.W + p.W}
		if w.X < c.AABBMin.X-1e-12 || w.X > c.AABBMax.X+1e-12 ||
			w.Y < c.AABBMin.Y-1e-12 || w.Y > c.AABBMax.Y+1e-12 ||
			w.Z < c.AABBMin.Z-1e-12 || w.Z > c.AABBMax.Z+1e-12 ||
			w.W < c.AABBMin.W-1e-12 || w.W > c.AABBMax.W+1e-12 {
			t.Fatalf("sample %d outside AABB: %+v not in [%+v .. %+v]", i, w, c.AABBMin, c.AABBMax)
		}
	}
}
