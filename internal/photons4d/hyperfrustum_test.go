package photons4d

import (
	"math"
	"testing"
)

func TestNewHyperFrustumAndIntersect(t *testing.T) {
	f, err := NewHyperFrustum(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.2, 0.2, 0.15},
		0.5,
		1.0,
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewHyperFrustum: %v", err)
	}

	// Off-axis so we hit the side, not the lower cap.
	O := Point4{0.15, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayHyperFrustum(O, D, f)
	if !ok {
		t.Fatalf("expected hyperfrustum hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.6)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.Len()-1)) > 1e-9 {
		t.Fatalf("normal not unit: %.12g", h.Nw.Len())
	}
	if h.Nw.X <= 0 || h.Nw.W >= 0 {
		t.Fatalf("unexpected side normal orientation: %+v", h.Nw)
	}

	// Start inside and exit through the upper cap.
	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayHyperFrustum(O2, D2, f)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.15)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W cap normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithHyperFrustum(t *testing.T) {
	f, err := NewHyperFrustum(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		0.6,
		1.1,
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(f.AABBMin.X <= f.AABBMax.X && f.AABBMin.Y <= f.AABBMax.Y &&
		f.AABBMin.Z <= f.AABBMax.Z && f.AABBMin.W <= f.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", f.AABBMin, f.AABBMax)
	}

	lower := f.LowerRadius
	upper := f.UpperRadius
	localSamples := []Vector4{
		{0, 0, 0, -f.Scale.W},
		{0, 0, 0, +f.Scale.W},
		{+lower * f.Scale.X, 0, 0, -f.Scale.W},
		{-lower * f.Scale.X, 0, 0, -f.Scale.W},
		{0, +lower * f.Scale.Y, 0, -f.Scale.W},
		{0, 0, +lower * f.Scale.Z, -f.Scale.W},
		{+upper * f.Scale.X, 0, 0, +f.Scale.W},
		{-upper * f.Scale.X, 0, 0, +f.Scale.W},
		{0, +upper * f.Scale.Y, 0, +f.Scale.W},
		{0, 0, +upper * f.Scale.Z, +f.Scale.W},
	}
	for i, v := range localSamples {
		p := f.R.MulVec(v)
		w := Point4{f.Center.X + p.X, f.Center.Y + p.Y, f.Center.Z + p.Z, f.Center.W + p.W}
		if w.X < f.AABBMin.X-1e-12 || w.X > f.AABBMax.X+1e-12 ||
			w.Y < f.AABBMin.Y-1e-12 || w.Y > f.AABBMax.Y+1e-12 ||
			w.Z < f.AABBMin.Z-1e-12 || w.Z > f.AABBMax.Z+1e-12 ||
			w.W < f.AABBMin.W-1e-12 || w.W > f.AABBMax.W+1e-12 {
			t.Fatalf("sample %d outside AABB: %+v not in [%+v .. %+v]", i, w, f.AABBMin, f.AABBMax)
		}
	}
}
