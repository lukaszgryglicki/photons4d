package photons4d

import (
	"math"
	"testing"
)

func TestNewCell5AndIntersect(t *testing.T) {
	sx, err := NewCell5(
		Point4{0, 0, 0, 0.25},       // slightly in front of scene plane
		Vector4{0.4, 0.4, 0.4, 0.4}, // no anisotropic scale
		Rot4{},                      // identity rotation
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewCell5: %v", err)
	}
	// Simple axis ray along +W from origin
	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayCell5(O, D, sx)
	if !ok {
		t.Fatalf("expected cell5 hit")
	}
	if h.t <= 0 {
		t.Fatalf("t must be positive, got %.12g", h.t)
	}
	// normal should be unit
	if math.Abs(float64(h.Nw.Len()-1)) > 1e-9 {
		t.Fatalf("normal not unit: %.12g", h.Nw.Len())
	}

	// Start inside: shoot from (near) center outward
	O2 := sx.Center
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayCell5(O2, D2, sx)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok2, h2.inv)
	}
}

func TestSceneAABBWithCell5(t *testing.T) {
	sx, err := NewCell5(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.3, 0.225, 0.275, 0.2},
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	// AABB must enclose all vertices
	for i := 0; i < 5; i++ {
		v := sx.Verts[i]
		if v.X < sx.AABBMin.X-1e-12 || v.X > sx.AABBMax.X+1e-12 ||
			v.Y < sx.AABBMin.Y-1e-12 || v.Y > sx.AABBMax.Y+1e-12 ||
			v.Z < sx.AABBMin.Z-1e-12 || v.Z > sx.AABBMax.Z+1e-12 ||
			v.W < sx.AABBMin.W-1e-12 || v.W > sx.AABBMax.W+1e-12 {
			t.Fatalf("vertex %d outside AABB", i)
		}
	}
}
