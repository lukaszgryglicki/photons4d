package photons4d

import (
	"math"
	"testing"
)

func TestNewCell16AndIntersect(t *testing.T) {
	c16, err := NewCell16(
		Point4{0, 0, 0, 0.6}, // clearly in front of W=0 plane
		0.3,                  // edge length
		Vector4{1, 1, 1, 1},  // no anisotropic scale
		Rot4{},               // identity rotation
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatalf("NewCell16: %v", err)
	}

	// Simple axis ray along +W from origin
	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayCell16(O, D, c16)
	if !ok {
		t.Fatalf("expected 16-cell hit")
	}
	if h.t <= 0 {
		t.Fatalf("t must be positive, got %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.Len()-1)) > 1e-9 {
		t.Fatalf("normal not unit: %.12g", h.Nw.Len())
	}

	// Start inside: shoot from center outward
	O2 := c16.Center
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayCell16(O2, D2, c16)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok2, h2.inv)
	}
}

func TestSceneAABBWithCell16(t *testing.T) {
	c16, err := NewCell16(
		Point4{0.1, -0.2, 0.05, 0.4},
		0.25,
		Vector4{1.1, 0.9, 1.0, 0.8},
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < len(c16.Verts); i++ {
		v := c16.Verts[i]
		if v.X < c16.AABBMin.X-1e-12 || v.X > c16.AABBMax.X+1e-12 ||
			v.Y < c16.AABBMin.Y-1e-12 || v.Y > c16.AABBMax.Y+1e-12 ||
			v.Z < c16.AABBMin.Z-1e-12 || v.Z > c16.AABBMax.Z+1e-12 ||
			v.W < c16.AABBMin.W-1e-12 || v.W > c16.AABBMax.W+1e-12 {
			t.Fatalf("vertex %d outside AABB", i)
		}
	}
}
