package photons4d

import (
	"math"
	"testing"
)

func TestNewTwentyFourCellAndIntersect(t *testing.T) {
	c24, err := NewTwentyFourCell(
		Point4{0, 0, 0, 0.6}, // clearly in front of W=0 plane
		0.3,                  // edge length
		Vector4{1, 1, 1, 1},  // no anisotropic scale
		Rot4{},               // identity rotation
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatalf("NewTwentyFourCell: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayTwentyFour(O, D, c24)
	if !ok {
		t.Fatalf("expected 24-cell hit")
	}
	if h.t <= 0 {
		t.Fatalf("t must be positive, got %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.Len()-1)) > 1e-9 {
		t.Fatalf("normal not unit: %.12g", h.Nw.Len())
	}

	// Start inside: shoot from center outward
	O2 := c24.Center
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayTwentyFour(O2, D2, c24)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok2, h2.inv)
	}
}

func TestSceneAABBWithTwentyFourCell(t *testing.T) {
	c24, err := NewTwentyFourCell(
		Point4{-0.15, 0.2, -0.05, 0.35},
		0.22,
		Vector4{1.0, 1.2, 0.9, 0.85},
		Rot4{XY: -0.3, XZ: 0.25, YW: -0.2},
		RGB{1, 1, 1}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.25, 1.3, 1.35},
	)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < len(c24.Verts); i++ {
		v := c24.Verts[i]
		if v.X < c24.AABBMin.X-1e-12 || v.X > c24.AABBMax.X+1e-12 ||
			v.Y < c24.AABBMin.Y-1e-12 || v.Y > c24.AABBMax.Y+1e-12 ||
			v.Z < c24.AABBMin.Z-1e-12 || v.Z > c24.AABBMax.Z+1e-12 ||
			v.W < c24.AABBMin.W-1e-12 || v.W > c24.AABBMax.W+1e-12 {
			t.Fatalf("vertex %d outside AABB", i)
		}
	}
}
