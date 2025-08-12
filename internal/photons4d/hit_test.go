package photons4d

import (
	"math"
	"testing"
)

func TestPlaneHit(t *testing.T) {
	scene := NewScene(Point4{0, 0, 0, 0.25}, 2, 2, 2, 2, 2, 1, 4)

	// Parallel to plane
	tp := planeHit(scene, Point4{0, 0, 0, 0}, Vector4{1, 0, 0, 0})
	if !math.IsInf(float64(tp), 1) {
		t.Fatalf("expected +Inf for parallel ray, got %.6g", tp)
	}

	// Towards plane: t = (Wc - Wo)/Dw
	tok := planeHit(scene, Point4{0, 0, 0, 1}, Vector4{0, 0, 0, -1})
	if math.Abs(float64(tok-0.75)) > 1e-12 {
		t.Fatalf("t wrong: %.12g", tok)
	}

	// Behind (negative or tiny t) => +Inf
	tneg := planeHit(scene, Point4{0, 0, 0, 0.2500000001}, Vector4{0, 0, 0, 1})
	if !math.IsInf(float64(tneg), 1) {
		t.Fatalf("expected +Inf for behind/tiny, got %.6g", tneg)
	}
}
