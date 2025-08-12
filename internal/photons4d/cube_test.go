package photons4d

import (
	"math"
	"testing"
)

func newUnitCubeAtW0p5() *HyperCube {
	h, err := NewHyperCube(
		Point4{0, 0, 0, 0.5},
		Vector4{0.2, 0.2, 0.2, 0.2}, // edges 0.2 ⇒ half 0.1
		Rot4{},                      // identity
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		panic(err)
	}
	return h
}

func TestIntersectRayHypercube_EnterExit(t *testing.T) {
	h := newUnitCubeAtW0p5()

	// Outside, coming along +W axis towards cube center
	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	hit, ok := intersectRayHypercube(O, D, h)
	if !ok {
		t.Fatal("expected hit")
	}
	if hit.inv {
		t.Fatal("expected entering, not inv")
	}
	// Cube spans W in [0.4, 0.6]. From W=0 along +W → first hit at t=0.4
	if math.Abs(float64(hit.t-0.4)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", hit.t)
	}
	// Entering face normal should be -W (outward normal of +W face flipped by sign)
	if math.Abs(float64(hit.Nw.W+1)) > 1e-12 {
		t.Fatalf("normal not -W: %+v", hit.Nw)
	}

	// Start inside the cube: O.W=0.5 (center), shoot +W, expect exit on +W face at t=+0.1
	O2 := Point4{0, 0, 0, 0.5}
	D2 := Vector4{0, 0, 0, 1}
	hit2, ok2 := intersectRayHypercube(O2, D2, h)
	if !ok2 {
		t.Fatal("expected hit from inside")
	}
	if !hit2.inv {
		t.Fatal("expected inv (exiting)")
	}
	if math.Abs(float64(hit2.t-0.1)) > 1e-12 {
		t.Fatalf("exit t wrong: %.12g", hit2.t)
	}
	// Exit face normal is +W
	if math.Abs(float64(hit2.Nw.W-1)) > 1e-12 {
		t.Fatalf("exit normal not +W: %+v", hit2.Nw)
	}
}
