package photons4d

import (
	"math"
	"testing"
)

func TestIntersectRayHyperSphere_AxisCase(t *testing.T) {
	// Unit sphere centered at W=2, identity rotation, equal radii
	hs, err := NewHyperSphere(
		Point4{0, 0, 0, 2},
		Vector4{1, 1, 1, 1}, Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.5, 1.5, 1.5},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	hit, ok := intersectRayHyperSphere(O, D, hs)
	if !ok {
		t.Fatal("expected sphere hit")
	}
	// Solutions along +W: (t-2)^2=1 => t in {1,3}. First positive is 1.
	if math.Abs(float64(hit.t-1)) > 1e-12 {
		t.Fatalf("t wrong: %.12g", hit.t)
	}
	if hit.inv {
		t.Fatal("should be entering (inv=false)")
	}
	// Normal at entry point (0,0,0,1): outward is from center to surface: (0,0,0,-1)
	if math.Abs(float64(hit.Nw.W+1)) > 1e-12 {
		t.Fatalf("normal not -W: %+v", hit.Nw)
	}

	// Start inside: O.W=2, shoot +W => exit at W=3 => t=1
	O2 := Point4{0, 0, 0, 2}
	D2 := Vector4{0, 0, 0, 1}
	hit2, ok2 := intersectRayHyperSphere(O2, D2, hs)
	if !ok2 || !hit2.inv || math.Abs(float64(hit2.t-1)) > 1e-12 {
		t.Fatalf("inside->exit wrong: ok=%v inv=%v t=%.12g", ok2, hit2.inv, hit2.t)
	}
	// Exit normal at (0,0,0,3): (0,0,0,+1)
	if math.Abs(float64(hit2.Nw.W-1)) > 1e-12 {
		t.Fatalf("exit normal not +W: %+v", hit2.Nw)
	}
}
