package photons4d

import (
	"math"
	"sync"
	"testing"
)

func TestDuotorus_BVHConsistency(t *testing.T) {
	bvhCache = sync.Map{}

	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 4, 4, 1, 8, false)
	obj, err := NewDuotorus(
		Point4{0, 0, 0, 0.6},
		Vector4{1, 1, 1, 1},
		0.3,
		0.2,
		0.08,
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatal(err)
	}
	s.AddDuotorus(obj)

	got := collectSceneObjects(s)
	if len(got) != 1 {
		t.Fatalf("expected 1 collected object, got %d", len(got))
	}

	O := Point4{0, 0, 0.2, 0.6}
	D := Vector4{1, 0, 0, 0}
	h1, ok1 := nearestHit(s, O, D, Real(math.Inf(1)))
	h2, ok2 := nearestHitBVH(s, O, D, Real(math.Inf(1)))
	if ok1 != ok2 {
		t.Fatalf("nearestHit and BVH disagree on hit existence: %v vs %v", ok1, ok2)
	}
	if !ok1 {
		t.Fatalf("expected a hit")
	}
	if math.Abs(float64(h1.t-h2.t)) > 1e-12 {
		t.Fatalf("t mismatch: nearestHit=%g BVH=%g", h1.t, h2.t)
	}
}
