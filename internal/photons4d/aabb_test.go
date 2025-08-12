package photons4d

import "testing"

func TestRayAABB_HitAndMiss(t *testing.T) {
	min := Point4{-1, -1, -1, -1}
	max := Point4{1, 1, 1, 1}

	O := Point4{-2, 0, 0, 0}
	D := Vector4{1, 0, 0, 0}
	rr := rayRecips{invX: 1 / D.X, invY: 0, invZ: 0, invW: 0, parY: true, parZ: true, parW: true}
	ok, tnear := rayAABB(O, min, max, rr)
	if !ok || tnear < 1 || tnear > 3 {
		t.Fatalf("expected hit with tnear in [1,3], got ok=%v tnear=%.6g", ok, tnear)
	}

	// Ray parallel outside
	O2 := Point4{-2, 2, 0, 0}
	D2 := Vector4{1, 0, 0, 0}
	rr2 := rayRecips{invX: 1 / D2.X, parY: true, parZ: true, parW: true}
	ok2, _ := rayAABB(O2, min, max, rr2)
	if ok2 {
		t.Fatal("expected no hit for parallel outside slab")
	}
}
