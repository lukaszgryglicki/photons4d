package photons4d

import (
	"math"
	"sync"
	"testing"
)

func r(v float64) Real              { return Real(v) }
func p4(x, y, z, w float64) Point4  { return Point4{X: r(x), Y: r(y), Z: r(z), W: r(w)} }
func v4(x, y, z, w float64) Vector4 { return Vector4{X: r(x), Y: r(y), Z: r(z), W: r(w)} }

// Build a synthetic leaf whose intersect always returns the provided t.
// NOTE: In tests that traverse, make sure t >= the leaf AABB's entry tmin.
func mkLeaf(min, max Point4, t Real) bvhLeaf {
	return bvhLeaf{
		min: min,
		max: max,
		intersect: func(O Point4, D Vector4) (objectHit, bool) {
			return objectHit{t: t}, true
		},
	}
}

func almostEq(a, b Real) bool { return math.Abs(float64(a-b)) < 1e-9 }

func TestAABBUnion(t *testing.T) {
	aMin, aMax := p4(0, 0, 0, 0), p4(1, 2, 3, 4)
	bMin, bMax := p4(-1, 1, 2, 5), p4(2, 1.5, 3.5, 6)

	uMin, uMax := aabbUnion(aMin, aMax, bMin, bMax)
	if uMin.X != r(-1) || uMin.Y != r(0) || uMin.Z != r(0) || uMin.W != r(0) {
		t.Fatalf("uMin wrong: %+v", uMin)
	}
	if uMax.X != r(2) || uMax.Y != r(2) || uMax.Z != r(3.5) || uMax.W != r(6) {
		t.Fatalf("uMax wrong: %+v", uMax)
	}
}

func TestGetCentroidAxis(t *testing.T) {
	l := bvhLeaf{min: p4(1, 2, 3, 4), max: p4(3, 6, 7, 8)}
	if !almostEq(getCentroidAxis(l, 0), r(2)) {
		t.Errorf("centroid X")
	}
	if !almostEq(getCentroidAxis(l, 1), r(4)) {
		t.Errorf("centroid Y")
	}
	if !almostEq(getCentroidAxis(l, 2), r(5)) {
		t.Errorf("centroid Z")
	}
	if !almostEq(getCentroidAxis(l, 3), r(6)) {
		t.Errorf("centroid W")
	}
}

func TestComputeRayRecips_Basics(t *testing.T) {
	rr := computeRayRecips(v4(0, 1, -2, 0))
	if !rr.parX || !rr.parW {
		t.Fatalf("expected X and W parallel: %+v", rr)
	}
	if rr.parY || rr.parZ {
		t.Fatalf("unexpected parallel flags: %+v", rr)
	}
	if !almostEq(rr.invY, r(1)) || !almostEq(rr.invZ, r(-0.5)) {
		t.Fatalf("wrong recip values: %+v", rr)
	}
}

func TestComputeRayRecips_AllZero(t *testing.T) {
	rr := computeRayRecips(v4(0, 0, 0, 0))
	if !rr.parX || !rr.parY || !rr.parZ || !rr.parW {
		t.Fatalf("all parallel expected: %+v", rr)
	}
	if rr.invX != 0 || rr.invY != 0 || rr.invZ != 0 || rr.invW != 0 {
		t.Fatalf("recips should be zero: %+v", rr)
	}
}

func TestBuildBVH_LeafAndInternalNodes(t *testing.T) {
	// ≤ AABBBVHMaxLeafSize ⇒ single leaf
	makeLeaf := func(i int) bvhLeaf {
		x0 := float64(1 + 3*i)
		// non-overlapping boxes along X so splits are deterministic enough
		return mkLeaf(
			p4(x0, 0, 0, 0),
			p4(x0+1, 1, 1, 1),
			r(float64(i+1)),
		)
	}

	small := make([]bvhLeaf, 0, AABBBVHMaxLeafSize)
	for i := 0; i < AABBBVHMaxLeafSize; i++ {
		small = append(small, makeLeaf(i))
	}
	leaf := buildBVH(small)
	if leaf == nil || leaf.leafObjs == nil || len(leaf.leafObjs) != AABBBVHMaxLeafSize {
		t.Fatalf("expected leaf node of size %d, got %+v", AABBBVHMaxLeafSize, leaf)
	}

	// > AABBBVHMaxLeafSize ⇒ internal split
	large := make([]bvhLeaf, 0, AABBBVHMaxLeafSize+1)
	for i := 0; i < AABBBVHMaxLeafSize+1; i++ {
		large = append(large, makeLeaf(i))
	}
	root := buildBVH(large)
	if root == nil || root.leafObjs != nil || root.left == nil || root.right == nil {
		t.Fatalf("expected internal root with children, got %+v", root)
	}
}

func TestTraverseNearest_ChoosesClosestHit(t *testing.T) {
	// Make hits consistent with each AABB's entry tmin:
	// [1,2] -> t≈1.5 (closest), [4,5] -> t≈4.1, [7,8] -> t≈7.2
	root := buildBVH([]bvhLeaf{
		mkLeaf(p4(1, 0, 0, 0), p4(2, 1, 1, 1), r(1.5)),
		mkLeaf(p4(4, 0, 0, 0), p4(5, 1, 1, 1), r(4.1)),
		mkLeaf(p4(7, 0, 0, 0), p4(8, 1, 1, 1), r(7.2)),
	})

	O := p4(0, 0, 0, 0)
	D := v4(1, 0, 0, 0)
	h, ok := traverseNearest(root, O, D, r(math.Inf(1)))
	if !ok {
		t.Fatalf("expected a hit")
	}
	if !almostEq(h.t, r(1.5)) {
		t.Fatalf("expected t=1.5, got t=%v", h.t)
	}
}

func TestTraverseNearest_RespectsTMax(t *testing.T) {
	root := buildBVH([]bvhLeaf{
		mkLeaf(p4(5, 0, 0, 0), p4(6, 1, 1, 1), r(10)), // beyond tMax
	})
	O := p4(0, 0, 0, 0)
	D := v4(1, 0, 0, 0)
	if _, ok := traverseNearest(root, O, D, r(2)); ok {
		t.Fatalf("expected no hit due to tMax pruning")
	}
}

func TestNearestHitBVH_EmptyScene(t *testing.T) {
	bvhCache = sync.Map{}
	s := &Scene{}
	O := p4(0, 0, 0, 0)
	D := v4(1, 0, 0, 0)
	if _, ok := nearestHitBVH(s, O, D, r(math.Inf(1))); ok {
		t.Fatalf("expected no hit on empty scene")
	}
}

func TestNearestHitBVH_UsesCacheAndReturnsHit(t *testing.T) {
	bvhCache = sync.Map{}
	s := &Scene{}
	root := buildBVH([]bvhLeaf{
		mkLeaf(p4(1, 0, 0, 0), p4(2, 1, 1, 1), r(3)),
		mkLeaf(p4(3, 0, 0, 0), p4(4, 1, 1, 1), r(3.5)),
	})
	bvhCache.Store(s, root)
	O := p4(0, 0, 0, 0)
	D := v4(1, 0, 0, 0)
	h, ok := nearestHitBVH(s, O, D, r(math.Inf(1)))
	if !ok {
		t.Fatalf("expected hit via cached BVH")
	}
	if !almostEq(h.t, r(3)) {
		t.Fatalf("expected t=3, got %v", h.t)
	}
}

func TestCollectSceneObjects_Empty(t *testing.T) {
	s := &Scene{}
	if got := collectSceneObjects(s); len(got) != 0 {
		t.Fatalf("expected 0 objects, got %d", len(got))
	}
}

func TestRayAABB_InsideOrigin_AllowsTraversal(t *testing.T) {
	O := Point4{0, 0, 0, 0}  // inside the box
	D := Vector4{1, 0, 0, 0} // shoot +X
	rr := computeRayRecips(D)
	ok, tmin := rayAABB(O,
		Point4{-1, -1, -1, -1},
		Point4{+1, +1, +1, +1},
		rr,
	)
	if !ok {
		t.Fatalf("expected ok for inside-origin ray")
	}
	if tmin >= 0 {
		t.Logf("rayAABB already clamps tmin>=0 (fine). Keeping or removing the tmin<0 cull is harmless.")
	} else {
		t.Logf("rayAABB returns tmin<0 for inside-origin (typical). You must NOT cull on tmin<0.")
	}
}
