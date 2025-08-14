package photons4d

import "math"

// nearestHitBVH returns the closest positive t hit among all scene objects.
// It uses the BVH (AABB tree) above. tMax can be +Inf to search everything.
func nearestHitBVH(scene *Scene, O Point4, D Vector4, tMax Real) (objectHit, bool) {
	root := getOrBuildBVH(scene)
	// DebugLog("nearestHitBVH: O=%+v, D=%+v, tMax=%f -> %+v", O, D, tMax, root)
	if root == nil {
		return objectHit{}, false
	}
	return traverseNearest(root, O, D, tMax)
}

// (Optional helper): a convenience that searches unbounded.
func nearestHitBVHAll(scene *Scene, O Point4, D Vector4) (objectHit, bool) {
	return nearestHitBVH(scene, O, D, Real(math.Inf(1)))
}
