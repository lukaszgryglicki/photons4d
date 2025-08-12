package photons4d

import "testing"

func TestRayLogCache(t *testing.T) {
	// reset
	cache = &RayLogCache{rays: make(map[string][]RayLog)}
	logRay("foo", Hit, Point4{}, Vector4{}, Point4{}, 0, 1)
	logRay("foo", Miss, Point4{}, Vector4{}, Point4{}, 1, 2)
	logRay("bar", Reflect, Point4{}, Vector4{}, Point4{}, 2, 3)
	if len(cache.rays["foo"]) != 2 || len(cache.rays["bar"]) != 1 {
		t.Fatalf("unexpected cache sizes: %+v", cache.rays)
	}
}
