package photons4d

import (
	"sync"
)

type Category uint8

const (
	Hit            Category = iota // ray hit
	Miss                           // ray missed
	Diffuse                        // ray scattered diffusely
	Absorb                         // ray absorbed
	Reflect                        // ray reflected
	Refract                        // ray refracted
	TIR                            // total internal reflection (ray did not exit object)
	Escape                         // ray escaped the scene (in hypersphere alpha, beta, gamma mode)
	RecurenceLimit                 // ray hit a recurrence limit (e.g. max bounces exceeded)
)

// RayLogCache accumulates per-event-name counters for DEBUG ray statistics.
// Only counts are retained: raysStats only ever reported counts, and storing a
// full record per ray event grows O(rays×bounces) and OOM-kills large renders.
type RayLogCache struct {
	mu   sync.Mutex
	rays map[string]int64 // map of ray event name to count
}

var cache = &RayLogCache{
	rays: make(map[string]int64),
}

func logRay(name string, _ Category, _ Point4, _ Vector4, _ Point4, _ int, _ Real) {
	cache.mu.Lock()
	cache.rays[name]++
	cache.mu.Unlock()
}

func raysStats() {
	var allRays int64
	for _, n := range cache.rays {
		allRays += n
	}
	if allRays == 0 {
		return
	}
	for k, n := range cache.rays {
		logf("Ray type %s: %d logs (%f%%)\n", k, n, Real(n)/Real(allRays)*100)
	}
}
