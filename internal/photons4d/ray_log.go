package photons4d

import (
	"fmt"
	"sync"
)

type Category uint8

const (
	Hit            Category = iota // ray hit
	Miss                           // ray missed
	Absorb                         // ray absorbed
	Reflect                        // ray reflected
	Refract                        // ray refracted
	TIR                            // total internal reflection (ray did not exit object)
	RecurenceLimit                 // ray hit a recurrence limit (e.g. max bounces exceeded)
)

type RayLog struct {
	Name      string
	Category  Category
	Origin    Point4
	Direction Vector4
	Point     Point4 // hit point, if any
	Bounce    int    // bounce number (0 for first ray)
	Distance  Real   // distance traveled by the ray
}

type RayLogCache struct {
	mu   sync.Mutex
	rays map[string][]RayLog // map of ray name to logs
}

var cache = &RayLogCache{
	rays: make(map[string][]RayLog),
}

func logRay(name string, category Category, origin Point4, direction Vector4, point Point4, bounce int, distance Real) {
	cache.mu.Lock()
	defer cache.mu.Unlock()
	cache.rays[name] = append(cache.rays[name], RayLog{
		Name:      name,
		Category:  category,
		Origin:    origin,
		Direction: direction,
		Point:     point,
		Bounce:    bounce,
		Distance:  distance,
	})
}

func raysStats() {
	for k, v := range cache.rays {
		fmt.Printf("Ray type %s: %d logs\n", k, len(v))
		//for _, log := range v {
		//	fmt.Printf("  Bounce %d: Category=%d, Origin=%+v, Direction=%+v, Point=%+v, Distance=%.6f\n",
		//		log.Bounce, log.Category, log.Origin, log.Direction, log.Point, log.Distance)
		//}
	}
}
