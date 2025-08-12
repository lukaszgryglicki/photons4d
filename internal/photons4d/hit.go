package photons4d

import "math"

type cubeHit struct {
	t   Real    // param distance along ray
	Nw  Vector4 // world-space unit normal at hit
	hc  *HyperCube
	inv bool // true if we were inside and are exiting
}

// plane distance to W = scene.Center.W (from O along D), +Inf if none
func planeHit(scene *Scene, O Point4, D Vector4) Real {
	if math.Abs(D.W) < 1e-12 {
		return math.Inf(1)
	}
	t := (scene.Center.W - O.W) / D.W
	if t <= 1e-9 {
		return math.Inf(1)
	}
	return t
}
