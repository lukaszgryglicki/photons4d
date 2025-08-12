package photons4d

import (
	"math"
)

type objectHit struct {
	t   Real
	Nw  Vector4
	hc  *HyperCube
	hs  *HyperSphere
	s5  *Simplex5
	inv bool
}

// ---- helpers to access material caches regardless of shape ----
func (h objectHit) pAbsCh(c int) Real {
	if h.hc != nil {
		return h.hc.pAbs[c]
	}
	if h.hs != nil {
		return h.hs.pAbs[c]
	}
	return h.s5.pAbs[c]
}
func (h objectHit) colorCh(c int) Real {
	if h.hc != nil {
		return h.hc.colorArr[c]
	}
	if h.hs != nil {
		return h.hs.colorArr[c]
	}
	return h.s5.colorArr[c]
}
func (h objectHit) f0Ch(c int) Real {
	if h.hc != nil {
		return h.hc.f0[c]
	}
	if h.hs != nil {
		return h.hs.f0[c]
	}
	return h.s5.f0[c]
}
func (h objectHit) reflCh(c int) Real {
	if h.hc != nil {
		return h.hc.refl[c]
	}
	if h.hs != nil {
		return h.hs.refl[c]
	}
	return h.s5.refl[c]
}
func (h objectHit) refrCh(c int) Real {
	if h.hc != nil {
		return h.hc.refr[c]
	}
	if h.hs != nil {
		return h.hs.refr[c]
	}
	return h.s5.refr[c]
}
func (h objectHit) iorCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorArr[c]
	}
	if h.hs != nil {
		return h.hs.iorArr[c]
	}
	return h.s5.iorArr[c]
}
func (h objectHit) iorInvCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorInv[c]
	}
	if h.hs != nil {
		return h.hs.iorInv[c]
	}
	return h.s5.iorInv[c]
}

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

func nearestHit(scene *Scene, O Point4, D Vector4, tMax Real) (objectHit, bool) {
	best := objectHit{}
	okAny := false
	bestT := tMax
	if !isFinite(bestT) {
		bestT = 1e300
	}

	const eps = 1e-12
	parX := math.Abs(D.X) < eps
	parY := math.Abs(D.Y) < eps
	parZ := math.Abs(D.Z) < eps
	parW := math.Abs(D.W) < eps
	rr := rayRecips{
		parX: parX, parY: parY, parZ: parZ, parW: parW,
	}
	if !parX {
		rr.invX = 1 / D.X
	}
	if !parY {
		rr.invY = 1 / D.Y
	}
	if !parZ {
		rr.invZ = 1 / D.Z
	}
	if !parW {
		rr.invW = 1 / D.W
	}

	// cubes
	for _, h := range scene.Hypercubes {
		if ok, tNear := rayAABB(O, h.AABBMin, h.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayHypercube(O, D, h); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// hyperspheres (ellipsoids)
	for _, s := range scene.Hyperspheres {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayHyperSphere(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// simplexes (5-cell)
	for _, s := range scene.Simplexes {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRaySimplex5(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	return best, okAny
}
