package photons4d

import (
	"math"
)

type objectHit struct {
	t   Real
	Nw  Vector4
	hc  *HyperCube
	hs  *HyperSphere
	inv bool
}

// ---- helpers to access material caches regardless of shape ----
func (h objectHit) pAbsCh(c int) Real {
	if h.hc != nil {
		return h.hc.pAbs[c]
	}
	return h.hs.pAbs[c]
}
func (h objectHit) colorCh(c int) Real {
	if h.hc != nil {
		return h.hc.colorArr[c]
	}
	return h.hs.colorArr[c]
}
func (h objectHit) f0Ch(c int) Real {
	if h.hc != nil {
		return h.hc.f0[c]
	}
	return h.hs.f0[c]
}
func (h objectHit) reflCh(c int) Real {
	if h.hc != nil {
		return h.hc.refl[c]
	}
	return h.hs.refl[c]
}
func (h objectHit) refrCh(c int) Real {
	if h.hc != nil {
		return h.hc.refr[c]
	}
	return h.hs.refr[c]
}
func (h objectHit) iorCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorArr[c]
	}
	return h.hs.iorArr[c]
}
func (h objectHit) iorInvCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorInv[c]
	}
	return h.hs.iorInv[c]
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
