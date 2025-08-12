package photons4d

import (
	"math"
)

type objectHit struct {
	t    Real
	Nw   Vector4
	hc   *Cell8
	hs   *HyperSphere
	s5   *Cell5
	c16  *Cell16
	c24  *TwentyFourCell
	poly *cellPoly // for 120/600 cells
	inv  bool
}

func (h objectHit) pAbsCh(c int) Real {
	if h.hc != nil {
		return h.hc.pAbs[c]
	}
	if h.hs != nil {
		return h.hs.pAbs[c]
	}
	if h.s5 != nil {
		return h.s5.pAbs[c]
	}
	if h.c16 != nil {
		return h.c16.pAbs[c]
	}
	if h.c24 != nil {
		return h.c24.pAbs[c]
	}
	return h.poly.pAbs[c]
}
func (h objectHit) colorCh(c int) Real {
	if h.hc != nil {
		return h.hc.colorArr[c]
	}
	if h.hs != nil {
		return h.hs.colorArr[c]
	}
	if h.s5 != nil {
		return h.s5.colorArr[c]
	}
	if h.c16 != nil {
		return h.c16.colorArr[c]
	}
	if h.c24 != nil {
		return h.c24.colorArr[c]
	}
	return h.poly.colorArr[c]
}
func (h objectHit) f0Ch(c int) Real {
	if h.hc != nil {
		return h.hc.f0[c]
	}
	if h.hs != nil {
		return h.hs.f0[c]
	}
	if h.s5 != nil {
		return h.s5.f0[c]
	}
	if h.c16 != nil {
		return h.c16.f0[c]
	}
	if h.c24 != nil {
		return h.c24.f0[c]
	}
	return h.poly.f0[c]
}
func (h objectHit) reflCh(c int) Real {
	if h.hc != nil {
		return h.hc.refl[c]
	}
	if h.hs != nil {
		return h.hs.refl[c]
	}
	if h.s5 != nil {
		return h.s5.refl[c]
	}
	if h.c16 != nil {
		return h.c16.refl[c]
	}
	if h.c24 != nil {
		return h.c24.refl[c]
	}
	return h.poly.refl[c]
}
func (h objectHit) refrCh(c int) Real {
	if h.hc != nil {
		return h.hc.refr[c]
	}
	if h.hs != nil {
		return h.hs.refr[c]
	}
	if h.s5 != nil {
		return h.s5.refr[c]
	}
	if h.c16 != nil {
		return h.c16.refr[c]
	}
	if h.c24 != nil {
		return h.c24.refr[c]
	}
	return h.poly.refr[c]
}
func (h objectHit) iorCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorArr[c]
	}
	if h.hs != nil {
		return h.hs.iorArr[c]
	}
	if h.s5 != nil {
		return h.s5.iorArr[c]
	}
	if h.c16 != nil {
		return h.c16.iorArr[c]
	}
	if h.c24 != nil {
		return h.c24.iorArr[c]
	}
	return h.poly.iorArr[c]
}
func (h objectHit) iorInvCh(c int) Real {
	if h.hc != nil {
		return h.hc.iorInv[c]
	}
	if h.hs != nil {
		return h.hs.iorInv[c]
	}
	if h.s5 != nil {
		return h.s5.iorInv[c]
	}
	if h.c16 != nil {
		return h.c16.iorInv[c]
	}
	if h.c24 != nil {
		return h.c24.iorInv[c]
	}
	return h.poly.iorInv[c]
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

	// cells8
	for _, h := range scene.Cells8 {
		if ok, tNear := rayAABB(O, h.AABBMin, h.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCell8(O, D, h); ok && hit.t > 1e-12 && hit.t < bestT {
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

	// (5-cell)
	for _, s := range scene.Cells5 {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCell5(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 16-cells
	for _, s := range scene.Cells16 {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCell16(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 24-cells
	for _, s := range scene.Cells24 {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayTwentyFour(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 120-cell
	for _, c := range scene.Cells120 {
		if ok, tNear := rayAABB(O, c.AABBMin, c.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCellPoly(O, D, &c.cellPoly); ok && hit.t > 1e-12 && hit.t < bestT {
			hit.poly = &c.cellPoly
			bestT, best, okAny = hit.t, hit, true
		}
	}
	// 600-cell
	for _, c := range scene.Cells600 {
		if ok, tNear := rayAABB(O, c.AABBMin, c.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCellPoly(O, D, &c.cellPoly); ok && hit.t > 1e-12 && hit.t < bestT {
			hit.poly = &c.cellPoly
			bestT, best, okAny = hit.t, hit, true
		}
	}

	return best, okAny
}
