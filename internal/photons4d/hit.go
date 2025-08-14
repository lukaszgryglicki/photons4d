package photons4d

import (
	"math"
)

type objectHit struct {
	t   Real
	Nw  Vector4
	mat material // the only thing you need for channels
	inv bool
}

func (h objectHit) pAbsCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.PAbsCh(c)
}
func (h objectHit) diffCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.DiffCh(c)
}
func (h objectHit) colorCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.ColorCh(c)
}
func (h objectHit) f0Ch(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.F0Ch(c)
}
func (h objectHit) reflCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.ReflCh(c)
}
func (h objectHit) refrCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.RefrCh(c)
}
func (h objectHit) iorCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.IORCh(c)
}
func (h objectHit) iorInvCh(c int) Real {
	if h.mat == nil {
		return 0
	}
	return h.mat.IORInvCh(c)
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

	// hyperspheres (ellipsoids)
	for _, s := range scene.Hyperspheres {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayHyperSphere(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 5-cells
	for _, s := range scene.Cells5 {
		if ok, tNear := rayAABB(O, s.AABBMin, s.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCell5(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 8-cells
	for _, h := range scene.Cells8 {
		if ok, tNear := rayAABB(O, h.AABBMin, h.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCell8(O, D, h); ok && hit.t > 1e-12 && hit.t < bestT {
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
		if hit, ok := intersectRayCell24(O, D, s); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	// 120-cell
	for _, c := range scene.Cells120 {
		if ok, tNear := rayAABB(O, c.AABBMin, c.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCellPoly(O, D, &c.cellPoly); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}
	// 600-cell
	for _, c := range scene.Cells600 {
		if ok, tNear := rayAABB(O, c.AABBMin, c.AABBMax, rr); !ok || tNear > bestT {
			continue
		}
		if hit, ok := intersectRayCellPoly(O, D, &c.cellPoly); ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}

	return best, okAny
}
