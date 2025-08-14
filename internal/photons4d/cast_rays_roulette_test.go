package photons4d

import "testing"

// helper mirrors the roulette buckets from castSingleRay using splitEventWeights
func decideEvent(refl, refr, diff, F, pAbs Real, u float64) string {
	pR, pT, pD := splitEventWeights(refl, refr, diff, F, pAbs)
	x := Real(u)
	if x < pAbs {
		return "absorb"
	}
	if x < pAbs+pR {
		return "reflect"
	}
	if x < pAbs+pR+pT {
		return "refract"
	}
	if x < pAbs+pR+pT+pD {
		return "diffuse"
	}
	return "none"
}

func mid(a, b Real) float64 { return float64((a + b) * 0.5) }

func TestRouletteBoundaries(t *testing.T) {
	refl, refr, diff := Real(0.3), Real(0.4), Real(0.3)
	F := Real(0.6)    // reasonably reflective
	pAbs := Real(0.1) // some absorption

	pR, pT, pD := splitEventWeights(refl, refr, diff, F, pAbs)
	if pR <= 0 || pT <= 0 || pD <= 0 {
		t.Fatalf("expected positive splits, got pR=%g pT=%g pD=%g", pR, pT, pD)
	}

	// Absorb bucket midpoint
	uAbs := float64(pAbs) * 0.5
	if decideEvent(refl, refr, diff, F, pAbs, uAbs) != "absorb" {
		t.Fatalf("u in absorb bucket misclassified")
	}

	// Reflect bucket midpoint
	uRefl := mid(pAbs, pAbs+pR)
	if decideEvent(refl, refr, diff, F, pAbs, uRefl) != "reflect" {
		t.Fatalf("u in reflect bucket misclassified")
	}

	// Refract bucket midpoint
	uRefr := mid(pAbs+pR, pAbs+pR+pT)
	if decideEvent(refl, refr, diff, F, pAbs, uRefr) != "refract" {
		t.Fatalf("u in refract bucket misclassified")
	}

	// Diffuse bucket midpoint
	uDiff := mid(pAbs+pR+pT, pAbs+pR+pT+pD)
	if decideEvent(refl, refr, diff, F, pAbs, uDiff) != "diffuse" {
		t.Fatalf("u in diffuse bucket misclassified")
	}
}

func TestRoulettePureCases(t *testing.T) {
	// Pure diffuse (all budget to diffuse)
	refl, refr, diff := Real(0), Real(0), Real(1)
	F, pAbs := Real(0.5), Real(0.2)
	pR, pT, pD := splitEventWeights(refl, refr, diff, F, pAbs)
	if pR != 0 || pT != 0 || !approxEqual(pD, 1-pAbs, 1e-12) {
		t.Fatalf("pure diffuse wrong: %g %g %g", pR, pT, pD)
	}
	if decideEvent(refl, refr, diff, F, pAbs, float64(pAbs)+float64(pD)*0.5) != "diffuse" {
		t.Fatalf("pure diffuse selection failed")
	}

	// Pure reflection (F=1, refl>0)
	refl, refr, diff = Real(1), Real(0), Real(0)
	F, pAbs = Real(1), Real(0.05)
	pR, pT, pD = splitEventWeights(refl, refr, diff, F, pAbs)
	if !approxEqual(pR, 1-pAbs, 1e-12) || pT != 0 || pD != 0 {
		t.Fatalf("pure reflection wrong: %g %g %g", pR, pT, pD)
	}
	if decideEvent(refl, refr, diff, F, pAbs, float64(pAbs)+float64(pR)*0.5) != "reflect" {
		t.Fatalf("pure reflection selection failed")
	}
}
