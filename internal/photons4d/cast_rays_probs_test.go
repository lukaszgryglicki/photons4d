package photons4d

import "testing"

func approxEqual(a, b Real, tol Real) bool {
	if a > b {
		return a-b <= tol
	}
	return b-a <= tol
}

func TestSplitEventWeights_SumsToAvail(t *testing.T) {
	// refl=0.2, refr=0.3, diff=0.5, Fresnel=0.5, pAbs=0.1
	// rW=0.2*0.5=0.1, tW=0.3*0.5=0.15, dW=0.5
	// sumW=0.75, avail=0.9
	// pR=0.9*(0.1/0.75)=0.12, pT=0.9*(0.15/0.75)=0.18, pD=0.9*(0.5/0.75)=0.6
	pR, pT, pD := splitEventWeights(0.2, 0.3, 0.5, 0.5, 0.1)
	if !approxEqual(pR, 0.12, 1e-12) || !approxEqual(pT, 0.18, 1e-12) || !approxEqual(pD, 0.6, 1e-12) {
		t.Fatalf("unexpected split: pR=%.15g pT=%.15g pD=%.15g", pR, pT, pD)
	}
	sum := pR + pT + pD
	if !approxEqual(sum, 0.9, 1e-12) {
		t.Fatalf("sum(p) != avail: got %.15g want 0.9", sum)
	}
}

func TestSplitEventWeights_EdgeCases(t *testing.T) {
	// No budget (fully absorbing)
	pR, pT, pD := splitEventWeights(0.7, 0.2, 0.1, 0.3, 1.0)
	if pR != 0 || pT != 0 || pD != 0 {
		t.Fatalf("expected zeros for full absorption, got %.15g %.15g %.15g", pR, pT, pD)
	}
	// No weights
	pR, pT, pD = splitEventWeights(0, 0, 0, 0.5, 0.25)
	if pR != 0 || pT != 0 || pD != 0 {
		t.Fatalf("expected zeros for zero weights, got %.15g %.15g %.15g", pR, pT, pD)
	}
	// Pure reflection (after Fresnel)
	pR, pT, pD = splitEventWeights(1, 0, 0, 1, 0.1)
	if !approxEqual(pR, 0.9, 1e-12) || pT != 0 || pD != 0 {
		t.Fatalf("pure reflection split wrong: %.15g %.15g %.15g", pR, pT, pD)
	}
	// Pure refraction (F=0)
	pR, pT, pD = splitEventWeights(1, 1, 0, 0, 0.2)
	if pR != 0 || !approxEqual(pT, 0.8, 1e-12) || pD != 0 {
		t.Fatalf("pure refraction split wrong: %.15g %.15g %.15g", pR, pT, pD)
	}
	// Pure diffuse
	pR, pT, pD = splitEventWeights(0, 0, 1, 0.5, 0.3)
	if pR != 0 || pT != 0 || !approxEqual(pD, 0.7, 1e-12) {
		t.Fatalf("pure diffuse split wrong: %.15g %.15g %.15g", pR, pT, pD)
	}
}

func TestSplitEventWeights_NumericStability(t *testing.T) {
	// Tiny budget, almost-zero weights
	pR, pT, pD := splitEventWeights(1e-12, 2e-12, 3e-12, 0.42, 0.999999999999)
	sum := pR + pT + pD
	avail := 1 - 0.999999999999
	if !approxEqual(sum, avail, 1e-15) {
		t.Fatalf("sum(p) not ~avail: sum=%.18g avail=%.18g", sum, avail)
	}
	// Ensure no negatives due to roundoff
	if pR < 0 || pT < 0 || pD < 0 {
		t.Fatalf("negative probability: pR=%.18g pT=%.18g pD=%.18g", pR, pT, pD)
	}
}
