package photons4d

import (
	"math"
)

func isFinite(x Real) bool { return !math.IsInf(x, 0) && !math.IsNaN(x) }

func imax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// splitEventWeights computes the dynamic (per-hit) probabilities for
// reflection, refraction and diffuse given:
//   - refl, refr, diff: per-channel material knobs (0..1)
//   - F: Fresnel factor in [0,1] (Schlick)
//   - pAbs: absorption probability in [0,1]
//
// The remaining budget "avail = 1 - pAbs" is split proportionally between:
//
//	rW = refl * F
//	tW = refr * (1 - F)
//	dW = diff
//
// Returns (pReflDyn, pRefrDyn, pDiffDyn), which sum to at most avail (with
// small numeric tolerance). If avail <= 0 or all weights are 0, returns zeros.
func splitEventWeights(refl, refr, diff, F, pAbs Real) (Real, Real, Real) {
	avail := 1 - pAbs
	if avail <= 0 {
		return 0, 0, 0
	}
	rW := refl * F
	tW := refr * (1 - F)
	dW := diff
	sumW := rW + tW + dW
	if sumW <= 0 {
		return 0, 0, 0
	}
	pRefl := avail * (rW / sumW)
	pRefr := avail * (tW / sumW)
	pDiff := avail - pRefl - pRefr // keep numerically consistent
	if pDiff < 0 {
		pDiff = 0
	}
	return pRefl, pRefr, pDiff
}
