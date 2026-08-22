package photons4d

import (
	crand "crypto/rand"
	"encoding/binary"
	"fmt"
	"math"
	"sync/atomic"
	"time"
)

func isFinite(x Real) bool { return !math.IsInf(x, 0) && !math.IsNaN(x) }

var rngSeedCounter int64

// freshSeed returns a unique, unpredictable RNG seed. Crypto entropy makes
// seeds independent across machines (two distributed clients starting the
// same nanosecond must never replay the same ray stream); the counter and
// clock guarantee uniqueness across workers even if entropy reads fail.
func freshSeed(wid int) int64 {
	var b [8]byte
	var e int64
	if _, err := crand.Read(b[:]); err == nil {
		e = int64(binary.LittleEndian.Uint64(b[:]))
	}
	c := atomic.AddInt64(&rngSeedCounter, 1)
	seed := e ^ time.Now().UnixNano() ^ int64(uint64(c)*0x9e3779b97f4a7c15) ^ int64(uint64(wid)*0xbf58476d1ce4e5b9)
	if seed == 0 {
		seed = 0x2545f4914f6cdd1d
	}
	return seed
}

func imax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func rmin(a, b Real) Real {
	if a < b {
		return a
	}
	return b
}
func rmax(a, b Real) Real {
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

// logf prints one log line prefixed with a UTC timestamp. All human-readable
// run logs (server, client, progress, output stages) go through it.
func logf(format string, a ...any) {
	fmt.Printf("%s "+format, append([]any{time.Now().UTC().Format("2006-01-02 15:04:05")}, a...)...)
}
