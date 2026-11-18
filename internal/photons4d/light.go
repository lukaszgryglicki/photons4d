package photons4d

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"time"
)

// Light is a 4D "spot" light.
type Light struct {
	Origin    Point4
	Direction Vector4 // unit
	Color     RGB     // clamped to [0,1]
	Angle     Real    // half-angle in radians, (0, π]
	Intensity Real

	// cached
	cosAngle    Real
	oneMinusCos Real
	colorSum    Real
	thrR        Real // cumulative threshold for R
	thrG        Real // cumulative threshold for R+G
	// orthonormal basis for the 3D subspace orthogonal to Direction
	U, V, W Vector4
	// S^3-cap sampling caches
	acceptProb Real // area ratio of cap (rejection acceptance)
	useReject  bool // true => reject uniform S^3 samples, false => sample the polar angle
}

// robust 3D orthonormal basis in the subspace orthogonal to 'a' (unit)
func orthonormal3(a Vector4) (u, v, w Vector4) {
	const eps = 1e-12
	// try up to a few random seeds to avoid degeneracy
	tryBuild := func(seed int64) (Vector4, Vector4, Vector4, bool) {
		rng := rand.New(rand.NewSource(seed))
		// three random candidates
		r1 := Vector4{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		r2 := Vector4{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		r3 := Vector4{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}

		// project each to ⟂ a and Gram–Schmidt
		proj := func(x Vector4) Vector4 { return x.Sub(a.Mul(x.Dot(a))) }

		u := proj(r1)
		lu := u.Len()
		if lu < eps {
			return Vector4{}, Vector4{}, Vector4{}, false
		}
		u = u.Mul(1 / lu)

		v := proj(r2).Sub(u.Mul(r2.Dot(u)))
		lv := v.Len()
		if lv < eps {
			return Vector4{}, Vector4{}, Vector4{}, false
		}
		v = v.Mul(1 / lv)

		w := proj(r3).Sub(u.Mul(r3.Dot(u))).Sub(v.Mul(r3.Dot(v)))
		lw := w.Len()
		if lw < eps {
			return Vector4{}, Vector4{}, Vector4{}, false
		}
		w = w.Mul(1 / lw)

		return u, v, w, true
	}

	seed := time.Now().UnixNano()
	for tries := 0; tries < 8; tries++ {
		if uu, vv, ww, ok := tryBuild(seed + int64(tries)*0x4f1bbcdcbfa53e0a); ok {
			return uu, vv, ww
		}
	}
	// ultra-conservative fallback: deterministic helpers
	h := Vector4{1, 0, 0, 0}
	if math.Abs(a.X) > 0.9 {
		h = Vector4{0, 1, 0, 0}
	}
	u = h.Sub(a.Mul(h.Dot(a))).Norm()

	h2 := Vector4{0, 0, 1, 0}
	v = h2.Sub(a.Mul(h2.Dot(a))).Sub(u.Mul(h2.Dot(u))).Norm()

	h3 := Vector4{0, 0, 0, 1}
	w = h3.Sub(a.Mul(h3.Dot(a))).Sub(u.Mul(h3.Dot(u))).Sub(v.Mul(h3.Dot(v))).Norm()
	return
}

// helper: integral for S^3 cap CDF piece
func jS3(t Real) Real { // J(t) = 1/2 * ( t*sqrt(1-t^2) + asin(t) )
	u := 1 - t*t
	if u < 0 {
		u = 0
	}
	return 0.5 * (t*math.Sqrt(u) + math.Asin(t))
}

// NewLight constructs a cone light and precomputes caches.
func NewLight(origin Point4, dir Vector4, color RGB, angle, intensity Real) (*Light, error) {
	if angle <= 0 || angle > math.Pi {
		return nil, errors.New("angle must be in (0, π]")
	}
	if math.Abs(Real(intensity)) < 1e-9 {
		DebugLog("Light intensity is zero (or very close to zero), setting to 1.0")
		intensity = 1.0
	}
	n := dir.Norm()
	if n.Len() == 0 {
		return nil, errors.New("direction must be non-zero")
	}
	c := color.clamp01()
	csum := c.R + c.G + c.B
	if csum <= 0 {
		return nil, errors.New("colorSum must be positive; got " + fmt.Sprintf("%.6g", csum))
	}
	cosA := math.Cos(angle)

	L := &Light{
		Origin:      origin,
		Direction:   n,
		Color:       c,
		Angle:       angle,
		Intensity:   intensity,
		cosAngle:    cosA,
		oneMinusCos: 1 - cosA,
		colorSum:    csum,
		thrR:        c.R,
		thrG:        c.R + c.G,
	}
	L.U, L.V, L.W = orthonormal3(n)

	Ja := jS3(L.cosAngle)
	L.acceptProb = (math.Pi/4 - Real(Ja)) / (math.Pi / 2)
	// Reject whole-sphere samples if the cap is reasonably wide; sample the
	// polar angle if it is narrow (fewer trials == faster).
	L.useReject = L.acceptProb >= CapRejectThreshold

	DebugLog("Created light %+v", L)
	return L, nil
}

// Fast uniform sample on S^3 without trig/logs (Marsaglia two-disc).
func unitS3(rng *rand.Rand) Vector4 {
	for {
		// first disc
		x1 := 2*rng.Float64() - 1
		x2 := 2*rng.Float64() - 1
		s1 := x1*x1 + x2*x2
		if s1 >= 1 {
			continue
		}
		// second disc
		x3 := 2*rng.Float64() - 1
		x4 := 2*rng.Float64() - 1
		s2 := x3*x3 + x4*x4
		if s2 <= 1e-18 || s2 >= 1 {
			continue
		}
		f := math.Sqrt((1 - s1) / s2)
		// return Vector4{2 * x1, 2 * x2, 2 * x3 * f, 2 * x4 * f}
		return Vector4{x1, x2, x3 * f, x4 * f}
	}
}

// Unit vector on S^2 (correct Marsaglia)
func sampleS2(rng *rand.Rand) (x, y, z Real) {
	for {
		u := 2*rng.Float64() - 1
		v := 2*rng.Float64() - 1
		s := u*u + v*v
		if s > 0 && s < 1 {
			f := 2 * math.Sqrt(1-s)
			return u * f, v * f, 1 - 2*s // already unit
		}
	}
}

// Exact uniform cap on S^3 via rejection.
// Accept if dot(v, axis) >= cos(theta).
func (l *Light) SampleDir(rng *rand.Rand) Vector4 {
	if l.oneMinusCos == 0 {
		return l.Direction
	}
	if l.useReject {
		// Exact uniform cap on S^3 via rejection.
		a := l.Direction
		c := l.cosAngle
		for {
			v := unitS3(rng)
			if a.Dot(v) >= c {
				return v
			}
		}
	}
	// Narrow cone: the polar angle φ has the density sin²φ on [0,θ]. Propose
	// φ = θ·u^(1/3) (density φ²) and accept with (sinφ/φ)²: exact.
	var sinPhi, cosPhi Real
	for {
		phi := l.Angle * math.Cbrt(rng.Float64())
		sinPhi, cosPhi = math.Sincos(phi)
		if phi == 0 {
			break
		}
		if q := sinPhi / phi; rng.Float64() < q*q {
			break
		}
	}
	// Uniform orientation in the 3D subspace orthogonal to axis
	x, y, z := sampleS2(rng)
	ortho := l.U.Mul(x).Add(l.V.Mul(y)).Add(l.W.Mul(z)) // unit, ⟂ axis
	return l.Direction.Mul(cosPhi).Add(ortho.Mul(sinPhi))
}
