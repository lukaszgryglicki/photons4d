package photons4d

import "math"

// reflection & refraction in 4D (assume unit I,N; periodic renorm in loop)
func reflect4(I, N Vector4) Vector4 {
	return I.Sub(N.Mul(2 * I.Dot(N)))
}

// Refraction with side awareness.
// Contract: eta must be n1/n2 for the *current* interface:
//   - outside → inside : eta = 1/ior
//   - inside  → outside: eta = ior
//
// I and N must be unit; N is the outward normal.
func refract4(I, N Vector4, eta Real) (Vector4, bool) {
	// If we're inside (I·N > 0), flip the normal so cosθ is positive.
	n := N
	cosi := I.Dot(N)
	if cosi > 0 {
		n = N.Mul(-1)
		cosi = -cosi
	} else {
		cosi = -cosi
	}
	// Numeric clamp to [0,1] to avoid tiny negatives/overs.
	if cosi < 0 {
		cosi = 0
	} else if cosi > 1 {
		cosi = 1
	}
	k := 1 - eta*eta*(1-cosi*cosi)
	if k < 0 {
		return Vector4{}, false // total internal reflection
	}
	T := I.Mul(eta).Add(n.Mul(eta*cosi - math.Sqrt(k)))
	return T, true
}
