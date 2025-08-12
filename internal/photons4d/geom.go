package photons4d

import "math"

// reflection & refraction in 4D (assume unit I,N; periodic renorm in loop)
func reflect4(I, N Vector4) Vector4 {
	return I.Sub(N.Mul(2 * I.Dot(N)))
}

func refract4(I, N Vector4, eta Real) (Vector4, bool) {
	// Use cosi = cos(angle between -I and N). Start with outward N.
	cosi := -(I.Dot(N))
	n := N
	// If cosi < 0, ray is inside w.r.t outward N -> flip N so cosi >= 0.
	if cosi < 0 {
		n = N.Mul(-1)
		cosi = -cosi
	}
	// Clamp just in case of tiny numeric drift.
	if cosi > 1 {
		cosi = 1
	}
	// Snell: k < 0 => total internal reflection.
	k := 1 - eta*eta*(1-cosi*cosi)
	if k < 0 {
		return Vector4{}, false
	}
	T := I.Mul(eta).Add(n.Mul(eta*cosi - math.Sqrt(k)))
	return T, true
}
