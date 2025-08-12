package photons4d

import "math"

// reflection & refraction in 4D (assume unit I,N; periodic renorm in loop)
func reflect4(I, N Vector4) Vector4 {
	return I.Sub(N.Mul(2 * I.Dot(N)))
}

func refract4(I, N Vector4, eta Real) (Vector4, bool) {
	cosi := -I.Dot(N)
	if cosi < -1 {
		cosi = -1
	} else if cosi > 1 {
		cosi = 1
	}
	k := 1 - eta*eta*(1-cosi*cosi)
	if k < 0 {
		return Vector4{}, false // TIR
	}
	T := I.Mul(eta).Add(N.Mul(eta*cosi - math.Sqrt(k)))
	return T, true
}
