package photons4d

/*
import (
	"math"
	"testing"
)

const eps = 1e-10

func nearly(a, b Real, tol Real) bool { return math.Abs(float64(a-b)) <= float64(tol) }

func vecLen(v Vector4) Real { return math.Sqrt(v.Dot(v)) }

func vecAlmostEq(a, b Vector4, tol Real) bool {
	d := a.Sub(b)
	return vecLen(d) <= tol
}

// tangent component of v w.r.t. unit N
func tangent(v, N Vector4) Vector4 {
	return v.Sub(N.Mul(v.Dot(N)))
}

// build a unit vector T orthogonal to unit N (simple one-step Gram-Schmidt)
func anyTangent(N Vector4) Vector4 {
	e := Vector4{1, 0, 0, 0}
	if math.Abs(N.X) > 0.9 {
		e = Vector4{0, 1, 0, 0}
	}
	T := e.Sub(N.Mul(e.Dot(N)))
	return T.Norm()
}

func TestReflect4_Properties(t *testing.T) {
	normals := []Vector4{
		{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
		{1, 2, 3, 4}.Norm(),
	}
	angles := []Real{math.Pi / 6, math.Pi / 3} // 30°, 60°

	for _, N := range normals {
		for _, a := range angles {
			Tt := anyTangent(N)
			I := N.Mul(-Real(math.Cos(a))).Add(Tt.Mul(Real(math.Sin(a)))) // towards surface

			R := reflect4(I, N)

			// |R| == 1 (inputs are unit)
			if !nearly(vecLen(R), 1, 1e-12) {
				t.Fatalf("reflect length != 1, got %.15g", vecLen(R))
			}
			// normal component flips sign: R·N == - I·N
			if !nearly(R.Dot(N), -I.Dot(N), eps) {
				t.Fatalf("R·N != -I·N: R·N=%.15g I·N=%.15g", R.Dot(N), I.Dot(N))
			}
			// tangential component preserved
			It := tangent(I, N)
			Rt := tangent(R, N)
			if !vecAlmostEq(It, Rt, 1e-12) {
				t.Fatalf("tangent not preserved: |It-Rt|=%.15g", vecLen(It.Sub(Rt)))
			}
		}
	}
}

func TestRefract4_Snell_Entering(t *testing.T) {
	// Entering: n1=1, n2=n>1 ⇒ eta=n1/n2 = 1/n < 1
	n := Real(1.5)
	eta := 1 / n
	N := Vector4{0, 0, 0, 1} // unit
	Tt := anyTangent(N)

	angles := []Real{10 * math.Pi / 180, 40 * math.Pi / 180}
	for _, a := range angles {
		I := N.Mul(-Real(math.Cos(a))).Add(Tt.Mul(Real(math.Sin(a)))) // towards surface
		T, ok := refract4(I, N, eta)
		if !ok {
			t.Fatalf("unexpected TIR on entering at angle %.3f rad", a)
		}
		if !nearly(vecLen(T), 1, 1e-12) {
			t.Fatalf("refract length != 1, got %.15g", vecLen(T))
		}
		// Snell: sin θ_t = η sin θ_i
		cosi := -I.Dot(N)
		sin2i := 1 - cosi*cosi
		sint := math.Sqrt(float64(max(0, eta*eta*sin2i)))
		cost := math.Sqrt(1 - sint*sint)
		// dot(T,N) should be -cos_t (see derivation in analysis)
		if !nearly(-T.Dot(N), Real(cost), 1e-9) {
			t.Fatalf("wrong normal component: -T·N=%.15g, cos_t=%.15g", -T.Dot(N), cost)
		}
		// tangential magnitude scales by η
		It := tangent(I, N)
		Tt2 := tangent(T, N)
		if !nearly(vecLen(Tt2), eta*vecLen(It), 1e-9) {
			t.Fatalf("|T_t| != η|I_t|: |T_t|=%.15g, η|I_t|=%.15g", vecLen(Tt2), eta*vecLen(It))
		}
	}
}

func TestRefract4_Snell_Exiting_NoTIR(t *testing.T) {
	// Exiting: n1=n>1, n2=1 ⇒ eta=n1/n2 = n > 1, choose small angle to avoid TIR.
	eta := Real(1.3)
	N := Vector4{1, 0, 0, 0}
	Tt := anyTangent(N)
	a := Real(10 * math.Pi / 180) // small
	I := N.Mul(-Real(math.Cos(a))).Add(Tt.Mul(Real(math.Sin(a))))
	T, ok := refract4(I, N, eta)
	if !ok {
		t.Fatalf("unexpected TIR when exiting at small angle")
	}
	if !nearly(vecLen(T), 1, 1e-12) {
		t.Fatalf("refract length != 1, got %.15g", vecLen(T))
	}
	cosi := -I.Dot(N)
	sin2i := 1 - cosi*cosi
	k := 1 - eta*eta*sin2i
	if k <= 0 {
		t.Fatalf("expected k>0, got %.15g", k)
	}
	cost := math.Sqrt(float64(k))
	if !nearly(-T.Dot(N), Real(cost), 1e-9) {
		t.Fatalf("wrong normal component (exit): -T·N=%.15g, cos_t=%.15g", -T.Dot(N), cost)
	}
	// tangential scales by η as well
	It := tangent(I, N)
	Tt2 := tangent(T, N)
	if !nearly(vecLen(Tt2), eta*vecLen(It), 1e-9) {
		t.Fatalf("|T_t| != η|I_t| (exit): |T_t|=%.15g, η|I_t|=%.15g", vecLen(Tt2), eta*vecLen(It))
	}
}

func TestRefract4_TIR(t *testing.T) {
	// Total internal reflection when eta>1 and sin θ_i > 1/eta
	eta := Real(1.5)
	N := Vector4{0, 1, 0, 0}
	Tt := anyTangent(N)
	theta := Real(60 * math.Pi / 180) // sin ~ 0.866 > 1/eta ~ 0.666
	I := N.Mul(-Real(math.Cos(theta))).Add(Tt.Mul(Real(math.Sin(theta))))
	_, ok := refract4(I, N, eta)
	if ok {
		t.Fatalf("expected TIR, got transmission")
	}
}

func max(a, b Real) Real {
	if a > b {
		return a
	}
	return b
}
*/
