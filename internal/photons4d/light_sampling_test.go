package photons4d

import (
	"math"
	"math/rand"
	"sort"
	"testing"
)

// KS statistic for a continuous target CDF F on sorted samples xs.
func ksD(xs []float64, F func(float64) float64) float64 {
	sort.Float64s(xs)
	n := len(xs)
	var d float64
	for i, x := range xs {
		Fi := F(x)
		empUpper := float64(i+1) / float64(n)
		empLower := float64(i) / float64(n)
		di := math.Max(Fi-empLower, empUpper-Fi)
		if di > d {
			d = di
		}
	}
	return d
}

// Build a non-axis-aligned unit direction for robustness.
func testAxis() Vector4 {
	return Vector4{0.3, -0.7, 0.2, 0.6}.Norm()
}

// CDF for uniform-on-cap in S^3, expressed in t = cos(phi).
func capCDF(cosTheta Real) func(float64) float64 {
	J1 := math.Pi / 4
	Ja := float64(jS3(cosTheta))
	den := J1 - Ja
	return func(t float64) float64 {
		if t <= float64(cosTheta) {
			return 0
		}
		if t >= 1 {
			return 1
		}
		Jt := float64(jS3(Real(t)))
		return (Jt - Ja) / den
	}
}

// Sanity checks: unit length and within the cap.
func checkGeom(t *testing.T, v Vector4, a Vector4, cosTheta Real) {
	// unit
	l2 := v.Dot(v)
	if math.Abs(float64(l2-1)) > 1e-9 {
		t.Fatalf("non-unit sample, |v|^2=%g", l2)
	}
	// in-cap
	d := a.Dot(v)
	if float64(d) < float64(cosTheta)-1e-12 {
		t.Fatalf("sample outside cap: a·v=%g, cosθ=%g", d, cosTheta)
	}
}

func runOneKSTest(t *testing.T, angle Real, n int, seed int64) {
	a := testAxis()
	L, err := NewLight(Point4{}, a, RGB{1, 0.7, 0.4}, angle, 1.0)
	if err != nil {
		t.Fatalf("NewLight failed: %v", err)
	}
	rng := rand.New(rand.NewSource(seed))
	cosTheta := math.Cos(angle)
	F := capCDF(cosTheta)

	ts := make([]float64, n)
	for i := 0; i < n; i++ {
		v := L.SampleDir(rng)
		checkGeom(t, v, a, cosTheta)
		ts[i] = float64(a.Dot(v)) // t = cos(phi)
	}

	D := ksD(ts, F)
	crit := 1.36 / math.Sqrt(float64(n)) // α≈0.05
	if D > crit {
		t.Fatalf("KS failed for angle=%.6g rad: D=%.6g > crit=%.6g (n=%d)",
			angle, D, crit, n)
	}
}

func TestLightSampleDir_UniformS3Cap_Wide(t *testing.T) {
	// Wide cone → exercises rejection path too
	runOneKSTest(t, 1.2, 50000, 12345)
}

func TestLightSampleDir_UniformS3Cap_Narrow(t *testing.T) {
	// Narrow cone → exercises LUT inverse-CDF path
	runOneKSTest(t, 0.05, 50000, 67890)
}
