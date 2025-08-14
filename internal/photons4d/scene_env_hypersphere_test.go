package photons4d

import (
	"math"
	"testing"
)

// expectedBins replicates the hyperspherical mapping used by AngleIndexOf.
// Parameterization on S^3 (unit 4D direction):
//
//	x = cos(α)
//	y = sin(α) cos(β)
//	z = sin(α) sin(β) cos(γ)
//	w = sin(α) sin(β) sin(γ)
//
// with α∈[0,π], β∈[0,π], γ∈[0,2π).
func expectedBins(D Vector4, nx, ny, nz int) (i, j, k int) {
	// Normalize
	n := math.Sqrt(D.X*D.X + D.Y*D.Y + D.Z*D.Z + D.W*D.W)
	if n == 0 {
		panic("expectedBins: zero direction")
	}
	x := D.X / n
	y := D.Y / n
	z := D.Z / n
	w := D.W / n

	// Recover angles (robust to edge cases).
	alpha := math.Acos(clamp(x, -1, 1)) // [0,π]
	s1 := math.Sin(alpha)

	var beta, gamma float64
	if math.Abs(s1) < 1e-15 {
		// α≈0 or α≈π ⇒ direction ~±X; set β=0, γ=0 canonically.
		beta = 0
		gamma = 0
	} else {
		yy := clamp(y/s1, -1, 1)
		beta = math.Acos(yy) // [0,π]
		s2 := math.Sin(beta)
		if math.Abs(s2) < 1e-15 {
			// β≈0 or β≈π ⇒ direction in XY-plane along ±Y; set γ=0 canonically.
			gamma = 0
		} else {
			// z = s1*s2*cos(γ), w = s1*s2*sin(γ)
			gamma = math.Atan2(w, z) // (-π,π]
			if gamma < 0 {
				gamma += 2 * math.Pi
			}
		}
	}

	// Map to bins.
	fi := int(alpha / math.Pi * float64(nx))
	fj := int(beta / math.Pi * float64(ny))
	fk := int(gamma / (2 * math.Pi) * float64(nz))

	// Upper-edge clamp (e.g., α==π → nx).
	if fi == nx {
		fi = nx - 1
	}
	if fj == ny {
		fj = ny - 1
	}
	if fk == nz {
		fk = nz - 1
	}
	return fi, fj, fk
}

func clamp(x, a, b float64) float64 {
	if x < a {
		return a
	}
	if x > b {
		return b
	}
	return x
}

func TestAngleIndexOf_CardinalAxes(t *testing.T) {
	// Angular grid 8x8x8 so midpoints & quarters land exactly.
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 8, 8, 8, 4, true)

	tests := []struct {
		name    string
		D       Vector4
		wantIJK [3]int
	}{
		// +X: α=0; expect (0,0,0).
		{"+X", Vector4{X: 1}, [3]int{0, 0, 0}},
		// -X: α=π; expect (nx-1,0,0) = (7,0,0).
		{"-X", Vector4{X: -1}, [3]int{7, 0, 0}},
		// +Y: α=π/2, β=0; expect (4,0,0).
		{"+Y", Vector4{Y: 1}, [3]int{4, 0, 0}},
		// -Y: α=π/2, β=π; expect (4,7,0).
		{"-Y", Vector4{Y: -1}, [3]int{4, 7, 0}},
		// +Z: α=π/2, β=π/2, γ=0; expect (4,4,0).
		{"+Z", Vector4{Z: 1}, [3]int{4, 4, 0}},
		// -Z: α=π/2, β=π/2, γ=π; expect (4,4,4).
		{"-Z", Vector4{Z: -1}, [3]int{4, 4, 4}},
		// +W: α=π/2, β=π/2, γ=π/2; expect (4,4,2).
		{"+W", Vector4{W: 1}, [3]int{4, 4, 2}},
		// -W: α=π/2, β=π/2, γ=3π/2; expect (4,4,6).
		{"-W", Vector4{W: -1}, [3]int{4, 4, 6}},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			i, j, k, _, _, _ := s.AngleIndexOf(tc.D)
			if got := [3]int{i, j, k}; got != tc.wantIJK {
				t.Fatalf("AngleIndexOf(%s) = %v, want %v", tc.name, got, tc.wantIJK)
			}
			// Cross-check against explicit math to ensure mapping consistency.
			ei, ej, ek := expectedBins(tc.D, s.Nx, s.Ny, s.Nz)
			if i != ei || j != ej || k != ek {
				t.Fatalf("explicit mapping disagrees: AngleIndexOf=%d,%d,%d expectedBins=%d,%d,%d",
					i, j, k, ei, ej, ek)
			}
		})
	}
}

func TestAngleIndexOf_ScaleInvariance(t *testing.T) {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 8, 8, 8, 4, true)

	vecs := []Vector4{
		{X: 1, Y: 2, Z: 3, W: 4},
		{X: -2, Y: 1, Z: 0.5, W: -0.25},
		{X: 0, Y: 1, Z: -1, W: 1},
	}
	scales := []float64{0.5, 1, 7.25}

	for vi, v := range vecs {
		i0, j0, k0, _, _, _ := s.AngleIndexOf(v)
		for si, a := range scales {
			i, j, k, _, _, _ := s.AngleIndexOf(Vector4{X: v.X * a, Y: v.Y * a, Z: v.Z * a, W: v.W * a})
			if i != i0 || j != j0 || k != k0 {
				t.Fatalf("scale invariance failed (vec #%d, scale #%d): got %d,%d,%d want %d,%d,%d",
					vi, si, i, j, k, i0, j0, k0)
			}
		}
	}
}

func TestAngleIndexOf_BoundsAndEdges(t *testing.T) {
	// Very coarse grid to exercise clamping on α≈π etc.
	s := NewScene(Point4{0, 0, 0, 0}, 1, 1, 1, 3, 3, 6, 1, true)

	// α→π (direction → -X) should clamp to i=nx-1.
	i, j, k, _, _, _ := s.AngleIndexOf(Vector4{X: -1, Y: 1e-16, Z: 0, W: 0})
	if i != s.Nx-1 {
		t.Fatalf("alpha upper-edge clamp failed: i=%d nx-1=%d", i, s.Nx-1)
	}
	// β→0 (direction → +Y) should map j=0 and γ canonical to 0.
	i, j, k, _, _, _ = s.AngleIndexOf(Vector4{Y: 1, Z: 1e-18, W: 1e-18})
	if j != 0 || k != 0 {
		t.Fatalf("beta≈0 canonicalization failed: j=%d k=%d (want 0,0)", j, k)
	}
	// γ wraparound: γ≈2π should map to k=0.
	// Achieve γ≈2π by picking (z>0, w≈0+) with α=β=π/2.
	i, j, k, _, _, _ = s.AngleIndexOf(Vector4{Z: 1, W: 1e-18})
	if k != 0 {
		t.Fatalf("gamma wraparound failed: k=%d (want 0)", k)
	}
}
