// internal/photons4d/diffuse_sampling_test.go
package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestOrthonormal3Fast_NearAxis(t *testing.T) {
	N := Vector4{0.9999999999995, 1e-6, 0, 0}.Norm()
	U, V, W := orthonormal3Fast(N)
	B := []Vector4{U, V, W}
	for i := 0; i < 3; i++ {
		if math.Abs(float64(B[i].Len()-1)) > 1e-9 {
			t.Fatalf("basis[%d] not unit: %.12g", i, B[i].Len())
		}
		if math.Abs(float64(B[i].Dot(N))) > 1e-9 {
			t.Fatalf("basis[%d] not orthogonal to N: %.12g", i, B[i].Dot(N))
		}
		for j := 0; j < i; j++ {
			if math.Abs(float64(B[i].Dot(B[j]))) > 1e-9 {
				t.Fatalf("basis[%d] not orthogonal to basis[%d]: %.12g", i, j, B[i].Dot(B[j]))
			}
		}
	}
}

func TestSampleDiffuseDir_IsCosineWeightedS3(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	N := Vector4{0, 0, 0, 1} // unit normal
	const M = 200000
	var sum Real
	for i := 0; i < M; i++ {
		d := sampleDiffuseDir(N, rng)
		if a := d.Dot(d); math.Abs(float64(a-1)) > 1e-9 {
			t.Fatalf("direction not unit: |d|^2=%g", a)
		}
		if dn := d.Dot(N); dn < -1e-12 {
			t.Fatalf("direction went below hemisphere: d·N=%g", dn)
		}
		sum += d.Dot(N)
	}
	mean := float64(sum) / M
	want := 3 * math.Pi / 16 // ≈ 0.5890486
	if math.Abs(mean-want) > 0.01 {
		t.Fatalf("mean(d·N)=%g, want≈%g", mean, want)
	}
}

func TestSampleDiffuseDir_IsCosineWeightedS3_NearAxis(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	N := Vector4{0.9999999999995, 1e-6, 0, 0}.Norm()
	const M = 200000
	var sum Real
	for i := 0; i < M; i++ {
		d := sampleDiffuseDir(N, rng)
		if a := d.Dot(d); math.Abs(float64(a-1)) > 1e-9 {
			t.Fatalf("direction not unit: |d|^2=%g", a)
		}
		if dn := d.Dot(N); dn < -1e-12 {
			t.Fatalf("direction went below hemisphere: d·N=%g", dn)
		}
		sum += d.Dot(N)
	}
	mean := float64(sum) / M
	want := 3 * math.Pi / 16
	if math.Abs(mean-want) > 0.01 {
		t.Fatalf("mean(d·N)=%g, want≈%g", mean, want)
	}
}
