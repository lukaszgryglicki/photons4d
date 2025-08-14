// internal/photons4d/diffuse_sampling_test.go
package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

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
