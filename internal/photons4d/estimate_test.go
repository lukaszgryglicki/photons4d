package photons4d

import (
	"math"
	"math/rand"
	"testing"
	"time"
)

func TestEstimateHitProbAndCastRays(t *testing.T) {
	// deterministic global flags for test
	Debug, PNG, UseLocks = false, false, false

	scene := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 4, 4, 1, 4, false)
	L, err := NewLight(Point4{0, 0, 0, 1}, Vector4{0, 0, 0, -1}, RGB{1, 1, 1}, math.Pi/8, 1.0)
	if err != nil {
		t.Fatal(err)
	}
	// Estimate should be close to 1 for this trivial case
	p := estimateHitProb(L, scene, 1000)
	if p <= 0.5 || p > 1.0000001 {
		t.Fatalf("unexpected hit prob: %.6f", p)
	}

	// castRays should deposit something
	before := 0.0
	for _, x := range scene.Buf {
		before += float64(x)
	}
	// Small number of rays to keep test snappy
	castRays([]*Light{L}, scene, []int{2000})
	after := 0.0
	for _, x := range scene.Buf {
		after += float64(x)
	}
	if after <= before {
		t.Fatalf("no energy deposited: before=%.6g after=%.6g", before, after)
	}

	// Additional sanity: castSingleRay returns deterministically true often.
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	hit := 0
	for i := 0; i < 100; i++ {
		if castSingleRay(L, scene, rng, nil, false) {
			hit++
		}
	}
	if hit < 50 {
		t.Fatalf("too few hits in direct cast: %d/100", hit)
	}
}
