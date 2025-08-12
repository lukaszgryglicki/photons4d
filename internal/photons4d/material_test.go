package photons4d

import (
	"math/rand"
	"testing"
)

func TestRGBClamp(t *testing.T) {
	c := RGB{-1, 0.5, 2}.clamp01()
	if c != (RGB{0, 0.5, 1}) {
		t.Fatalf("clamp01 wrong: %+v", c)
	}
}

func TestPickChannelThresholds(t *testing.T) {
	// colorSum=3, thrR=1, thrG=2
	rng := rand.New(rand.NewSource(123))
	const N = 60000
	count := [3]int{}
	for i := 0; i < N; i++ {
		c := pickChannel(3, 1, 2, rng)
		count[c]++
	}
	// Expect ~1/3 each; allow 2% tolerance
	for ch := 0; ch < 3; ch++ {
		if frac := float64(count[ch]) / N; frac < 0.31 || frac > 0.36 {
			t.Fatalf("channel %d fraction %.3f out of bounds; counts=%v", ch, frac, count)
		}
	}
}
