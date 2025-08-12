package photons4d

import (
	"math"
	"testing"
)

func TestIsFinite(t *testing.T) {
	if !isFinite(1) || isFinite(math.Inf(1)) || isFinite(math.NaN()) {
		t.Fatal("isFinite failed")
	}
}

func TestIMax(t *testing.T) {
	if imax(3, 5) != 5 || imax(5, 3) != 5 {
		t.Fatal("imax failed")
	}
}
