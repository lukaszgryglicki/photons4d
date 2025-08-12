package photons4d

import (
	"math"
)

func isFinite(x Real) bool { return !math.IsInf(x, 0) && !math.IsNaN(x) }
