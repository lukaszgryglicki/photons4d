package photons4d

import (
	"math/rand"
)

// RGB stores color components; each should be in [0,1].
type RGB struct {
	R, G, B Real
}

type material interface {
	PAbsCh(int) Real
	DiffCh(int) Real
	ColorCh(int) Real
	F0Ch(int) Real
	ReflCh(int) Real
	RefrCh(int) Real
	IORCh(int) Real
	IORInvCh(int) Real
}

// clamp01 clamps each channel to [0,1].
func (c RGB) clamp01() RGB {
	cl := func(x Real) Real {
		if x < 0 {
			return 0
		}
		if x > 1 {
			return 1
		}
		return x
	}
	return RGB{cl(c.R), cl(c.G), cl(c.B)}
}

// fast channel pick with precomputed thresholds
func pickChannel(colorSum, thrR, thrG Real, rng *rand.Rand) int {
	u := rng.Float64() * colorSum
	if u < thrR {
		return ChR
	}
	if u < thrG {
		return ChG
	}
	return ChB
}
