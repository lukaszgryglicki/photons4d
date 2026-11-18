package photons4d

import (
	"math"
	"testing"
)

// Unit ball around the origin met by the ray O + t*D (D of length speed).
func ballDist(O, D Vector4, exact bool) func(Real) rayDist {
	return func(t Real) rayDist {
		p := O.Add(D.Mul(t))
		r := p.Len()
		if !exact {
			// only a bound of the distance, without derivatives
			return rayDist{d: 0.3 * (r - 1)}
		}
		return rayDist{d: r - 1, rate: p.Dot(D) / r, smooth: r}
	}
}

func TestTraceFirstCrossing_Ball(t *testing.T) {
	D := Vector4{0, 0, 0, 2}
	for _, exact := range []bool{true, false} {
		// from outside: enters at w=-1
		tt, inside, ok := traceFirstCrossing(ballDist(Vector4{0, 0, 0, -3}, D, exact), 2, 0, 10)
		if !ok || inside || math.Abs(tt-1) > 1e-12 {
			t.Fatalf("exact=%v entry: ok=%v inside=%v t=%.15g, want t=1", exact, ok, inside, tt)
		}
		// from inside: leaves at w=+1
		tt, inside, ok = traceFirstCrossing(ballDist(Vector4{0.6, 0, 0, 0}, D, exact), 2, 0, 10)
		if !ok || !inside || math.Abs(tt-0.4) > 1e-12 {
			t.Fatalf("exact=%v exit: ok=%v inside=%v t=%.15g, want t=0.4", exact, ok, inside, tt)
		}
		// passes the ball at the distance 1e-9: no crossing
		if _, _, ok = traceFirstCrossing(ballDist(Vector4{1 + 1e-9, 0, 0, -3}, D, exact), 2, 0, 10); ok {
			t.Fatalf("exact=%v: crossing reported for a ray that passes the ball", exact)
		}
		// cuts a chord of the length 2e-4 near the rim: its entry is not passed over
		h := math.Sqrt(1 - 1e-8)
		tt, inside, ok = traceFirstCrossing(ballDist(Vector4{h, 0, 0, -3}, D, exact), 2, 0, 10)
		if want := (3 - 1e-4) / 2; !ok || inside || math.Abs(tt-want) > 1e-9 {
			t.Fatalf("exact=%v thin chord: ok=%v inside=%v t=%.15g, want %.15g", exact, ok, inside, tt, want)
		}
	}
}

func TestTraceSolid_StartOnTheSurface(t *testing.T) {
	D := Vector4{0, 0, 0, 1}
	// 1e-13 outside the ball and going in: the solid is left at the far side
	tt, inside, ok := traceSolid(ballDist(Vector4{0, 0, 0, -1 - 1e-13}, D, true), 1, -5, 5)
	if !ok || !inside || math.Abs(tt-2) > 1e-9 {
		t.Fatalf("going in: ok=%v inside=%v t=%.15g, want the exit at t=2", ok, inside, tt)
	}
	// 1e-13 inside and going out: nothing ahead
	if tt, _, ok = traceSolid(ballDist(Vector4{0, 0, 0, 1 - 1e-13}, D, true), 1, -5, 5); ok {
		t.Fatalf("going out: crossing reported at t=%.15g", tt)
	}
	// the solid is behind the ray
	if _, _, ok = traceSolid(ballDist(Vector4{0, 0, 0, 3}, D, true), 1, -4, -2); ok {
		t.Fatalf("solid behind the ray: crossing reported")
	}
}
