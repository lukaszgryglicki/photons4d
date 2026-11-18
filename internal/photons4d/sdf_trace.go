package photons4d

import "math"

// rayDist is one sample of a signed distance along a ray, taken in a space
// where the ray moves with a constant speed:
//   - d is a lower bound of the distance to the surface with the sign of the
//     side (negative inside); it is zero only on the surface,
//   - rate is its derivative along the ray parameter,
//   - smooth is the radius of the ball around the sample in which the second
//     derivatives of d stay below 1/(distance to the border of that ball);
//     zero when unknown, which leaves only the first-order step |d|/speed.
type rayDist struct {
	d, rate, smooth Real
}

const traceMaxSteps = 1 << 17

// traceFirstCrossing walks from t0 to t1 and returns the first parameter at
// which the distance changes its sign, together with the side of t0. A step
// never passes a point of the surface: it is the larger of |d|/speed and the
// root of the quadratic bound |d| + m*s - k*s*s/2 with k = 2*speed²/smooth,
// which holds for s <= smooth/(2*speed). Steps below the floor are replaced by
// the floor, so features thinner than 1e-9 of the interval may be passed; the
// floor grows on a ray that creeps along the surface for thousands of steps.
func traceFirstCrossing(dist func(Real) rayDist, speed, t0, t1 Real) (Real, bool, bool) {
	cur := dist(t0)
	inside := cur.d <= 0
	if !(t1 > t0) || !(speed > 0) {
		return 0, inside, false
	}
	floor := 1e-9*(t1-t0) + 1e-13*math.Max(math.Abs(t0), math.Abs(t1))
	t := t0
	for i := 0; i < traceMaxSteps; i++ {
		if i&4095 == 4095 {
			floor *= 4
		}
		ad := math.Abs(cur.d)
		step := ad / speed
		if cur.smooth > 0 {
			k := 2 * speed * speed / cur.smooth
			m := cur.rate
			if cur.d < 0 {
				m = -m
			}
			s := (m + math.Sqrt(m*m+2*k*ad)) / k
			if lim := cur.smooth / (2 * speed); s > lim {
				s = lim
			}
			if s > step {
				step = s
			}
		}
		if !(step >= floor) {
			step = floor
		}
		tn := t + step
		last := tn >= t1
		if last {
			tn = t1
		}
		next := dist(tn)
		if (next.d <= 0) != inside {
			return refineCrossing(dist, t, cur.d, tn, next.d, inside), inside, true
		}
		if last {
			return 0, inside, false
		}
		t, cur = tn, next
	}
	return 0, inside, false
}

// refineCrossing narrows a bracket whose ends lie on different sides (a on
// the side given by inside) and returns its end beyond the surface.
func refineCrossing(dist func(Real) rayDist, a, fa, b, fb Real, inside bool) Real {
	side := 0
	for i := 0; i < 200 && b-a > 1e-15*math.Max(math.Abs(a), math.Abs(b)); i++ {
		m := (a*fb - b*fa) / (fb - fa)
		if !(m > a && m < b) {
			m = 0.5 * (a + b)
			if !(m > a && m < b) {
				break
			}
		}
		fm := dist(m).d
		if (fm <= 0) == inside {
			a, fa = m, fm
			if side < 0 {
				fb *= 0.5
			}
			side = -1
		} else {
			b, fb = m, fm
			if side > 0 {
				fa *= 0.5
			}
			side = 1
		}
	}
	return b
}

// traceSolid returns the first crossing of the surface at a positive ray
// parameter and whether the ray starts inside the solid. [t0, t1] is the part
// of the ray inside a bound of the solid; it is widened a little, so that a
// surface which touches the bound is not met exactly at an end.
func traceSolid(dist func(Real) rayDist, speed, t0, t1 Real) (Real, bool, bool) {
	const eps = 1e-12
	if !(t1 >= t0) {
		return 0, false, false
	}
	pad := 1e-6*(t1-t0) + eps
	t0 -= pad
	t1 += pad
	if t0 < 0 {
		t0 = 0
	}
	for i := 0; i < 4; i++ {
		t, inside, ok := traceFirstCrossing(dist, speed, t0, t1)
		if !ok {
			return 0, false, false
		}
		if t > eps {
			return t, inside, true
		}
		t0 = t
	}
	return 0, false, false
}
