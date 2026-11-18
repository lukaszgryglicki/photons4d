package photons4d

import (
	"math/rand"
	"testing"
)

// copy of intersectRayKleinSurface with counters
func zzKlein(O Point4, D Vector4, k *KleinSurface, cheap, exact, lower *int) (objectHit, bool) {
	Ol := k.RT.MulVec(Vector4{O.X - k.Center.X, O.Y - k.Center.Y, O.Z - k.Center.Z, O.W - k.Center.W})
	Dl := k.RT.MulVec(D)
	A := Vector4{Ol.X / k.Scale.X, Ol.Y / k.Scale.Y, Ol.Z / k.Scale.Z, Ol.W / k.Scale.W}
	B := Vector4{Dl.X / k.Scale.X, Dl.Y / k.Scale.Y, Dl.Z / k.Scale.Z, Dl.W / k.Scale.W}
	hxy, hzw := k.MajorRadius+k.tube, k.MinorRadius+k.tube
	if !k.Band {
		hxy += k.MinorRadius
	}
	iv, ok := boxIntervalLocal(A, B, hxy, hxy, hzw, hzw)
	if !ok || !(iv.t1 > 0) {
		return objectHit{}, false
	}
	speed := B.Len()
	dist := func(t Real) rayDist {
		q := A.Add(B.Mul(t))
		if far := k.aside(q) - k.tube; far > k.tube {
			*cheap++
			return rayDist{d: far}
		}
		d, foot, ex := k.nearest(q, true)
		if ex {
			*exact++
		} else {
			*lower++
		}
		rd := rayDist{d: d - k.tube}
		if ex && rd.d < 0 && d > 0 {
			rd.rate = q.Sub(foot).Dot(B) / d
			rd.smooth = d
		}
		return rd
	}
	t, inside, ok := traceSolid(dist, speed, iv.t0, iv.t1)
	if !ok {
		return objectHit{}, false
	}
	return objectHit{t: t, inv: inside}, true
}

func TestZZCount(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	for _, k := range kleinSamples(t) {
		type tally struct{ rays, cheap, exact, lower int }
		var hitT, missT, inT tally
		for n := 0; n < 4000; n++ {
			O := Point4{}.Add(unitS3(rng).Mul(4))
			T := Vector4{(rng.Float64()*2 - 1) * 1.4, (rng.Float64()*2 - 1) * 1.4, (rng.Float64()*2 - 1) * 0.5, (rng.Float64()*2 - 1) * 0.5}
			D := T.Sub(Vector4{O.X, O.Y, O.Z, O.W}).Norm()
			var c, e, l int
			h, ok := zzKlein(O, D, k, &c, &e, &l)
			tl := &missT
			if ok {
				tl = &hitT
			}
			tl.rays++
			tl.cheap += c
			tl.exact += e
			tl.lower += l
			if ok {
				c, e, l = 0, 0, 0
				P := O.Add(D.Mul(h.t + 1e-6))
				if _, ok2 := zzKlein(P, D, k, &c, &e, &l); ok2 {
					inT.rays++
					inT.cheap += c
					inT.exact += e
					inT.lower += l
				}
			}
		}
		pr := func(name string, x tally) {
			if x.rays == 0 {
				return
			}
			t.Logf("%-18s tw%d %-5s rays=%5d cheap=%6.1f lower=%5.1f exact=%5.1f per ray", k.Name, k.Twist, name, x.rays,
				float64(x.cheap)/float64(x.rays), float64(x.lower)/float64(x.rays), float64(x.exact)/float64(x.rays))
		}
		pr("hit", hitT)
		pr("miss", missT)
		pr("inner", inT)
	}
}
