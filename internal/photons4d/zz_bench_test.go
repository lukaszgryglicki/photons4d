package photons4d

import (
	"math/rand"
	"os"
	"regexp"
	"testing"
)

func BenchmarkZZRays(b *testing.B) {
	t := &testing.T{}
	rng := rand.New(rand.NewSource(7))
	solids := oracleSolids(t, rng, 1)
	only := regexp.MustCompile(os.Getenv("ZZ_ONLY"))
	for i := range solids {
		s := &solids[i]
		if !only.MatchString(s.name) {
			continue
		}
		mid := Point4{(s.min.X + s.max.X) / 2, (s.min.Y + s.max.Y) / 2, (s.min.Z + s.max.Z) / 2, (s.min.W + s.max.W) / 2}
		within := func(k Real) Point4 {
			return Point4{
				mid.X + (rng.Float64()-0.5)*(s.max.X-s.min.X)*k, mid.Y + (rng.Float64()-0.5)*(s.max.Y-s.min.Y)*k,
				mid.Z + (rng.Float64()-0.5)*(s.max.Z-s.min.Z)*k, mid.W + (rng.Float64()-0.5)*(s.max.W-s.min.W)*k,
			}
		}
		type ray struct {
			O Point4
			D Vector4
		}
		var outer, inner []ray
		hits := 0
		for n := 0; n < 2048; n++ {
			O := mid.Add(unitS3(rng).Mul(-(s.size + 1)))
			T := within(1)
			D := Vector4{T.X - O.X, T.Y - O.Y, T.Z - O.Z, T.W - O.W}.Norm()
			outer = append(outer, ray{O, D})
			if h, ok := s.hit(O, D); ok {
				hits++
				// continue inside: refracted direction
				if T, ok := refract4(D, h.Nw, 1/1.5); ok {
					P := O.Add(D.Mul(h.t))
					inner = append(inner, ray{P.Add(T.Norm().Mul(bumpShift)), T.Norm()})
				}
			}
		}
		b.Run(s.name+"/outer", func(b *testing.B) {
			b.ReportMetric(float64(hits)/2048, "hit-share")
			for i := 0; i < b.N; i++ {
				r := outer[i%len(outer)]
				s.hit(r.O, r.D)
			}
		})
		if len(inner) > 0 {
			b.Run(s.name+"/inner", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					r := inner[i%len(inner)]
					s.hit(r.O, r.D)
				}
			})
		}
	}
}
