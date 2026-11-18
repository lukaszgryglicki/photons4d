package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func BenchmarkZZNearest(b *testing.B) {
	t := &testing.T{}
	rng := rand.New(rand.NewSource(7))
	for _, k := range kleinSamples(t) {
		reach := k.MajorRadius + k.MinorRadius
		var near, far, in []Vector4
		for n := 0; n < 4096; n++ {
			s := (rng.Float64()*2 - 1) * math.Pi
			if k.Band {
				s = (rng.Float64()*2 - 1) * k.MinorRadius
			}
			p := kleinPoint(k, rng.Float64()*2*math.Pi, s)
			near = append(near, p.Add(unitS3(rng).Mul(k.tube*(1+0.2*rng.Float64()))))
			in = append(in, p.Add(unitS3(rng).Mul(k.tube*rng.Float64())))
			far = append(far, Vector4{(rng.Float64()*2 - 1) * reach, (rng.Float64()*2 - 1) * reach, (rng.Float64()*2 - 1) * k.MinorRadius, (rng.Float64()*2 - 1) * k.MinorRadius})
		}
		for _, c := range []struct {
			name string
			pts  []Vector4
		}{{"near", near}, {"in", in}, {"far", far}} {
			pts := c.pts
			exact := 0
			for _, q := range pts {
				if _, _, e := k.nearest(q, true); e {
					exact++
				}
			}
			b.Run(k.Name+"/"+c.name, func(b *testing.B) {
				b.ReportMetric(float64(exact)/float64(len(pts)), "exact-share")
				for i := 0; i < b.N; i++ {
					k.nearest(pts[i%len(pts)], true)
				}
			})
		}
	}
}
