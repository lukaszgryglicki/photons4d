package photons4d

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"regexp"
	"strconv"
	"strings"
	"testing"
	"time"
)

func TestZZSweep(t *testing.T) {
	sizes := []Real{1, 0.02, 0.2, 10, 100}
	raysOut := 160
	seed := int64(1)
	if v := os.Getenv("ZZ_SEED"); v != "" {
		seed, _ = strconv.ParseInt(v, 10, 64)
	}
	if v := os.Getenv("ZZ_RAYS"); v != "" {
		n, _ := strconv.Atoi(v)
		raysOut = n
	}
	if v := os.Getenv("ZZ_SIZES"); v != "" {
		sizes = nil
		for _, f := range strings.Split(v, ",") {
			x, _ := strconv.ParseFloat(f, 64)
			sizes = append(sizes, x)
		}
	}
	only := os.Getenv("ZZ_ONLY")
	spent := map[string]time.Duration{}
	rng := rand.New(rand.NewSource(seed))
	tally := map[string]*oracleTally{}
	var names []string
	for _, size := range sizes {
		solids := oracleSolids(t, rng, size)
		for i := range solids {
			s := &solids[i]
			if only != "" && !regexp.MustCompile(only).MatchString(s.name) {
				continue
			}
			began := time.Now()
			c := tally[s.name]
			if c == nil {
				c = &oracleTally{}
				tally[s.name] = c
				names = append(names, s.name)
			}
			mid := Point4{(s.min.X + s.max.X) / 2, (s.min.Y + s.max.Y) / 2, (s.min.Z + s.max.Z) / 2, (s.min.W + s.max.W) / 2}
			within := func(k Real) Point4 {
				return Point4{
					mid.X + (rng.Float64()-0.5)*(s.max.X-s.min.X)*k, mid.Y + (rng.Float64()-0.5)*(s.max.Y-s.min.Y)*k,
					mid.Z + (rng.Float64()-0.5)*(s.max.Z-s.min.Z)*k, mid.W + (rng.Float64()-0.5)*(s.max.W-s.min.W)*k,
				}
			}

			// the AABB holds the whole solid
			for n := 0; n < 4000; n++ {
				x := within(1.6)
				if s.f(x) <= 0 && (x.X < s.min.X || x.X > s.max.X || x.Y < s.min.Y || x.Y > s.max.Y ||
					x.Z < s.min.Z || x.Z > s.max.Z || x.W < s.min.W || x.W > s.max.W) {
					c.boxed++
					c.note("point of the solid outside the AABB: %v", x)
				}
			}

			type surf struct {
				P Point4
				N Vector4
				D Vector4
			}
			var surfs []surf
			// from outside towards the solid
			for n := 0; n < raysOut; n++ {
				O := mid.Add(unitS3(rng).Mul(-(s.size + 1)))
				T := within(1)
				D := Vector4{T.X - O.X, T.Y - O.Y, T.Z - O.Z, T.W - O.W}.Norm()
				if h, ok := oracleCheck(s, O, D, c, "outside"); ok {
					surfs = append(surfs, surf{O.Add(D.Mul(h.t)), h.Nw, D})
				}
			}
			// from inside the solid
			for n, tries := 0, 0; n < raysOut/2 && tries < raysOut*200; tries++ {
				if O := within(1); s.f(O) < 0 {
					n++
					oracleCheck(s, O, unitS3(rng), c, "inside")
				}
			}
			// the way the tracer continues from a hit: reflected, refracted and scattered rays
			for _, sf := range surfs {
				dirs := []Vector4{reflect4(sf.D, sf.N).Norm(), sampleDiffuseDir(sf.N, rng)}
				if T, ok := refract4(sf.D, sf.N, 1/1.5); ok {
					dirs = append(dirs, T.Norm())
				}
				for _, d := range dirs {
					oracleCheck(s, sf.P.Add(d.Mul(bumpShift)), d, c, "continued")
				}
			}
			// along the surface, a little above or below it
			for n, sf := range surfs {
				if n >= raysOut/2 {
					break
				}
				U, V, W := orthonormal3Fast(sf.N)
				a, b, w := sampleS2(rng)
				T := U.Mul(a).Add(V.Mul(b)).Add(W.Mul(w)).Norm()
				off := (rng.Float64()*2 - 1) * math.Pow(10, -1-5*rng.Float64()) * (1 + s.size)
				oracleCheck(s, sf.P.Add(sf.N.Mul(off)).Add(T.Mul(-(s.size + 1))), T, c, "grazing")
			}
			spent[s.name] += time.Since(began)
		}
	}
	for _, name := range names {
		c := tally[name]
		t.Logf("%-40s rays=%6d hits=%6d bad=%d time=%v", name, c.rays, c.hits, c.bad(), spent[name].Round(time.Millisecond))
		if c.hits == 0 {
			t.Errorf("%s: no hit in %d rays", name, c.rays)
		}
		if c.bad() > 0 {
			t.Errorf("%s: %d rays, %d hits: missed=%d late=%d phantom=%d side=%d normal=%d off-surface=%d aabb=%d\n  %s",
				name, c.rays, c.hits, c.missed, c.late, c.phantom, c.side, c.normal, c.away, c.boxed, fmt.Sprint(c.notes))
		}
	}
}
