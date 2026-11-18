package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestNewRevolvedQuadricsAndIntersect(t *testing.T) {
	one := Vector4{1, 1, 1, 1}
	must := func(s *RevolvedQuadric, err error) *RevolvedQuadric {
		t.Helper()
		if err != nil {
			t.Fatal(err)
		}
		return s
	}
	bowl := must(NewParaboloid(Point4{}, one, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	glass := must(NewHyperboloid(Point4{}, one, 0.5, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	cones := must(NewHyperboloid(Point4{}, one, 0, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	dome := must(NewSphericalSegment(Point4{}, one, 0, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	s := math.Sqrt(0.5)
	// the height at which the hourglass has the radius 0.8, and its normal there
	h := math.Sqrt(0.52)
	side := Vector4{1.6, 0, 0, 1.5 * h}.Norm()
	for _, c := range []struct {
		name string
		q    *RevolvedQuadric
		O    Point4
		D    Vector4
		ok   bool
		t    Real
		N    Vector4
		inv  bool
	}{
		{"paraboloid: onto the tip", bowl, Point4{0, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 4, Vector4{0, 0, 0, -1}, false},
		{"paraboloid: out of the base", bowl, Point4{0, 0, 0, 0}, Vector4{0, 0, 0, 1}, true, 1, Vector4{0, 0, 0, 1}, true},
		{"paraboloid: into the base", bowl, Point4{0.5, 0, 0, 5}, Vector4{0, 0, 0, -1}, true, 4, Vector4{0, 0, 0, 1}, false},
		{"paraboloid: onto its side", bowl, Point4{5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, true, 5 - s, Vector4{2 * s, 0, 0, -0.5}.Norm(), false},
		{"paraboloid: out of its side", bowl, Point4{0, 0, 0, 0}, Vector4{0, -1, 0, 0}, true, s, Vector4{0, -2 * s, 0, -0.5}.Norm(), true},
		{"paraboloid: past the base", bowl, Point4{2, 0, 0, -5}, Vector4{0, 0, 0, 1}, false, 0, Vector4{}, false},
		{"paraboloid: below the tip", bowl, Point4{5, 0, 0, -1.5}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"paraboloid: above the base", bowl, Point4{5, 0, 0, 1.1}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"hyperboloid: onto the waist", glass, Point4{5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, true, 4.5, Vector4{1, 0, 0, 0}, false},
		{"hyperboloid: out of the waist", glass, Point4{0, 0, 0, 0}, Vector4{0, 0, 1, 0}, true, 0.5, Vector4{0, 0, 1, 0}, true},
		{"hyperboloid: into an end", glass, Point4{0.8, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 4, Vector4{0, 0, 0, -1}, false},
		{"hyperboloid: out of the lower half by its side", glass, Point4{0.8, 0, 0, -0.9}, Vector4{0, 0, 0, 1}, true, 0.9 - h, Vector4{side.X, 0, 0, side.W}, true},
		{"hyperboloid: back into the upper half", glass, Point4{0.8, 0, 0, 0}, Vector4{0, 0, 0, 1}, true, h, Vector4{side.X, 0, 0, -side.W}, false},
		{"hyperboloid: out of the upper end", glass, Point4{0.8, 0, 0, 0.9}, Vector4{0, 0, 0, 1}, true, 0.1, Vector4{0, 0, 0, 1}, true},
		{"hyperboloid: past the waist", glass, Point4{5, 0.6, 0, 0}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"two cones: onto the upper one", cones, Point4{5, 0, 0, 0.5}, Vector4{-1, 0, 0, 0}, true, 4.5, Vector4{s, 0, 0, -s}, false},
		{"two cones: between them", cones, Point4{5, 0.01, 0, 0}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"two cones: from one into the other", cones, Point4{0.2, 0, 0, -0.5}, Vector4{0, 0, 0, 1}, true, 0.3, Vector4{s, 0, 0, s}, true},
		{"half ball: onto the flat side", dome, Point4{0.3, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 5, Vector4{0, 0, 0, -1}, false},
		{"half ball: onto the round side", dome, Point4{0, 0, 0, 5}, Vector4{0, 0, 0, -1}, true, 4, Vector4{0, 0, 0, 1}, false},
		{"half ball: out of the round side", dome, Point4{0, 0, 0, 0.5}, Vector4{0, 0.6, 0, 0.8}, true, 0.5 * (math.Sqrt(0.64+3) - 0.8), Vector4{0, 0.6 * 0.5 * (math.Sqrt(0.64+3) - 0.8), 0, 0.5 + 0.8*0.5*(math.Sqrt(0.64+3)-0.8)}, true},
		{"half ball: where the ball is cut away", dome, Point4{5, 0, 0, -0.1}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"half ball: onto its side", dome, Point4{5, 0, 0, 0.6}, Vector4{-1, 0, 0, 0}, true, 4.2, Vector4{0.8, 0, 0, 0.6}, false},
	} {
		got, ok := intersectRayRevolvedQuadric(c.O, c.D, c.q)
		if ok != c.ok || (ok && (math.Abs(got.t-c.t) > 1e-12 || got.Nw.Sub(c.N).Len() > 1e-12 || got.inv != c.inv)) {
			t.Errorf("%s: ok=%v t=%.15g N=%v inv=%v, want ok=%v t=%.15g N=%v inv=%v", c.name, ok, got.t, got.Nw, got.inv, c.ok, c.t, c.N, c.inv)
		}
	}

	if bowl.AABBMin != (Point4{-1, -1, -1, -1}) || bowl.AABBMax != (Point4{1, 1, 1, 1}) ||
		dome.AABBMin != (Point4{-1, -1, -1, 0}) || dome.AABBMax != (Point4{1, 1, 1, 1}) {
		t.Errorf("AABB of the paraboloid %v..%v, of the half ball %v..%v", bowl.AABBMin, bowl.AABBMax, dome.AABBMin, dome.AABBMax)
	}
	for _, c := range []struct {
		name string
		err  error
	}{
		{"a waist above 1", func() error {
			_, err := NewHyperboloid(Point4{}, one, 1.2, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a negative waist", func() error {
			_, err := NewHyperboloid(Point4{}, one, -0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a segment upside down", func() error {
			_, err := NewSphericalSegment(Point4{}, one, 0.5, 0.2, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a segment without a height", func() error {
			_, err := NewSphericalSegment(Point4{}, one, 0.5, 0.5, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a segment beyond the ball", func() error {
			_, err := NewSphericalSegment(Point4{}, one, 0, 1.5, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a segment that is not a number", func() error {
			_, err := NewSphericalSegment(Point4{}, one, math.NaN(), 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a zero scale", func() error {
			_, err := NewParaboloid(Point4{}, Vector4{1, 1, 1, 0}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a material above 1", func() error {
			_, err := NewParaboloid(Point4{}, one, Rot4{}, testWhite, testHalf, testHalf, testHalf, testIOR)
			return err
		}()},
	} {
		if c.err == nil {
			t.Errorf("%s must be rejected", c.name)
		}
	}
}

// The whole ball as a segment is the hypersphere, the hyperboloid with the
// waist 1 the spherinder, and the two cones are two hypercones.
func TestRevolvedQuadricsEqualKnownPrimitives(t *testing.T) {
	rng := rand.New(rand.NewSource(43))
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	must := func(err error) {
		t.Helper()
		if err != nil {
			t.Fatal(err)
		}
	}
	ball, err := NewSphericalSegment(Point4{}, sc, -1, 1, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	sphere, err := NewHyperSphere(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	sameHits(t, "the whole ball and the hypersphere", rng, 1.5,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayRevolvedQuadric(O, D, ball) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, sphere) })

	tube, err := NewHyperboloid(Point4{}, sc, 1, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	spherinder, err := NewSpherinder(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	sameHits(t, "the hyperboloid without a waist and the spherinder", rng, 1.5,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayRevolvedQuadric(O, D, tube) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRaySpherinder(O, D, spherinder) })

	// The upper cone is a hypercone with its tip at the centre, the lower one
	// such a hypercone turned over (the half-turn in ZW keeps X and Y, and Z
	// does not matter to a round base).
	cones, err := NewHyperboloid(Point4{}, Vector4{1.5, 1.5, 1.5, 1}, 0, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	upper, err := NewHyperCone(Point4{0, 0, 0, 0.5}, Vector4{1.5, 1.5, 1.5, 0.5}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	lower, err := NewHyperCone(Point4{0, 0, 0, -0.5}, Vector4{1.5, 1.5, 1.5, 0.5}, Rot4{ZW: math.Pi}, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	sameHits(t, "the hyperboloid with the waist 0 and two hypercones", rng, 1,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayRevolvedQuadric(O, D, cones) },
		func(O Point4, D Vector4) (objectHit, bool) {
			a, oka := intersectRayHyperCone(O, D, upper)
			b, okb := intersectRayHyperCone(O, D, lower)
			if okb && (!oka || b.t < a.t) {
				return b, true
			}
			return a, oka
		})
}

// The box is the exact one: the largest and the least coordinate over the
// balls of the solid at many heights.
func TestRevolvedQuadricAABBIsExact(t *testing.T) {
	rng := rand.New(rand.NewSource(44))
	for n := 0; n < 40; n++ {
		rot := Rot4{XY: rng.Float64() * 6, XZ: rng.Float64() * 6, XW: rng.Float64() * 6, YZ: rng.Float64() * 6, YW: rng.Float64() * 6, ZW: rng.Float64() * 6}
		ctr := Point4{rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5}
		sc := Vector4{0.5 + rng.Float64(), 0.5 + rng.Float64(), 0.5 + rng.Float64(), 0.5 + rng.Float64()}
		var q *RevolvedQuadric
		var err error
		switch n % 4 {
		case 0:
			q, err = NewParaboloid(ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		case 1:
			q, err = NewHyperboloid(ctr, sc, rng.Float64(), rot, testWhite, testNone, testHalf, testHalf, testIOR)
		case 2:
			lower := -1 + 1.5*rng.Float64()
			q, err = NewSphericalSegment(ctr, sc, lower, lower+(1-lower)*(0.1+0.9*rng.Float64()), rot, testWhite, testNone, testHalf, testHalf, testIOR)
		default:
			q, err = NewHyperboloid(ctr, sc, 0, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		}
		if err != nil {
			t.Fatal(err)
		}
		lo, hi := [4]Real{}, [4]Real{}
		for axis := 0; axis < 4; axis++ {
			a := Vector4{sc.X * q.R.M[axis][0], sc.Y * q.R.M[axis][1], sc.Z * q.R.M[axis][2], sc.W * q.R.M[axis][3]}
			lat := math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z)
			least, most := math.Inf(1), math.Inf(-1)
			const steps = 200000
			for i := 0; i <= steps; i++ {
				w := q.Lower + (q.Upper-q.Lower)*Real(i)/steps
				r := math.Sqrt(math.Max(0, q.K0+q.K1*w+q.K2*w*w))
				least, most = math.Min(least, a.W*w-lat*r), math.Max(most, a.W*w+lat*r)
			}
			c := [4]Real{ctr.X, ctr.Y, ctr.Z, ctr.W}[axis]
			lo[axis], hi[axis] = c+least, c+most
		}
		got := [8]Real{-q.AABBMin.X, -q.AABBMin.Y, -q.AABBMin.Z, -q.AABBMin.W, q.AABBMax.X, q.AABBMax.Y, q.AABBMax.Z, q.AABBMax.W}
		seen := [8]Real{-lo[0], -lo[1], -lo[2], -lo[3], hi[0], hi[1], hi[2], hi[3]}
		for i := range got {
			if got[i] < seen[i]-1e-12 || got[i] > seen[i]+1e-8 {
				t.Errorf("%s (K %g %g %g, w %g..%g): AABB %v..%v, the solid reaches %v..%v", q.Name, q.K0, q.K1, q.K2, q.Lower, q.Upper, q.AABBMin, q.AABBMax, lo, hi)
				break
			}
		}
	}
}
