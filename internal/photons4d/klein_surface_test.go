package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

// kleinPoint is a point of the surface by its definition: the turn u and the
// place s on the profile (an angle of the circle, an offset on the segment).
func kleinPoint(k *KleinSurface, u, s Real) Vector4 {
	a, b := k.MajorRadius, s
	if !k.Band {
		a, b = k.MajorRadius+k.MinorRadius*math.Cos(s), k.MinorRadius*math.Sin(s)
	}
	c := Real(k.Twist) / 2
	return Vector4{a * math.Cos(u), a * math.Sin(u), b * math.Cos(c*u), b * math.Sin(c*u)}
}

// kleinGridDistance returns the least distance from q to the surface found
// on a grid of the two parameters and by a pattern search from its lowest
// cells: it shares nothing with the search of the solid over the turn alone.
func kleinGridDistance(k *KleinSurface, q Vector4) (Real, Vector4) {
	nu, ns := 120*(1+k.Twist/2), 60
	s0, s1 := Real(0), 2*math.Pi
	if k.Band {
		s0, s1 = -k.MinorRadius, k.MinorRadius
	}
	hu, hs := 2*math.Pi/Real(nu), (s1-s0)/Real(ns)
	at := func(u, s Real) Real {
		if k.Band {
			s = math.Max(s0, math.Min(s1, s))
		}
		return q.Sub(kleinPoint(k, u, s)).Len()
	}
	type cell struct{ d, u, s Real }
	var low [6]cell
	for i := range low {
		low[i].d = math.Inf(1)
	}
	for i := 0; i < nu; i++ {
		for j := 0; j <= ns; j++ {
			c := cell{at(Real(i)*hu, s0+Real(j)*hs), Real(i) * hu, s0 + Real(j)*hs}
			for n := range low {
				if c.d < low[n].d {
					c, low[n] = low[n], c
				}
			}
		}
	}
	best, bu, bs := math.Inf(1), Real(0), Real(0)
	for _, c := range low {
		d, u, s := c.d, c.u, c.s
		for du, ds := hu, hs; du > 1e-11; {
			moved := false
			for _, m := range [][2]Real{{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}} {
				if e := at(u+m[0]*du, s+m[1]*ds); e < d {
					d, u, s, moved = e, u+m[0]*du, s+m[1]*ds, true
				}
			}
			if !moved {
				du, ds = du/2, ds/2
			}
		}
		if d < best {
			best, bu, bs = d, u, s
		}
	}
	if k.Band {
		bs = math.Max(s0, math.Min(s1, bs))
	}
	return best, kleinPoint(k, bu, bs)
}

func kleinSamples(t *testing.T) []*KleinSurface {
	t.Helper()
	one := Vector4{1, 1, 1, 1}
	var out []*KleinSurface
	add := func(k *KleinSurface, err error) {
		t.Helper()
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, k)
	}
	add(NewKleinBottle(Point4{}, one, 1, 0.4, 0.16, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewKleinBottle(Point4{}, one, 0.7, 0.6, 0.1, 3, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewKleinBottle(Point4{}, one, 1.2, 0.5, 0.3, 2, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewKleinBottle(Point4{}, one, 1, 0.3, 0.1, 0, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewProjectivePlane(Point4{}, one, 0.5, 0.12, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewMobiusBand(Point4{}, one, 1, 0.45, 0.14, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewMobiusBand(Point4{}, one, 0.8, 0.5, 0.1, 3, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	add(NewMobiusBand(Point4{}, one, 0.9, 0.3, 0.2, 2, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	return out
}

// The distance to the surface, which the solid finds by a search over the
// turn of the profile alone, against a search over the whole surface.
func TestKleinSurfaceNearestAgainstGridSearch(t *testing.T) {
	rng := rand.New(rand.NewSource(41))
	points := 200
	if testing.Short() {
		points = 75
	}
	for _, k := range kleinSamples(t) {
		reach := k.MajorRadius + k.MinorRadius
		cheap := 0
		for n := 0; n < points; n++ {
			var q Vector4
			onSurface := func() Vector4 {
				s := (rng.Float64()*2 - 1) * math.Pi
				if k.Band {
					s = (rng.Float64()*2 - 1) * k.MinorRadius
				}
				return kleinPoint(k, rng.Float64()*2*math.Pi, s)
			}
			switch n % 5 {
			case 0, 1:
				// next to the surface, inside the wall and outside
				q = onSurface().Add(unitS3(rng).Mul(k.Thickness * 1.5 * rng.Float64()))
			case 2:
				q = unitS3(rng).Mul(reach * 1.3 * rng.Float64())
			case 3:
				// around the plane of the twist, where every turn is as near as any other
				q = Vector4{(rng.Float64() - 0.5) * 0.02, (rng.Float64() - 0.5) * 0.02, (rng.Float64() - 0.5) * 2 * reach, (rng.Float64() - 0.5) * 2 * reach}
			default:
				// between two parts of the surface, a little nearer to one of them:
				// the turn of the point itself is no guide then
				a, b := onSurface(), onSurface()
				q = a.Add(b.Sub(a).Mul(0.35 + 0.3*rng.Float64()))
			}
			want, wantFoot := kleinGridDistance(k, q)
			d, foot, exact := k.nearest(q, false)
			if !exact || math.Abs(d-want) > 1e-9 || math.Abs(q.Sub(foot).Len()-d) > 1e-12 {
				t.Fatalf("%s twist %d: distance of %v is %.12g (exact=%v, foot %v at %.12g), the grid search gives %.12g at %v",
					k.Name, k.Twist, q, d, exact, foot, q.Sub(foot).Len(), want, wantFoot)
			}
			if far, _ := kleinGridDistance(k, foot); far > 1e-9 {
				t.Fatalf("%s twist %d: the nearest point %v is off the surface by %.3g", k.Name, k.Twist, foot, far)
			}
			// the quick answer for a far point is a lower bound that keeps at least half of the way to the wall
			lo, _, exact := k.nearest(q, true)
			if exact {
				if lo != d {
					t.Fatalf("%s twist %d: two exact answers %.15g and %.15g for %v", k.Name, k.Twist, lo, d, q)
				}
			} else {
				cheap++
				if lo > want+1e-12 || lo-k.tube < (want-k.tube)/2-1e-12 {
					t.Fatalf("%s twist %d: bound %.12g of the distance %.12g of %v (wall %.3g)", k.Name, k.Twist, lo, want, q, k.tube)
				}
			}
		}
		if cheap == 0 || cheap > points*3/4 {
			t.Errorf("%s twist %d: %d of %d points took the quick answer", k.Name, k.Twist, cheap, points)
		}
	}
}

func TestNewKleinSurfacesAndIntersect(t *testing.T) {
	one := Vector4{1, 1, 1, 1}
	must := func(k *KleinSurface, err error) *KleinSurface {
		t.Helper()
		if err != nil {
			t.Fatal(err)
		}
		return k
	}
	// At the turn 0 the profile lies in the plane XZ around (1,0,0,0), at the
	// turn pi in the plane XW around (-1,0,0,0).
	bottle := must(NewKleinBottle(Point4{}, one, 1, 0.4, 0.2, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	band := must(NewMobiusBand(Point4{}, one, 1, 0.4, 0.2, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	// all the circles of the projective plane pass through the centre, where the surface lies along ZW
	plane := must(NewProjectivePlane(Point4{}, one, 0.5, 0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	s := math.Sqrt(0.5)
	for _, c := range []struct {
		name string
		k    *KleinSurface
		O    Point4
		D    Vector4
		ok   bool
		t    Real
		N    Vector4
		inv  bool
	}{
		{"bottle: into the outer rim", bottle, Point4{-5, 0, 0, 0}, Vector4{1, 0, 0, 0}, true, 3.5, Vector4{-1, 0, 0, 0}, false},
		{"bottle: out of the wall inwards", bottle, Point4{-1.4, 0, 0, 0}, Vector4{1, 0, 0, 0}, true, 0.1, Vector4{1, 0, 0, 0}, true},
		{"bottle: from the middle to the inner rim", bottle, Point4{0, 0, 0, 0}, Vector4{1, 0, 0, 0}, true, 0.5, Vector4{-1, 0, 0, 0}, false},
		{"bottle: onto the top of the circle at the turn 0", bottle, Point4{1, 0, 5, 0}, Vector4{0, 0, -1, 0}, true, 4.5, Vector4{0, 0, 1, 0}, false},
		{"bottle: the circle at the turn pi stands along W", bottle, Point4{-1, 0, 0, 5}, Vector4{0, 0, 0, -1}, true, 4.5, Vector4{0, 0, 0, 1}, false},
		{"bottle: and not along Z", bottle, Point4{-1, 0, 5, 0}, Vector4{0, 0, -1, 0}, false, 0, Vector4{}, false},
		{"bottle: through the hole along W", bottle, Point4{0, 0, 0, -5}, Vector4{0, 0, 0, 1}, false, 0, Vector4{}, false},
		{"bottle: away", bottle, Point4{-5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"band: onto its face", band, Point4{5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, true, 3.9, Vector4{1, 0, 0, 0}, false},
		{"band: out of its face", band, Point4{1, 0, 0.2, 0}, Vector4{-1, 0, 0, 0}, true, 0.1, Vector4{-1, 0, 0, 0}, true},
		{"band: onto its round edge", band, Point4{1, 0, 5, 0}, Vector4{0, 0, -1, 0}, true, 4.5, Vector4{0, 0, 1, 0}, false},
		{"band: onto the edge at the turn pi", band, Point4{-1, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 4.5, Vector4{0, 0, 0, -1}, false},
		{"band: past the edge", band, Point4{1, 0, 5, 0.11}, Vector4{0, 0, -1, 0}, false, 0, Vector4{}, false},
		{"band: across the edge of the wall", band, Point4{1 + 0.1*s, 0, 5, 0}, Vector4{0, 0, -1, 0}, true, 4.6 - 0.1*s, Vector4{s, 0, s, 0}, false},
		{"plane: onto the far side of a circle", plane, Point4{5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, true, 3.95, Vector4{1, 0, 0, 0}, false},
		{"plane: onto the sheet through the centre", plane, Point4{0.5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, true, 0.45, Vector4{1, 0, 0, 0}, false},
		{"plane: out of that sheet", plane, Point4{0, 0, 0, 0}, Vector4{0, 1, 0, 0}, true, 0.05, Vector4{0, 1, 0, 0}, true},
		// along W the ray stays in the plane of the circle at the turn pi, around (-0.5,0,0,0)
		{"plane: along the sheet, which bends away", plane, Point4{0, 0, 0, 0}, Vector4{0, 0, 0, 1}, true, math.Sqrt(0.0525), Vector4{0.5 / 0.55, 0, 0, math.Sqrt(0.0525) / 0.55}, true},
	} {
		h, ok := intersectRayKleinSurface(c.O, c.D, c.k)
		if ok != c.ok || (ok && (math.Abs(h.t-c.t) > 1e-9 || h.Nw.Sub(c.N).Len() > 1e-7 || h.inv != c.inv)) {
			t.Errorf("%s: ok=%v t=%.15g N=%v inv=%v, want ok=%v t=%g N=%v inv=%v", c.name, ok, h.t, h.Nw, h.inv, c.ok, c.t, c.N, c.inv)
		}
	}

	for _, c := range []struct {
		name string
		err  error
	}{
		{"equal radii of a bottle", func() error {
			_, err := NewKleinBottle(Point4{}, one, 0.5, 0.5, 0.1, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a minor radius above the major one", func() error {
			_, err := NewKleinBottle(Point4{}, one, 0.5, 0.6, 0.1, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"no thickness", func() error {
			_, err := NewKleinBottle(Point4{}, one, 1, 0.4, 0, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a negative twist", func() error {
			_, err := NewKleinBottle(Point4{}, one, 1, 0.4, 0.1, -1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"too many half-turns", func() error {
			_, err := NewMobiusBand(Point4{}, one, 1, 0.4, 0.1, maxKleinTwist+1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a band without a width", func() error {
			_, err := NewMobiusBand(Point4{}, one, 1, 0, 0.1, 1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a projective plane without a radius", func() error {
			_, err := NewProjectivePlane(Point4{}, one, 0, 0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a radius that is not a number", func() error {
			_, err := NewProjectivePlane(Point4{}, one, math.NaN(), 0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a zero scale", func() error {
			_, err := NewProjectivePlane(Point4{}, Vector4{1, 1, 0, 1}, 0.5, 0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
			return err
		}()},
		{"a material above 1", func() error {
			_, err := NewMobiusBand(Point4{}, one, 1, 0.4, 0.1, 1, Rot4{}, testWhite, testHalf, testHalf, testHalf, testIOR)
			return err
		}()},
	} {
		if c.err == nil {
			t.Errorf("%s must be rejected", c.name)
		}
	}
	for _, c := range []struct {
		k    *KleinSurface
		name string
	}{
		{bottle, "klein bottle"}, {band, "mobius band"}, {plane, "projective plane"},
		{must(NewKleinBottle(Point4{}, one, 1, 0.4, 0.2, 2, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)), "twisted torus"},
		{must(NewMobiusBand(Point4{}, one, 1, 0.4, 0.2, 4, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)), "twisted ring"},
	} {
		if c.k.Name != c.name {
			t.Errorf("%q is named %q", c.name, c.k.Name)
		}
	}
}

// Without a twist the circle sweeps a torus of XYZ, and the wall around it
// is the ditorus.
func TestKleinSurfaceWithoutTwistIsTheDitorus(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	k, err := NewKleinBottle(Point4{}, sc, 1, 0.45, 0.4, 0, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	d, err := NewDitorus(Point4{}, sc, 1, 0.45, 0.2, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	sameHits(t, "klein surface without a twist and ditorus", rng, 1.5,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayKleinSurface(O, D, k) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayDitorus(O, D, d) })
	// the box of the ditorus is exact; the one of the wall may be a little wider
	got := [8]Real{-k.AABBMin.X, -k.AABBMin.Y, -k.AABBMin.Z, -k.AABBMin.W, k.AABBMax.X, k.AABBMax.Y, k.AABBMax.Z, k.AABBMax.W}
	want := [8]Real{-d.AABBMin.X, -d.AABBMin.Y, -d.AABBMin.Z, -d.AABBMin.W, d.AABBMax.X, d.AABBMax.Y, d.AABBMax.Z, d.AABBMax.W}
	for i := range got {
		if got[i] < want[i]-1e-12 || got[i] > want[i]+0.02 {
			t.Fatalf("AABB %v..%v, the ditorus has %v..%v", k.AABBMin, k.AABBMax, d.AABBMin, d.AABBMax)
		}
	}
}

// The box holds the wall and comes close to it on all eight sides, whatever
// the twist; the wall reaches half of the thickness beyond the surface.
func TestKleinSurfaceAABB(t *testing.T) {
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	ctr := Point4{0.3, -0.2, 0.5, 0.1}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	for _, base := range kleinSamples(t) {
		k, err := newKleinSurface(base.Name, base.Band, base.MajorRadius, base.MinorRadius, base.Thickness, base.Twist, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		lo, hi := [4]Real{}, [4]Real{}
		for axis := 0; axis < 4; axis++ {
			a := Vector4{sc.X * k.R.M[axis][0], sc.Y * k.R.M[axis][1], sc.Z * k.R.M[axis][2], sc.W * k.R.M[axis][3]}
			least, most := math.Inf(1), math.Inf(-1)
			for i := 0; i < 2880; i++ {
				for j := 0; j <= 96; j++ {
					s := Real(j) / 96 * 2 * math.Pi
					if k.Band {
						s = (Real(j)/48 - 1) * k.MinorRadius
					}
					v := a.Dot(kleinPoint(k, Real(i)/2880*2*math.Pi, s))
					least, most = math.Min(least, v), math.Max(most, v)
				}
			}
			c := [4]Real{ctr.X, ctr.Y, ctr.Z, ctr.W}[axis]
			lo[axis], hi[axis] = c+least-k.Thickness/2*a.Len(), c+most+k.Thickness/2*a.Len()
		}
		got := [8]Real{-k.AABBMin.X, -k.AABBMin.Y, -k.AABBMin.Z, -k.AABBMin.W, k.AABBMax.X, k.AABBMax.Y, k.AABBMax.Z, k.AABBMax.W}
		seen := [8]Real{-lo[0], -lo[1], -lo[2], -lo[3], hi[0], hi[1], hi[2], hi[3]}
		for i := range got {
			if got[i] < seen[i]-1e-12 || got[i] > seen[i]+0.03 {
				t.Errorf("%s twist %d: AABB %v..%v, the wall reaches %v..%v", k.Name, k.Twist, k.AABBMin, k.AABBMax, lo, hi)
				break
			}
		}
	}
}
