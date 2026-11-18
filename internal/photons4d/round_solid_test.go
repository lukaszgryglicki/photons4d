package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestLens(t *testing.T) {
	rng := rand.New(rand.NewSource(45))
	one := Vector4{1, 1, 1, 1}
	lens, err := NewLens(Point4{}, one, 0.6, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	// the faces belong to balls of the radius rho around (0,0,0,±k)
	k := (1 - 0.09) / 0.6
	rho := 0.3 + k
	low := k - math.Sqrt(rho*rho-0.25)
	for _, c := range []struct {
		name string
		O    Point4
		D    Vector4
		ok   bool
		t    Real
		N    Vector4
		inv  bool
	}{
		{"along the axis", Point4{0, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 4.7, Vector4{0, 0, 0, -1}, false},
		{"out along the axis", Point4{0, 0, 0, 0}, Vector4{0, 0, 0, 1}, true, 0.3, Vector4{0, 0, 0, 1}, true},
		{"off the axis", Point4{0.5, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 5 + low, Vector4{0.5 / rho, 0, 0, (low - k) / rho}, false},
		{"out off the axis", Point4{0, 0, -0.5, 0}, Vector4{0, 0, 0, 1}, true, -low, Vector4{0, 0, -0.5 / rho, (k - low) / rho}, true},
		{"past the rim", Point4{1.01, 0, 0, -5}, Vector4{0, 0, 0, 1}, false, 0, Vector4{}, false},
		{"beside a face", Point4{5, 0, 0, 0.31}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
	} {
		h, ok := intersectRayConvexSolid(c.O, c.D, lens)
		if ok != c.ok || (ok && (math.Abs(h.t-c.t) > 1e-12 || h.Nw.Sub(c.N).Len() > 1e-12 || h.inv != c.inv)) {
			t.Errorf("%s: ok=%v t=%.15g N=%v inv=%v, want ok=%v t=%.15g N=%v inv=%v", c.name, ok, h.t, h.Nw, h.inv, c.ok, c.t, c.N, c.inv)
		}
	}
	if lo, hi := lens.AABBMin, lens.AABBMax; (Vector4{lo.X + 1, lo.Y + 1, lo.Z + 1, lo.W + 0.3}).Len() > 1e-12 || (Vector4{hi.X - 1, hi.Y - 1, hi.Z - 1, hi.W - 0.3}).Len() > 1e-12 {
		t.Errorf("AABB %v..%v", lo, hi)
	}

	// Rays next to the axis meet behind a lens of glass at the distance that
	// the lensmaker's equation for a thick lens gives from its back.
	const n = 1.5
	d := 0.6
	f := 1 / ((n - 1) * (2/rho - (n-1)*d/(n*rho*rho)))
	back := f * (1 - (n-1)*d/(n*rho))
	for _, off := range []Vector4{{0.01, 0, 0, 0}, {0, -0.006, 0.008, 0}, {0.004, 0.004, -0.004, 0}} {
		O, D := Point4{off.X, off.Y, off.Z, -5}, Vector4{0, 0, 0, 1}
		h, ok := intersectRayConvexSolid(O, D, lens)
		T, okT := refract4(D, h.Nw, 1/n)
		if !ok || !okT {
			t.Fatalf("no way into the lens at %v", off)
		}
		O, D = O.Add(D.Mul(h.t)).Add(T.Norm().Mul(bumpShift)), T.Norm()
		h, ok = intersectRayConvexSolid(O, D, lens)
		T, okT = refract4(D, h.Nw, n)
		if !ok || !h.inv || !okT {
			t.Fatalf("no way out of the lens at %v", off)
		}
		O, D = O.Add(D.Mul(h.t)), T.Norm()
		// nearest to the axis where the part of the ray across it vanishes
		s := -(O.X*D.X + O.Y*D.Y + O.Z*D.Z) / (D.X*D.X + D.Y*D.Y + D.Z*D.Z)
		P := O.Add(D.Mul(s))
		if miss := math.Sqrt(P.X*P.X + P.Y*P.Y + P.Z*P.Z); miss > 1e-9 || math.Abs(P.W-0.3-back) > 1e-3*back {
			t.Errorf("ray at %v crosses the axis at w=%.9g (%.3g beside it), the focus is at %.9g", off, P.W, miss, 0.3+back)
		}
	}

	// the thickest lens is the ball
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	ball, err := NewLens(Point4{}, sc, 2, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	sphere, err := NewHyperSphere(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	sameHits(t, "the lens of the thickness 2 and the hypersphere", rng, 1.5,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayConvexSolid(O, D, ball) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, sphere) })

	for _, bad := range []Real{0, -0.5, 2.01, math.NaN()} {
		if _, err := NewLens(Point4{}, one, bad, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
			t.Errorf("the thickness %g must be rejected", bad)
		}
	}
}

func TestReuleauxSolids(t *testing.T) {
	rng := rand.New(rand.NewSource(46))
	one := Vector4{1, 1, 1, 1}
	rotor, err := NewReuleauxSolid("rotor", Point4{}, one, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	simplex, err := NewReuleauxSolid("Simplex", Point4{}, one, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	top, base := 2/math.Sqrt(3), -1/math.Sqrt(3)
	// at the height 0 the side is a piece of the sphere around the far corner (-1, base) of the triangle
	x := math.Sqrt(4-base*base) - 1
	corner := canonicalCell5()[0].Norm()
	edge := math.Sqrt(2.5)
	for _, c := range []struct {
		name string
		s    *ConvexSolid
		O    Point4
		D    Vector4
		t    Real
		N    Vector4
		inv  bool
	}{
		{"rotor: onto the bottom", rotor, Point4{0, 0, 0, -5}, Vector4{0, 0, 0, 1}, 5 + top - 2, Vector4{0, 0, 0, -1}, false},
		{"rotor: onto the side", rotor, Point4{5, 0, 0, 0}, Vector4{-1, 0, 0, 0}, 5 - x, Vector4{(x + 1) / 2, 0, 0, -base / 2}, false},
		{"rotor: out of the side", rotor, Point4{0, 0, 0, 0}, Vector4{0, 0, -1, 0}, x, Vector4{0, 0, -(x + 1) / 2, -base / 2}, true},
		{"simplex: onto the face across a corner", simplex, Point4{}.Add(corner.Mul(-5)), corner, 5 - (edge - 1), corner.Mul(-1), false},
		{"simplex: out of that face", simplex, Point4{}, corner.Mul(-1), edge - 1, corner.Mul(-1), true},
	} {
		h, ok := intersectRayConvexSolid(c.O, c.D, c.s)
		if !ok || math.Abs(h.t-c.t) > 1e-12 || h.Nw.Sub(c.N).Len() > 1e-12 || h.inv != c.inv {
			t.Errorf("%s: ok=%v t=%.15g N=%v inv=%v, want t=%.15g N=%v inv=%v", c.name, ok, h.t, h.Nw, h.inv, c.t, c.N, c.inv)
		}
	}
	// the corners of the simplex are its points furthest from the centre
	if h, ok := intersectRayConvexSolid(Point4{}, corner, simplex); !ok || math.Abs(h.t-1) > 1e-12 {
		t.Errorf("the corner of the simplex is at %.15g", h.t)
	}

	// The rotor has the width 2 in every direction: turned in any way, its
	// box is the box of the unit ball. The simplex has no single width.
	uneven := 0
	for n := 0; n < 50; n++ {
		rot := Rot4{XY: rng.Float64() * 6, XZ: rng.Float64() * 6, XW: rng.Float64() * 6, YZ: rng.Float64() * 6, YW: rng.Float64() * 6, ZW: rng.Float64() * 6}
		ctr := Point4{rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5}
		size := 0.3 + rng.Float64()
		r, err := NewReuleauxSolid("", ctr, one.Mul(size), rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		s, err := NewReuleauxSolid("simplex", ctr, one.Mul(size), rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		wide := Vector4{r.AABBMax.X - r.AABBMin.X, r.AABBMax.Y - r.AABBMin.Y, r.AABBMax.Z - r.AABBMin.Z, r.AABBMax.W - r.AABBMin.W}
		if wide.Sub(one.Mul(2*size)).Len() > 1e-10 {
			t.Fatalf("the rotor of the size %g is %v wide", size, wide)
		}
		wide = Vector4{s.AABBMax.X - s.AABBMin.X, s.AABBMax.Y - s.AABBMin.Y, s.AABBMax.Z - s.AABBMin.Z, s.AABBMax.W - s.AABBMin.W}
		if math.Max(math.Max(wide.X, wide.Y), math.Max(wide.Z, wide.W))-math.Min(math.Min(wide.X, wide.Y), math.Min(wide.Z, wide.W)) > 1e-3*size {
			uneven++
		}
	}
	if uneven < 25 {
		t.Errorf("the simplex has one width in %d of 50 positions", 50-uneven)
	}
	// The same width, measured on the surface itself: the furthest points of
	// the rotor along a direction and against it.
	for n := 0; n < 6; n++ {
		dir := unitS3(rng)
		reach := func(sign Real) Real {
			D, turn := dir.Mul(sign), Real(0.5)
			at := func(D Vector4) Real {
				h, ok := intersectRayConvexSolid(Point4{}, D, rotor)
				if !ok {
					t.Fatalf("a ray from the centre must leave the rotor")
				}
				return sign * D.Mul(h.t).Dot(dir)
			}
			best := at(D)
			for try := 0; try < 6000; try++ {
				C := D.Add(unitS3(rng).Mul(turn)).Norm()
				if r := at(C); r > best {
					D, best = C, r
				} else {
					turn = math.Max(turn*0.995, 1e-9)
				}
			}
			return best
		}
		if w := reach(1) + reach(-1); w > 2+1e-9 || w < 2-2e-3 {
			t.Errorf("the rotor is %.9g wide along %v", w, dir)
		}
	}

	if _, err := NewReuleauxSolid("tetrahedron", Point4{}, one, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("an unknown kind must be rejected")
	}
	if _, err := NewReuleauxSolid("rotor", Point4{}, Vector4{1, -1, 1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("a negative scale must be rejected")
	}
}

// The support of the part common to balls against the furthest of many
// points of its surface.
func TestBallsSupport(t *testing.T) {
	rng := rand.New(rand.NewSource(47))
	canon := canonicalCell5()
	var corners []Vector4
	for _, v := range canon {
		corners = append(corners, v.Norm())
	}
	for _, c := range []struct {
		name    string
		centers []Vector4
		rho     Real
	}{
		{"one ball", []Vector4{{0.3, -0.2, 0.1, 0.4}}, 0.7},
		{"lens", []Vector4{{0, 0, 0, -1.5}, {0, 0, 0, 1.5}}, 1.8},
		{"three balls", []Vector4{{0.5, 0, 0, 0}, {-0.3, 0.4, 0, 0.1}, {0, -0.3, 0.5, -0.2}}, 0.9},
		{"simplex", corners, math.Sqrt(2.5)},
		{"four of the five", corners[:4], 1.7},
	} {
		// points of the surface: from a point inside, out to the nearest of the spheres
		var mid Vector4
		for _, p := range c.centers {
			mid = mid.Add(p.Mul(1 / Real(len(c.centers))))
		}
		along := func(d Vector4) Vector4 {
			out := math.Inf(1)
			for _, p := range c.centers {
				// |mid + s·d - p| = rho
				m := mid.Sub(p)
				b := m.Dot(d)
				out = math.Min(out, -b+math.Sqrt(b*b-m.Dot(m)+c.rho*c.rho))
			}
			return mid.Add(d.Mul(out))
		}
		var surface []Vector4
		for n := 0; n < 20000; n++ {
			surface = append(surface, along(unitS3(rng)))
		}
		for n := 0; n < 40; n++ {
			a := unitS3(rng).Mul(0.5 + rng.Float64())
			if n < 4 {
				a = vecWith(Vector4{}, n, 1)
			}
			seen, D := math.Inf(-1), Vector4{}
			for _, p := range surface {
				if v := a.Dot(p); v > seen {
					seen, D = v, p.Sub(mid).Norm()
				}
			}
			// and closer to the furthest point by ever smaller turns
			for try, turn := 0, Real(0.05); try < 3000; try++ {
				C := D.Add(unitS3(rng).Mul(turn)).Norm()
				if v := a.Dot(along(C)); v > seen {
					seen, D = v, C
				} else {
					turn = math.Max(turn*0.99, 1e-9)
				}
			}
			if got := ballsSupport(c.centers, c.rho, a); got < seen-1e-12 || got > seen+5e-3 {
				t.Errorf("%s: support %.9g along %v, the surface reaches %.9g", c.name, got, a, seen)
			}
		}
	}
}

func TestBracketopeNamedAxes(t *testing.T) {
	rng := rand.New(rand.NewSource(48))
	for _, bad := range []string{"(xyz)", "[(xy)(xz)]", "[(xy)(zw)I]", "[(II)(zw)]", "(xyzq)", "x", "(xyzw", "[(xy)(zw)](xy)"} {
		if _, _, err := parseBracketope(bad); err == nil {
			t.Errorf("%q must not parse", bad)
		}
	}
	q := Vector4{0.3, -0.4, 1.2, -0.5}
	for _, c := range []struct {
		shape  string
		shared bool
		gauge  Real
	}{
		{"(xyzw)", false, q.Len()},
		{"[(zw)(yx)]", false, 1.3},
		{"<(XW)[yz]>", false, math.Hypot(0.3, 0.5) + 1.2},
		{"[(xyz)(xyw)]", true, math.Sqrt(0.09 + 0.16 + 1.44)},
		{"[(xy)(xz)(xw)(yz)(yw)(zw)]", true, 1.3},
		{"<(xy)(yz)(zw)>", true, 0.5 + math.Hypot(0.4, 1.2) + 1.3},
		{"([xy](yzw))", true, math.Sqrt(0.16 + 0.16 + 1.44 + 0.25)},
		{"[(xy)(xz)w]", true, math.Hypot(0.3, 1.2)},
	} {
		n, shared, err := parseBracketope(c.shape)
		if err != nil {
			t.Fatalf("%s: %v", c.shape, err)
		}
		g, grad := n.eval(q)
		if shared != c.shared || math.Abs(g-c.gauge) > 1e-12 || math.Abs(grad.Dot(q)-g) > 1e-12 {
			t.Errorf("%s: shared=%v gauge %.12g gradient %v, want %v and %.12g", c.shape, shared, g, grad, c.shared, c.gauge)
		}
		// The bound of the support holds for the points of the surface and is
		// not far above them; it is the polar gauge when no axis is shared.
		var surface []Vector4
		for i := 0; i < 40000; i++ {
			d := unitS3(rng)
			v, _ := n.eval(d)
			surface = append(surface, d.Mul(1/v))
		}
		for i := 0; i < 30; i++ {
			a := unitS3(rng)
			seen := math.Inf(-1)
			for _, p := range surface {
				seen = math.Max(seen, a.Dot(p))
			}
			got := n.cover(a)
			if got < seen-1e-12 || got > 1.25*seen {
				t.Errorf("%s: bound %.9g of the support along %v, the surface reaches %.9g", c.shape, got, a, seen)
			}
			if !shared && math.Abs(got-n.dual(a)) > 1e-12 {
				t.Errorf("%s: bound %.12g, polar gauge %.12g", c.shape, got, n.dual(a))
			}
		}
	}

	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	for _, pair := range [][2]string{{"[(xy)(zw)]", "[(II)(II)]"}, {"[(ZW)(YX)]", "[(II)(II)]"}, {"<(xyz)w>", "<(III)I>"}, {"(x[yz]w)", "(I[II]I)"}} {
		a, errA := NewBracketope(pair[0], Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		b, errB := NewBracketope(pair[1], Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if errA != nil || errB != nil {
			t.Fatalf("%v: %v, %v", pair, errA, errB)
		}
		sameHits(t, pair[0]+" and "+pair[1], rng, 1.5,
			func(O Point4, D Vector4) (objectHit, bool) { return intersectRayConvexSolid(O, D, a) },
			func(O Point4, D Vector4) (objectHit, bool) { return intersectRayConvexSolid(O, D, b) })
		lo := Vector4{a.AABBMin.X - b.AABBMin.X, a.AABBMin.Y - b.AABBMin.Y, a.AABBMin.Z - b.AABBMin.Z, a.AABBMin.W - b.AABBMin.W}
		hi := Vector4{a.AABBMax.X - b.AABBMax.X, a.AABBMax.Y - b.AABBMax.Y, a.AABBMax.Z - b.AABBMax.Z, a.AABBMax.W - b.AABBMax.W}
		if lo.Len() > 1e-12 || hi.Len() > 1e-12 {
			t.Errorf("%v: AABB %v..%v and %v..%v", pair, a.AABBMin, a.AABBMax, b.AABBMin, b.AABBMax)
		}
	}
}
