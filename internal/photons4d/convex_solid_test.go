package photons4d

import (
	"math"
	"math/rand"
	"strings"
	"testing"
)

func TestBracketopeNotation(t *testing.T) {
	for _, bad := range []string{"", "IIII", "(III)", "(IIIII)", "(II)(II)", "((II)(II)", "(II]II)", "[()IIII]", "(IIxI)", "I", "<(II)(II)>>"} {
		if _, _, err := parseBracketope(bad); err == nil {
			t.Errorf("%q must not parse", bad)
		}
	}
	q := Vector4{0.3, -0.4, 1.2, -0.5}
	for _, c := range []struct {
		shape string
		gauge Real
		dual  Real
	}{
		{"(IIII)", q.Len(), q.Len()},
		{"[IIII]", 1.2, 2.4},
		{"<IIII>", 2.4, 1.2},
		{"[(II)(II)]", 1.3, 1.8},
		{"<(II)(II)>", 1.8, 1.3},
		{"[(III)I]", math.Sqrt(0.09 + 0.16 + 1.44), math.Sqrt(0.09+0.16+1.44) + 0.5},
		{"([II]<II>)", math.Hypot(0.4, 1.7), math.Hypot(0.7, 1.2)},
		{"<[I(II)]I>", math.Hypot(0.4, 1.2) + 0.5, math.Max(0.3+math.Hypot(0.4, 1.2), 0.5)},
	} {
		n, _, err := parseBracketope(c.shape)
		if err != nil {
			t.Fatalf("%s: %v", c.shape, err)
		}
		g, grad := n.eval(q)
		if math.Abs(g-c.gauge) > 1e-12 || math.Abs(n.dual(q)-c.dual) > 1e-12 {
			t.Errorf("%s: gauge %.12g and polar gauge %.12g, want %.12g and %.12g", c.shape, g, n.dual(q), c.gauge, c.dual)
		}
		// the gauge is its gradient times the point, and the gradient has the polar gauge 1
		if math.Abs(grad.Dot(q)-g) > 1e-12 || math.Abs(n.dual(grad)-1) > 1e-12 {
			t.Errorf("%s: gradient %v", c.shape, grad)
		}
	}
	for name, notation := range bracketopeNames {
		a, errA := NewBracketope(name, Point4{}, Vector4{1, 1, 1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
		b, errB := NewBracketope(notation, Point4{}, Vector4{1, 1, 1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
		if errA != nil || errB != nil || a.Name != b.Name || !strings.HasSuffix(a.Name, notation) {
			t.Errorf("%q is not %q (%v, %v)", name, notation, errA, errB)
		}
	}
	if _, err := NewBracketope("(IIII)", Point4{}, Vector4{1, 1, -1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("a negative scale must be rejected")
	}
}

// The bracketopes that are primitives of their own must be the same solids.
func TestBracketopesEqualKnownPrimitives(t *testing.T) {
	rng := rand.New(rand.NewSource(21))
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	must := func(err error) {
		t.Helper()
		if err != nil {
			t.Fatal(err)
		}
	}
	sphere, err := NewHyperSphere(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	cube, err := NewCell8(Point4{}, sc.Mul(2), rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	cross, err := NewCell16(Point4{}, sc.Mul(math.Sqrt2), rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	spherinder, err := NewSpherinder(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	cubinder, err := NewCubinder(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	duocylinder, err := NewDuocylinder(Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	must(err)
	for _, c := range []struct {
		shape string
		hit   func(Point4, Vector4) (objectHit, bool)
	}{
		{"(IIII)", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, sphere) }},
		{"((II)(II))", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, sphere) }},
		{"[IIII]", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell8(O, D, cube) }},
		{"[[II][II]]", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell8(O, D, cube) }},
		{"<IIII>", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell16(O, D, cross) }},
		{"<I<II>I>", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell16(O, D, cross) }},
		{"[(III)I]", func(O Point4, D Vector4) (objectHit, bool) { return intersectRaySpherinder(O, D, spherinder) }},
		{"[(II)II]", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCubinder(O, D, cubinder) }},
		{"[(II)(II)]", func(O Point4, D Vector4) (objectHit, bool) { return intersectRayDuocylinder(O, D, duocylinder) }},
	} {
		b, err := NewBracketope(c.shape, Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		must(err)
		sameHits(t, "bracketope "+c.shape, rng, 1.5, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayConvexSolid(O, D, b) }, c.hit)
	}
}

// In the hyperplane z = 0 the cylindrone and the dicone are the cone that the
// hypercone is there, and so is the coninder in w = 0 after Z and W change
// places.
func TestConeSolidsEqualTheHyperConeInASlice(t *testing.T) {
	rng := rand.New(rand.NewSource(22))
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	cone, err := NewHyperCone(Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	cylindrone, err := NewCylindrone(Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	dicone, err := NewDicone(Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	coninder, err := NewConinder(Point4{}, Vector4{sc.X, sc.Y, sc.W, 5}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	flat := func(v Vector4) Vector4 { return Vector4{v.X, v.Y, 0, v.W} }
	swap := func(v Vector4) Vector4 { return Vector4{v.X, v.Y, v.W, v.Z} }
	hits := 0
	for i := 0; i < 2000; i++ {
		o := flat(unitS3(rng).Mul(0.05 + 2.5*rng.Float64()))
		d := flat(unitS3(rng).Mul(1.2 * rng.Float64()).Sub(o)).Norm()
		O := Point4{o.X, o.Y, o.Z, o.W}
		want, ok := intersectRayHyperCone(O, d, cone)
		for name, got := range map[string]func() (objectHit, bool){
			"cylindrone": func() (objectHit, bool) { return intersectRayConvexSolid(O, d, cylindrone) },
			"dicone":     func() (objectHit, bool) { return intersectRayConvexSolid(O, d, dicone) },
			"coninder": func() (objectHit, bool) {
				h, ok := intersectRayConvexSolid(Point4{o.X, o.Y, o.W, o.Z}, swap(d), coninder)
				h.Nw = swap(h.Nw)
				return h, ok
			},
		} {
			h, okh := got()
			if ok != okh || (ok && (math.Abs(h.t-want.t) > 1e-9 || h.inv != want.inv)) {
				t.Fatalf("%s: ray %d O=%v D=%v: ok=%v t=%.12g inv=%v, the hypercone gives ok=%v t=%.12g inv=%v", name, i, O, d, okh, h.t, h.inv, ok, want.t, want.inv)
			}
			if ok && h.Nw.Sub(want.Nw).Len() > 1e-6 {
				// at the rim of the base either normal is right
				P := o.Add(d.Mul(h.t))
				if math.Abs(P.W-sc.W) > 1e-6 {
					t.Fatalf("%s: ray %d O=%v D=%v: normal %v, the hypercone gives %v", name, i, O, d, h.Nw, want.Nw)
				}
			}
		}
		if ok {
			hits++
		}
	}
	if hits < 500 {
		t.Fatalf("only %d hits", hits)
	}
}

func TestPrismicCylinder(t *testing.T) {
	rng := rand.New(rand.NewSource(23))
	sc := Vector4{0.7, 1.1, 0.9, 0.9}
	// with 4 sides it is the cubinder with its square turned by 45 degrees
	p, err := NewPrismicCylinder(4, Point4{}, sc, Rot4{XY: 0.4, YW: 0.3}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	// the last turn of rotFromAngles is the one in ZW, so it acts on the local axes first
	c, err := NewCubinder(Point4{}, Vector4{sc.X, sc.Y, sc.Z / math.Sqrt2, sc.W / math.Sqrt2}, Rot4{XY: 0.4, YW: 0.3, ZW: math.Pi / 4}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	sameHits(t, "4-gonal prismic cylinder and cubinder", rng, 1.3,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayConvexSolid(O, D, p) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCubinder(O, D, c) })

	// with many sides it tends to the duocylinder
	p, err = NewPrismicCylinder(360, Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	d, err := NewDuocylinder(Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 500; i++ {
		O := Point4{}.Add(unitS3(rng).Mul(3))
		D := unitS3(rng).Mul(0.5).Sub(Vector4{O.X, O.Y, O.Z, O.W}).Norm()
		hp, okp := intersectRayConvexSolid(O, D, p)
		hd, okd := intersectRayDuocylinder(O, D, d)
		if okp && okd && math.Abs(hp.t-hd.t) > 1e-3 {
			t.Fatalf("ray %d: t=%.9g, the duocylinder gives %.9g", i, hp.t, hd.t)
		}
		if okp && !okd {
			t.Fatalf("ray %d: the prism lies inside the duocylinder", i)
		}
	}
	for _, n := range []int{-1, 0, 2, 361} {
		if _, err := NewPrismicCylinder(n, Point4{}, sc, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
			t.Errorf("%d sides must be rejected", n)
		}
	}
}

// The AABB of a convex solid holds it and comes close to it on all eight
// sides (the sharp corners are found by the climb only roughly).
func TestConvexSolidAABBIsTight(t *testing.T) {
	rng := rand.New(rand.NewSource(24))
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	ctr := Point4{0.3, -0.2, 0.5, 0.1}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	var solids []*ConvexSolid
	for _, shape := range []string{"<(II)(II)>", "[(II)<II>]", "(<II>[II])", "<(III)I>"} {
		s, err := NewBracketope(shape, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		solids = append(solids, s)
	}
	for _, build := range []func(Point4, Vector4, Rot4, RGB, RGB, RGB, RGB, RGB) (*ConvexSolid, error){NewConinder, NewCylindrone, NewDicone} {
		s, err := build(ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		solids = append(solids, s)
	}
	s, err := NewPrismicCylinder(5, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	solids = append(solids, s)
	for _, thick := range []Real{0.2, 1.3} {
		s, err = NewLens(ctr, sc, thick, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		solids = append(solids, s)
	}
	for _, kind := range []string{"rotor", "simplex"} {
		s, err = NewReuleauxSolid(kind, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		solids = append(solids, s)
	}
	for _, s := range solids {
		lo, hi := [4]Real{}, [4]Real{}
		// A coordinate has no other top on a convex surface than the highest
		// one: climb from the point above the centre by ever smaller random
		// turns of the ray from the centre.
		for axis := 0; axis < 4; axis++ {
			for _, sign := range []Real{1, -1} {
				reach := func(D Vector4) Real {
					h, ok := intersectRayConvexSolid(ctr, D, s)
					if !ok || !h.inv {
						t.Fatalf("%s: a ray from the centre must leave the solid", s.Name)
					}
					return sign * vecAt(D.Mul(h.t), axis)
				}
				D := vecWith(Vector4{}, axis, sign)
				best, turn := reach(D), Real(0.5)
				for try := 0; try < 6000; try++ {
					C := D.Add(unitS3(rng).Mul(turn)).Norm()
					if r := reach(C); r > best {
						D, best = C, r
					} else {
						turn = math.Max(turn*0.995, 1e-9)
					}
				}
				c := [4]Real{ctr.X, ctr.Y, ctr.Z, ctr.W}[axis]
				if sign > 0 {
					hi[axis] = c + best
				} else {
					lo[axis] = c - best
				}
			}
		}
		got := [8]Real{s.AABBMin.X, s.AABBMin.Y, s.AABBMin.Z, s.AABBMin.W, s.AABBMax.X, s.AABBMax.Y, s.AABBMax.Z, s.AABBMax.W}
		seen := [8]Real{lo[0], lo[1], lo[2], lo[3], hi[0], hi[1], hi[2], hi[3]}
		for i := range got {
			out := got[i] - seen[i]
			if i < 4 {
				out = -out
			}
			if out < -1e-9 || out > 0.05 {
				t.Errorf("%s: AABB %v, the surface reaches %v", s.Name, got, seen)
				break
			}
		}
	}
}
