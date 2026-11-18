package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestNewTorinderAndIntersect(t *testing.T) {
	mk := func(sw, R, r Real) *Torinder {
		t.Helper()
		tor, err := NewTorinder(Point4{}, Vector4{1, 1, 1, sw}, R, r, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		return tor
	}
	tor := mk(2, 1, 0.25)
	for _, c := range []struct {
		name string
		O    Point4
		D    Vector4
		ok   bool
		t    Real
		N    Vector4
		inv  bool
	}{
		{"into an end", Point4{1, 0, 0, -5}, Vector4{0, 0, 0, 1}, true, 3, Vector4{0, 0, 0, -1}, false},
		{"into the other end", Point4{0, -1.1, 0.1, 6}, Vector4{0, 0, 0, -1}, true, 4, Vector4{0, 0, 0, 1}, false},
		{"out of an end", Point4{1, 0, 0, 0.5}, Vector4{0, 0, 0, 1}, true, 1.5, Vector4{0, 0, 0, 1}, true},
		{"through the hole", Point4{0, 0, 0, -5}, Vector4{0, 0, 0, 1}, false, 0, Vector4{}, false},
		{"past the tube", Point4{1.3, 0, 0, -5}, Vector4{0, 0, 0, 1}, false, 0, Vector4{}, false},
		{"into the tube", Point4{5, 0, 0, 1}, Vector4{-1, 0, 0, 0}, true, 3.75, Vector4{1, 0, 0, 0}, false},
		{"out of the tube into the hole", Point4{1, 0, 0, -1.5}, Vector4{-1, 0, 0, 0}, true, 0.25, Vector4{-1, 0, 0, 0}, true},
		{"across the hole", Point4{0, 0, 0, 1.9}, Vector4{0, 1, 0, 0}, true, 0.75, Vector4{0, -1, 0, 0}, false},
		{"out of the tube upwards", Point4{0, 1, 0, 0}, Vector4{0, 0, 1, 0}, true, 0.25, Vector4{0, 0, 1, 0}, true},
		{"beyond the end", Point4{5, 0, 0, 2.1}, Vector4{-1, 0, 0, 0}, false, 0, Vector4{}, false},
		{"away", Point4{5, 0, 0, 0}, Vector4{1, 0, 0, 0}, false, 0, Vector4{}, false},
	} {
		h, ok := intersectRayTorinder(c.O, c.D, tor)
		if ok != c.ok || (ok && (math.Abs(h.t-c.t) > 1e-12 || h.Nw.Sub(c.N).Len() > 1e-12 || h.inv != c.inv)) {
			t.Errorf("%s: ok=%v t=%.15g N=%v inv=%v, want ok=%v t=%g N=%v inv=%v", c.name, ok, h.t, h.Nw, h.inv, c.ok, c.t, c.N, c.inv)
		}
	}
	// slanted through an end: enters the slab outside the tube and meets it later
	O, D := Point4{2.5, 0, 0, -3}, Vector4{-1, 0, 0, 1}.Norm()
	h, ok := intersectRayTorinder(O, D, tor)
	P := O.Add(D.Mul(h.t))
	if !ok || h.inv || math.Abs(P.X-1.25) > 1e-12 || math.Abs(P.W+1.75) > 1e-12 || h.Nw.Sub(Vector4{1, 0, 0, 0}).Len() > 1e-12 {
		t.Errorf("slanted ray: ok=%v P=%v N=%v inv=%v", ok, P, h.Nw, h.inv)
	}

	if tor.AABBMin != (Point4{-1.25, -1.25, -0.25, -2}) || tor.AABBMax != (Point4{1.25, 1.25, 0.25, 2}) {
		t.Errorf("AABB %v..%v", tor.AABBMin, tor.AABBMax)
	}
	for _, c := range [][2]Real{{1, 0}, {1, -0.1}, {0.3, 0.3}, {0.2, 0.5}} {
		if _, err := NewTorinder(Point4{}, Vector4{1, 1, 1, 1}, c[0], c[1], Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
			t.Errorf("radii %v must be rejected", c)
		}
	}
	if _, err := NewTorinder(Point4{}, Vector4{1, 1, 1, 0}, 1, 0.25, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("a zero scale must be rejected")
	}
}

// In the middle hyperplane of a long torinder the solid is a torus of XYZ;
// so is the torisphere in its hyperplane z = 0, with W in the place of Z.
func TestTorinderEqualsTheTorisphereInASlice(t *testing.T) {
	rng := rand.New(rand.NewSource(31))
	R, r := 0.9, 0.3
	tor, err := NewTorinder(Point4{}, Vector4{0.7, 1.1, 0.9, 50}, R, r, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	ts, err := NewTorisphere(Point4{}, Vector4{0.7, 1.1, 1.7, 0.9}, R, r, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	swap := func(v Vector4) Vector4 { return Vector4{v.X, v.Y, v.W, v.Z} }
	hits, inner := 0, 0
	for i := 0; i < 3000; i++ {
		o := unitS3(rng).Mul(0.05 + 2.5*rng.Float64())
		o.W = 0
		d := unitS3(rng).Mul(1.2 * rng.Float64())
		d.W = 0
		d = d.Sub(o).Norm()
		h, ok := intersectRayTorinder(Point4{o.X, o.Y, o.Z, 0}, d, tor)
		want, okw := intersectRayTorisphere(Point4{o.X, o.Y, 0, o.Z}, swap(d), ts)
		if ok != okw || (ok && (math.Abs(h.t-want.t) > 1e-9 || h.inv != want.inv || h.Nw.Sub(swap(want.Nw)).Len() > 1e-7)) {
			t.Fatalf("ray %d O=%v D=%v: ok=%v t=%.12g inv=%v N=%v, the torisphere gives ok=%v t=%.12g inv=%v N=%v",
				i, o, d, ok, h.t, h.inv, h.Nw, okw, want.t, want.inv, swap(want.Nw))
		}
		if ok {
			hits++
			if h.inv {
				inner++
			}
		}
	}
	if hits < 1000 || inner < 100 {
		t.Fatalf("%d hits, %d from inside", hits, inner)
	}
}
