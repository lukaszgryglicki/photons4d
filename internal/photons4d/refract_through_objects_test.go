package photons4d

import (
	"math"
	"testing"
)

// tiny helper
func bumpPoint(p Point4, d Vector4, eps Real) Point4 {
	return Point4{p.X + d.X*eps, p.Y + d.Y*eps, p.Z + d.Z*eps, p.W + d.W*eps}
}

func TestRefractionThrough_Hypercube(t *testing.T) {
	ior := Real(1.1) // mild to avoid TIR in weird angles
	h, err := NewHyperCube(
		Point4{0, 0, 0, 0.6},
		Vector4{0.3, 0.3, 0.3, 0.3},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1} // shoot along +W

	// 1) enter
	h1, ok := intersectRayHypercube(O, D, h)
	if !ok || h1.inv {
		t.Fatalf("expected entering hit, ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))

	T1, ok := refract4(D, h1.Nw, 1/ior) // air->material
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	// 2) inside -> move a bit and hit far side
	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayHypercube(O2, T1, h)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	// 3) exit refraction
	T2, ok := refract4(T1, h2.Nw, ior) // material->air
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	// sanity: exiting means outward normal has positive dot with outgoing dir
	if h2.Nw.Dot(T2) <= 1e-9 {
		t.Fatalf("not exiting outward: N·T2=%.12g", h2.Nw.Dot(T2))
	}
	// and we shouldn't immediately re-enter a convex object
	O3 := bumpPoint(P2, T2, bumpShift)
	if h3, ok := intersectRayHypercube(O3, T2, h); ok && h3.t < 1e-6 {
		t.Fatalf("ray re-entered immediately after exit")
	}
}

func TestRefractionThrough_HyperSphere(t *testing.T) {
	ior := Real(1.1)
	s, err := NewHyperSphere(
		Point4{0, 0, 0, 0.6},
		Vector4{0.25, 0.25, 0.25, 0.25},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRayHyperSphere(O, D, s)
	if !ok || h1.inv {
		t.Fatalf("expected entering sphere hit, ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayHyperSphere(O2, T1, s)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit sphere hit, ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	if h2.Nw.Dot(T2) <= 1e-9 {
		t.Fatalf("not exiting outward: N·T2=%.12g", h2.Nw.Dot(T2))
	}
	O3 := bumpPoint(P2, T2, bumpShift)
	if h3, ok := intersectRayHyperSphere(O3, T2, s); ok && h3.t < 1e-6 {
		t.Fatalf("ray re-entered immediately after exit")
	}

	// axis symmetry sanity: for axis case, final dir should still be "forward"
	if D.Dot(T2) <= 0 {
		t.Fatalf("emerged going backwards: D·T2=%.12g", D.Dot(T2))
	}
}

func TestRefractionThrough_Cell5(t *testing.T) {
	ior := Real(1.1)
	sx, err := NewCell5(
		Point4{0, 0, 0, 0.6},
		0.30,
		Vector4{1, 1, 1, 1},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRayCell5(O, D, sx)
	if !ok || h1.inv {
		t.Fatalf("expected entering 5-cell hit, ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayCell5(O2, T1, sx)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit 5-cell hit, ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	if h2.Nw.Dot(T2) <= 1e-9 {
		t.Fatalf("not exiting outward: N·T2=%.12g", h2.Nw.Dot(T2))
	}
	O3 := bumpPoint(P2, T2, bumpShift)
	if h3, ok := intersectRayCell5(O3, T2, sx); ok && h3.t < 1e-6 {
		t.Fatalf("ray re-entered immediately after exit")
	}
}

func TestRefractionThrough_SixteenCell(t *testing.T) {
	// Requires NewSixteenCell + intersectRaySixteen to be implemented.
	ior := Real(1.1)
	c16, err := NewSixteenCell(
		Point4{0, 0, 0, 0.6},
		0.30,
		Vector4{1, 1, 1, 1},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRaySixteen(O, D, c16)
	if !ok || h1.inv {
		t.Fatalf("expected entering 16-cell hit, ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRaySixteen(O2, T1, c16)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit 16-cell hit, ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	if h2.Nw.Dot(T2) <= 1e-9 {
		t.Fatalf("not exiting outward: N·T2=%.12g", h2.Nw.Dot(T2))
	}
	O3 := bumpPoint(P2, T2, bumpShift)
	if h3, ok := intersectRaySixteen(O3, T2, c16); ok && h3.t < 1e-6 {
		t.Fatalf("ray re-entered immediately after exit")
	}
}

func TestRefractionThrough_TwentyFourCell(t *testing.T) {
	// Requires NewTwentyFourCell + intersectRayTwentyFour to be implemented.
	ior := Real(1.1)
	c24, err := NewTwentyFourCell(
		Point4{0, 0, 0, 0.6},
		0.30,
		Vector4{1, 1, 1, 1},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatal(err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRayTwentyFour(O, D, c24)
	if !ok || h1.inv {
		t.Fatalf("expected entering 24-cell hit, ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayTwentyFour(O2, T1, c24)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit 24-cell hit, ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	if h2.Nw.Dot(T2) <= 1e-9 {
		t.Fatalf("not exiting outward: N·T2=%.12g", h2.Nw.Dot(T2))
	}
	O3 := bumpPoint(P2, T2, bumpShift)
	if h3, ok := intersectRayTwentyFour(O3, T2, c24); ok && h3.t < 1e-6 {
		t.Fatalf("ray re-entered immediately after exit")
	}
}

// Optional tiny numeric sanity for unit-length directions
func almostUnit(v Vector4) bool {
	l := v.Len()
	return math.Abs(float64(l-1)) < 1e-9
}

func TestRefractionThrough_Cell120(t *testing.T) {
	ior := Real(1.05)
	obj, err := NewCell120(
		Point4{0, 0, 0, 0.6},
		0.20,
		Vector4{1, 1, 1, 1},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatalf("NewCell120: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRayCellPoly(O, D, &obj.cellPoly)
	if !ok || h1.inv {
		t.Fatalf("expected entry hit; ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayCellPoly(O2, T1, &obj.cellPoly)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok, h2.inv)
	}
	P2 := O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()

	// should continue generally forward in +W
	if T2.W <= 0 {
		t.Fatalf("expected +W after exit, got %+v", T2)
	}
	_ = P2
}

func TestRefractionThrough_Cell600(t *testing.T) {
	ior := Real(1.05)
	obj, err := NewCell600(
		Point4{0, 0, 0, 0.6},
		0.20,
		Vector4{1, 1, 1, 1},
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{ior, ior, ior},
	)
	if err != nil {
		t.Fatalf("NewCell600: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}

	h1, ok := intersectRayCellPoly(O, D, &obj.cellPoly)
	if !ok || h1.inv {
		t.Fatalf("expected entry hit; ok=%v inv=%v", ok, h1.inv)
	}
	P1 := O.Add(D.Mul(h1.t))
	T1, ok := refract4(D, h1.Nw, 1/ior)
	if !ok {
		t.Fatalf("unexpected TIR at entry")
	}
	T1 = T1.Norm()

	O2 := bumpPoint(P1, T1, bumpShift)
	h2, ok := intersectRayCellPoly(O2, T1, &obj.cellPoly)
	if !ok || !h2.inv {
		t.Fatalf("expected inside->exit hit; ok=%v inv=%v", ok, h2.inv)
	}
	_ = O2.Add(T1.Mul(h2.t))

	T2, ok := refract4(T1, h2.Nw, ior)
	if !ok {
		t.Fatalf("unexpected TIR at exit")
	}
	T2 = T2.Norm()
	if T2.W <= 0 {
		t.Fatalf("expected +W after exit, got %+v", T2)
	}
}
