package photons4d

import (
	"math"
	"testing"
)

func TestNewSuperquadricAndIntersect(t *testing.T) {
	obj, err := NewSuperquadric(
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.3, 0.25, 0.15},
		4,
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewSuperquadric: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRaySuperquadric(O, D, obj)
	if !ok {
		t.Fatalf("expected superquadric hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.45)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.W+1)) > 1e-9 {
		t.Fatalf("expected -W normal, got %+v", h.Nw)
	}

	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRaySuperquadric(O2, D2, obj)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.15)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithSuperquadric(t *testing.T) {
	power := Real(4)
	obj, err := NewSuperquadric(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		power,
		Rot4{XY: 0.2, XZ: -0.15, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(obj.AABBMin.X <= obj.AABBMax.X && obj.AABBMin.Y <= obj.AABBMax.Y &&
		obj.AABBMin.Z <= obj.AABBMax.Z && obj.AABBMin.W <= obj.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", obj.AABBMin, obj.AABBMax)
	}

	c := math.Pow(0.25, 1/float64(power))
	localSamples := []Vector4{
		{+obj.Scale.X, 0, 0, 0},
		{-obj.Scale.X, 0, 0, 0},
		{0, +obj.Scale.Y, 0, 0},
		{0, -obj.Scale.Y, 0, 0},
		{0, 0, +obj.Scale.Z, 0},
		{0, 0, -obj.Scale.Z, 0},
		{0, 0, 0, +obj.Scale.W},
		{0, 0, 0, -obj.Scale.W},
		{obj.Scale.X * c, obj.Scale.Y * c, obj.Scale.Z * c, obj.Scale.W * c},
		{-obj.Scale.X * c, obj.Scale.Y * c, -obj.Scale.Z * c, obj.Scale.W * c},
	}
	for i, v := range localSamples {
		p := obj.R.MulVec(v)
		w := Point4{obj.Center.X + p.X, obj.Center.Y + p.Y, obj.Center.Z + p.Z, obj.Center.W + p.W}
		if w.X < obj.AABBMin.X-1e-12 || w.X > obj.AABBMax.X+1e-12 ||
			w.Y < obj.AABBMin.Y-1e-12 || w.Y > obj.AABBMax.Y+1e-12 ||
			w.Z < obj.AABBMin.Z-1e-12 || w.Z > obj.AABBMax.Z+1e-12 ||
			w.W < obj.AABBMin.W-1e-12 || w.W > obj.AABBMax.W+1e-12 {
			t.Fatalf("sample %d outside AABB: %+v not in [%+v .. %+v]", i, w, obj.AABBMin, obj.AABBMax)
		}
	}
}

// A ray that starts on the surface (within rounding of it) and leaves must not
// hit anything: the bounding box of the shape is not a part of it.
func TestSuperquadricNoHitWhenLeavingFromTheSurface(t *testing.T) {
	sq, err := NewSuperquadric(Point4{0.1, -0.2, 0.3, 0.5}, Vector4{0.5, 0.4, 0.6, 0.3}, 4, Rot4{XY: 0.3, XW: -0.4, ZW: 0.7},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.5, 0.5, 0.5}, RGB{0.5, 0.5, 0.5}, RGB{1.5, 1.5, 1.5})
	if err != nil {
		t.Fatal(err)
	}
	O := Point4{2, 1.5, -1, 2.5}
	D := Vector4{sq.Center.X - O.X + 0.05, sq.Center.Y - O.Y - 0.03, sq.Center.Z - O.Z + 0.02, sq.Center.W - O.W + 0.04}.Norm()
	h, ok := intersectRaySuperquadric(O, D, sq)
	if !ok || h.inv {
		t.Fatalf("expected entering hit, ok=%v inv=%v", ok, h.inv)
	}
	P := O.Add(D.Mul(h.t))
	for _, off := range []Real{1e-9, 1e-10, 1e-11, 3e-12, 1e-12, 1e-13, 0} {
		for _, dir := range []Vector4{h.Nw, reflect4(D, h.Nw).Norm()} {
			if h2, ok2 := intersectRaySuperquadric(P.Add(h.Nw.Mul(off)), dir, sq); ok2 {
				t.Fatalf("offset %g: ray leaving the surface hit at t=%.12g (inv=%v)", off, h2.t, h2.inv)
			}
		}
	}
	// the same point, going in: the exit on the far side
	h3, ok3 := intersectRaySuperquadric(P.Add(h.Nw.Mul(1e-13)), D, sq)
	if !ok3 || !h3.inv || h3.t < 0.1 {
		t.Fatalf("expected the far exit, ok=%v inv=%v t=%.12g", ok3, h3.inv, h3.t)
	}
}
