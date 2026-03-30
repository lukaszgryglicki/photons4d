package photons4d

import (
	"math"
	"testing"
)

func TestNewCell600AndIntersect(t *testing.T) {
	obj, err := NewCell600(
		Point4{0, 0, 0, 0.25},
		Vector4{.2, .2, .2, .2},
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.05, 1.05, 1.05},
	)
	if err != nil {
		t.Fatalf("NewCell600: %v", err)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayCellPoly(O, D, &obj.cellPoly)
	if !ok || h.t <= 0 {
		t.Fatalf("expected 600-cell hit, ok=%v t=%.12g", ok, h.t)
	}
	if math.Abs(float64(h.Nw.Len()-1)) > 1e-9 {
		t.Fatalf("normal not unit: %.12g", h.Nw.Len())
	}
}

func TestSceneAABBWithCell600(t *testing.T) {
	obj, err := NewCell600(
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{1.05 * .15, 0.95 * .15, 1.1 * .15, 0.85 * .15},
		Rot4{XY: 0.15, XZ: -0.2, ZW: 0.25},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatal(err)
	}

	if !(obj.AABBMin.X <= obj.AABBMax.X && obj.AABBMin.Y <= obj.AABBMax.Y &&
		obj.AABBMin.Z <= obj.AABBMax.Z && obj.AABBMin.W <= obj.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", obj.AABBMin, obj.AABBMax)
	}

	for i, v := range verts600Unit() {
		p := Vector4{v.X * obj.Scale.X, v.Y * obj.Scale.Y, v.Z * obj.Scale.Z, v.W * obj.Scale.W}
		p = obj.R.MulVec(p)
		w := Point4{obj.Center.X + p.X, obj.Center.Y + p.Y, obj.Center.Z + p.Z, obj.Center.W + p.W}
		if w.X < obj.AABBMin.X-1e-12 || w.X > obj.AABBMax.X+1e-12 ||
			w.Y < obj.AABBMin.Y-1e-12 || w.Y > obj.AABBMax.Y+1e-12 ||
			w.Z < obj.AABBMin.Z-1e-12 || w.Z > obj.AABBMax.Z+1e-12 ||
			w.W < obj.AABBMin.W-1e-12 || w.W > obj.AABBMax.W+1e-12 {
			t.Fatalf("transformed 600-cell vertex %d outside AABB: %+v not in [%+v .. %+v]", i, w, obj.AABBMin, obj.AABBMax)
		}
	}

	c := obj.Center
	if c.X < obj.AABBMin.X-1e-12 || c.X > obj.AABBMax.X+1e-12 ||
		c.Y < obj.AABBMin.Y-1e-12 || c.Y > obj.AABBMax.Y+1e-12 ||
		c.Z < obj.AABBMin.Z-1e-12 || c.Z > obj.AABBMax.Z+1e-12 ||
		c.W < obj.AABBMin.W-1e-12 || c.W > obj.AABBMax.W+1e-12 {
		t.Fatalf("center outside AABB: center=%+v aabb=[%+v..%+v]", c, obj.AABBMin, obj.AABBMax)
	}

	O := Point4{obj.AABBMin.X - 10, obj.AABBMin.Y - 10, obj.AABBMin.Z - 10, obj.AABBMin.W - 10}
	D := Vector4{1, 1, 1, 1}.Norm()
	rr := rayRecips{
		invX: 1 / D.X, invY: 1 / D.Y, invZ: 1 / D.Z, invW: 1 / D.W,
		parX: false, parY: false, parZ: false, parW: false,
	}
	ok, _ := rayAABB(O, obj.AABBMin, obj.AABBMax, rr)
	if !ok {
		t.Fatalf("rayAABB should intersect 600-cell AABB")
	}
}
