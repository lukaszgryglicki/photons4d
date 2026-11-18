package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestNewDuoprismAndIntersect(t *testing.T) {
	// Two squares with the corners on the axes: |x|+|y| <= 0.3 and |z|+|w-0.6| <= 0.2.
	obj, err := NewDuoprism(
		Point4{0, 0, 0, 0.6},
		Vector4{0.3, 0.3, 0.2, 0.2},
		4, 4,
		Rot4{},
		RGB{1, 1, 1},
		RGB{0, 0, 0},
		RGB{0, 0, 0},
		RGB{1, 1, 1},
		RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewDuoprism: %v", err)
	}
	if len(obj.N) != 8 {
		t.Fatalf("4-4 duoprism has %d cells, want 8", len(obj.N))
	}

	s := Real(math.Sqrt(0.5))
	h, ok := intersectRayCellPoly(Point4{0, 0, 0.05, 0}, Vector4{0, 0, 0, 1}, &obj.cellPoly)
	if !ok || h.inv {
		t.Fatalf("expected entering hit, ok=%v inv=%v", ok, h.inv)
	}
	if math.Abs(float64(h.t-0.45)) > 1e-9 || h.Nw.Sub(Vector4{0, 0, s, -s}).Len() > 1e-9 {
		t.Fatalf("entering: t=%.12g N=%+v", h.t, h.Nw)
	}
	h2, ok2 := intersectRayCellPoly(Point4{0, 0, 0.05, 0.6}, Vector4{0, 0, 0, 1}, &obj.cellPoly)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.15)) > 1e-9 || h2.Nw.Sub(Vector4{0, 0, s, s}).Len() > 1e-9 {
		t.Fatalf("leaving: t=%.12g N=%+v", h2.t, h2.Nw)
	}
	h3, ok3 := intersectRayCellPoly(Point4{-1, 0.1, 0, 0.6}, Vector4{1, 0, 0, 0}, &obj.cellPoly)
	if !ok3 || h3.inv || math.Abs(float64(h3.t-0.8)) > 1e-9 || h3.Nw.Sub(Vector4{-s, s, 0, 0}).Len() > 1e-9 {
		t.Fatalf("XY side: ok=%v inv=%v t=%.12g N=%+v", ok3, h3.inv, h3.t, h3.Nw)
	}

	for _, pq := range [][2]int{{2, 4}, {4, 2}, {0, 0}, {361, 3}} {
		if _, err := NewDuoprism(Point4{}, Vector4{1, 1, 1, 1}, pq[0], pq[1], Rot4{}, RGB{1, 1, 1}, RGB{}, RGB{}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2}); err == nil {
			t.Fatalf("p=%d q=%d accepted", pq[0], pq[1])
		}
	}
}

// Every vertex lies on two sides of each polygon and inside all other cells,
// and the AABB is the box of the vertices.
func TestDuoprismVerticesCellsAndAABB(t *testing.T) {
	for _, pq := range [][2]int{{3, 3}, {3, 5}, {6, 4}, {7, 12}} {
		p, q := pq[0], pq[1]
		obj, err := NewDuoprism(
			Point4{0.1, -0.2, 0.05, 0.3},
			Vector4{0.24, 0.18, 0.21, 0.16},
			p, q,
			Rot4{XY: 0.2, XZ: -0.15, XW: 0.4, YW: -0.5, ZW: 0.3},
			RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
		)
		if err != nil {
			t.Fatal(err)
		}
		if len(obj.N) != p+q {
			t.Fatalf("%d-%d: %d cells, want %d", p, q, len(obj.N), p+q)
		}
		lo := [4]Real{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
		hi := [4]Real{math.Inf(-1), math.Inf(-1), math.Inf(-1), math.Inf(-1)}
		for i := 0; i < p; i++ {
			a := 2 * math.Pi * Real(i) / Real(p)
			for j := 0; j < q; j++ {
				b := 2 * math.Pi * Real(j) / Real(q)
				v := obj.R.MulVec(Vector4{math.Cos(a) * obj.Scale.X, math.Sin(a) * obj.Scale.Y, math.Cos(b) * obj.Scale.Z, math.Sin(b) * obj.Scale.W})
				w := Vector4{obj.Center.X + v.X, obj.Center.Y + v.Y, obj.Center.Z + v.Z, obj.Center.W + v.W}
				on := 0
				for k := range obj.N {
					d := (obj.N[k].Dot(w) - obj.D[k]) / obj.N[k].Len()
					if d > 1e-12 {
						t.Fatalf("%d-%d: vertex (%d,%d) is outside cell %d by %g", p, q, i, j, k, d)
					}
					if d > -1e-12 {
						on++
					}
				}
				if on != 4 {
					t.Fatalf("%d-%d: vertex (%d,%d) lies on %d cells, want 4", p, q, i, j, on)
				}
				for c, x := range [4]Real{w.X, w.Y, w.Z, w.W} {
					lo[c] = math.Min(lo[c], x)
					hi[c] = math.Max(hi[c], x)
				}
			}
		}
		min := [4]Real{obj.AABBMin.X, obj.AABBMin.Y, obj.AABBMin.Z, obj.AABBMin.W}
		max := [4]Real{obj.AABBMax.X, obj.AABBMax.Y, obj.AABBMax.Z, obj.AABBMax.W}
		for c := 0; c < 4; c++ {
			if math.Abs(lo[c]-min[c]) > 1e-12 || math.Abs(hi[c]-max[c]) > 1e-12 {
				t.Fatalf("%d-%d axis %d: AABB [%g, %g], vertices span [%g, %g]", p, q, c, min[c], max[c], lo[c], hi[c])
			}
		}
	}
}

// The 4-4 duoprism is an 8-cell turned by 45° in both planes, and a duoprism
// of two 360-gons is within 4e-5 of its scale from the duocylinder.
func TestDuoprismMatchesCell8AndDuocylinder(t *testing.T) {
	ctr := Point4{0.1, -0.2, 0.05, 0.3}
	white, none, half, ior := RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.5, 0.5, 0.5}, RGB{1.5, 1.5, 1.5}
	rot := Rot4{XZ: 0.3, XW: -0.2, YZ: 0.5, YW: 0.1}
	sq, err := NewDuoprism(ctr, Vector4{0.3, 0.3, 0.2, 0.2}, 4, 4, rot, white, none, half, half, ior)
	if err != nil {
		t.Fatal(err)
	}
	box, err := NewCell8(ctr, Vector4{0.3 * math.Sqrt2, 0.3 * math.Sqrt2, 0.2 * math.Sqrt2, 0.2 * math.Sqrt2}, rot, white, none, half, half, ior)
	if err != nil {
		t.Fatal(err)
	}
	// the 8-cell gets the rotation of the duoprism after an eighth of a turn in each of the two planes
	box.R = sq.R.Mul(rotXY(math.Pi / 4).Mul(rotZW(math.Pi / 4)))
	box.RT = box.R.Transpose()
	for k := 0; k < 4; k++ {
		box.Normals[k] = Vector4{box.R.M[0][k], box.R.M[1][k], box.R.M[2][k], box.R.M[3][k]}
	}

	round, err := NewDuoprism(ctr, Vector4{0.3, 0.25, 0.2, 0.35}, 360, 360, rot, white, none, half, half, ior)
	if err != nil {
		t.Fatal(err)
	}
	cyl, err := NewDuocylinder(ctr, Vector4{0.3, 0.25, 0.2, 0.35}, rot, white, none, half, half, ior)
	if err != nil {
		t.Fatal(err)
	}

	rng := rand.New(rand.NewSource(5))
	hits := 0
	for i := 0; i < 4000; i++ {
		O := ctr.Add(unitS3(rng).Mul(2))
		T := ctr.Add(unitS3(rng).Mul(0.3 * rng.Float64()))
		D := Vector4{T.X - O.X, T.Y - O.Y, T.Z - O.Z, T.W - O.W}.Norm()

		h1, ok1 := intersectRayCellPoly(O, D, &sq.cellPoly)
		h2, ok2 := intersectRayCell8(O, D, box)
		if ok1 != ok2 {
			t.Fatalf("ray %d: 4-4 duoprism hit=%v, 8-cell hit=%v", i, ok1, ok2)
		}
		if ok1 {
			hits++
			if math.Abs(float64(h1.t-h2.t)) > 1e-12 || h1.inv != h2.inv || h1.Nw.Sub(h2.Nw).Len() > 1e-12 {
				t.Fatalf("ray %d: 4-4 duoprism t=%.15g N=%+v, 8-cell t=%.15g N=%+v", i, h1.t, h1.Nw, h2.t, h2.Nw)
			}
		}

		h3, ok3 := intersectRayCellPoly(O, D, &round.cellPoly)
		h4, ok4 := intersectRayDuocylinder(O, D, cyl)
		if ok3 && ok4 && math.Abs(float64(h3.t-h4.t)) > 1e-3 {
			t.Fatalf("ray %d: 360-360 duoprism t=%.9g, duocylinder t=%.9g", i, h3.t, h4.t)
		}
	}
	if hits < 1000 {
		t.Fatalf("only %d of 4000 rays hit", hits)
	}
}
