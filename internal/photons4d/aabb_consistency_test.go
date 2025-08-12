package photons4d

import "testing"

// For a diagonal ray fired from well outside AABBMin, rayAABB should report an intersection.
func TestAABB_RayConsistency_CellPolys(t *testing.T) {
	makeRR := func(d Vector4) rayRecips {
		rr := rayRecips{}
		const eps = 1e-18
		if abs := d.X; abs > eps || abs < -eps {
			rr.invX = 1 / d.X
		} else {
			rr.parX = true
		}
		if abs := d.Y; abs > eps || abs < -eps {
			rr.invY = 1 / d.Y
		} else {
			rr.parY = true
		}
		if abs := d.Z; abs > eps || abs < -eps {
			rr.invZ = 1 / d.Z
		} else {
			rr.parZ = true
		}
		if abs := d.W; abs > eps || abs < -eps {
			rr.invW = 1 / d.W
		} else {
			rr.parW = true
		}
		return rr
	}

	check := func(name string, minP, maxP Point4) {
		O := Point4{minP.X - 10, minP.Y - 10, minP.Z - 10, minP.W - 10}
		D := Vector4{1, 1, 1, 1}.Norm()
		ok, _ := rayAABB(O, minP, maxP, makeRR(D))
		if !ok {
			t.Fatalf("%s: rayAABB should intersect the AABB", name)
		}
	}

	// 120-cell
	cell120, err := NewCell120(Point4{0.12, -0.2, 0.05, 0.33}, 0.18, Vector4{1, 1, 1, 1}, Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.1, 1.1, 1.1})
	if err != nil {
		t.Fatal(err)
	}
	check("120-cell", cell120.AABBMin, cell120.AABBMax)

	// 600-cell
	cell600, err := NewCell600(Point4{-0.15, 0.22, -0.07, 0.25}, 0.18, Vector4{1, 1, 1, 1}, Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.1, 1.1, 1.1})
	if err != nil {
		t.Fatal(err)
	}
	check("600-cell", cell600.AABBMin, cell600.AABBMax)
}
