package photons4d

import "testing"

// The 120-cell and 600-cell vertex/normal sets must be exact polar duals:
// scaling each vertex direction by 1/max_i(n_i·v̂) must use a CONSTANT scale
// (all vertices equidistant from the center), every scaled vertex must lie
// on or inside every plane, and each facet must touch exactly 20 vertices
// (dodecahedron, 120-cell) or 4 vertices (tetrahedron, 600-cell).
func TestCellDualPolytopeRegularity(t *testing.T) {
	v600 := verts600Unit()
	v120 := verts120Unit()
	if len(v600) != 120 {
		t.Fatalf("verts600Unit: want 120 vertices, got %d", len(v600))
	}
	if len(v120) != 600 {
		t.Fatalf("verts120Unit: want 600 vertices, got %d", len(v120))
	}

	check := func(name string, normals, dirs []Vector4, wantFacetVerts int) {
		lmin, lmax := 1e300, -1e300
		scaled := make([]Vector4, len(dirs))
		for i, v := range dirs {
			sup := -1e300
			for _, n := range normals {
				if d := n.Dot(v); d > sup {
					sup = d
				}
			}
			l := 1.0 / sup
			if l < lmin {
				lmin = l
			}
			if l > lmax {
				lmax = l
			}
			scaled[i] = v.Mul(l)
		}
		if lmax-lmin > 1e-9 {
			t.Errorf("%s: vertex support scale not constant: [%.12f .. %.12f]", name, lmin, lmax)
		}
		for _, n := range normals {
			c := 0
			for _, v := range scaled {
				d := n.Dot(v)
				if d > 1+1e-9 {
					t.Fatalf("%s: vertex outside facet plane by %.3g", name, d-1)
				}
				if d > 1-1e-9 {
					c++
				}
			}
			if c != wantFacetVerts {
				t.Fatalf("%s: facet touches %d vertices, want %d", name, c, wantFacetVerts)
			}
		}
	}
	check("120-cell", v600, v120, 20)
	check("600-cell", v120, v600, 4)
}
