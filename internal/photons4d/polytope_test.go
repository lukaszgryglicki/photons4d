package photons4d

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"strings"
	"testing"
)

var testWhite, testNone, testHalf, testIOR = RGB{1, 1, 1}, RGB{}, RGB{0.5, 0.5, 0.5}, RGB{1.5, 1.5, 1.5}

// Order of the reflection group of a connected linear diagram.
func coxeterOrder(marks []Real) int {
	switch len(marks) {
	case 0:
		return 2
	case 1:
		return 2 * int(marks[0])
	}
	key := ""
	for _, m := range marks {
		key += fmt.Sprint(int(m))
	}
	if key > reverseString(key) {
		key = reverseString(key)
	}
	return map[string]int{"33": 24, "34": 48, "35": 120, "333": 120, "334": 384, "343": 1152, "335": 14400}[key]
}

func reverseString(s string) string {
	b := []byte(s)
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
	return string(b)
}

// Product of the orders of the groups of the runs of kept nodes; ringed tells
// whether every run holds a ringed node.
func runsOrder(marks []Real, rings []bool, keep func(int) bool) (order int, ringed bool) {
	order, ringed = 1, true
	for i := 0; i < len(rings); {
		if !keep(i) {
			i++
			continue
		}
		j, any := i, false
		for j < len(rings) && keep(j) {
			any = any || rings[j]
			j++
		}
		order *= coxeterOrder(marks[i : j-1])
		ringed = ringed && any
		i = j
	}
	return order, ringed
}

func diagramString(marks []Real, rings []bool) string {
	s := ""
	for i, r := range rings {
		if i > 0 {
			s += fmt.Sprint(int(marks[i-1]))
		}
		if r {
			s += "x"
		} else {
			s += "o"
		}
	}
	return s
}

// The vertices are the images of a point whose stabilizer is the group of
// the nodes without a ring; the cells of one kind are the images of the
// polyhedron of the diagram without one node, which is solid when each of
// its parts has a ring.
func TestWythoffPolychoraHaveTheCountsOfTheirSubgroups(t *testing.T) {
	for _, marks := range [][]Real{{3, 3, 3}, {4, 3, 3}, {3, 4, 3}, {5, 3, 3}} {
		order := coxeterOrder(marks)
		for pattern := 1; pattern < 16; pattern++ {
			rings := []bool{pattern&1 != 0, pattern&2 != 0, pattern&4 != 0, pattern&8 != 0}
			name := diagramString(marks, rings)
			verts, planes, err := uniformPolychoron(name)
			if err != nil {
				t.Fatalf("%s: %v", name, err)
			}
			stab, _ := runsOrder(marks, rings, func(i int) bool { return !rings[i] })
			cells := 0
			for gone := 0; gone < 4; gone++ {
				if sub, solid := runsOrder(marks, rings, func(i int) bool { return i != gone }); solid {
					cells += order / sub
				}
			}
			if len(verts) != order/stab || len(planes) != cells {
				t.Errorf("%s: %d vertices and %d cells, want %d and %d", name, len(verts), len(planes), order/stab, cells)
			}
			for _, v := range verts {
				if math.Abs(v.Len()-1) > 1e-12 {
					t.Fatalf("%s: vertex %v is not on the unit sphere", name, v)
				}
			}
		}
	}
}

// The cells of every Wythoffian polychoron and of its dual against the
// kaleidoscope: a point folded into the chamber of the mirrors sees one
// plane per mirror only.
func TestUniformPolychoraMatchTheKaleidoscope(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	for _, marks := range [][]Real{{3, 3, 3}, {4, 3, 3}, {3, 4, 3}, {5, 3, 3}} {
		for pattern := 1; pattern < 16; pattern++ {
			rings := []bool{pattern&1 != 0, pattern&2 != 0, pattern&4 != 0, pattern&8 != 0}
			name := diagramString(marks, rings)
			verts, planes, err := uniformPolychoron(name)
			if err != nil {
				t.Fatalf("%s: %v", name, err)
			}
			k := newKaleidoscope(t, name)
			dverts, dplanes := dualShape(verts, planes)
			far := Real(0)
			for _, v := range dverts {
				far = math.Max(far, v.Len())
			}
			if math.Abs(far-1) > 1e-12 {
				t.Errorf("dual %s: farthest vertex at %.15g, want 1", name, far)
			}
			gauge := func(ns []Vector4, q Vector4) Real {
				g := math.Inf(-1)
				for _, n := range ns {
					g = math.Max(g, n.Dot(q))
				}
				return g
			}
			n := 40
			if len(verts) > 2000 {
				n = 8
			}
			for i := 0; i < n; i++ {
				q := unitS3(rng).Mul(0.2 + rng.Float64())
				if g, want := gauge(planes, q), k.gauge(q); math.Abs(g-want) > 1e-9 {
					t.Fatalf("%s: gauge %.12g at %v, the kaleidoscope gives %.12g", name, g, q, want)
				}
				if g, want := gauge(dplanes, q), k.dualGauge(q); math.Abs(g-want) > 1e-9 {
					t.Fatalf("dual %s: gauge %.12g at %v, the kaleidoscope gives %.12g", name, g, q, want)
				}
			}
		}
	}
}

// Vertices and cells of the 47 convex uniform polychora that are not prisms,
// under their usual names (three of them have two names).
func TestNamedUniformPolychora(t *testing.T) {
	for _, c := range []struct {
		name  string
		v, c  int
		edges int // edges at a vertex
	}{
		{"5-cell", 5, 5, 4}, {"rectified 5-cell", 10, 10, 6}, {"truncated 5-cell", 20, 10, 4}, {"cantellated 5-cell", 30, 20, 6},
		{"runcinated 5-cell", 20, 30, 6}, {"bitruncated 5-cell", 30, 10, 4}, {"cantitruncated 5-cell", 60, 20, 4},
		{"runcitruncated 5-cell", 60, 30, 5}, {"omnitruncated 5-cell", 120, 30, 4},

		{"tesseract", 16, 8, 4}, {"rectified tesseract", 32, 24, 6}, {"truncated tesseract", 64, 24, 4},
		{"cantellated tesseract", 96, 56, 6}, {"runcinated tesseract", 64, 80, 6}, {"bitruncated tesseract", 96, 24, 4},
		{"cantitruncated tesseract", 192, 56, 4}, {"runcitruncated tesseract", 192, 80, 5}, {"omnitruncated tesseract", 384, 80, 4},
		{"16-cell", 8, 16, 6}, {"rectified 16-cell", 24, 24, 8}, {"truncated 16-cell", 48, 24, 5}, {"cantellated 16-cell", 96, 48, 6},
		{"cantitruncated 16-cell", 192, 48, 4}, {"runcitruncated 16-cell", 192, 80, 5},

		{"24-cell", 24, 24, 8}, {"rectified 24-cell", 96, 48, 6}, {"truncated 24-cell", 192, 48, 4}, {"cantellated 24-cell", 288, 144, 6},
		{"runcinated 24-cell", 144, 240, 8}, {"bitruncated 24-cell", 288, 48, 4}, {"cantitruncated 24-cell", 576, 144, 4},
		{"runcitruncated 24-cell", 576, 240, 5}, {"omnitruncated 24-cell", 1152, 240, 4}, {"snub 24-cell", 96, 144, 9},

		{"120-cell", 600, 120, 4}, {"rectified 120-cell", 1200, 720, 6}, {"truncated 120-cell", 2400, 720, 4},
		{"cantellated 120-cell", 3600, 1920, 6}, {"runcinated 120-cell", 2400, 2640, 6}, {"bitruncated 120-cell", 3600, 720, 4},
		{"cantitruncated 120-cell", 7200, 1920, 4}, {"runcitruncated 120-cell", 7200, 2640, 5}, {"omnitruncated 120-cell", 14400, 2640, 4},
		{"600-cell", 120, 600, 12}, {"rectified 600-cell", 720, 720, 10}, {"truncated 600-cell", 1440, 720, 6},
		{"cantellated 600-cell", 3600, 1440, 6}, {"cantitruncated 600-cell", 7200, 1440, 4}, {"runcitruncated 600-cell", 7200, 2640, 5},

		{"grand antiprism", 100, 320, 10},
	} {
		verts, planes, err := uniformPolychoron(c.name)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if len(verts) != c.v || len(planes) != c.c {
			t.Errorf("%s: %d vertices and %d cells, want %d and %d", c.name, len(verts), len(planes), c.v, c.c)
		}
		// uniform: all vertices on one sphere and alike, all edges equal
		if len(verts) <= 1500 {
			edge, degree := uniformEdges(verts)
			for i, d := range degree {
				if d != c.edges {
					t.Errorf("%s: vertex %d has %d edges of the length %.6g, want %d", c.name, i, d, edge, c.edges)
					break
				}
			}
		}
	}
	// names written in other ways, and names that do not exist
	for _, pair := range [][2]string{
		{"  Truncated   120-Cell ", "x5x3o3o"}, {"pentachoron", "x3o3o3o"}, {"hypercube", "x4o3o3o"},
		{"runcicantitruncated 24-cell", "x3x4x3x"}, {"birectified 600-cell", "o3o3x5o"}, {"runcicantellated tesseract", "x4o3x3x"},
		{"truncated hexacosichoron", "x3x3o5o"}, {"rectified icositetrachoron", "o3x4o3o"}, {"hecatonicosachoron", "x5o3o3o"},
		{"cantellated hexadecachoron", "x3o3x4o"}, {"octachoron", "x4o3o3o"},
	} {
		a, _, errA := uniformPolychoron(pair[0])
		b, _, errB := uniformPolychoron(pair[1])
		if errA != nil || errB != nil || !oracleSameSet(a, b) {
			t.Errorf("%q and %q differ (%v, %v)", pair[0], pair[1], errA, errB)
		}
	}
	for _, bad := range []string{"", "truncated", "snub 120-cell", "rectified cube", "x5/2o5o3o", "x4o3o4o", "x5o3o5o", "x3o3o", "x3o3o3o3o", "x1o3o3o", "y3o3o3o"} {
		if _, _, err := uniformPolychoron(bad); err == nil {
			t.Errorf("%q must not name a polychoron", bad)
		}
	}
}

// uniformEdges returns the least distance between two points and, for each
// point, the number of points at that distance.
func uniformEdges(verts []Vector4) (Real, []int) {
	edge := math.Inf(1)
	for i := range verts {
		for j := i + 1; j < len(verts); j++ {
			edge = math.Min(edge, verts[i].Sub(verts[j]).Len())
		}
	}
	degree := make([]int, len(verts))
	for i := range verts {
		for j := range verts {
			if i != j && verts[i].Sub(verts[j]).Len() < edge*(1+1e-9) {
				degree[i]++
			}
		}
	}
	return edge, degree
}

func cellSizes(verts, planes []Vector4) map[int]int {
	sizes := map[int]int{}
	for _, n := range planes {
		on := 0
		for _, v := range verts {
			if math.Abs(n.Dot(v)-1) < 1e-9 {
				on++
			} else if n.Dot(v) > 1 {
				return nil
			}
		}
		sizes[on]++
	}
	return sizes
}

func TestSnub24CellAndGrandAntiprismCells(t *testing.T) {
	verts, planes, err := uniformPolychoron("snub 24-cell")
	if err != nil {
		t.Fatal(err)
	}
	if s := cellSizes(verts, planes); len(s) != 2 || s[4] != 120 || s[12] != 24 {
		t.Errorf("snub 24-cell: cells by the number of vertices %v, want 120 tetrahedra and 24 icosahedra", s)
	}
	verts, planes, err = uniformPolychoron("grand antiprism")
	if err != nil {
		t.Fatal(err)
	}
	if s := cellSizes(verts, planes); len(s) != 2 || s[4] != 300 || s[10] != 20 {
		t.Errorf("grand antiprism: cells by the number of vertices %v, want 300 tetrahedra and 20 pentagonal antiprisms", s)
	}
	// the vertices left out are those of the independent constructions
	for _, c := range []struct {
		name string
		kept []Vector4
		gone []Vector4
	}{
		{"snub 24-cell", snub24Vertices(), oracleHurwitzUnits()},
		{"grand antiprism", grandAntiprismVertices(), oracleDecagonRings()},
	} {
		all := append(append([]Vector4{}, c.kept...), c.gone...)
		if !oracleSameSet(all, oracleIcosians()) {
			t.Errorf("%s: its %d vertices and the %d left out are not the 600-cell", c.name, len(c.kept), len(c.gone))
		}
	}
}

func TestUniformPolyhedra(t *testing.T) {
	for _, c := range []struct {
		name         string
		v, f, degree int
		radius       Real // circumradius over edge
	}{
		{"tetrahedron", 4, 4, 3, math.Sqrt(6) / 4}, {"cube", 8, 6, 3, math.Sqrt(3) / 2}, {"octahedron", 6, 8, 4, math.Sqrt2 / 2},
		{"dodecahedron", 20, 12, 3, math.Sqrt(3) * (1 + math.Sqrt(5)) / 4}, {"icosahedron", 12, 20, 5, math.Sqrt(10+2*math.Sqrt(5)) / 4},
		{"truncated tetrahedron", 12, 8, 3, math.Sqrt(22) / 4}, {"cuboctahedron", 12, 14, 4, 1},
		{"truncated cube", 24, 14, 3, math.Sqrt(7+4*math.Sqrt2) / 2}, {"truncated octahedron", 24, 14, 3, math.Sqrt(10) / 2},
		{"rhombicuboctahedron", 24, 26, 4, math.Sqrt(5+2*math.Sqrt2) / 2}, {"truncated cuboctahedron", 48, 26, 3, math.Sqrt(13+6*math.Sqrt2) / 2},
		{"snub cube", 24, 38, 5, 1.3437133737446}, {"icosidodecahedron", 30, 32, 4, (1 + math.Sqrt(5)) / 2},
		{"truncated dodecahedron", 60, 32, 3, math.Sqrt(74+30*math.Sqrt(5)) / 4}, {"truncated icosahedron", 60, 32, 3, math.Sqrt(58+18*math.Sqrt(5)) / 4},
		{"rhombicosidodecahedron", 60, 62, 4, math.Sqrt(11+4*math.Sqrt(5)) / 2}, {"truncated icosidodecahedron", 120, 62, 3, math.Sqrt(31+12*math.Sqrt(5)) / 2},
		{"snub dodecahedron", 60, 92, 5, 2.155837375115},
		{"hexahedron", 8, 6, 3, 0}, {"x5x3o", 60, 32, 3, 0}, {"x4o2x", 8, 6, 3, 0},
		{"triangular prism", 6, 5, 3, 0}, {"square prism", 8, 6, 3, math.Sqrt(3) / 2}, {"7-gonal prism", 14, 9, 3, 0}, {"12-prism", 24, 14, 3, 0},
		{"3-gonal antiprism", 6, 8, 4, math.Sqrt2 / 2}, {"square antiprism", 8, 10, 4, 0}, {"pentagonal antiprism", 10, 12, 4, 0},
		{"17-antiprism", 34, 36, 4, 0}, {"Dodecagonal  Antiprism", 24, 26, 4, 0},
	} {
		verts, planes, err := uniformPolyhedron(c.name)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if len(verts) != c.v || len(planes) != c.f {
			t.Errorf("%s: %d vertices and %d faces, want %d and %d", c.name, len(verts), len(planes), c.v, c.f)
		}
		edge, degree := uniformEdges(verts)
		for i, v := range verts {
			if math.Abs(v.Len()-1) > 1e-12 || v.W != 0 || degree[i] != c.degree {
				t.Errorf("%s: vertex %v at the distance %.15g with %d edges, want 1 and %d", c.name, v, v.Len(), degree[i], c.degree)
				break
			}
		}
		if c.radius != 0 && math.Abs(1/edge-c.radius) > 1e-10 {
			t.Errorf("%s: circumradius over edge %.13g, want %.13g", c.name, 1/edge, c.radius)
		}
		// the faces found from the edges are those of the hull
		if hull := oracleHull3(verts); !oracleSameSet(hull, planes) {
			t.Errorf("%s: %d faces, the hull has %d", c.name, len(planes), len(hull))
		}
	}
	for _, bad := range []string{"", "prism", "2-gonal prism", "1000-gonal prism", "snub tetrahedron", "x5/2o5o", "x4o4o", "x5o3o3o", "great icosahedron"} {
		if _, _, err := uniformPolyhedron(bad); err == nil {
			t.Errorf("%q must not name a polyhedron", bad)
		}
	}
}

// The prisms over the uniform polyhedra are uniform polychora: their
// vertices lie on a sphere and all their edges are equal. The pyramids and
// bipyramids stand on the base of the circumradius 1.
func TestPolyhedralSolidsShape(t *testing.T) {
	for _, base := range []string{"tetrahedron", "cube", "truncated icosahedron", "snub cube", "octagonal prism", "heptagonal antiprism"} {
		bverts, bplanes, err := uniformPolyhedron(base)
		if err != nil {
			t.Fatal(err)
		}
		_, bdeg := uniformEdges(bverts)
		verts, planes, err := liftPolyhedron("prism", bverts, bplanes)
		if err != nil {
			t.Fatal(err)
		}
		_, degree := uniformEdges(verts)
		radius := verts[0].Len()
		for i, v := range verts {
			if math.Abs(v.Len()-radius) > 1e-12 || degree[i] != bdeg[0]+1 {
				t.Errorf("%s prism: vertex %v at %.15g with %d edges, want %.15g and %d", base, v, v.Len(), degree[i], radius, bdeg[0]+1)
				break
			}
		}
		if len(verts) != 2*len(bverts) || len(planes) != len(bplanes)+2 || cellSizes(verts, planes) == nil {
			t.Errorf("%s prism: %d vertices, %d cells", base, len(verts), len(planes))
		}

		verts, planes, _ = liftPolyhedron("pyramid", bverts, bplanes)
		sizes := cellSizes(verts, planes)
		if len(verts) != len(bverts)+1 || len(planes) != len(bplanes)+1 || sizes == nil || sizes[len(bverts)] < 1 {
			t.Errorf("%s pyramid: %d vertices, %d cells, cells by vertices %v", base, len(verts), len(planes), sizes)
		}
		verts, planes, _ = liftPolyhedron("bipyramid", bverts, bplanes)
		if len(verts) != len(bverts)+2 || len(planes) != 2*len(bplanes) || cellSizes(verts, planes) == nil {
			t.Errorf("%s bipyramid: %d vertices, %d cells", base, len(verts), len(planes))
		}
	}
	if _, _, err := liftPolyhedron("cone", nil, nil); err == nil {
		t.Errorf("an unknown lift must be rejected")
	}
	if _, err := NewPolyhedralSolid("prism", "sphere", Point4{}, Vector4{1, 1, 1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("an unknown base must be rejected")
	}
	if _, err := NewPolyhedralSolid("prism", "cube", Point4{}, Vector4{1, 0, 1, 1}, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("a zero scale must be rejected")
	}
	if _, err := NewPolychoron("24-cell", false, Point4{}, Vector4{1, 1, 1, 1}, Rot4{}, testWhite, testHalf, testHalf, testHalf, testIOR); err == nil {
		t.Errorf("a material that returns more light than it gets must be rejected")
	}
}

// sameHits sends rays at two objects that must be the same solid.
func sameHits(t *testing.T, name string, rng *rand.Rand, reach Real, a, b func(Point4, Vector4) (objectHit, bool)) {
	t.Helper()
	hits := 0
	for i := 0; i < 600; i++ {
		O := Point4{}.Add(unitS3(rng).Mul(reach * (0.05 + 1.6*rng.Float64())))
		T := Point4{}.Add(unitS3(rng).Mul(reach * 0.7 * rng.Float64()))
		D := Vector4{T.X - O.X, T.Y - O.Y, T.Z - O.Z, T.W - O.W}.Norm()
		ha, oka := a(O, D)
		hb, okb := b(O, D)
		if oka != okb {
			t.Fatalf("%s: ray %d O=%v D=%v: hit %v against %v", name, i, O, D, oka, okb)
		}
		if !oka {
			continue
		}
		hits++
		if math.Abs(ha.t-hb.t) > 1e-9*(1+ha.t) || ha.inv != hb.inv {
			t.Fatalf("%s: ray %d O=%v D=%v: t=%.12g inv=%v against t=%.12g inv=%v", name, i, O, D, ha.t, ha.inv, hb.t, hb.inv)
		}
		// at an edge either normal is right
		if ha.Nw.Sub(hb.Nw).Len() > 1e-6 {
			if ga, _ := a(O.Add(D.Mul(ha.t-1e-4*reach)), D); math.Abs(ga.t-1e-4*reach) < 1e-9 {
				t.Fatalf("%s: ray %d O=%v D=%v: normal %v against %v", name, i, O, D, ha.Nw, hb.Nw)
			}
		}
	}
	if hits < 200 {
		t.Fatalf("%s: only %d hits", name, hits)
	}
}

func TestDuopyramid(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	rot := Rot4{XY: 0.3, XW: -0.7, YZ: 1.1, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	// the 4-4 duopyramid is the 16-cell |x|+|y|+|z|+|w| <= 1, which Cell16 writes with 1/sqrt(2)
	d, err := NewDuopyramid(Point4{}, sc, 4, 4, rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	c, err := NewCell16(Point4{}, sc.Mul(math.Sqrt2), rot, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	sameHits(t, "4-4 duopyramid and 16-cell", rng, 1.3,
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &d.cellPoly) },
		func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell16(O, D, c) })

	// the dual of the duoprism: a vertex for each of its P + Q cells and a
	// tetrahedron, between an edge of each polygon, for each of its P·Q vertices
	d, err = NewDuopyramid(Point4{}, Vector4{1, 1, 1, 1}, 5, 7, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR)
	if err != nil {
		t.Fatal(err)
	}
	var verts []Vector4
	for k := 0; k < 5; k++ {
		verts = append(verts, Vector4{math.Cos(2 * math.Pi * Real(k) / 5), math.Sin(2 * math.Pi * Real(k) / 5), 0, 0})
	}
	for k := 0; k < 7; k++ {
		verts = append(verts, Vector4{0, 0, math.Cos(2 * math.Pi * Real(k) / 7), math.Sin(2 * math.Pi * Real(k) / 7)})
	}
	if s := cellSizes(verts, d.N); len(d.N) != 35 || len(s) != 1 || s[4] != 35 {
		t.Fatalf("5-7 duopyramid: %d cells, by the number of vertices %v, want 35 tetrahedra", len(d.N), s)
	}
	for _, pq := range [][2]int{{2, 5}, {5, 2}, {200, 200}} {
		if _, err := NewDuopyramid(Point4{}, Vector4{1, 1, 1, 1}, pq[0], pq[1], Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR); err == nil {
			t.Errorf("a %d-%d duopyramid must be rejected", pq[0], pq[1])
		}
	}
}

// The AABB of a polytope is that of its vertices.
func TestPolytopeAABBIsTight(t *testing.T) {
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	ctr := Point4{0.3, -0.2, 0.5, 0.1}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	for _, kind := range []string{"cantellated 5-cell", "snub 24-cell", "runcinated 120-cell"} {
		for _, dual := range []bool{false, true} {
			o, err := NewPolychoron(kind, dual, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
			if err != nil {
				t.Fatal(err)
			}
			verts, planes, _ := uniformPolychoron(kind)
			if dual {
				verts, _ = dualShape(verts, planes)
			}
			lo, hi := [4]Real{1e9, 1e9, 1e9, 1e9}, [4]Real{-1e9, -1e9, -1e9, -1e9}
			for _, v := range verts {
				p := o.R.MulVec(Vector4{v.X * sc.X, v.Y * sc.Y, v.Z * sc.Z, v.W * sc.W})
				for i, x := range []Real{ctr.X + p.X, ctr.Y + p.Y, ctr.Z + p.Z, ctr.W + p.W} {
					lo[i], hi[i] = math.Min(lo[i], x), math.Max(hi[i], x)
				}
			}
			got := [8]Real{o.AABBMin.X, o.AABBMin.Y, o.AABBMin.Z, o.AABBMin.W, o.AABBMax.X, o.AABBMax.Y, o.AABBMax.Z, o.AABBMax.W}
			want := [8]Real{lo[0], lo[1], lo[2], lo[3], hi[0], hi[1], hi[2], hi[3]}
			for i := range got {
				if math.Abs(got[i]-want[i]) > 1e-12 {
					t.Errorf("%s dual=%v: AABB %v, the vertices span %v", kind, dual, got, want)
					break
				}
			}
		}
	}
}

func TestPolychoronNamesAreListed(t *testing.T) {
	// every operation on every regular polychoron builds
	var names []string
	for base := range polychoronBases {
		for op := range polychoronOps {
			names = append(names, strings.TrimSpace(op+" "+base))
		}
	}
	sort.Strings(names)
	for _, n := range names {
		if _, _, err := uniformPolychoron(n); err != nil {
			t.Errorf("%s: %v", n, err)
		}
	}
}
