package photons4d

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"sync"
	"testing"
)

// The ten regular star polychora with their numbers of vertices, edges,
// faces and cells and their densities (how many times they wind around the
// centre), and the six convex regular polychora.
var regularPolychoraTable = []struct {
	name, symbol        string
	v, e, f, c, density int
	cellDensity         int
}{
	{"icosahedral 120-cell", "{3,5,5/2}", 120, 720, 1200, 120, 4, 1},
	{"small stellated 120-cell", "{5/2,5,3}", 120, 1200, 720, 120, 4, 3},
	{"great 120-cell", "{5,5/2,5}", 120, 720, 720, 120, 6, 3},
	{"grand 120-cell", "{5,3,5/2}", 120, 720, 720, 120, 20, 1},
	{"great stellated 120-cell", "{5/2,3,5}", 120, 720, 720, 120, 20, 7},
	{"grand stellated 120-cell", "{5/2,5,5/2}", 120, 720, 720, 120, 66, 3},
	{"great grand 120-cell", "{5,5/2,3}", 120, 1200, 720, 120, 76, 3},
	{"great icosahedral 120-cell", "{3,5/2,5}", 120, 720, 1200, 120, 76, 7},
	{"grand 600-cell", "{3,3,5/2}", 120, 720, 1200, 600, 191, 1},
	{"great grand stellated 120-cell", "{5/2,3,3}", 600, 1200, 720, 120, 191, 7},
	{"", "{3,3,3}", 5, 10, 10, 5, 1, 1},
	{"", "{4,3,3}", 16, 32, 24, 8, 1, 1},
	{"", "{3,3,4}", 8, 24, 32, 16, 1, 1},
	{"", "{3,4,3}", 24, 96, 96, 24, 1, 1},
	{"", "{5,3,3}", 600, 1200, 720, 120, 1, 1},
	{"", "{3,3,5}", 120, 720, 1200, 600, 1, 1},
}

func TestRegularPolychoraMatchThePublishedTables(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	for _, c := range regularPolychoraTable {
		s, err := regularStarShape(c.symbol)
		if err != nil {
			t.Fatalf("%s: %v", c.symbol, err)
		}
		if len(s.verts) != c.v || s.edgeCount != c.e || s.faceCount != c.f || len(s.normals) != c.c {
			t.Errorf("%s: %d vertices, %d edges, %d faces, %d cells; want %d, %d, %d, %d",
				c.symbol, len(s.verts), s.edgeCount, s.faceCount, len(s.normals), c.v, c.e, c.f, c.c)
		}
		// Euler: the four numbers of a convex polychoron cancel out (those of
		// four of the star polychora do not: they give 480 or -480)
		if c.density == 1 && len(s.verts)-s.edgeCount+s.faceCount-len(s.normals) != 0 {
			t.Errorf("%s: V - E + F - C is not 0", c.symbol)
		}
		for i := 0; i < 50; i++ {
			if d := s.density(unitS3(rng).Mul(1e-3)); d != c.density {
				t.Fatalf("%s: density %d around the centre, want %d", c.symbol, d, c.density)
			}
			a, b, w := sampleS2(rng)
			if d := s.cellDensity([3]Real{a * 1e-4, b * 1e-4, w * 1e-4}); d != c.cellDensity {
				t.Fatalf("%s: density %d around the centre of a cell, want %d", c.symbol, d, c.cellDensity)
			}
			if d := s.density(unitS3(rng).Mul(1 + 1e-9)); d != 0 {
				t.Fatalf("%s: density %d outside the sphere of the vertices", c.symbol, d)
			}
		}
		if c.name != "" {
			named, err := regularStarShape(c.name)
			if err != nil || named != s {
				t.Errorf("%q does not name %s (%v)", c.name, c.symbol, err)
			}
		}
	}
	for _, bad := range []string{"", "great 600-cell", "{5/2,5}", "{3,5/2,3}", "{4,3,5/2}", "{5/2,3,5/2}", "{5,3,5}", "{4,3,4}", "{2,3,3}", "{3,3,3,3}", "5/2,5,3"} {
		if _, err := regularStarShape(bad); err == nil {
			t.Errorf("%q must not name a regular polychoron", bad)
		}
	}
}

// The density at a point is the number of the simplices (centre, cell
// centre, face centre, edge) that hold it; here they are listed one by one.
func TestStarPolychoronDensityAgainstFlagSimplices(t *testing.T) {
	rng := rand.New(rand.NewSource(12))
	for _, c := range regularPolychoraTable {
		s, err := regularStarShape(c.symbol)
		if err != nil {
			t.Fatal(err)
		}
		flags := newOracleFlags(t, c.symbol)
		// two cells at each face, and a simplex for each edge of the face
		marks, _ := parseSchlafli(c.symbol)
		sides := map[Real]int{3: 3, 4: 4, 5: 5, 2.5: 5}[marks[0]]
		n := 0
		for _, w := range flags.walls {
			n += len(w)
		}
		if len(flags.cells) != c.c || n != 2*c.f*sides || len(s.flags)*c.c != n {
			t.Fatalf("%s: %d simplices in %d cells, want %d in %d", c.symbol, n, len(flags.cells), 2*c.f*sides, c.c)
		}
		// the tree over the mirror planes of a cell gives the count of the
		// simplices that hold a point of the cell
		for i := 0; i < 2000; i++ {
			a, b, w := sampleS2(rng)
			r := 1.1 * math.Sqrt(1-s.height*s.height) * math.Pow(rng.Float64(), 0.4)
			y := [3]Real{a * r, b * r, w * r}
			var word uint32
			for j, n := range s.walls {
				if dot3(n, y) >= 0 {
					word |= 1 << uint(j)
				}
			}
			want := 0
			for _, f := range s.flags {
				if word&f.mask == f.want && dot3(s.faces[f.face], y) <= s.faceDist {
					want++
				}
			}
			if got := s.cellDensity(y); got != want {
				t.Fatalf("%s: cell density %d at %v, %d simplices hold the point", c.symbol, got, y, want)
			}
		}
		inside := 0
		for i := 0; i < 300; i++ {
			q := unitS3(rng).Mul(math.Pow(rng.Float64(), 0.4) * 1.02)
			got, want := s.density(q), flags.count(q)
			if got != want {
				t.Fatalf("%s: density %d at %v, %d simplices hold the point", c.symbol, got, q, want)
			}
			if want > 0 {
				inside++
			}
		}
		if inside == 0 || inside == 300 {
			t.Fatalf("%s: %d of 300 points inside", c.symbol, inside)
		}
	}
}

type flagSpan struct {
	t0, t1 Real
	c0, c1 int // cell whose hyperplane ends the span, -1 for a side shared with another simplex
}

// firstCrossing returns the first point at t > eps where the ray A + t·B
// enters or leaves the union of the simplices: each simplex is cut by its
// five sides into a span of the ray, and the spans are joined.
func (o *oracleFlags) firstCrossing(A, B Vector4, eps Real) (t Real, cell int, inside, ok bool) {
	var spans []flagSpan
	for ci, n := range o.cells {
		for _, w := range o.walls[ci] {
			s := flagSpan{math.Inf(-1), math.Inf(1), -1, -1}
			clip := func(a, b Real, cell int) bool { // a + t·b <= 0
				if b == 0 {
					return a <= 0
				}
				if x := -a / b; b > 0 {
					if x < s.t1 {
						s.t1, s.c1 = x, cell
					}
				} else if x > s.t0 {
					s.t0, s.c0 = x, cell
				}
				return s.t0 < s.t1
			}
			if clip(n.Dot(A)-o.height, n.Dot(B), ci) && clip(w[0].Dot(A), w[0].Dot(B), -1) && clip(w[1].Dot(A), w[1].Dot(B), -1) &&
				clip(w[2].Dot(A), w[2].Dot(B), -1) && clip(w[3].Dot(A), w[3].Dot(B), -1) && s.t1 > eps {
				spans = append(spans, s)
			}
		}
	}
	sort.Slice(spans, func(i, j int) bool { return spans[i].t0 < spans[j].t0 })
	const join = 1e-9
	end, cell := math.Inf(-1), -1
	for _, s := range spans {
		if s.t0 > eps && end < eps {
			// nothing holds the start: the first span is the way in
			return s.t0, s.c0, false, s.c0 >= 0
		}
		if s.t0 > end+join && end >= eps {
			break
		}
		if s.t1 > end {
			end, cell = s.t1, s.c1
		}
	}
	if end < eps {
		return 0, -1, false, false
	}
	return end, cell, true, cell >= 0
}

// Rays against the list of simplices, for all sixteen regular polychora at
// several sizes: from outside, from inside, and the way the tracer goes on
// from a hit.
func TestStarPolychoraAgainstFlagSimplices(t *testing.T) {
	rays := 120
	if testing.Short() {
		rays = 40
	}
	rng := rand.New(rand.NewSource(13))
	for _, c := range regularPolychoraTable {
		flags := newOracleFlags(t, c.symbol)
		for _, size := range []Real{1, 0.03, 30} {
			rnd := func() Real { return (rng.Float64() - 0.5) * 2 * math.Pi }
			rot := Rot4{XY: rnd(), XZ: rnd(), XW: rnd(), YZ: rnd(), YW: rnd(), ZW: rnd()}
			ctr := Point4{rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5, rng.Float64() - 0.5}
			rs := func() Real { return (0.5 + rng.Float64()) * size }
			sc := Vector4{rs(), rs(), rs(), rs()}
			o, err := NewStarPolychoron(c.symbol, ctr, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
			if err != nil {
				t.Fatal(err)
			}
			world := func(q Vector4) Point4 {
				return ctr.Add(o.R.MulVec(Vector4{q.X * sc.X, q.Y * sc.Y, q.Z * sc.Z, q.W * sc.W}))
			}
			hits, inner := 0, 0
			check := func(O Point4, D Vector4, kind string) (objectHit, bool) {
				l := o.RT.MulVec(Vector4{O.X - ctr.X, O.Y - ctr.Y, O.Z - ctr.Z, O.W - ctr.W})
				d := o.RT.MulVec(D)
				A := Vector4{l.X / sc.X, l.Y / sc.Y, l.Z / sc.Z, l.W / sc.W}
				B := Vector4{d.X / sc.X, d.Y / sc.Y, d.Z / sc.Z, d.W / sc.W}
				tr, cell, inside, okr := flags.firstCrossing(A, B, 1e-12)
				h, ok := intersectRayStarPolychoron(O, D, o)
				if ok != okr {
					t.Fatalf("%s size %g, %s ray O=%v D=%v: hit %v (t=%.12g), the simplices give %v (t=%.12g)", c.symbol, size, kind, O, D, ok, h.t, okr, tr)
				}
				if !ok {
					return h, false
				}
				hits++
				if inside {
					inner++
				}
				n := flags.cells[cell]
				want := o.R.MulVec(Vector4{n.X / sc.X, n.Y / sc.Y, n.Z / sc.Z, n.W / sc.W}).Norm()
				if math.Abs(h.t-tr) > 1e-9*(1+tr) || h.inv != inside || h.Nw.Sub(want).Len() > 1e-9 {
					t.Fatalf("%s size %g, %s ray O=%v D=%v: t=%.12g inv=%v N=%v, the simplices give t=%.12g inv=%v N=%v",
						c.symbol, size, kind, O, D, h.t, h.inv, h.Nw, tr, inside, want)
				}
				P := O.Add(D.Mul(h.t))
				if P.X < o.AABBMin.X-1e-9*size || P.X > o.AABBMax.X+1e-9*size || P.Y < o.AABBMin.Y-1e-9*size || P.Y > o.AABBMax.Y+1e-9*size ||
					P.Z < o.AABBMin.Z-1e-9*size || P.Z > o.AABBMax.Z+1e-9*size || P.W < o.AABBMin.W-1e-9*size || P.W > o.AABBMax.W+1e-9*size {
					t.Fatalf("%s size %g: hit point %v outside the AABB", c.symbol, size, P)
				}
				return h, true
			}
			for i := 0; i < rays; i++ {
				// from outside towards a point of the unit ball, then on from the hit
				O := world(unitS3(rng).Mul(1.2 + 2*rng.Float64()))
				T := world(unitS3(rng).Mul(math.Pow(rng.Float64(), 0.3)))
				D := Vector4{T.X - O.X, T.Y - O.Y, T.Z - O.Z, T.W - O.W}.Norm()
				for bounce := 0; bounce < 6; bounce++ {
					h, ok := check(O, D, "outer")
					if !ok {
						break
					}
					P := O.Add(D.Mul(h.t))
					N := h.Nw
					if h.inv {
						N = N.Mul(-1)
					}
					switch rng.Intn(3) {
					case 0:
						D = reflect4(D, N).Norm()
					case 1:
						eta := Real(1 / 1.5)
						if h.inv {
							eta = 1.5
						}
						if T, ok := refract4(D, N, eta); ok {
							D = T.Norm()
						} else {
							D = reflect4(D, N).Norm()
						}
					default:
						D = sampleDiffuseDir(N, rng)
					}
					O = P.Add(D.Mul(bumpShift * size))
				}
				// from a point of the unit ball
				check(world(unitS3(rng).Mul(math.Pow(rng.Float64(), 0.3))), unitS3(rng), "inner")
			}
			if hits < rays/2 || inner < rays/8 {
				t.Fatalf("%s size %g: %d hits, %d of them from inside", c.symbol, size, hits, inner)
			}
		}
	}
}

// The convex regular polychora built as stars are the polytopes of the same
// mirrors.
func TestConvexRegularStarsEqualPolytopes(t *testing.T) {
	rng := rand.New(rand.NewSource(14))
	rot := Rot4{XY: 0.3, XZ: 0.9, XW: -0.7, YZ: 1.1, YW: 2.2, ZW: 0.4}
	sc := Vector4{0.7, 1.1, 0.9, 1.3}
	for _, pair := range [][2]string{{"{3,3,3}", "x3o3o3o"}, {"{4,3,3}", "x4o3o3o"}, {"{3,3,4}", "x3o3o4o"}, {"{3,4,3}", "x3o4o3o"}, {"{5,3,3}", "x5o3o3o"}, {"{3,3,5}", "x3o3o5o"}} {
		s, err := NewStarPolychoron(pair[0], Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		p, err := NewPolychoron(pair[1], false, Point4{}, sc, rot, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			t.Fatal(err)
		}
		sameHits(t, pair[0]+" and "+pair[1], rng, 1.3,
			func(O Point4, D Vector4) (objectHit, bool) { return intersectRayStarPolychoron(O, D, s) },
			func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &p.cellPoly) })
		for _, d := range []Real{
			s.AABBMin.X - p.AABBMin.X, s.AABBMin.Y - p.AABBMin.Y, s.AABBMin.Z - p.AABBMin.Z, s.AABBMin.W - p.AABBMin.W,
			s.AABBMax.X - p.AABBMax.X, s.AABBMax.Y - p.AABBMax.Y, s.AABBMax.Z - p.AABBMax.Z, s.AABBMax.W - p.AABBMax.W,
		} {
			if math.Abs(d) > 1e-12 {
				t.Errorf("%s: AABB %v..%v, the polytope has %v..%v", pair[0], s.AABBMin, s.AABBMax, p.AABBMin, p.AABBMax)
				break
			}
		}
	}
}

// The shapes are built on first use and shared between objects, and rays are
// traced from many goroutines: all of them must see what one sees alone.
func TestStarPolychoraFromManyGoroutines(t *testing.T) {
	trace := func() (string, error) {
		s, err := NewStarPolychoron("grand 600-cell", Point4{}, Vector4{1, 1, 1, 1}, Rot4{XY: 0.3, ZW: 1.1}, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			return "", err
		}
		p, err := NewPolychoron("snub 24-cell", false, Point4{}, Vector4{1, 1, 1, 1}, Rot4{XY: 0.3, ZW: 1.1}, testWhite, testNone, testHalf, testHalf, testIOR)
		if err != nil {
			return "", err
		}
		rng := rand.New(rand.NewSource(15))
		hits, sum := 0, Real(0)
		for i := 0; i < 400; i++ {
			O := Point4{}.Add(unitS3(rng).Mul(0.2 + 2*rng.Float64()))
			D := unitS3(rng)
			if h, ok := intersectRayStarPolychoron(O, D, s); ok {
				hits++
				sum += h.t
			}
			if h, ok := intersectRayCellPoly(O, D, &p.cellPoly); ok {
				hits++
				sum += h.t
			}
		}
		return fmt.Sprintf("%d hits, sum of t %.17g", hits, sum), nil
	}
	want, err := trace()
	if err != nil {
		t.Fatal(err)
	}
	// build the shapes again, this time from all goroutines at once
	starShapes.Delete("{3,3,5/2}")
	polytopeShapes.Delete("polychoron snub 24-cell")
	var wg sync.WaitGroup
	got := make([]string, 8)
	errs := make([]error, len(got))
	for g := range got {
		wg.Add(1)
		go func(g int) {
			defer wg.Done()
			got[g], errs[g] = trace()
		}(g)
	}
	wg.Wait()
	for g := range got {
		if errs[g] != nil || got[g] != want {
			t.Errorf("goroutine %d: %s (%v), alone: %s", g, got[g], errs[g], want)
		}
	}
}
