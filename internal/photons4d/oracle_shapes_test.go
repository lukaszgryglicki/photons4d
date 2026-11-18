package photons4d

import (
	"math"
	"testing"
)

// Descriptions of the newer solids that share no code with the way the
// objects are built, used by the tests of the geometry.

// kaleidoscope folds a point into the chamber between the mirrors of a
// Coxeter diagram. There the polytope of the diagram is cut by one plane per
// mirror only: the one square to the corner of the chamber off that mirror.
type kaleidoscope struct {
	mirrors []Vector4
	seed    Vector4   // the vertex inside the chamber, on the unit sphere
	corners []Vector4 // unit vectors along the corners of the chamber
	heights []Real    // their products with the seed
}

func newKaleidoscope(t *testing.T, diagram string) *kaleidoscope {
	t.Helper()
	marks, rings, ok := parseCoxeterDiagram(diagram)
	if !ok {
		t.Fatalf("bad diagram %q", diagram)
	}
	mirrors, dual, err := coxeterMirrors(marks)
	if err != nil {
		t.Fatalf("%s: %v", diagram, err)
	}
	k := &kaleidoscope{mirrors: mirrors}
	for i, on := range rings {
		if on {
			k.seed = k.seed.Add(dual[i])
		}
	}
	k.seed = k.seed.Norm()
	for _, d := range dual {
		k.corners = append(k.corners, d.Norm())
		k.heights = append(k.heights, d.Norm().Dot(k.seed))
	}
	return k
}

func (k *kaleidoscope) fold(q Vector4) Vector4 {
	for moved, n := true, 0; moved && n < 10000; n++ {
		moved = false
		for _, m := range k.mirrors {
			if d := q.Dot(m); d < 0 {
				q = q.Sub(m.Mul(2 * d))
				moved = true
			}
		}
	}
	return q
}

// gauge is 1 on the surface of the polytope with the unit circumradius.
func (k *kaleidoscope) gauge(q Vector4) Real {
	q = k.fold(q)
	g := Real(0)
	for i, c := range k.corners {
		g = math.Max(g, c.Dot(q)/k.heights[i])
	}
	return g
}

// dualGauge is 1 on the surface of the dual with the unit circumradius: the
// planes of the dual are square to the vertices, the nearest vertex to a
// point of the chamber is the seed, and the farthest vertex of the dual
// belongs to the nearest cell.
func (k *kaleidoscope) dualGauge(q Vector4) Real {
	least := math.Inf(1)
	for _, h := range k.heights {
		least = math.Min(least, h)
	}
	return k.seed.Dot(k.fold(q)) / least
}

func oracleQuatMul(a, b Vector4) Vector4 {
	return Vector4{
		a.X*b.X - a.Y*b.Y - a.Z*b.Z - a.W*b.W,
		a.X*b.Y + a.Y*b.X + a.Z*b.W - a.W*b.Z,
		a.X*b.Z - a.Y*b.W + a.Z*b.X + a.W*b.Y,
		a.X*b.W + a.Y*b.Z - a.Z*b.Y + a.W*b.X,
	}
}

// oracleHurwitzUnits returns the 24 vertices of a 24-cell inscribed in the
// 600-cell of the icosians.
func oracleHurwitzUnits() []Vector4 {
	var out []Vector4
	for i := 0; i < 4; i++ {
		for _, s := range []Real{1, -1} {
			out = append(out, vecWith(Vector4{}, i, s))
		}
	}
	for s := 0; s < 16; s++ {
		v := Vector4{0.5, 0.5, 0.5, 0.5}
		for i := 0; i < 4; i++ {
			if s>>i&1 == 1 {
				v = vecWith(v, i, -0.5)
			}
		}
		out = append(out, v)
	}
	return out
}

// oracleDecagonRings returns two rings of ten vertices of the 600-cell in
// planes square to each other: the powers of an icosian of order 10, and
// their products with a unit quaternion square to all of them.
func oracleDecagonRings() []Vector4 {
	phi := (1 + math.Sqrt(5)) / 2
	q := Vector4{phi / 2, 0.5, 0.5 / phi, 0}
	p := Vector4{1, 0, 0, 0}
	var out []Vector4
	for i := 0; i < 10; i++ {
		out = append(out, p, oracleQuatMul(p, Vector4{0, 0, 0, 1}))
		p = oracleQuatMul(p, q)
	}
	return out
}

// oracleHull3 returns the face planes n·x <= 1 of the hull of points of XYZ
// by trying every plane through three of them.
func oracleHull3(verts []Vector4) []Vector4 {
	var out []Vector4
	for i := range verts {
		for j := i + 1; j < len(verts); j++ {
		next:
			for k := j + 1; k < len(verts); k++ {
				n := cross3_4D(verts[j].Sub(verts[i]), verts[k].Sub(verts[i]), Vector4{0, 0, 0, 1})
				h := n.Dot(verts[i])
				if n.Len() < 1e-9 || math.Abs(h) < 1e-12 {
					continue
				}
				n = n.Mul(1 / h)
				for _, w := range verts {
					if n.Dot(w) > 1+1e-9 {
						continue next
					}
				}
				for _, m := range out {
					if m.Sub(n).Len() < 1e-9 {
						continue next
					}
				}
				out = append(out, n)
			}
		}
	}
	return out
}

// oracleFlags holds a regular polychoron {p,q,r} as the list of the
// 4-simplices spanned by its centre, the centre of a cell, the centre of a
// face of the cell and an edge of the face. They are found one by one: the
// simplex of the base flag and its images under the mirrors.
type oracleFlags struct {
	height Real
	cells  []Vector4      // unit normals of the cell hyperplanes
	walls  [][][4]Vector4 // per cell, per simplex: unit normals of its four sides through the centre
	corner [][][4]Vector4 // per cell, per simplex: cell centre, face centre and the ends of the edge
}

func newOracleFlags(t *testing.T, symbol string) *oracleFlags {
	t.Helper()
	marks, ok := parseSchlafli(symbol)
	if !ok {
		t.Fatalf("bad symbol %q", symbol)
	}
	mirrors, dual, err := coxeterMirrors(marks)
	if err != nil {
		t.Fatalf("%s: %v", symbol, err)
	}
	// The base flag from the definitions: the vertex lies on the mirrors 1,
	// 2, 3; the edge joins it to its image in the mirror 0; the centre of the
	// face is the nearest point to the vertex that lies on the mirrors 0 and
	// 1; the centre of the cell the nearest one on the mirrors 0, 1 and 2.
	v := dual[0].Norm()
	w := mirrorImage(v, mirrors[0])
	cdir := dual[3].Norm()
	c := cdir.Mul(cdir.Dot(v))
	// the points on the mirrors 0 and 1 are spanned by dual[2] and dual[3]
	g22, g23, g33 := dual[2].Dot(dual[2]), dual[2].Dot(dual[3]), dual[3].Dot(dual[3])
	r2, r3 := dual[2].Dot(v), dual[3].Dot(v)
	det := g22*g33 - g23*g23
	f := dual[2].Mul((r2*g33 - r3*g23) / det).Add(dual[3].Mul((r3*g22 - r2*g23) / det))

	type simplex [4]Vector4
	key := func(s simplex) [12]int64 {
		m := s[2].Add(s[3])
		var k [12]int64
		for i := 0; i < 4; i++ {
			k[i] = int64(math.Round(vecAt(s[0], i) * 1e7))
			k[4+i] = int64(math.Round(vecAt(s[1], i) * 1e7))
			k[8+i] = int64(math.Round(vecAt(m, i) * 1e7))
		}
		return k
	}
	base := simplex{c, f, v, w}
	seen := map[[12]int64]bool{key(base): true}
	all := []simplex{base}
	for i := 0; i < len(all); i++ {
		for _, m := range mirrors {
			var s simplex
			for j := range s {
				s[j] = mirrorImage(all[i][j], m)
			}
			if k := key(s); !seen[k] {
				seen[k] = true
				all = append(all, s)
			}
		}
		if len(all) > 20000 {
			t.Fatalf("%s: too many simplices", symbol)
		}
	}

	o := &oracleFlags{height: c.Len()}
	index := map[[4]int64]int{}
	for _, s := range all {
		n := s[0].Norm()
		k := [4]int64{int64(math.Round(n.X * 1e7)), int64(math.Round(n.Y * 1e7)), int64(math.Round(n.Z * 1e7)), int64(math.Round(n.W * 1e7))}
		ci, ok := index[k]
		if !ok {
			ci = len(o.cells)
			index[k] = ci
			o.cells = append(o.cells, n)
			o.walls = append(o.walls, nil)
			o.corner = append(o.corner, nil)
		}
		var walls [4]Vector4
		for j := 0; j < 4; j++ {
			var rest []Vector4
			for i := 0; i < 4; i++ {
				if i != j {
					rest = append(rest, s[i])
				}
			}
			n := cross3_4D(rest[0], rest[1], rest[2]).Norm()
			if n.Dot(s[j]) > 0 {
				n = n.Mul(-1)
			}
			walls[j] = n
		}
		o.walls[ci] = append(o.walls[ci], walls)
		o.corner[ci] = append(o.corner[ci], [4]Vector4(s))
	}
	return o
}

// distance is negative inside the union of the simplices and positive
// outside; on its surface it is the distance to the hyperplane of the cell.
func (o *oracleFlags) distance(q Vector4) Real {
	best := math.Inf(1)
	for ci, n := range o.cells {
		top := n.Dot(q) - o.height
		if top >= best {
			continue
		}
	simplices:
		for _, w := range o.walls[ci] {
			d := top
			for j := 0; j < 4; j++ {
				if v := w[j].Dot(q); v > d {
					if v >= best {
						continue simplices
					}
					d = v
				}
			}
			best = d
		}
	}
	return best
}

// count returns the number of simplices that hold the point.
func (o *oracleFlags) count(q Vector4) int {
	n := 0
	for ci, c := range o.cells {
		if c.Dot(q) > o.height {
			continue
		}
		for _, w := range o.walls[ci] {
			if w[0].Dot(q) <= 0 && w[1].Dot(q) <= 0 && w[2].Dot(q) <= 0 && w[3].Dot(q) <= 0 {
				n++
			}
		}
	}
	return n
}

// oracleSweptDistance is the distance from a point to the surface that a
// flat curve sweeps while its plane turns once around XY and by twist
// half-turns in ZW: a circle of the radius m around a point at the distance
// r from the centre, or the segment of the half-length m through that point
// that is square to the radius. Every position of the curve is measured with
// vectors of 4D; the best turns of a grid are narrowed by a golden section.
func oracleSweptDistance(q Vector4, r, m Real, twist int, band bool) Real {
	at := func(u Real) Real {
		e := Vector4{math.Cos(u), math.Sin(u), 0, 0}
		n := Vector4{0, 0, math.Cos(u * Real(twist) / 2), math.Sin(u * Real(twist) / 2)}
		d := q.Sub(e.Mul(r))
		x, z := d.Dot(e), d.Dot(n)
		off := d.Sub(e.Mul(x)).Sub(n.Mul(z)).Len()
		if band {
			return math.Hypot(off, math.Hypot(x, math.Max(0, math.Abs(z)-m)))
		}
		return math.Hypot(off, math.Hypot(x, z)-m)
	}
	var grid [48 * 9]Real
	steps := 48 * (1 + (twist+1)/2)
	h := 2 * math.Pi / Real(steps)
	for i := 0; i < steps; i++ {
		grid[i] = at(Real(i) * h)
	}
	const golden = 0.6180339887498949
	best := math.Inf(1)
	for i := 0; i < steps; i++ {
		if grid[i] > grid[(i+steps-1)%steps] || grid[i] > grid[(i+1)%steps] {
			continue
		}
		lo, hi := Real(i-1)*h, Real(i+1)*h
		a, b := hi-golden*(hi-lo), lo+golden*(hi-lo)
		fa, fb := at(a), at(b)
		for k := 0; k < 34; k++ {
			if fa < fb {
				hi, b, fb = b, a, fa
				a = hi - golden*(hi-lo)
				fa = at(a)
			} else {
				lo, a, fa = a, b, fb
				b = lo + golden*(hi-lo)
				fb = at(b)
			}
		}
		best = math.Min(best, math.Min(math.Min(fa, fb), grid[i]))
	}
	return best
}
