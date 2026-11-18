package photons4d

import (
	"fmt"
	"math"
	"strconv"
	"strings"
	"sync"
)

// StarPolychoron is a regular polychoron {p,q,r} as a solid. It serves the
// ten regular star polychora of Schläfli and Hess, whose cells pass through
// each other; the six convex ones are accepted as well.
//
// The solid is the set of points that the polychoron winds around: the union
// of the 4-simplices spanned by the centre of the polychoron, the centre of a
// cell, the centre of a face of that cell and an edge of that face. The
// number of those simplices that hold a point is the density at that point;
// along a ray it only changes on the hyperplane of a cell, by the density of
// the cell (a star polyhedron of its own) at the crossing. Before Scale is
// applied the vertices lie on the unit sphere.
type StarPolychoron struct {
	Center Point4
	Scale  Vector4
	R      Mat4
	RT     Mat4
	Name   string

	Color   RGB
	Diffuse RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	AABBMin Point4
	AABBMax Point4

	shape *starShape

	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

// Schläfli symbols of the regular star polychora.
var starPolychoronKinds = map[string]string{
	"icosahedral 120-cell":           "{3,5,5/2}",
	"faceted 600-cell":               "{3,5,5/2}",
	"small stellated 120-cell":       "{5/2,5,3}",
	"great 120-cell":                 "{5,5/2,5}",
	"grand 120-cell":                 "{5,3,5/2}",
	"great stellated 120-cell":       "{5/2,3,5}",
	"grand stellated 120-cell":       "{5/2,5,5/2}",
	"great grand 120-cell":           "{5,5/2,3}",
	"great icosahedral 120-cell":     "{3,5/2,5}",
	"great faceted 600-cell":         "{3,5/2,5}",
	"grand 600-cell":                 "{3,3,5/2}",
	"great grand stellated 120-cell": "{5/2,3,3}",
}

// One simplex of a cell: the cell centre, the centre of the face and an edge
// of that face. A point of the cell hyperplane lies in it when it is below
// the face plane and on the wanted side of three mirror planes of the cell.
type starFlag struct {
	face       int
	mask, want uint32 // which mirror planes bound the simplex, and its side of each
}

// A node of the tree that finds the region of a point of a cell between the
// mirror planes: an inner node tests the plane wall and goes on at pos (the
// point is not below the plane) or at neg; a leaf (wall < 0) holds the faces
// cover[pos:neg], whose simplices cover the whole region.
type starNode struct {
	wall     int32
	pos, neg int32
}

type starCover struct {
	face, times int32
}

type starShape struct {
	symbol string
	verts  []Vector4 // on the unit sphere
	height Real      // distance of the cell hyperplanes from the centre

	// per cell: the unit normal of its hyperplane and the directions which
	// give the coordinates of a point of that hyperplane in the frame of the
	// base cell (origin at the cell centre)
	normals []Vector4
	frames  [][3]Vector4

	// the base cell in its frame
	walls    [][3]Real // unit normals of its mirror planes
	faces    [][3]Real // unit normals of its face planes
	faceDist Real      // distance of the face planes from the cell centre
	flags    []starFlag
	nodes    []starNode
	cover    []starCover

	faceCount, edgeCount int
}

func dot3(a, b [3]Real) Real { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }

func cross3(a, b [3]Real) [3]Real {
	return [3]Real{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

// cellDensity returns how many times the base cell covers the point y of its
// hyperplane: the number of its simplices that hold y. The mirror planes cut
// the hyperplane into regions, and a simplex that reaches a region holds all
// of its points below the face of the simplex.
func (s *starShape) cellDensity(y [3]Real) int {
	nd := &s.nodes[0]
	for nd.wall >= 0 {
		if dot3(s.walls[nd.wall], y) >= 0 {
			nd = &s.nodes[nd.pos]
		} else {
			nd = &s.nodes[nd.neg]
		}
	}
	n := 0
	for _, c := range s.cover[nd.pos:nd.neg] {
		if dot3(s.faces[c.face], y) <= s.faceDist {
			n += int(c.times)
		}
	}
	return n
}

// buildRegions makes the tree of cellDensity from the simplices of the base
// cell. points holds one point inside each chamber of the mirrors of the cell.
func (s *starShape) buildRegions(points [][3]Real) {
	seen := map[uint32]bool{}
	var words []uint32
	for _, y := range points {
		var w uint32
		for j, n := range s.walls {
			if dot3(n, y) >= 0 {
				w |= 1 << uint(j)
			}
		}
		if !seen[w] {
			seen[w] = true
			words = append(words, w)
		}
	}
	var build func(words []uint32) int32
	build = func(words []uint32) int32 {
		at := int32(len(s.nodes))
		s.nodes = append(s.nodes, starNode{wall: -1})
		if len(words) == 1 {
			times := make([]int32, len(s.faces))
			for _, f := range s.flags {
				if words[0]&f.mask == f.want {
					times[f.face]++
				}
			}
			s.nodes[at].pos = int32(len(s.cover))
			for face, n := range times {
				if n > 0 {
					s.cover = append(s.cover, starCover{face: int32(face), times: n})
				}
			}
			s.nodes[at].neg = int32(len(s.cover))
			return at
		}
		// the plane that splits the regions most evenly
		best, even := -1, len(words)
		for j := range s.walls {
			on := 0
			for _, w := range words {
				on += int(w >> uint(j) & 1)
			}
			off := len(words) - on
			if larger := imax(on, off); on > 0 && off > 0 && larger < even {
				best, even = j, larger
			}
		}
		var pos, neg []uint32
		for _, w := range words {
			if w>>uint(best)&1 == 1 {
				pos = append(pos, w)
			} else {
				neg = append(neg, w)
			}
		}
		p := build(pos)
		n := build(neg)
		s.nodes[at] = starNode{wall: int32(best), pos: p, neg: n}
		return at
	}
	build(words)
}

// density returns how many times the polychoron covers the point x.
func (s *starShape) density(x Vector4) int {
	n := 0
	for k, c := range s.normals {
		d := c.Dot(x)
		if d <= 0 || d > s.height {
			continue
		}
		p := x.Mul(s.height / d)
		if p.Dot(p) > 1 {
			continue
		}
		f := &s.frames[k]
		n += s.cellDensity([3]Real{f[0].Dot(p), f[1].Dot(p), f[2].Dot(p)})
	}
	return n
}

func parseSchlafli(s string) ([]Real, bool) {
	if len(s) < 2 || s[0] != '{' || s[len(s)-1] != '}' {
		return nil, false
	}
	var marks []Real
	for _, part := range strings.Split(s[1:len(s)-1], ",") {
		num, den := strings.TrimSpace(part), "1"
		if k := strings.IndexByte(num, '/'); k >= 0 {
			num, den = num[:k], num[k+1:]
		}
		a, errA := strconv.Atoi(num)
		b, errB := strconv.Atoi(den)
		if errA != nil || errB != nil || b < 1 || a <= 2*b {
			return nil, false
		}
		marks = append(marks, Real(a)/Real(b))
	}
	return marks, true
}

var starShapes sync.Map // symbol -> *starShape or error

// regularStarShape builds the cells of the regular polychoron with the given
// name or Schläfli symbol such as "{5/2,5,3}".
func regularStarShape(kind string) (*starShape, error) {
	key := shapeKey(kind)
	if symbol, ok := starPolychoronKinds[key]; ok {
		key = symbol
	}
	key = strings.ReplaceAll(key, " ", "")
	if v, ok := starShapes.Load(key); ok {
		if s, ok := v.(*starShape); ok {
			return s, nil
		}
		return nil, v.(error)
	}
	marks, ok := parseSchlafli(key)
	if !ok || len(marks) != 3 {
		err := fmt.Errorf("unknown regular polychoron %q", kind)
		starShapes.Store(key, err)
		return nil, err
	}
	s, err := buildStarShape(marks)
	if err != nil {
		err = fmt.Errorf("regular polychoron %q: %w", kind, err)
		starShapes.Store(key, err)
		return nil, err
	}
	s.symbol = key
	starShapes.Store(key, s)
	return s, nil
}

func buildStarShape(marks []Real) (*starShape, error) {
	mirrors, dual, err := coxeterMirrors(marks)
	if err != nil {
		return nil, err
	}
	// The base flag: a vertex, the edge to its image in the first mirror, the
	// face kept by the last two mirrors and the cell kept by the last one.
	v0 := dual[0].Norm()
	v1 := mirrorImage(v0, mirrors[0])
	n0 := dual[3].Norm()
	s := &starShape{height: n0.Dot(v0)}
	if !(s.height > 1e-9 && s.height < 1-1e-9) {
		return nil, fmt.Errorf("the cells do not surround the centre")
	}
	c0 := n0.Mul(s.height)
	e2 := dual[2].Sub(n0.Mul(dual[2].Dot(n0))).Norm()
	f0 := n0.Mul(v0.Dot(n0)).Add(e2.Mul(v0.Dot(e2)))

	if s.verts, err = reflectionOrbit(mirrors, v0); err != nil {
		return nil, err
	}
	faces, err := reflectionOrbit(mirrors, f0)
	if err != nil {
		return nil, err
	}
	edges, err := reflectionOrbit(mirrors, v0.Add(v1).Mul(0.5))
	if err != nil {
		return nil, err
	}
	s.faceCount, s.edgeCount = len(faces), len(edges)

	// frame of the base cell
	var b [3]Vector4
	for i, d := range []Vector4{f0.Sub(c0), v0.Sub(c0), v1.Sub(c0)} {
		for j := 0; j < i; j++ {
			d = d.Sub(b[j].Mul(d.Dot(b[j])))
		}
		if d.Len() < 1e-9 {
			return nil, fmt.Errorf("the cells are flat")
		}
		b[i] = d.Norm()
	}
	local := func(p Vector4) [3]Real {
		d := p.Sub(c0)
		return [3]Real{d.Dot(b[0]), d.Dot(b[1]), d.Dot(b[2])}
	}

	// all cells, each with the frame carried to it by the mirrors
	set := newPointSet()
	set.add(n0)
	s.normals = []Vector4{n0}
	s.frames = [][3]Vector4{b}
	for i := 0; i < len(s.normals); i++ {
		for _, m := range mirrors {
			if _, fresh := set.add(mirrorImage(s.normals[i], m)); fresh {
				f := s.frames[i]
				s.normals = append(s.normals, mirrorImage(s.normals[i], m))
				s.frames = append(s.frames, [3]Vector4{mirrorImage(f[0], m), mirrorImage(f[1], m), mirrorImage(f[2], m)})
				if len(s.normals) > maxOrbit {
					return nil, fmt.Errorf("the reflection group is not finite")
				}
			}
		}
	}

	// the simplices of the base cell: the images of (face, edge) under its mirrors
	type flag struct{ f, a, b Vector4 }
	fset, vset := newPointSet(), newPointSet()
	key := func(g flag) [3]int {
		fi, _ := fset.add(g.f)
		ai, _ := vset.add(g.a)
		bi, _ := vset.add(g.b)
		if ai > bi {
			ai, bi = bi, ai
		}
		return [3]int{fi, ai, bi}
	}
	seen := map[[3]int]bool{key(flag{f0, v0, v1}): true}
	queue := []flag{{f0, v0, v1}}
	for i := 0; i < len(queue); i++ {
		for _, m := range mirrors[:3] {
			g := flag{mirrorImage(queue[i].f, m), mirrorImage(queue[i].a, m), mirrorImage(queue[i].b, m)}
			if k := key(g); !seen[k] {
				seen[k] = true
				queue = append(queue, g)
				if len(queue) > maxOrbit {
					return nil, fmt.Errorf("the reflection group is not finite")
				}
			}
		}
	}
	s.faceDist = f0.Sub(c0).Len()
	for _, f := range fset.pts {
		l := local(f)
		s.faces = append(s.faces, [3]Real{l[0] / s.faceDist, l[1] / s.faceDist, l[2] / s.faceDist})
	}
	wall := func(p, q, inner [3]Real) (int, uint32, error) {
		n := cross3(p, q)
		l := math.Sqrt(dot3(n, n))
		if l < 1e-9 {
			return 0, 0, fmt.Errorf("a face is flat")
		}
		n = [3]Real{n[0] / l, n[1] / l, n[2] / l}
		for i := 0; i < 3; i++ {
			if math.Abs(n[i]) > 1e-9 {
				if n[i] < 0 {
					n = [3]Real{-n[0], -n[1], -n[2]}
				}
				break
			}
		}
		side := uint32(0)
		if dot3(n, inner) > 0 {
			side = 1
		}
		for j, w := range s.walls {
			if math.Abs(w[0]-n[0])+math.Abs(w[1]-n[1])+math.Abs(w[2]-n[2]) < 1e-8 {
				return j, side, nil
			}
		}
		if len(s.walls) == 32 {
			return 0, 0, fmt.Errorf("a cell with more than 32 mirrors is not supported")
		}
		s.walls = append(s.walls, n)
		return len(s.walls) - 1, side, nil
	}
	for _, g := range queue {
		f, a, c := local(g.f), local(g.a), local(g.b)
		fl := starFlag{}
		fl.face, _ = fset.add(g.f)
		for _, w := range [][3][3]Real{{f, a, c}, {f, c, a}, {a, c, f}} {
			j, side, err := wall(w[0], w[1], w[2])
			if err != nil {
				return nil, err
			}
			fl.mask |= 1 << uint(j)
			fl.want |= side << uint(j)
		}
		s.flags = append(s.flags, fl)
	}

	// A point of the base simplex off its planes of symmetry, and its images
	// under the mirrors of the cell: one point in each of their chambers.
	chambers := newPointSet()
	chambers.add(c0.Add(f0.Sub(c0).Mul(0.45)).Add(v0.Sub(c0).Mul(0.27)).Add(v1.Sub(c0).Mul(0.18)))
	for i := 0; i < len(chambers.pts) && len(chambers.pts) <= maxOrbit; i++ {
		for _, m := range mirrors[:3] {
			chambers.add(mirrorImage(chambers.pts[i], m))
		}
	}
	points := make([][3]Real, len(chambers.pts))
	for i, p := range chambers.pts {
		points[i] = local(p)
	}
	s.buildRegions(points)
	return s, nil
}

// NewStarPolychoron builds a regular polychoron as a solid. kind is one of
// the names of the ten regular star polychora ("small stellated 120-cell",
// "great grand stellated 120-cell", ...) or a Schläfli symbol ("{5/2,5,3}").
func NewStarPolychoron(kind string, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*StarPolychoron, error) {
	if err := checkSolid("star polychoron", scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	shape, err := regularStarShape(kind)
	if err != nil {
		return nil, err
	}
	R := rotFromAngles(angles)
	s := &StarPolychoron{
		Center:  center,
		Scale:   scale,
		R:       R,
		RT:      R.Transpose(),
		Name:    shape.symbol,
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
		shape:   shape,
	}
	// the solid lies in the hull of its vertices
	s.AABBMin = Point4{math.Inf(1), math.Inf(1), math.Inf(1), math.Inf(1)}
	s.AABBMax = Point4{math.Inf(-1), math.Inf(-1), math.Inf(-1), math.Inf(-1)}
	for _, v := range shape.verts {
		p := R.MulVec(Vector4{v.X * scale.X, v.Y * scale.Y, v.Z * scale.Z, v.W * scale.W})
		w := Point4{center.X + p.X, center.Y + p.Y, center.Z + p.Z, center.W + p.W}
		s.AABBMin = Point4{math.Min(s.AABBMin.X, w.X), math.Min(s.AABBMin.Y, w.Y), math.Min(s.AABBMin.Z, w.Z), math.Min(s.AABBMin.W, w.W)}
		s.AABBMax = Point4{math.Max(s.AABBMax.X, w.X), math.Max(s.AABBMax.Y, w.Y), math.Max(s.AABBMax.Z, w.Z), math.Max(s.AABBMax.W, w.W)}
	}

	s.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	s.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	s.colorArr = [3]Real{color.R, color.G, color.B}
	s.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	s.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - s.refl[i] - s.refr[i] - s.diff[i]
		if p < 0 {
			p = 0
		}
		s.pAbs[i] = p
		s.iorInv[i] = 1 / s.iorArr[i]
		n := s.iorArr[i]
		r0 := (n - 1) / (n + 1)
		s.f0[i] = r0 * r0
	}
	DebugLog("Created regular polychoron %s: %d cells, %d vertices", shape.symbol, len(shape.normals), len(shape.verts))
	return s, nil
}

type starEvent struct {
	t    Real
	jump int
	cell int
}

// Room for the events of the shapes with more cells than fit on the stack.
var starEventPool = sync.Pool{New: func() any { return new([]starEvent) }}

func starEventDown(ev []starEvent, i int) {
	for {
		l := 2*i + 1
		if l >= len(ev) {
			return
		}
		if r := l + 1; r < len(ev) && ev[r].t < ev[l].t {
			l = r
		}
		if ev[i].t <= ev[l].t {
			return
		}
		ev[i], ev[l] = ev[l], ev[i]
		i = l
	}
}

// firstCrossing returns the first parameter in (t0, t1) at which the ray
// A + t·B enters or leaves the solid, the cell it meets there and whether it
// starts inside. outside tells that the ray is known to start outside the
// solid; events is an empty buffer with room for one event per cell.
func (s *starShape) firstCrossing(A, B Vector4, t0, t1 Real, outside bool, events []starEvent) (Real, int, bool, bool) {
	// Crossings of the cells: the density falls by the density of the cell
	// where the ray leaves the side of the centre, and rises where it enters.
	density := 0
	for k := range s.normals {
		n := &s.normals[k]
		dn := n.X*B.X + n.Y*B.Y + n.Z*B.Z + n.W*B.W
		if dn == 0 || (outside && dn > 0) {
			continue
		}
		t := (s.height - (n.X*A.X + n.Y*A.Y + n.Z*A.Z + n.W*A.W)) / dn
		if !(t > t0 && t < t1) {
			continue
		}
		c := 1
		if !outside {
			if c = s.crossingDensity(A, B, t, k); c == 0 {
				continue
			}
			if dn > 0 {
				c = -c
			}
			// beyond all crossings the density is 0
			density -= c
		}
		events = append(events, starEvent{t: t, jump: c, cell: k})
	}
	// in the order of t, as far as the first change between nothing and something
	for i := len(events)/2 - 1; i >= 0; i-- {
		starEventDown(events, i)
	}
	for len(events) > 0 {
		e := events[0]
		if outside {
			// the first cell that is entered where it is solid
			if s.crossingDensity(A, B, e.t, e.cell) > 0 {
				return e.t, e.cell, false, true
			}
		} else if next := density + e.jump; (density > 0) != (next > 0) {
			return e.t, e.cell, density > 0, true
		} else {
			density = next
		}
		events[0] = events[len(events)-1]
		events = events[:len(events)-1]
		starEventDown(events, 0)
	}
	return 0, 0, false, false
}

// crossingDensity returns the density of the cell at the point A + t·B of
// its hyperplane.
func (s *starShape) crossingDensity(A, B Vector4, t Real, cell int) int {
	X := Vector4{A.X + t*B.X, A.Y + t*B.Y, A.Z + t*B.Z, A.W + t*B.W}
	f := &s.frames[cell]
	return s.cellDensity([3]Real{f[0].Dot(X), f[1].Dot(X), f[2].Dot(X)})
}

func intersectRayStarPolychoron(O Point4, D Vector4, sp *StarPolychoron) (objectHit, bool) {
	const eps = 1e-12
	Ol := sp.RT.MulVec(Vector4{O.X - sp.Center.X, O.Y - sp.Center.Y, O.Z - sp.Center.Z, O.W - sp.Center.W})
	Dl := sp.RT.MulVec(D)
	A := Vector4{Ol.X / sp.Scale.X, Ol.Y / sp.Scale.Y, Ol.Z / sp.Scale.Z, Ol.W / sp.Scale.W}
	B := Vector4{Dl.X / sp.Scale.X, Dl.Y / sp.Scale.Y, Dl.Z / sp.Scale.Z, Dl.W / sp.Scale.W}

	// the part of the ray inside the unit ball, which holds the solid
	bb, ab := B.Dot(B), A.Dot(B)
	disc := ab*ab - bb*(A.Dot(A)-1)
	if !(bb > 0) || disc <= 0 {
		return objectHit{}, false
	}
	root := math.Sqrt(disc)
	tIn, tOut := (-ab-root)/bb, (-ab+root)/bb
	if tOut <= eps {
		return objectHit{}, false
	}
	// a ray that comes from outside the ball starts outside the solid
	outside := tIn >= eps
	if !outside {
		tIn = eps
	}

	s := sp.shape
	var t Real
	var cell int
	var inside, ok bool
	var store [128]starEvent
	if len(s.normals) <= len(store) {
		t, cell, inside, ok = s.firstCrossing(A, B, tIn, tOut, outside, store[:0])
	} else {
		buf := starEventPool.Get().(*[]starEvent)
		if cap(*buf) < len(s.normals) {
			*buf = make([]starEvent, 0, len(s.normals))
		}
		t, cell, inside, ok = s.firstCrossing(A, B, tIn, tOut, outside, (*buf)[:0])
		starEventPool.Put(buf)
	}
	if !ok {
		return objectHit{}, false
	}
	n := s.normals[cell]
	Nw := sp.R.MulVec(Vector4{n.X / sp.Scale.X, n.Y / sp.Scale.Y, n.Z / sp.Scale.Z, n.W / sp.Scale.W}).Norm()
	return objectHit{t: t, Nw: Nw, mat: sp, inv: inside}, true
}

func (s *StarPolychoron) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *StarPolychoron) DiffCh(i int) Real   { return s.diff[i] }
func (s *StarPolychoron) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *StarPolychoron) F0Ch(i int) Real     { return s.f0[i] }
func (s *StarPolychoron) ReflCh(i int) Real   { return s.refl[i] }
func (s *StarPolychoron) RefrCh(i int) Real   { return s.refr[i] }
func (s *StarPolychoron) IORCh(i int) Real    { return s.iorArr[i] }
func (s *StarPolychoron) IORInvCh(i int) Real { return s.iorInv[i] }
