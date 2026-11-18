package photons4d

import (
	"fmt"
	"math"
	"strings"
	"sync"
)

// Polytope is a convex polychoron given by its cells. It serves the uniform
// polychora and their duals, the prisms, pyramids and bipyramids over the
// uniform polyhedra and the duopyramids.
type Polytope struct {
	cellPoly
	Name string
}

// Branch marks of the Coxeter diagram of a regular polychoron; its own node
// is the first one.
var polychoronBases = map[string][]Real{
	"5-cell":             {3, 3, 3},
	"pentachoron":        {3, 3, 3},
	"8-cell":             {4, 3, 3},
	"tesseract":          {4, 3, 3},
	"hypercube":          {4, 3, 3},
	"octachoron":         {4, 3, 3},
	"16-cell":            {3, 3, 4},
	"hexadecachoron":     {3, 3, 4},
	"24-cell":            {3, 4, 3},
	"icositetrachoron":   {3, 4, 3},
	"120-cell":           {5, 3, 3},
	"hecatonicosachoron": {5, 3, 3},
	"600-cell":           {3, 3, 5},
	"hexacosichoron":     {3, 3, 5},
}

// Ringed nodes of the operations that turn a regular polychoron into a
// uniform one.
var polychoronOps = map[string][]bool{
	"":                    {true, false, false, false},
	"rectified":           {false, true, false, false},
	"birectified":         {false, false, true, false},
	"truncated":           {true, true, false, false},
	"bitruncated":         {false, true, true, false},
	"cantellated":         {true, false, true, false},
	"runcinated":          {true, false, false, true},
	"cantitruncated":      {true, true, true, false},
	"runcitruncated":      {true, true, false, true},
	"runcicantellated":    {true, false, true, true},
	"omnitruncated":       {true, true, true, true},
	"runcicantitruncated": {true, true, true, true},
}

// snub24Vertices returns the 96 vertices of the snub 24-cell: those of the
// 600-cell without the 24 of an inscribed 24-cell (the points with the
// coordinates 0, ±1/2 and ±1).
func snub24Vertices() []Vector4 {
	var out []Vector4
	for _, v := range verts600Unit() {
		a := [4]Real{math.Abs(v.X), math.Abs(v.Y), math.Abs(v.Z), math.Abs(v.W)}
		halves := 0
		for _, c := range a {
			if math.Abs(c-0.5) < 1e-9 {
				halves++
			}
		}
		if halves == 4 || math.Max(math.Max(a[0], a[1]), math.Max(a[2], a[3])) > 1-1e-9 {
			continue
		}
		out = append(out, v)
	}
	return out
}

// grandAntiprismVertices returns the 100 vertices of the grand antiprism:
// those of the 600-cell without two rings of ten that lie in two planes
// square to each other.
func grandAntiprismVertices() []Vector4 {
	phi := (1 + math.Sqrt(5)) / 2
	e := Vector4{1, 0, 0, 0}
	q := Vector4{phi / 2, 0.5, 0.5 / phi, 0}
	u := q.Sub(e.Mul(q.Dot(e))).Norm()
	var out []Vector4
	for _, v := range verts600Unit() {
		a, b := v.Dot(e), v.Dot(u)
		if in := a*a + b*b; in > 1e-9 && in < 1-1e-9 {
			out = append(out, v)
		}
	}
	return out
}

type polytopeShape struct {
	verts, planes []Vector4
	err           error
}

var polytopeShapes sync.Map // name -> *polytopeShape

func cachedShape(key string, build func() ([]Vector4, []Vector4, error)) ([]Vector4, []Vector4, error) {
	if v, ok := polytopeShapes.Load(key); ok {
		s := v.(*polytopeShape)
		return s.verts, s.planes, s.err
	}
	s := &polytopeShape{}
	s.verts, s.planes, s.err = build()
	polytopeShapes.Store(key, s)
	return s.verts, s.planes, s.err
}

// uniformPolychoron returns the vertices (unit circumradius) and the cell
// planes n·x <= 1 of a convex uniform polychoron given by
//   - an operation and a regular polychoron: "truncated 120-cell",
//     "rectified tesseract", "omnitruncated 24-cell", "600-cell", ...
//   - "snub 24-cell" or "grand antiprism",
//   - a Coxeter diagram: "x5o3o3x" (a mark 2 gives prisms and duoprisms).
func uniformPolychoron(kind string) (verts, planes []Vector4, err error) {
	key := shapeKey(kind)
	return cachedShape("polychoron "+key, func() ([]Vector4, []Vector4, error) {
		switch key {
		case "snub 24-cell":
			v := snub24Vertices()
			return v, hullPlanes(v, 4), nil
		case "grand antiprism":
			v := grandAntiprismVertices()
			return v, hullPlanes(v, 4), nil
		}
		marks, rings, ok := parseCoxeterDiagram(key)
		if !ok {
			op, base := "", key
			if i := strings.LastIndexByte(key, ' '); i >= 0 {
				op, base = key[:i], key[i+1:]
			}
			if marks, ok = polychoronBases[base]; ok {
				rings, ok = polychoronOps[op]
			}
		}
		if !ok || len(rings) != 4 || !wholeMarks(marks) {
			return nil, nil, fmt.Errorf("unknown polychoron %q", kind)
		}
		v, p, err := wythoffPolytope(marks, rings)
		if err != nil {
			return nil, nil, fmt.Errorf("polychoron %q: %w", kind, err)
		}
		return v, p, nil
	})
}

// dualShape returns the polar dual of a polytope around the origin, scaled to
// the unit circumradius: its vertices are the plane vectors of the polytope
// and its planes are the vertices.
func dualShape(verts, planes []Vector4) (dverts, dplanes []Vector4) {
	far := Real(0)
	for _, n := range planes {
		far = math.Max(far, n.Len())
	}
	dverts = make([]Vector4, len(planes))
	for i, n := range planes {
		dverts[i] = n.Mul(1 / far)
	}
	dplanes = make([]Vector4, len(verts))
	for i, v := range verts {
		dplanes[i] = v.Mul(far)
	}
	return dverts, dplanes
}

// liftPolyhedron turns a polyhedron of XYZ (planes n·x <= 1) into a
// polychoron along W:
//   - "prism": the polyhedron times the segment |w| <= half of its edge, so
//     that the prism over a uniform polyhedron is uniform,
//   - "pyramid": the polyhedron at w = +1 joined to the apex at w = -1,
//   - "bipyramid": the polyhedron at w = 0 joined to the apexes at w = ±1.
func liftPolyhedron(shape string, verts, planes []Vector4) (lverts, lplanes []Vector4, err error) {
	switch shape {
	case "prism":
		h := math.Inf(1)
		for i := range verts {
			for j := i + 1; j < len(verts); j++ {
				h = math.Min(h, verts[i].Sub(verts[j]).Len()/2)
			}
		}
		for _, v := range verts {
			lverts = append(lverts, Vector4{v.X, v.Y, v.Z, h}, Vector4{v.X, v.Y, v.Z, -h})
		}
		lplanes = append(lplanes, planes...)
		lplanes = append(lplanes, Vector4{0, 0, 0, 1 / h}, Vector4{0, 0, 0, -1 / h})
	case "pyramid":
		for _, v := range verts {
			lverts = append(lverts, Vector4{v.X, v.Y, v.Z, 1})
		}
		lverts = append(lverts, Vector4{0, 0, 0, -1})
		for _, n := range planes {
			lplanes = append(lplanes, Vector4{2 * n.X, 2 * n.Y, 2 * n.Z, -1})
		}
		lplanes = append(lplanes, Vector4{0, 0, 0, 1})
	case "bipyramid":
		lverts = append(lverts, verts...)
		lverts = append(lverts, Vector4{0, 0, 0, 1}, Vector4{0, 0, 0, -1})
		for _, n := range planes {
			lplanes = append(lplanes, Vector4{n.X, n.Y, n.Z, 1}, Vector4{n.X, n.Y, n.Z, -1})
		}
	default:
		return nil, nil, fmt.Errorf("unknown polyhedral solid %q (want prism, pyramid or bipyramid)", shape)
	}
	return lverts, lplanes, nil
}

// checkSolid validates the scale and the material of an object.
func checkSolid(what string, scale Vector4, diffuse, reflectivity, refractivity, ior RGB) error {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return fmt.Errorf("%s scale must be > 0 on all axes, got %+v", what, scale)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	for _, c := range []struct {
		n       string
		r, t, d Real
	}{
		{"R", reflectivity.R, refractivity.R, diffuse.R},
		{"G", reflectivity.G, refractivity.G, diffuse.G},
		{"B", reflectivity.B, refractivity.B, diffuse.B},
	} {
		if !in01(c.r) || !in01(c.t) || !in01(c.d) {
			return fmt.Errorf("reflect/refract/diffuse must be in [0,1]; channel %s got reflect=%.6g refract=%.6g diffuse=%.6g", c.n, c.r, c.t, c.d)
		}
		if c.r+c.t+c.d > 1+1e-12 {
			return fmt.Errorf("per-channel reflect+refract+diffuse must be ≤1; channel %s got %.6g", c.n, c.r+c.t+c.d)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}
	return nil
}

func newPolytope(name string, verts, planes []Vector4, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Polytope, error) {
	if err := checkSolid(name, scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	R := rotFromAngles(angles)
	cp := cellPoly{
		Center: center,
		Scale:  scale,
		R:      R,
		RT:     R.Transpose(),
	}
	cp.materialFrom(color, diffuse, reflectivity, refractivity, ior)
	cp.buildPlanes(planes, 1.0)
	cp.computeAABBFromLocalVerts(verts)
	DebugLog("Created %s: %d vertices, %d cells, AABB=[%+v..%+v]", name, len(verts), len(planes), cp.AABBMin, cp.AABBMax)
	return &Polytope{cellPoly: cp, Name: name}, nil
}

// NewPolychoron builds a convex uniform polychoron (see uniformPolychoron for
// the names) or, with dual set, its dual. Before Scale is applied the
// vertices lie on the unit sphere; the farthest ones in the case of a dual.
func NewPolychoron(kind string, dual bool, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Polytope, error) {
	verts, planes, err := uniformPolychoron(kind)
	if err != nil {
		return nil, err
	}
	name := shapeKey(kind)
	if dual {
		verts, planes = dualShape(verts, planes)
		name = "dual " + name
	}
	return newPolytope(name, verts, planes, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewPolyhedralSolid builds the prism, the pyramid or the bipyramid (shape)
// over a uniform polyhedron (base, see uniformPolyhedron) that lies in the
// local XYZ space with the circumradius 1; liftPolyhedron tells where the
// apexes and the ends are before Scale is applied.
func NewPolyhedralSolid(shape, base string, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Polytope, error) {
	key := shapeKey(base)
	verts, planes, err := cachedShape(shape+" "+key, func() ([]Vector4, []Vector4, error) {
		v, p, err := uniformPolyhedron(base)
		if err != nil {
			return nil, nil, err
		}
		return liftPolyhedron(shape, v, p)
	})
	if err != nil {
		return nil, err
	}
	return newPolytope(key+" "+shape, verts, planes, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// Largest number of cells accepted for a duopyramid.
const maxDuopyramidCells = 16384

// NewDuopyramid builds the dual of the duoprism: the hull of a regular P-gon
// in the local XY plane and a regular Q-gon in the local ZW plane, both with
// the circumradius 1 before Scale is applied. It has P + Q vertices and P·Q
// tetrahedral cells; the 4-4 duopyramid is the 16-cell.
func NewDuopyramid(center Point4, scale Vector4, p, q int, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Polytope, error) {
	if p < 3 || q < 3 || p*q > maxDuopyramidCells {
		return nil, fmt.Errorf("duopyramid needs polygons with at least 3 sides and at most %d cells, got p=%d q=%d", maxDuopyramidCells, p, q)
	}
	polygon := func(n int) (corners, sides [][2]Real) {
		l := 1 / math.Cos(math.Pi/Real(n))
		for k := 0; k < n; k++ {
			a := 2 * math.Pi * Real(k) / Real(n)
			b := a + math.Pi/Real(n)
			corners = append(corners, [2]Real{math.Cos(a), math.Sin(a)})
			sides = append(sides, [2]Real{l * math.Cos(b), l * math.Sin(b)})
		}
		return corners, sides
	}
	cp, sp := polygon(p)
	cq, sq := polygon(q)
	verts := make([]Vector4, 0, p+q)
	for _, c := range cp {
		verts = append(verts, Vector4{c[0], c[1], 0, 0})
	}
	for _, c := range cq {
		verts = append(verts, Vector4{0, 0, c[0], c[1]})
	}
	planes := make([]Vector4, 0, p*q)
	for _, a := range sp {
		for _, b := range sq {
			planes = append(planes, Vector4{a[0], a[1], b[0], b[1]})
		}
	}
	return newPolytope(fmt.Sprintf("%d-%d duopyramid", p, q), verts, planes, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}
