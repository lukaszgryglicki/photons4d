package photons4d

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Coxeter diagrams of the uniform polyhedra with a reflection group.
var polyhedronDiagrams = map[string]string{
	"tetrahedron":                 "x3o3o",
	"cube":                        "x4o3o",
	"hexahedron":                  "x4o3o",
	"octahedron":                  "o4o3x",
	"dodecahedron":                "x5o3o",
	"icosahedron":                 "o5o3x",
	"truncated tetrahedron":       "x3x3o",
	"cuboctahedron":               "o4x3o",
	"truncated cube":              "x4x3o",
	"truncated octahedron":        "o4x3x",
	"rhombicuboctahedron":         "x4o3x",
	"truncated cuboctahedron":     "x4x3x",
	"icosidodecahedron":           "o5x3o",
	"truncated dodecahedron":      "x5x3o",
	"truncated icosahedron":       "o5x3x",
	"rhombicosidodecahedron":      "x5o3x",
	"truncated icosidodecahedron": "x5x3x",
}

var polygonAdjectives = map[string]int{
	"triangular": 3, "square": 4, "pentagonal": 5, "hexagonal": 6, "heptagonal": 7, "octagonal": 8,
	"enneagonal": 9, "nonagonal": 9, "decagonal": 10, "hendecagonal": 11, "dodecagonal": 12,
}

// shapeKey lowers a name and joins its words with single spaces.
func shapeKey(name string) string {
	return strings.Join(strings.Fields(strings.ToLower(name)), " ")
}

// parseCoxeterDiagram reads a linear Coxeter diagram such as "x5o3o3x": the
// nodes are x (ringed) or o, the numbers between them are the branch marks
// (2 for no branch, fractions such as 5/2 for star figures).
func parseCoxeterDiagram(s string) (marks []Real, rings []bool, ok bool) {
	for i := 0; ; {
		if i >= len(s) || (s[i] != 'x' && s[i] != 'o') {
			return nil, nil, false
		}
		rings = append(rings, s[i] == 'x')
		i++
		if i == len(s) {
			return marks, rings, len(rings) >= 2
		}
		j := i
		for j < len(s) && (s[j] >= '0' && s[j] <= '9' || s[j] == '/') {
			j++
		}
		num, den := s[i:j], "1"
		if k := strings.IndexByte(num, '/'); k >= 0 {
			num, den = num[:k], num[k+1:]
		}
		a, errA := strconv.Atoi(num)
		b, errB := strconv.Atoi(den)
		if errA != nil || errB != nil || b < 1 || a < 2*b {
			return nil, nil, false
		}
		marks = append(marks, Real(a)/Real(b))
		i = j
	}
}

func wholeMarks(marks []Real) bool {
	for _, m := range marks {
		if m != math.Trunc(m) {
			return false
		}
	}
	return true
}

// parsePrismName reads "7-gonal prism", "7-prism", "heptagonal antiprism", ...
func parsePrismName(key string) (n int, anti, ok bool) {
	i := strings.LastIndexAny(key, " -")
	if i < 0 {
		return 0, false, false
	}
	switch key[i+1:] {
	case "prism":
	case "antiprism":
		anti = true
	default:
		return 0, false, false
	}
	rest := key[:i]
	if v, found := polygonAdjectives[rest]; found {
		return v, anti, true
	}
	n, err := strconv.Atoi(strings.TrimSuffix(rest, "-gonal"))
	return n, anti, err == nil
}

// prismVertices returns the vertices of the uniform prism or antiprism over a
// regular polygon with n sides: all its edges are equal and its axis is Z.
func prismVertices(n int, anti bool) []Vector4 {
	a := math.Pi / Real(n)
	h := math.Sin(a)
	if anti {
		h = math.Sqrt((math.Cos(a) - math.Cos(2*a)) / 2)
	}
	s := 1 / math.Hypot(1, h)
	verts := make([]Vector4, 0, 2*n)
	for k := 0; k < n; k++ {
		top := 2 * a * Real(k)
		bottom := top
		if anti {
			bottom += a
		}
		verts = append(verts,
			Vector4{s * math.Cos(top), s * math.Sin(top), s * h, 0},
			Vector4{s * math.Cos(bottom), s * math.Sin(bottom), -s * h, 0})
	}
	return verts
}

// snubVertices returns the vertices of the snub cube (p = 4) or of the snub
// dodecahedron (p = 5): the images, under the rotations of the group of
// {p,3}, of the point that is equally far from its images under the three
// rotations about the corners of a mirror triangle.
func snubVertices(p Real) ([]Vector4, error) {
	mirrors, dual, err := coxeterMirrors([]Real{p, 3})
	if err != nil {
		return nil, err
	}
	turn := func(i, j int, v Vector4) Vector4 {
		return mirrorImage(mirrorImage(v, mirrors[j]), mirrors[i])
	}
	f := func(x Vector4) [3]Real {
		a, b, c := x.Sub(turn(0, 1, x)), x.Sub(turn(1, 2, x)), x.Sub(turn(0, 2, x))
		return [3]Real{a.Dot(a) - c.Dot(c), b.Dot(b) - c.Dot(c), x.Dot(x) - 1}
	}
	x := dual[0].Norm().Add(dual[1].Norm()).Add(dual[2].Norm()).Norm()
	for it := 0; it < 100; it++ {
		r := f(x)
		if math.Abs(r[0])+math.Abs(r[1])+math.Abs(r[2]) < 1e-15 {
			break
		}
		// the equations are quadratic: central differences give their exact slopes
		var j [3][3]Real
		for c := 0; c < 3; c++ {
			const h = 1e-3
			up, down := f(vecWith(x, c, vecAt(x, c)+h)), f(vecWith(x, c, vecAt(x, c)-h))
			for k := 0; k < 3; k++ {
				j[k][c] = (up[k] - down[k]) / (2 * h)
			}
		}
		det := j[0][0]*(j[1][1]*j[2][2]-j[1][2]*j[2][1]) - j[0][1]*(j[1][0]*j[2][2]-j[1][2]*j[2][0]) + j[0][2]*(j[1][0]*j[2][1]-j[1][1]*j[2][0])
		if det == 0 {
			break
		}
		var step [3]Real
		for c := 0; c < 3; c++ {
			m := j
			for k := 0; k < 3; k++ {
				m[k][c] = r[k]
			}
			step[c] = (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])) / det
		}
		x = Vector4{x.X - step[0], x.Y - step[1], x.Z - step[2], 0}
	}
	if r := f(x); math.Abs(r[0])+math.Abs(r[1])+math.Abs(r[2]) > 1e-12 {
		return nil, fmt.Errorf("no snub vertex found for the group of {%g,3}", p)
	}
	set := newPointSet()
	set.add(x)
	for i := 0; i < len(set.pts) && len(set.pts) <= maxOrbit; i++ {
		set.add(turn(0, 1, set.pts[i]))
		set.add(turn(1, 2, set.pts[i]))
	}
	return set.pts, nil
}

// uniformPolyhedron returns the vertices (unit circumradius, in XYZ) and the
// face planes n·x <= 1 of a convex uniform polyhedron: one of the 5 Platonic
// and 13 Archimedean solids by name, a prism or an antiprism ("7-gonal
// prism", "pentagonal antiprism"), or a Coxeter diagram such as "x5x3o".
func uniformPolyhedron(name string) (verts, planes []Vector4, err error) {
	key := shapeKey(name)
	if n, anti, ok := parsePrismName(key); ok {
		if n < 3 || n > maxDuoprismSides {
			return nil, nil, fmt.Errorf("polyhedron %q: the polygon needs 3..%d sides", name, maxDuoprismSides)
		}
		verts = prismVertices(n, anti)
	} else if key == "snub cube" {
		verts, err = snubVertices(4)
	} else if key == "snub dodecahedron" {
		verts, err = snubVertices(5)
	} else {
		diagram, known := polyhedronDiagrams[key]
		if !known {
			diagram = key
		}
		marks, rings, ok := parseCoxeterDiagram(diagram)
		if !ok || len(rings) != 3 || !wholeMarks(marks) {
			return nil, nil, fmt.Errorf("unknown polyhedron %q", name)
		}
		if verts, planes, err = wythoffPolytope(marks, rings); err != nil {
			return nil, nil, fmt.Errorf("polyhedron %q: %w", name, err)
		}
		return verts, planes, nil
	}
	if err != nil {
		return nil, nil, err
	}
	if planes = hullPlanes(verts, 3); len(planes) < 4 {
		return nil, nil, fmt.Errorf("polyhedron %q has no faces", name)
	}
	return verts, planes, nil
}
