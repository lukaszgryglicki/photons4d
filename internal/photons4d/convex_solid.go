package photons4d

import (
	"fmt"
	"math"
)

// ConvexSolid is a convex body given by a convex function of the local
// coordinates divided by Scale: negative inside, positive outside. It serves
// the bracketopes (bodies of the nested norms of the bracket notation) and
// the cone-like solids: coninder, cylindrone, dicone and the prisms over a
// disk and a polygon.
type ConvexSolid struct {
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

	f     func(q Vector4) (Real, Vector4) // the function and its gradient
	bound Real                            // radius of a ball around the local origin that holds the solid

	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

// newConvexSolid builds the solid {f <= 0}. support(a) is the largest value
// of a·q over the solid, bound the largest length of q.
func newConvexSolid(name string, f func(Vector4) (Real, Vector4), support func(Vector4) Real, bound Real, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	if err := checkSolid(name, scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	R := rotFromAngles(angles)
	s := &ConvexSolid{
		Center:  center,
		Scale:   scale,
		R:       R,
		RT:      R.Transpose(),
		Name:    name,
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
		f:       f,
		bound:   bound,
	}
	// A world axis sees the local point q as c + a·q with a = Scale times the row of R.
	extent := func(row int, c Real) (Real, Real) {
		a := Vector4{scale.X * R.M[row][0], scale.Y * R.M[row][1], scale.Z * R.M[row][2], scale.W * R.M[row][3]}
		return c - support(a.Mul(-1)), c + support(a)
	}
	minX, maxX := extent(0, center.X)
	minY, maxY := extent(1, center.Y)
	minZ, maxZ := extent(2, center.Z)
	minW, maxW := extent(3, center.W)
	s.AABBMin = Point4{minX, minY, minZ, minW}
	s.AABBMax = Point4{maxX, maxY, maxZ, maxW}

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
	DebugLog("Created %s: %+v", name, s)
	return s, nil
}

// traceConvex returns the first parameter in (t0, t1] at which a function
// that is convex along the ray changes its sign, and whether it is negative
// at t0. at returns the function and its slope along the ray; the function
// must be positive at t1 if it is not positive at t0.
func traceConvex(at func(Real) (Real, Real), t0, t1 Real) (Real, bool, bool) {
	value := func(t Real) rayDist {
		v, _ := at(t)
		return rayDist{d: v}
	}
	v0, _ := at(t0)
	if v0 <= 0 {
		v1, _ := at(t1)
		if v1 <= 0 {
			return 0, true, false
		}
		return refineCrossing(value, t0, v0, t1, v1, true), true, true
	}
	// outside: go down to the least value until a point of the solid is met
	lo, vlo, hi := t0, v0, t1
	for i := 0; i < 200 && hi-lo > 1e-15*math.Max(math.Abs(lo), math.Abs(hi)); i++ {
		m := 0.5 * (lo + hi)
		v, slope := at(m)
		if v <= 0 {
			return refineCrossing(value, lo, vlo, m, v, false), false, true
		}
		if slope > 0 {
			hi = m
		} else {
			lo, vlo = m, v
		}
	}
	return 0, false, false
}

func intersectRayConvexSolid(O Point4, D Vector4, s *ConvexSolid) (objectHit, bool) {
	const eps = 1e-12
	Ol := s.RT.MulVec(Vector4{O.X - s.Center.X, O.Y - s.Center.Y, O.Z - s.Center.Z, O.W - s.Center.W})
	Dl := s.RT.MulVec(D)
	A := Vector4{Ol.X / s.Scale.X, Ol.Y / s.Scale.Y, Ol.Z / s.Scale.Z, Ol.W / s.Scale.W}
	B := Vector4{Dl.X / s.Scale.X, Dl.Y / s.Scale.Y, Dl.Z / s.Scale.Z, Dl.W / s.Scale.W}

	// the part of the ray inside a ball that holds the solid with a margin
	rb := s.bound * (1 + 1e-6)
	bb, ab := B.Dot(B), A.Dot(B)
	disc := ab*ab - bb*(A.Dot(A)-rb*rb)
	if !(bb > 0) || disc <= 0 {
		return objectHit{}, false
	}
	root := math.Sqrt(disc)
	t0, t1 := (-ab-root)/bb, (-ab+root)/bb
	if t1 <= eps {
		return objectHit{}, false
	}
	if t0 < eps {
		t0 = eps
	}
	t, inside, ok := traceConvex(func(t Real) (Real, Real) {
		v, g := s.f(A.Add(B.Mul(t)))
		return v, g.Dot(B)
	}, t0, t1)
	if !ok {
		return objectHit{}, false
	}
	_, g := s.f(A.Add(B.Mul(t)))
	Nl := Vector4{g.X / s.Scale.X, g.Y / s.Scale.Y, g.Z / s.Scale.Z, g.W / s.Scale.W}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	return objectHit{t: t, Nw: s.R.MulVec(Nl).Norm(), mat: s, inv: inside}, true
}

func (s *ConvexSolid) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *ConvexSolid) DiffCh(i int) Real   { return s.diff[i] }
func (s *ConvexSolid) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *ConvexSolid) F0Ch(i int) Real     { return s.f0[i] }
func (s *ConvexSolid) ReflCh(i int) Real   { return s.refl[i] }
func (s *ConvexSolid) RefrCh(i int) Real   { return s.refr[i] }
func (s *ConvexSolid) IORCh(i int) Real    { return s.iorArr[i] }
func (s *ConvexSolid) IORInvCh(i int) Real { return s.iorInv[i] }

// ---- bracketopes ----

// gaugeNode is a node of the bracket notation: an axis (I), or a group that
// joins the values of its members by the root of the sum of squares (round
// brackets), by the largest (square brackets) or by the sum (angle brackets).
type gaugeNode struct {
	op   byte
	axis int
	kids []*gaugeNode
}

// Names of some bracketopes.
var bracketopeNames = map[string]string{
	"glome":               "(IIII)",
	"tesseract":           "[IIII]",
	"16-cell":             "<IIII>",
	"spherinder":          "[(III)I]",
	"cubinder":            "[(II)II]",
	"duocylinder":         "[(II)(II)]",
	"spherical bipyramid": "<(III)I>",
	"duospindle":          "<(II)(II)>",
	"bicircular tegum":    "<(II)(II)>",
}

// parseBracketope reads the bracket notation of a 4D shape: I stands for an
// axis (X, Y, Z, W in this order), ( ) for a round join, [ ] for a square one
// and < > for a diamond one. "[(II)(II)]" is the duocylinder, "<(II)(II)>"
// its dual, "(IIII)" the ball.
func parseBracketope(s string) (*gaugeNode, error) {
	closing := map[byte]byte{'(': ')', '[': ']', '<': '>'}
	axes := 0
	pos := 0
	var parse func() (*gaugeNode, error)
	parse = func() (*gaugeNode, error) {
		if pos >= len(s) {
			return nil, fmt.Errorf("the notation ends too early")
		}
		c := s[pos]
		pos++
		if c == 'I' {
			axes++
			return &gaugeNode{op: 'I', axis: axes - 1}, nil
		}
		end, ok := closing[c]
		if !ok {
			return nil, fmt.Errorf("unexpected %q", string(c))
		}
		n := &gaugeNode{op: c}
		for pos < len(s) && s[pos] != end {
			k, err := parse()
			if err != nil {
				return nil, err
			}
			n.kids = append(n.kids, k)
		}
		if pos >= len(s) {
			return nil, fmt.Errorf("missing %q", string(end))
		}
		pos++
		if len(n.kids) == 0 {
			return nil, fmt.Errorf("empty brackets")
		}
		return n, nil
	}
	root, err := parse()
	if err == nil && pos != len(s) {
		err = fmt.Errorf("more than one group")
	}
	if err == nil && (axes != 4 || root.op == 'I') {
		err = fmt.Errorf("a 4D shape needs brackets around 4 axes, got %d", axes)
	}
	if err != nil {
		return nil, fmt.Errorf("bracketope %q: %w", s, err)
	}
	return root, nil
}

// eval returns the gauge of the point and its gradient.
func (n *gaugeNode) eval(q Vector4) (Real, Vector4) {
	switch n.op {
	case 'I':
		if x := vecAt(q, n.axis); x < 0 {
			return -x, vecWith(Vector4{}, n.axis, -1)
		} else {
			return x, vecWith(Vector4{}, n.axis, 1)
		}
	case '[':
		best, grad := math.Inf(-1), Vector4{}
		for _, k := range n.kids {
			if v, g := k.eval(q); v > best {
				best, grad = v, g
			}
		}
		return best, grad
	case '<':
		sum, grad := Real(0), Vector4{}
		for _, k := range n.kids {
			v, g := k.eval(q)
			sum += v
			grad = grad.Add(g)
		}
		return sum, grad
	}
	var vs [4]Real
	var gs [4]Vector4
	sq := Real(0)
	for i, k := range n.kids {
		vs[i], gs[i] = k.eval(q)
		sq += vs[i] * vs[i]
	}
	r := math.Sqrt(sq)
	grad := Vector4{}
	if r > 0 {
		for i := range n.kids {
			grad = grad.Add(gs[i].Mul(vs[i] / r))
		}
	}
	return r, grad
}

// dual returns the gauge of the polar body: the largest value of a·q over
// the body of this gauge.
func (n *gaugeNode) dual(a Vector4) Real {
	switch n.op {
	case 'I':
		return math.Abs(vecAt(a, n.axis))
	case '[':
		sum := Real(0)
		for _, k := range n.kids {
			sum += k.dual(a)
		}
		return sum
	case '<':
		best := Real(0)
		for _, k := range n.kids {
			best = math.Max(best, k.dual(a))
		}
		return best
	}
	sq := Real(0)
	for _, k := range n.kids {
		v := k.dual(a)
		sq += v * v
	}
	return math.Sqrt(sq)
}

// NewBracketope builds the body of a bracket notation (see parseBracketope)
// or of one of the names in bracketopeNames: all points whose nested norm is
// at most 1 after the division of the local coordinates by Scale.
func NewBracketope(shape string, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	notation := shape
	if v, ok := bracketopeNames[shapeKey(shape)]; ok {
		notation = v
	}
	root, err := parseBracketope(notation)
	if err != nil {
		return nil, err
	}
	f := func(q Vector4) (Real, Vector4) {
		v, g := root.eval(q)
		return v - 1, g
	}
	// no coordinate of a point of the body is above 1
	return newConvexSolid("bracketope "+notation, f, root.dual, 2, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// ---- cone-like solids ----

// conicPart is the set where the length of the coordinates on the axes of
// mask is at most k + b·q: a ball, a cylinder or a slab for b = 0, a cone
// otherwise, a half-space for an empty mask.
type conicPart struct {
	mask [4]bool
	k    Real
	b    Vector4
}

// conicCap is a ball over the axes of mask around a point; the solid is the
// hull of its caps.
type conicCap struct {
	at   Vector4
	mask [4]bool
}

func masked(v Vector4, mask [4]bool) Vector4 {
	var m Vector4
	if mask[0] {
		m.X = v.X
	}
	if mask[1] {
		m.Y = v.Y
	}
	if mask[2] {
		m.Z = v.Z
	}
	if mask[3] {
		m.W = v.W
	}
	return m
}

func newConicSolid(name string, parts []conicPart, caps []conicCap, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	f := func(q Vector4) (Real, Vector4) {
		best, grad := math.Inf(-1), Vector4{}
		for i := range parts {
			p := &parts[i]
			m := masked(q, p.mask)
			l := m.Len()
			if v := l - p.k - p.b.Dot(q); v > best {
				best, grad = v, p.b.Mul(-1)
				if l > 0 {
					grad = grad.Add(m.Mul(1 / l))
				}
			}
		}
		return best, grad
	}
	support := func(a Vector4) Real {
		best := math.Inf(-1)
		for i := range caps {
			best = math.Max(best, a.Dot(caps[i].at)+masked(a, caps[i].mask).Len())
		}
		return best
	}
	bound := Real(0)
	for i := range caps {
		reach := caps[i].at.Len()
		if caps[i].mask != [4]bool{} {
			reach++
		}
		bound = math.Max(bound, reach)
	}
	return newConvexSolid(name, f, support, bound, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

var (
	axesXY = [4]bool{true, true, false, false}
	axisZ  = [4]bool{false, false, true, false}
	axisW  = [4]bool{false, false, false, true}
)

// NewConinder builds the prism over a cone: in local coordinates divided by
// Scale = (rx, ry, h, hw)
//
//	sqrt(x² + y²) <= (z + 1)/2,  z <= 1,  |w| <= 1
//
// The cone has its apex at z = -1 and its base disk at z = +1.
func NewConinder(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	parts := []conicPart{
		{mask: axesXY, k: 0.5, b: Vector4{0, 0, 0.5, 0}},
		{k: 1, b: Vector4{0, 0, -1, 0}},
		{mask: axisW, k: 1},
	}
	caps := []conicCap{
		{at: Vector4{0, 0, 1, 1}, mask: axesXY}, {at: Vector4{0, 0, 1, -1}, mask: axesXY},
		{at: Vector4{0, 0, -1, 1}}, {at: Vector4{0, 0, -1, -1}},
	}
	return newConicSolid("coninder", parts, caps, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewCylindrone builds the pyramid over a cylinder: in local coordinates
// divided by Scale = (rx, ry, hz, h), with u = (w + 1)/2,
//
//	sqrt(x² + y²) <= u,  |z| <= u,  w <= 1
//
// The apex is at w = -1 and the base cylinder at w = +1.
func NewCylindrone(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	up := Vector4{0, 0, 0, 0.5}
	parts := []conicPart{
		{mask: axesXY, k: 0.5, b: up},
		{mask: axisZ, k: 0.5, b: up},
		{k: 1, b: Vector4{0, 0, 0, -1}},
	}
	caps := []conicCap{
		{at: Vector4{0, 0, 1, 1}, mask: axesXY}, {at: Vector4{0, 0, -1, 1}, mask: axesXY},
		{at: Vector4{0, 0, 0, -1}},
	}
	return newConicSolid("cylindrone", parts, caps, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewDicone builds the pyramid over a cone, which is also the hull of a disk
// and a segment: in local coordinates divided by Scale = (rx, ry, hz, h),
// with u = (w + 1)/2,
//
//	sqrt(x² + y²) <= u,  |z| <= 1 - u
//
// The disk is at w = +1 and the segment (along Z) at w = -1.
func NewDicone(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	parts := []conicPart{
		{mask: axesXY, k: 0.5, b: Vector4{0, 0, 0, 0.5}},
		{mask: axisZ, k: 0.5, b: Vector4{0, 0, 0, -0.5}},
	}
	caps := []conicCap{
		{at: Vector4{0, 0, 0, 1}, mask: axesXY},
		{at: Vector4{0, 0, 1, -1}}, {at: Vector4{0, 0, -1, -1}},
	}
	return newConicSolid("dicone", parts, caps, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewPrismicCylinder builds the product of a disk in the local XY plane and
// a regular polygon in the local ZW plane, both with the circumradius 1
// before Scale is applied: the duoprism with a circle in place of its first
// polygon. With 3 sides it is the triangular prismic cylinder.
func NewPrismicCylinder(sides int, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	if sides < 3 || sides > maxDuoprismSides {
		return nil, fmt.Errorf("prismic cylinder needs a polygon with 3..%d sides, got %d", maxDuoprismSides, sides)
	}
	parts := []conicPart{{mask: axesXY, k: 1}}
	var caps []conicCap
	for k := 0; k < sides; k++ {
		a := 2 * math.Pi * Real(k) / Real(sides)
		m := a + math.Pi/Real(sides)
		parts = append(parts, conicPart{k: math.Cos(math.Pi / Real(sides)), b: Vector4{0, 0, -math.Cos(m), -math.Sin(m)}})
		caps = append(caps, conicCap{at: Vector4{0, 0, math.Cos(a), math.Sin(a)}, mask: axesXY})
	}
	return newConicSolid(fmt.Sprintf("%d-gonal prismic cylinder", sides), parts, caps, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}
