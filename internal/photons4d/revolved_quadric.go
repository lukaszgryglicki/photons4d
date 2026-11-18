package photons4d

import (
	"fmt"
	"math"
)

// RevolvedQuadric is a solid of revolution around the local W axis whose
// radius in XYZ follows a quadric: in local coordinates divided by Scale
//
//	x² + y² + z² <= K0 + K1·w + K2·w²,  Lower <= w <= Upper
//
// It serves the paraboloid (K0 = K1 = 1/2: the tip at w = -1 and the radius
// 1 at w = +1), the hyperboloid of one sheet (K0 = waist², K2 = 1 - waist²:
// an hourglass with the radius 1 at w = ±1, two cones tip to tip when the
// waist is 0) and the part of the unit ball between two parallel
// hyperplanes (K0 = 1, K2 = -1). The local object is then rotated by R and
// translated to Center.
type RevolvedQuadric struct {
	Center       Point4
	Scale        Vector4
	K0, K1, K2   Real
	Lower, Upper Real
	R            Mat4
	RT           Mat4
	Name         string

	Color   RGB
	Diffuse RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	AABBMin Point4
	AABBMax Point4

	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

// NewParaboloid builds the cap of a paraboloid of revolution: the tip at
// w = -1 and the base ball of the radius 1 at w = +1 before Scale.
func NewParaboloid(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*RevolvedQuadric, error) {
	return newRevolvedQuadric("paraboloid", 0.5, 0.5, 0, -1, 1, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewHyperboloid builds the hyperboloid of one sheet between w = -1 and
// w = +1, where its radius is 1; waist is its radius in the middle.
func NewHyperboloid(center Point4, scale Vector4, waist Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*RevolvedQuadric, error) {
	if !(waist >= 0 && waist <= 1) {
		return nil, fmt.Errorf("hyperboloid waist must be in [0,1], got %.6g", waist)
	}
	return newRevolvedQuadric("hyperboloid", waist*waist, 0, 1-waist*waist, -1, 1, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewSphericalSegment builds the part of the unit ball (an ellipsoid after
// Scale) between the hyperplanes w = lower and w = upper: a half of it, a
// cap or a zone.
func NewSphericalSegment(center Point4, scale Vector4, lower, upper Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*RevolvedQuadric, error) {
	if !(lower >= -1 && upper <= 1) {
		return nil, fmt.Errorf("spherical segment must lie within the ball, -1 <= lower < upper <= 1, got lower=%.6g upper=%.6g", lower, upper)
	}
	return newRevolvedQuadric("spherical segment", 1, 0, -1, lower, upper, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

func newRevolvedQuadric(name string, k0, k1, k2, lower, upper Real, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*RevolvedQuadric, error) {
	if err := checkSolid(name, scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	if !(lower < upper) || !isFinite(lower+upper+k0+k1+k2) {
		return nil, fmt.Errorf("%s needs lower < upper, got lower=%.6g upper=%.6g", name, lower, upper)
	}
	R := rotFromAngles(angles)
	s := &RevolvedQuadric{
		Center:  center,
		Scale:   scale,
		K0:      k0,
		K1:      k1,
		K2:      k2,
		Lower:   lower,
		Upper:   upper,
		R:       R,
		RT:      R.Transpose(),
		Name:    name,
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}
	// the squared radius must not be negative between the two ends
	least := math.Min(s.radius2(lower), s.radius2(upper))
	if v := -k1 / (2 * k2); k2 > 0 && v > lower && v < upper {
		least = math.Min(least, s.radius2(v))
	}
	if least < -1e-12 || !(math.Max(least, s.radius2(0.5*(lower+upper))) > 0 || s.radius2(lower) > 0 || s.radius2(upper) > 0) {
		return nil, fmt.Errorf("%s has no radius somewhere between w=%.6g and w=%.6g", name, lower, upper)
	}

	// A world axis sees the local point q as c + a·q with a = Scale times the row of R.
	extent := func(row int, c Real) (Real, Real) {
		a := Vector4{scale.X * R.M[row][0], scale.Y * R.M[row][1], scale.Z * R.M[row][2], scale.W * R.M[row][3]}
		return c - s.support(a.Mul(-1)), c + s.support(a)
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

func (s *RevolvedQuadric) radius2(w Real) Real {
	return s.K0 + w*(s.K1+w*s.K2)
}

// support returns the largest a·q over the solid: the best of
// a.W·w + |a.XYZ|·radius(w) over the ends and over the heights where it is
// level. The root of a quadric is convex or concave all the way, so there is
// at most one such height that counts.
func (s *RevolvedQuadric) support(a Vector4) Real {
	lat := math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z)
	at := func(w Real) Real {
		return a.W*w + lat*math.Sqrt(math.Max(0, s.radius2(w)))
	}
	best := math.Max(at(s.Lower), at(s.Upper))
	// a.W + lat·(K1 + 2·K2·w) / (2·radius) = 0, squared
	qa := 4 * s.K2 * (a.W*a.W - lat*lat*s.K2)
	qb := 4 * s.K1 * (a.W*a.W - lat*lat*s.K2)
	qc := 4*a.W*a.W*s.K0 - lat*lat*s.K1*s.K1
	r0, r1, n := solveQuadratic(qa, qb, qc)
	for i, w := range [2]Real{r0, r1} {
		if i < n && w > s.Lower && w < s.Upper {
			best = math.Max(best, at(w))
		}
	}
	return best
}

// solveQuadratic returns the real roots of a·t² + b·t + c in rising order
// and their number.
func solveQuadratic(a, b, c Real) (r0, r1 Real, n int) {
	if a == 0 {
		if b == 0 {
			return 0, 0, 0
		}
		return -c / b, 0, 1
	}
	disc := b*b - 4*a*c
	if disc < 0 {
		return 0, 0, 0
	}
	q := -0.5 * (b + math.Copysign(math.Sqrt(disc), b))
	if q == 0 {
		return 0, 0, 2
	}
	r0, r1 = q/a, c/q
	if r0 > r1 {
		r0, r1 = r1, r0
	}
	return r0, r1, 2
}

func intersectRayRevolvedQuadric(O Point4, D Vector4, s *RevolvedQuadric) (objectHit, bool) {
	const eps = 1e-12
	Ol := s.RT.MulVec(Vector4{O.X - s.Center.X, O.Y - s.Center.Y, O.Z - s.Center.Z, O.W - s.Center.W})
	Dl := s.RT.MulVec(D)
	sx, sy, sz, sw := s.Scale.X, s.Scale.Y, s.Scale.Z, s.Scale.W
	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	// the events of the ray: the two ends and the roots of
	// F(t) = x² + y² + z² - K0 - K1·w - K2·w²
	var events [4]Real
	n := 0
	tLower, tUpper := math.Inf(-1), math.Inf(1)
	if math.Abs(B.W) < eps {
		if A.W < s.Lower || A.W > s.Upper {
			return objectHit{}, false
		}
	} else {
		tLower, tUpper = (s.Lower-A.W)/B.W, (s.Upper-A.W)/B.W
		if math.Max(tLower, tUpper) <= eps {
			return objectHit{}, false
		}
		events[0], events[1] = tLower, tUpper
		n = 2
	}
	qa := B.X*B.X + B.Y*B.Y + B.Z*B.Z - s.K2*B.W*B.W
	qb := 2*(A.X*B.X+A.Y*B.Y+A.Z*B.Z) - s.K1*B.W - 2*s.K2*A.W*B.W
	qc := A.X*A.X + A.Y*A.Y + A.Z*A.Z - s.radius2(A.W)
	r0, r1, roots := solveQuadratic(qa, qb, qc)
	if roots > 0 {
		events[n] = r0
		n++
	}
	if roots > 1 {
		events[n] = r1
		n++
	}
	for i := 1; i < n; i++ {
		for j := i; j > 0 && events[j-1] > events[j]; j-- {
			events[j-1], events[j] = events[j], events[j-1]
		}
	}

	within := func(t Real) bool {
		w := A.W + t*B.W
		return w >= s.Lower && w <= s.Upper && qc+t*(qb+t*qa) <= 0
	}
	// the side of the ray between two events is taken in the middle
	from := eps
	i := 0
	for i < n && events[i] <= eps {
		i++
	}
	if i == n {
		return objectHit{}, false
	}
	inside := within(0.5 * (from + events[i]))
	for ; i < n; i++ {
		t := events[i]
		next := t + 1
		if i+1 < n {
			next = events[i+1]
		}
		if within(0.5*(t+next)) == inside {
			continue
		}
		var Nl Vector4
		switch {
		case t == tLower:
			Nl = Vector4{0, 0, 0, -1 / sw}
		case t == tUpper:
			Nl = Vector4{0, 0, 0, 1 / sw}
		default:
			q := A.Add(B.Mul(t))
			Nl = Vector4{2 * q.X / sx, 2 * q.Y / sy, 2 * q.Z / sz, -(s.K1 + 2*s.K2*q.W) / sw}
		}
		if Nl.Len() == 0 {
			return objectHit{}, false
		}
		return objectHit{t: t, Nw: s.R.MulVec(Nl).Norm(), mat: s, inv: inside}, true
	}
	return objectHit{}, false
}

func (s *RevolvedQuadric) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *RevolvedQuadric) DiffCh(i int) Real   { return s.diff[i] }
func (s *RevolvedQuadric) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *RevolvedQuadric) F0Ch(i int) Real     { return s.f0[i] }
func (s *RevolvedQuadric) ReflCh(i int) Real   { return s.refl[i] }
func (s *RevolvedQuadric) RefrCh(i int) Real   { return s.refr[i] }
func (s *RevolvedQuadric) IORCh(i int) Real    { return s.iorArr[i] }
func (s *RevolvedQuadric) IORInvCh(i int) Real { return s.iorInv[i] }
