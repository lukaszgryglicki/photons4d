package photons4d

import (
	"fmt"
	"math"
)

const maxKleinTwist = 16

// KleinSurface is a wall of even thickness around a surface that has no
// inside: a flat curve, the profile, goes once around the local XY plane and
// turns by Twist half-turns in the local ZW plane on its way. In local
// coordinates divided by Scale the surface is
//
//	(a cos u, a sin u, b cos(Twist·u/2), b sin(Twist·u/2))
//
// where (a, b) runs over the profile: the circle of the radius MinorRadius
// around (MajorRadius, 0) or, for a band, the segment a = MajorRadius,
// |b| <= MinorRadius. The solid holds the points within Thickness/2 of the
// surface.
//
// With one half-turn the circle sweeps a Klein bottle, or the real
// projective plane when the two radii are equal (all the circles then pass
// through the centre), and the segment sweeps a Möbius band. An even Twist
// gives a torus or a ring that is twisted in 4D, no twist at all the ditorus.
// The local object is then rotated by R and translated to Center.
type KleinSurface struct {
	Center      Point4
	Scale       Vector4
	MajorRadius Real
	MinorRadius Real
	Thickness   Real
	Twist       int
	Band        bool
	R           Mat4
	RT          Mat4
	Name        string

	Color   RGB
	Diffuse RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	AABBMin Point4
	AABBMax Point4

	tube  Real // half of the thickness
	half  Real // half of the twist: the turn in ZW per unit of u
	sweep Real // bound of the speed of a point of the profile as u grows
	du    Real // step of the search over u

	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

// NewKleinBottle builds the wall around a Klein bottle (an odd twist) or
// around a twisted torus (an even one).
func NewKleinBottle(center Point4, scale Vector4, majorRadius, minorRadius, thickness Real, twist int, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*KleinSurface, error) {
	if !(majorRadius > minorRadius) {
		return nil, fmt.Errorf("klein bottle requires majorRadius > minorRadius, got major=%.6g minor=%.6g", majorRadius, minorRadius)
	}
	name := "klein bottle"
	if twist%2 == 0 {
		name = "twisted torus"
	}
	return newKleinSurface(name, false, majorRadius, minorRadius, thickness, twist, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewProjectivePlane builds the wall around the real projective plane: the
// circles of the radius given, which all pass through the centre.
func NewProjectivePlane(center Point4, scale Vector4, radius, thickness Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*KleinSurface, error) {
	return newKleinSurface("projective plane", false, radius, radius, thickness, 1, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewMobiusBand builds the wall around a Möbius band (an odd twist) or
// around a twisted ring (an even one) of the width 2·halfWidth.
func NewMobiusBand(center Point4, scale Vector4, majorRadius, halfWidth, thickness Real, twist int, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*KleinSurface, error) {
	name := "mobius band"
	if twist%2 == 0 {
		name = "twisted ring"
	}
	return newKleinSurface(name, true, majorRadius, halfWidth, thickness, twist, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

func newKleinSurface(name string, band bool, major, minor, thickness Real, twist int, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*KleinSurface, error) {
	if err := checkSolid(name, scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	if !(major > 0 && minor > 0 && thickness > 0) || !isFinite(major+minor+thickness) {
		return nil, fmt.Errorf("%s radii and thickness must be > 0, got major=%.6g minor=%.6g thickness=%.6g", name, major, minor, thickness)
	}
	if twist < 0 || twist > maxKleinTwist {
		return nil, fmt.Errorf("%s twist must be 0..%d half-turns, got %d", name, maxKleinTwist, twist)
	}
	if !band && !(major > minor) && (major != minor || twist != 1) {
		return nil, fmt.Errorf("%s requires majorRadius > minorRadius, or equal radii with one half-turn, got major=%.6g minor=%.6g twist=%d", name, major, minor, twist)
	}
	R := rotFromAngles(angles)
	k := &KleinSurface{
		Center:      center,
		Scale:       scale,
		MajorRadius: major,
		MinorRadius: minor,
		Thickness:   thickness,
		Twist:       twist,
		Band:        band,
		R:           R,
		RT:          R.Transpose(),
		Name:        name,
		Color:       color,
		Diffuse:     diffuse,
		Reflect:     reflectivity,
		Refract:     refractivity,
		IOR:         ior,
		tube:        thickness / 2,
		half:        Real(twist) / 2,
	}
	far := major
	if !band {
		far += minor
	}
	k.sweep = math.Hypot(far, k.half*minor)
	// A quarter of the widest step that still found the nearest point of
	// every one of millions of points, for 0 to 16 half-turns.
	k.du = 0.3 / math.Max(1, k.half)

	// A world axis sees the local point q as c + a·q with a = Scale times the row of R.
	extent := func(row int, c Real) (Real, Real) {
		a := Vector4{scale.X * R.M[row][0], scale.Y * R.M[row][1], scale.Z * R.M[row][2], scale.W * R.M[row][3]}
		wall := k.tube * a.Len()
		return c - k.support(a.Mul(-1)) - wall, c + k.support(a) + wall
	}
	minX, maxX := extent(0, center.X)
	minY, maxY := extent(1, center.Y)
	minZ, maxZ := extent(2, center.Z)
	minW, maxW := extent(3, center.W)
	k.AABBMin = Point4{minX, minY, minZ, minW}
	k.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	k.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	k.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	k.colorArr = [3]Real{color.R, color.G, color.B}
	k.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	k.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - k.refl[i] - k.refr[i] - k.diff[i]
		if p < 0 {
			p = 0
		}
		k.pAbs[i] = p
		k.iorInv[i] = 1 / k.iorArr[i]
		n := k.iorArr[i]
		r0 := (n - 1) / (n + 1)
		k.f0[i] = r0 * r0
	}
	DebugLog("Created %s: %+v", name, k)
	return k, nil
}

// support returns a bound of the largest a·q over the surface that is at most
// a little above it. For a turn u the profile reaches
// MajorRadius·α + MinorRadius·hypot(α, β) (the circle) or
// MajorRadius·α + MinorRadius·|β| (the segment), where α and β are the parts
// of a along the two axes of the profile.
func (k *KleinSurface) support(a Vector4) Real {
	axy, azw := math.Hypot(a.X, a.Y), math.Hypot(a.Z, a.W)
	reach := func(alpha, beta Real) Real {
		if k.Band {
			return k.MajorRadius*alpha + k.MinorRadius*math.Abs(beta)
		}
		return k.MajorRadius*alpha + k.MinorRadius*math.Hypot(alpha, beta)
	}
	// both parts at their best at once
	bound := reach(axy, azw)
	// the best turn of a fine search, plus what the surface can gain between two turns
	n := 720 * (1 + k.Twist)
	h := 2 * math.Pi / Real(n)
	best := math.Inf(-1)
	for i := 0; i < n; i++ {
		u := Real(i) * h
		su, cu := math.Sincos(u)
		sh, ch := math.Sincos(k.half * u)
		best = math.Max(best, reach(a.X*cu+a.Y*su, a.Z*ch+a.W*sh))
	}
	return math.Min(bound, best+a.Len()*k.sweep*h/2)
}

// profile returns the point of the profile that is nearest to (x, z) in its
// plane, the unit tangent of the profile there (none at an end of the
// segment) and the second derivative of half of the squared distance along
// the profile.
func (k *KleinSurface) profile(x, z Real) (a, b, ta, tb, curv Real) {
	if k.Band {
		if z > k.MinorRadius {
			return k.MajorRadius, k.MinorRadius, 0, 0, 1
		}
		if z < -k.MinorRadius {
			return k.MajorRadius, -k.MinorRadius, 0, 0, 1
		}
		return k.MajorRadius, z, 0, 1, 1
	}
	dx := x - k.MajorRadius
	h := math.Hypot(dx, z)
	if h == 0 {
		return k.MajorRadius + k.MinorRadius, 0, 0, 0, 1
	}
	c, s := dx/h, z/h
	return k.MajorRadius + k.MinorRadius*c, k.MinorRadius * s, -s, c, h / k.MinorRadius
}

// sq returns the squared distance from q to the profile turned by u, given
// the sine and the cosine of u and of the twist angle.
func (k *KleinSurface) sq(q Vector4, su, cu, sh, ch Real) Real {
	x, y := q.X*cu+q.Y*su, q.Y*cu-q.X*su
	z, w := q.Z*ch+q.W*sh, q.W*ch-q.Z*sh
	dx := x - k.MajorRadius
	if k.Band {
		dz := math.Abs(z) - k.MinorRadius
		if dz < 0 {
			dz = 0
		}
		return y*y + w*w + dx*dx + dz*dz
	}
	e := math.Sqrt(dx*dx+z*z) - k.MinorRadius
	return y*y + w*w + e*e
}

// aside returns the distance in the plane of the profile from the lengths
// of q in XY and in ZW to the profile, which the distance from q to the
// surface is not below.
func (k *KleinSurface) aside(q Vector4) Real {
	dx, z := math.Hypot(q.X, q.Y)-k.MajorRadius, math.Hypot(q.Z, q.W)
	if k.Band {
		return math.Hypot(dx, math.Max(0, z-k.MinorRadius))
	}
	return math.Abs(math.Hypot(dx, z) - k.MinorRadius)
}

// kleinFoot is a point of the surface: its turn and its place on the profile.
type kleinFoot struct {
	u              Real
	su, cu, sh, ch Real
	a, b           Real
}

func (f kleinFoot) point() Vector4 {
	return Vector4{f.a * f.cu, f.a * f.su, f.b * f.ch, f.b * f.sh}
}

// at returns the squared distance from q to the profile turned by u, its
// first two derivatives over u and the nearest point of that profile.
func (k *KleinSurface) at(q Vector4, u Real) (g, g1, g2 Real, f kleinFoot) {
	c := k.half
	f.u = u
	f.su, f.cu = math.Sincos(u)
	f.sh, f.ch = math.Sincos(c * u)
	x, y := q.X*f.cu+q.Y*f.su, q.Y*f.cu-q.X*f.su
	z, w := q.Z*f.ch+q.W*f.sh, q.W*f.ch-q.Z*f.sh
	var ta, tb, curv Real
	f.a, f.b, ta, tb, curv = k.profile(x, z)
	g = y*y + w*w + (x-f.a)*(x-f.a) + (z-f.b)*(z-f.b)
	g1 = -2 * (y*f.a + c*w*f.b)
	m := y*ta + c*w*tb
	g2 = 2*(x*f.a+c*c*z*f.b) - 2*m*m/curv
	return
}

// descend goes down from the turn u to a least squared distance by Newton
// steps of at most lim, each of which must not go up.
func (k *KleinSurface) descend(q Vector4, u, lim Real) (Real, kleinFoot) {
	g, g1, g2, f := k.at(q, u)
	return k.descendFrom(q, u, lim, g, g1, g2, f)
}

// descendFrom is descend with the values at u at hand.
func (k *KleinSurface) descendFrom(q Vector4, u, lim, g, g1, g2 Real, f kleinFoot) (Real, kleinFoot) {
	for it := 0; it < 64; it++ {
		if g2 > 0 && math.Abs(g1) < 1e-7*g2 {
			// the last step: what is left after it is its square, below the rounding
			g, _, _, f = k.at(q, u-g1/g2)
			break
		}
		step := -lim
		if g2 > 0 && math.Abs(g1) < lim*g2 {
			step = -g1 / g2
		} else if g1 < 0 {
			step = lim
		} else if g1 == 0 && g2 >= 0 {
			// level and not on a top: nowhere to go
			break
		}
		moved := false
		for try := 0; try < 20 && !moved; try++ {
			if gn, g1n, g2n, fn := k.at(q, u+step); gn <= g {
				u, g, g1, g2, f = u+step, gn, g1n, g2n, fn
				moved = true
			} else {
				step /= 2
			}
		}
		if !moved {
			// the values no longer tell the two sides apart: trust the slope
			if g2 > 0 && math.Abs(g1) < 1e-6*g2 {
				g, _, _, f = k.at(q, u-g1/g2)
			}
			break
		}
	}
	return g, f
}

// nearest returns the distance from q to the surface and the point of the
// surface at that distance. With far set, a point that is clearly further
// from the surface than the wall gets only a lower bound of its distance,
// which is told by exact = false. last, if given, holds the turn of the
// nearest point of a point not far from q (NaN when there is none yet) and
// takes the turn found.
func (k *KleinSurface) nearest(q Vector4, far bool, last *Real) (dist Real, foot Vector4, exact bool) {
	const most = 1024
	rho := math.Hypot(q.X, q.Y)
	mid := math.Atan2(q.Y, q.X)

	// A bound of the distance from above: the dip next to the turn that was
	// the nearest a moment ago, or else the profile in the half-plane of q.
	var g, g1, g2 Real
	var f kleinFoot
	known := last != nil && !math.IsNaN(*last)
	if known {
		g, f = k.descend(q, *last, k.du)
	} else {
		g, g1, g2, f = k.at(q, mid)
	}
	// A profile turned away from q by the angle at which the line of its
	// half-plane is as far from q as that bound, or by more, cannot be nearer.
	width := Real(math.Pi)
	if bound := math.Sqrt(g); rho > bound {
		width = math.Min(width, 1.01*math.Asin(bound/rho)+1e-9)
	}
	if width > k.du && k.Twist == 1 {
		// With one half-turn there is one turn at which the plane of the
		// profile holds the part of q in ZW, and a plane turned away from it
		// leaves q at the distance of that part times the sine of half the
		// turn: a second bound of the same kind, around another turn.
		sigma := math.Hypot(q.Z, q.W)
		turn := 2 * math.Atan2(q.W, q.Z)
		if !known {
			if ga, g1a, g2a, fa := k.at(q, turn); ga < g {
				g, g1, g2, f = ga, g1a, g2a, fa
			}
		}
		bound := math.Sqrt(g)
		if rho > bound {
			width = math.Min(width, 1.01*math.Asin(bound/rho)+1e-9)
		}
		if sigma > bound {
			if wide := 2.02*math.Asin(bound/sigma) + 1e-9; wide < width {
				mid, width = turn, wide
			}
		}
	}
	if width <= k.du {
		// too narrow for more than one dip, which the bound belongs to
		if !known {
			g, f = k.descendFrom(q, f.u, width, g, g1, g2, f)
		}
		if last != nil {
			*last = f.u
		}
		return math.Sqrt(g), f.point(), true
	}
	n := 2 * int(math.Ceil(width/k.du))
	if n < 4 {
		n = 4
	}
	if n > most {
		n = most
	}
	h := 2 * width / Real(n)

	var gs [most + 1]Real
	su, cu := math.Sincos(mid - width)
	sh, ch := math.Sincos(k.half * (mid - width))
	sd, cd := math.Sincos(h)
	shd, chd := math.Sincos(k.half * h)
	low := math.Inf(1)
	for i := 0; i <= n; i++ {
		gs[i] = k.sq(q, su, cu, sh, ch)
		low = math.Min(low, gs[i])
		su, cu = su*cd+cu*sd, cu*cd-su*sd
		sh, ch = sh*chd+ch*shd, ch*chd-sh*shd
	}
	low = math.Sqrt(low)
	// between two turns of the search the distance changes by at most this
	slack := k.sweep * h / 2
	if far && low-k.tube >= 2*slack {
		return low - slack, Vector4{}, false
	}

	// the turns of the search that are not above their neighbours, the lowest first
	var cand [16]int
	nc := 0
	for i := 0; i <= n; i++ {
		if (i > 0 && gs[i] > gs[i-1]) || (i < n && gs[i] >= gs[i+1]) {
			continue
		}
		j := nc
		if nc < len(cand) {
			nc++
		} else if j--; gs[i] >= gs[cand[j]] {
			continue
		}
		for ; j > 0 && gs[cand[j-1]] > gs[i]; j-- {
			cand[j] = cand[j-1]
		}
		cand[j] = i
	}
	dist = math.Inf(1)
	if known {
		dist = math.Sqrt(g)
	}
	for _, i := range cand[:nc] {
		// nothing next to this turn is lower than its distance less two slacks
		if math.Sqrt(gs[i])-2*slack >= dist {
			break
		}
		gi, fi := k.descend(q, mid-width+Real(i)*h, h)
		if d := math.Sqrt(gi); d < dist {
			dist, f = d, fi
		}
	}
	if last != nil {
		*last = f.u
	}
	return dist, f.point(), true
}

func intersectRayKleinSurface(O Point4, D Vector4, k *KleinSurface) (objectHit, bool) {
	Ol := k.RT.MulVec(Vector4{O.X - k.Center.X, O.Y - k.Center.Y, O.Z - k.Center.Z, O.W - k.Center.W})
	Dl := k.RT.MulVec(D)
	A := Vector4{Ol.X / k.Scale.X, Ol.Y / k.Scale.Y, Ol.Z / k.Scale.Z, Ol.W / k.Scale.W}
	B := Vector4{Dl.X / k.Scale.X, Dl.Y / k.Scale.Y, Dl.Z / k.Scale.Z, Dl.W / k.Scale.W}

	// the part of the ray inside a box around the solid
	hxy, hzw := k.MajorRadius+k.tube, k.MinorRadius+k.tube
	if !k.Band {
		hxy += k.MinorRadius
	}
	iv, ok := boxIntervalLocal(A, B, hxy, hxy, hzw, hzw)
	if !ok || !(iv.t1 > 0) {
		return objectHit{}, false
	}

	// The distance to the surface less the half of the thickness is an exact
	// signed distance of the solid. Inside, the ball of that half around the
	// nearest point of the surface belongs to the solid, which bounds the
	// second derivative along the ray by speed²/distance. Far outside, the
	// lengths of the point in XY and in ZW do: taken as a point of the plane
	// of the profile, they are no further from it than the point is from the
	// surface.
	speed := B.Len()
	last := math.NaN()
	dist := func(t Real) rayDist {
		q := A.Add(B.Mul(t))
		if far := k.aside(q) - k.tube; far > k.tube {
			return rayDist{d: far}
		}
		d, foot, exact := k.nearest(q, true, &last)
		rd := rayDist{d: d - k.tube}
		if exact && rd.d < 0 && d > 0 {
			rd.rate = q.Sub(foot).Dot(B) / d
			rd.smooth = d
		}
		return rd
	}
	t, inside, ok := traceSolid(dist, speed, iv.t0, iv.t1)
	if !ok {
		return objectHit{}, false
	}
	q := A.Add(B.Mul(t))
	_, foot, _ := k.nearest(q, false, &last)
	n := q.Sub(foot)
	Nl := Vector4{n.X / k.Scale.X, n.Y / k.Scale.Y, n.Z / k.Scale.Z, n.W / k.Scale.W}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	return objectHit{t: t, Nw: k.R.MulVec(Nl).Norm(), mat: k, inv: inside}, true
}

func (k *KleinSurface) PAbsCh(i int) Real   { return k.pAbs[i] }
func (k *KleinSurface) DiffCh(i int) Real   { return k.diff[i] }
func (k *KleinSurface) ColorCh(i int) Real  { return k.colorArr[i] }
func (k *KleinSurface) F0Ch(i int) Real     { return k.f0[i] }
func (k *KleinSurface) ReflCh(i int) Real   { return k.refl[i] }
func (k *KleinSurface) RefrCh(i int) Real   { return k.refr[i] }
func (k *KleinSurface) IORCh(i int) Real    { return k.iorArr[i] }
func (k *KleinSurface) IORInvCh(i int) Real { return k.iorInv[i] }
