package photons4d

import (
	"fmt"
	"math"
)

// Torinder (local space): the prism over a torus,
//
//	(sqrt(u² + v²) - R)² + z² <= r²,  |s| <= 1
//
// where the scaled local coordinates are
//
//	u = x / sx,  v = y / sy,  z = z / sz,  s = w / sw
//
// with Scale = (sx, sy, sz, sw), MajorRadius = R and MinorRadius = r. The
// torus lies in the XYZ subspace with its ring in the XY plane, and Scale.W
// is half of the length of the prism. The local object is then rotated by R
// and translated to Center.
type Torinder struct {
	Center      Point4
	Scale       Vector4
	MajorRadius Real
	MinorRadius Real
	R           Mat4
	RT          Mat4

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

func NewTorinder(center Point4, scale Vector4, majorRadius, minorRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Torinder, error) {
	if err := checkSolid("torinder", scale, diffuse, reflectivity, refractivity, ior); err != nil {
		return nil, err
	}
	if !(minorRadius > 0) {
		return nil, fmt.Errorf("torinder minorRadius must be > 0, got %.6g", minorRadius)
	}
	if !(majorRadius > minorRadius) {
		return nil, fmt.Errorf("torinder requires majorRadius > minorRadius for a proper ring, got major=%.6g minor=%.6g", majorRadius, minorRadius)
	}
	R := rotFromAngles(angles)
	t := &Torinder{
		Center:      center,
		Scale:       scale,
		MajorRadius: majorRadius,
		MinorRadius: minorRadius,
		R:           R,
		RT:          R.Transpose(),
		Color:       color,
		Diffuse:     diffuse,
		Reflect:     reflectivity,
		Refract:     refractivity,
		IOR:         ior,
	}

	// Exact support-based AABB: in scaled local space the solid is the sum
	// of a circle (radius R in XY), a 3D ball (radius r in XYZ) and the
	// segment |w| <= 1.
	axisExtent := func(row int, c Real) (Real, Real) {
		ax, ay, az, aw := scale.X*R.M[row][0], scale.Y*R.M[row][1], scale.Z*R.M[row][2], scale.W*R.M[row][3]
		off := majorRadius*math.Hypot(ax, ay) + minorRadius*math.Sqrt(ax*ax+ay*ay+az*az) + math.Abs(aw)
		return c - off, c + off
	}
	minX, maxX := axisExtent(0, center.X)
	minY, maxY := axisExtent(1, center.Y)
	minZ, maxZ := axisExtent(2, center.Z)
	minW, maxW := axisExtent(3, center.W)
	t.AABBMin = Point4{minX, minY, minZ, minW}
	t.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	t.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	t.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	t.colorArr = [3]Real{color.R, color.G, color.B}
	t.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	t.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - t.refl[i] - t.refr[i] - t.diff[i]
		if p < 0 {
			p = 0
		}
		t.pAbs[i] = p
		t.iorInv[i] = 1 / t.iorArr[i]
		n := t.iorArr[i]
		r0 := (n - 1) / (n + 1)
		t.f0[i] = r0 * r0
	}

	DebugLog("Created torinder: %+v", t)
	return t, nil
}

func intersectRayTorinder(O Point4, D Vector4, tor *Torinder) (objectHit, bool) {
	const eps = 1e-12
	Ol := tor.RT.MulVec(Vector4{O.X - tor.Center.X, O.Y - tor.Center.Y, O.Z - tor.Center.Z, O.W - tor.Center.W})
	Dl := tor.RT.MulVec(D)
	sx, sy, sz, sw := tor.Scale.X, tor.Scale.Y, tor.Scale.Z, tor.Scale.W
	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}
	R, r := tor.MajorRadius, tor.MinorRadius

	// the part of the ray between the two ends of the prism
	ta, tb := math.Inf(-1), math.Inf(1)
	if math.Abs(B.W) < eps {
		if math.Abs(A.W) > 1 {
			return objectHit{}, false
		}
	} else {
		ta, tb = (-1-A.W)/B.W, (1-A.W)/B.W
		if ta > tb {
			ta, tb = tb, ta
		}
	}
	// and inside a box around the torus, with a margin
	t0, t1 := ta, tb
	for _, c := range [3][3]Real{{A.X, B.X, R + r}, {A.Y, B.Y, R + r}, {A.Z, B.Z, r}} {
		h := c[2] * (1 + 1e-6)
		if math.Abs(c[1]) < eps {
			if math.Abs(c[0]) > h {
				return objectHit{}, false
			}
			continue
		}
		u, v := (-h-c[0])/c[1], (h-c[0])/c[1]
		if u > v {
			u, v = v, u
		}
		t0, t1 = math.Max(t0, u), math.Min(t1, v)
	}
	if t0 < eps {
		t0 = eps
	}
	if !(t1 > t0) {
		return objectHit{}, false
	}

	// In XYZ the distance to the ring minus r is an exact signed distance
	// of the torus; the speed of the ray there is that of its XYZ part.
	speed := math.Sqrt(B.X*B.X + B.Y*B.Y + B.Z*B.Z)
	dist := func(tt Real) rayDist {
		x, y, z := A.X+tt*B.X, A.Y+tt*B.Y, A.Z+tt*B.Z
		rho := math.Hypot(x, y)
		p := rho - R
		delta := math.Hypot(p, z)
		rd := rayDist{d: delta - r}
		if rho > 0 && delta > 0 {
			rd.rate = (p/rho*(x*B.X+y*B.Y) + z*B.Z) / delta
			rd.smooth = math.Min(delta, rho)
		}
		return rd
	}
	end := func(t Real, inv bool) (objectHit, bool) {
		Nl := Vector4{0, 0, 0, math.Copysign(1, A.W+t*B.W) / sw}
		return objectHit{t: t, Nw: tor.R.MulVec(Nl).Norm(), mat: tor, inv: inv}, true
	}

	if ta > eps && t0 == ta && dist(ta).d <= 0 {
		// enters through an end of the prism
		return end(ta, false)
	}
	var tHit Real
	var inside, ok bool
	if speed > 0 {
		tHit, inside, ok = traceFirstCrossing(dist, speed, t0, t1)
	} else {
		inside = dist(t0).d <= 0
	}
	if !ok {
		if inside && t1 == tb {
			// leaves through an end of the prism
			return end(tb, true)
		}
		return objectHit{}, false
	}

	x, y, z := A.X+tHit*B.X, A.Y+tHit*B.Y, A.Z+tHit*B.Z
	rho := math.Hypot(x, y)
	k := Real(0)
	if rho > 0 {
		k = (rho - R) / rho
	}
	Nl := Vector4{k * x / sx, k * y / sy, z / sz, 0}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	return objectHit{t: tHit, Nw: tor.R.MulVec(Nl).Norm(), mat: tor, inv: inside}, true
}

func (t *Torinder) PAbsCh(i int) Real   { return t.pAbs[i] }
func (t *Torinder) DiffCh(i int) Real   { return t.diff[i] }
func (t *Torinder) ColorCh(i int) Real  { return t.colorArr[i] }
func (t *Torinder) F0Ch(i int) Real     { return t.f0[i] }
func (t *Torinder) ReflCh(i int) Real   { return t.refl[i] }
func (t *Torinder) RefrCh(i int) Real   { return t.refr[i] }
func (t *Torinder) IORCh(i int) Real    { return t.iorArr[i] }
func (t *Torinder) IORInvCh(i int) Real { return t.iorInv[i] }
