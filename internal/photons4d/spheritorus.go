package photons4d

import (
	"fmt"
	"math"
)

const SpheritorusConservativeAABB = false

// Spheritorus (local space, core circle in the local XW plane):
//
//	( sqrt(u² + s²) - R )² + v² + z² <= r²
//
// where scaled local coordinates are:
//
//	u = x / sx
//	v = y / sy
//	z = z / sz
//	s = w / sw
//
// with:
//
//	Scale       = (sx, sy, sz, sw)
//	MajorRadius = R
//	MinorRadius = r
//
// This is the natural 4D analog of a solid torus with topology S¹ × B³:
// a 3-ball tube around a circle in the XW plane.
//
// The local object is then rotated by R and translated to Center.
type Spheritorus struct {
	Center      Point4
	Scale       Vector4 // per-axis scale applied before torus evaluation
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

	// material caches
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func NewSpheritorus(center Point4, scale Vector4, majorRadius, minorRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Spheritorus, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("spheritorus scale must be >0 on all axes, got %+v", scale)
	}
	if !(majorRadius > 0) {
		return nil, fmt.Errorf("spheritorus majorRadius must be > 0, got %.6g", majorRadius)
	}
	if !(minorRadius > 0) {
		return nil, fmt.Errorf("spheritorus minorRadius must be > 0, got %.6g", minorRadius)
	}
	if !(majorRadius > minorRadius) {
		return nil, fmt.Errorf("spheritorus requires majorRadius > minorRadius for a proper ring, got major=%.6g minor=%.6g", majorRadius, minorRadius)
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
			return nil, fmt.Errorf("reflect/refract/diffuse must be in [0,1]; channel %s got reflect=%.6g refract=%.6g diffuse=%.6g", c.n, c.r, c.t, c.d)
		}
		if c.r+c.t+c.d > 1+1e-12 {
			return nil, fmt.Errorf("per-channel reflect+refract+diffuse must be ≤1; channel %s got %.6g", c.n, c.r+c.t+c.d)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	R := rotFromAngles(angles)
	t := &Spheritorus{
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

	if SpheritorusConservativeAABB {
		// Conservative analytic AABB from the rotated bounding box of the
		// axis-scaled canonical torus.
		hx := scale.X * (majorRadius + minorRadius)
		hy := scale.Y * minorRadius
		hz := scale.Z * minorRadius
		hw := scale.W * (majorRadius + minorRadius)

		abs := func(x Real) Real {
			if x < 0 {
				return -x
			}
			return x
		}
		extent := func(row int) (minV, maxV Real) {
			off := abs(R.M[row][0])*hx + abs(R.M[row][1])*hy + abs(R.M[row][2])*hz + abs(R.M[row][3])*hw
			var c Real
			switch row {
			case 0:
				c = center.X
			case 1:
				c = center.Y
			case 2:
				c = center.Z
			default:
				c = center.W
			}
			return c - off, c + off
		}
		minX, maxX := extent(0)
		minY, maxY := extent(1)
		minZ, maxZ := extent(2)
		minW, maxW := extent(3)
		t.AABBMin = Point4{minX, minY, minZ, minW}
		t.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	} else {
		// Exact support-based AABB.
		// In scaled local space this shape is:
		//   circle(radius=R in XW) Minkowski-summed with a 4D ball(radius=r).
		// So support along local direction a is:
		//   R * sqrt((sx*ax)^2 + (sw*aw)^2)
		// + r * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2 + (sw*aw)^2)
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			Ax := scale.X * ax
			Ay := scale.Y * ay
			Az := scale.Z * az
			Aw := scale.W * aw
			major := majorRadius * math.Sqrt(Ax*Ax+Aw*Aw)
			minor := minorRadius * math.Sqrt(Ax*Ax+Ay*Ay+Az*Az+Aw*Aw)
			off := major + minor
			var c Real
			switch row {
			case 0:
				c = center.X
			case 1:
				c = center.Y
			case 2:
				c = center.Z
			default:
				c = center.W
			}
			return c - off, c + off
		}
		minX, maxX := axisExtent(0)
		minY, maxY := axisExtent(1)
		minZ, maxZ := axisExtent(2)
		minW, maxW := axisExtent(3)
		t.AABBMin = Point4{minX, minY, minZ, minW}
		t.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	}

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

	DebugLog("Created spheritorus: %+v", t)
	return t, nil
}

func boxIntervalLocal(Ol, Dl Vector4, hx, hy, hz, hw Real) (tInterval, bool) {
	const eps = 1e-12
	tmin, tmax := Real(math.Inf(-1)), Real(math.Inf(1))

	update := func(o, d, h Real) bool {
		if math.Abs(float64(d)) < eps {
			return !(o < -h || o > h)
		}
		t0 := (-h - o) / d
		t1 := (h - o) / d
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		if t0 > tmin {
			tmin = t0
		}
		if t1 < tmax {
			tmax = t1
		}
		return tmin <= tmax
	}

	if !update(Ol.X, Dl.X, hx) {
		return tInterval{}, false
	}
	if !update(Ol.Y, Dl.Y, hy) {
		return tInterval{}, false
	}
	if !update(Ol.Z, Dl.Z, hz) {
		return tInterval{}, false
	}
	if !update(Ol.W, Dl.W, hw) {
		return tInterval{}, false
	}
	return tInterval{tmin, tmax}, true
}

func intersectRaySpheritorus(O Point4, D Vector4, t *Spheritorus) (hit objectHit, ok bool) {
	Op := Vector4{O.X - t.Center.X, O.Y - t.Center.Y, O.Z - t.Center.Z, O.W - t.Center.W}
	Ol := t.RT.MulVec(Op)
	Dl := t.RT.MulVec(D)

	R := t.MajorRadius
	r := t.MinorRadius
	sx, sy, sz, sw := t.Scale.X, t.Scale.Y, t.Scale.Z, t.Scale.W

	box, ok := boxIntervalLocal(Ol, Dl, sx*(R+r), sy*r, sz*r, sw*(R+r))
	if !ok {
		return objectHit{}, false
	}

	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	// In the scaled space the solid is the set of points within r of the core
	// circle: the distance to that circle minus r is an exact signed distance.
	dist := func(tt Real) rayDist {
		q := A.Add(B.Mul(tt))
		rho := math.Hypot(q.X, q.W)
		p := rho - R
		delta := math.Sqrt(p*p + q.Y*q.Y + q.Z*q.Z)
		rd := rayDist{d: delta - r}
		if rho > 0 && delta > 0 {
			rd.rate = (p/rho*(q.X*B.X+q.W*B.W) + q.Y*B.Y + q.Z*B.Z) / delta
			rd.smooth = math.Min(delta, rho)
		}
		return rd
	}

	tHit, inv, ok := traceSolid(dist, B.Len(), box.t0, box.t1)
	if !ok {
		return objectHit{}, false
	}

	q := A.Add(B.Mul(tHit))
	k := Real(0)
	if rho := math.Hypot(q.X, q.W); rho > 0 {
		k = (rho - R) / rho
	}
	gradL := Vector4{k * q.X / sx, q.Y / sy, q.Z / sz, k * q.W / sw}
	if gradL.Len() == 0 {
		return objectHit{}, false
	}
	Nw := t.R.MulVec(gradL).Norm()
	return objectHit{t: tHit, Nw: Nw, mat: t, inv: inv}, true
}

func (t *Spheritorus) PAbsCh(i int) Real   { return t.pAbs[i] }
func (t *Spheritorus) DiffCh(i int) Real   { return t.diff[i] }
func (t *Spheritorus) ColorCh(i int) Real  { return t.colorArr[i] }
func (t *Spheritorus) F0Ch(i int) Real     { return t.f0[i] }
func (t *Spheritorus) ReflCh(i int) Real   { return t.refl[i] }
func (t *Spheritorus) RefrCh(i int) Real   { return t.refr[i] }
func (t *Spheritorus) IORCh(i int) Real    { return t.iorArr[i] }
func (t *Spheritorus) IORInvCh(i int) Real { return t.iorInv[i] }
