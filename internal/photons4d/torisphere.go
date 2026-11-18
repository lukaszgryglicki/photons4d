package photons4d

import (
	"fmt"
	"math"
)

const TorisphereConservativeAABB = false

// Torisphere / Sphere-tube (local space, core 2-sphere in the local XYZ subspace):
//
//	( sqrt(u² + v² + z²) - R )² + s² <= r²
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
// This is a tube around a 2-sphere of radius R embedded in the XYZ subspace.
// Topologically the solid is S² × D².
//
// The local object is then rotated by R and translated to Center.
type Torisphere struct {
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

func NewTorisphere(center Point4, scale Vector4, majorRadius, minorRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Torisphere, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("torisphere scale must be >0 on all axes, got %+v", scale)
	}
	if !(majorRadius > 0) {
		return nil, fmt.Errorf("torisphere majorRadius must be > 0, got %.6g", majorRadius)
	}
	if !(minorRadius > 0) {
		return nil, fmt.Errorf("torisphere minorRadius must be > 0, got %.6g", minorRadius)
	}
	if !(majorRadius > minorRadius) {
		return nil, fmt.Errorf("torisphere requires majorRadius > minorRadius for a proper ring, got major=%.6g minor=%.6g", majorRadius, minorRadius)
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
	t := &Torisphere{
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

	if TorisphereConservativeAABB {
		// Conservative rotated bounding box.
		hx := scale.X * (majorRadius + minorRadius)
		hy := scale.Y * (majorRadius + minorRadius)
		hz := scale.Z * (majorRadius + minorRadius)
		hw := scale.W * minorRadius

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
		// In scaled local space:
		//   core support = R * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2)
		//   tube support = r * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2 + (sw*aw)^2)
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			Ax := scale.X * ax
			Ay := scale.Y * ay
			Az := scale.Z * az
			Aw := scale.W * aw
			core := math.Sqrt(Ax*Ax + Ay*Ay + Az*Az)
			off := majorRadius*core + minorRadius*math.Sqrt(Ax*Ax+Ay*Ay+Az*Az+Aw*Aw)
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

	DebugLog("Created torisphere: %+v", t)
	return t, nil
}

func intersectRayTorisphere(O Point4, D Vector4, tor *Torisphere) (hit objectHit, ok bool) {
	Op := Vector4{O.X - tor.Center.X, O.Y - tor.Center.Y, O.Z - tor.Center.Z, O.W - tor.Center.W}
	Ol := tor.RT.MulVec(Op)
	Dl := tor.RT.MulVec(D)

	R := tor.MajorRadius
	r := tor.MinorRadius
	sx, sy, sz, sw := tor.Scale.X, tor.Scale.Y, tor.Scale.Z, tor.Scale.W

	box, ok := boxIntervalLocal(Ol, Dl, sx*(R+r), sy*(R+r), sz*(R+r), sw*r)
	if !ok {
		return objectHit{}, false
	}

	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	// In the scaled space the solid is the set of points within r of the core
	// 2-sphere: the distance to that sphere minus r is an exact signed distance.
	dist := func(tt Real) rayDist {
		q := A.Add(B.Mul(tt))
		a := math.Sqrt(q.X*q.X + q.Y*q.Y + q.Z*q.Z)
		p := a - R
		delta := math.Hypot(p, q.W)
		rd := rayDist{d: delta - r}
		if a > 0 && delta > 0 {
			rd.rate = (p/a*(q.X*B.X+q.Y*B.Y+q.Z*B.Z) + q.W*B.W) / delta
			rd.smooth = math.Min(delta, a)
		}
		return rd
	}

	tHit, inv, ok := traceSolid(dist, B.Len(), box.t0, box.t1)
	if !ok {
		return objectHit{}, false
	}

	q := A.Add(B.Mul(tHit))
	k := Real(0)
	if a := math.Sqrt(q.X*q.X + q.Y*q.Y + q.Z*q.Z); a > 0 {
		k = (a - R) / a
	}
	Nl := Vector4{k * q.X / sx, k * q.Y / sy, k * q.Z / sz, q.W / sw}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	Nw := tor.R.MulVec(Nl).Norm()
	return objectHit{t: tHit, Nw: Nw, mat: tor, inv: inv}, true
}

func (t *Torisphere) PAbsCh(i int) Real   { return t.pAbs[i] }
func (t *Torisphere) DiffCh(i int) Real   { return t.diff[i] }
func (t *Torisphere) ColorCh(i int) Real  { return t.colorArr[i] }
func (t *Torisphere) F0Ch(i int) Real     { return t.f0[i] }
func (t *Torisphere) ReflCh(i int) Real   { return t.refl[i] }
func (t *Torisphere) RefrCh(i int) Real   { return t.refr[i] }
func (t *Torisphere) IORCh(i int) Real    { return t.iorArr[i] }
func (t *Torisphere) IORInvCh(i int) Real { return t.iorInv[i] }
