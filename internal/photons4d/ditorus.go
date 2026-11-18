package photons4d

import (
	"fmt"
	"math"
)

// Ditorus (local space):
//
//	( sqrt( (sqrt(u² + v²) - R)² + z² ) - M )² + s² <= r²
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
//	Scale        = (sx, sy, sz, sw)
//	MajorRadius  = R
//	MiddleRadius = M
//	MinorRadius  = r
//
// This is a tube of radius r around the surface of an ordinary torus (ring
// radius R, tube radius M) that lies in the XYZ subspace. Like the duotorus
// it is a T² × D² solid with a 3-torus surface, but its core is the round
// torus and not the flat one.
//
// The local object is then rotated by R and translated to Center.
type Ditorus struct {
	Center       Point4
	Scale        Vector4
	MajorRadius  Real
	MiddleRadius Real
	MinorRadius  Real
	R            Mat4
	RT           Mat4

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

func NewDitorus(center Point4, scale Vector4, majorRadius, middleRadius, minorRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Ditorus, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("ditorus scale must be >0 on all axes, got %+v", scale)
	}
	if !(minorRadius > 0) {
		return nil, fmt.Errorf("ditorus minorRadius must be > 0, got %.6g", minorRadius)
	}
	if !(middleRadius > minorRadius) {
		return nil, fmt.Errorf("ditorus requires middleRadius > minorRadius, got middle=%.6g minor=%.6g", middleRadius, minorRadius)
	}
	if !(majorRadius > middleRadius+minorRadius) {
		return nil, fmt.Errorf("ditorus requires majorRadius > middleRadius + minorRadius for a proper ring, got major=%.6g middle=%.6g minor=%.6g", majorRadius, middleRadius, minorRadius)
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
	d := &Ditorus{
		Center:       center,
		Scale:        scale,
		MajorRadius:  majorRadius,
		MiddleRadius: middleRadius,
		MinorRadius:  minorRadius,
		R:            R,
		RT:           R.Transpose(),
		Color:        color,
		Diffuse:      diffuse,
		Reflect:      reflectivity,
		Refract:      refractivity,
		IOR:          ior,
	}

	// Exact support-based AABB.
	// In scaled local space this shape is:
	//   circle(radius=R in XY) + 3D ball(radius=M in XYZ) + 4D ball(radius=r).
	// So support along local direction a is:
	//   R * sqrt((sx*ax)^2 + (sy*ay)^2)
	// + M * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2)
	// + r * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2 + (sw*aw)^2)
	axisExtent := func(row int) (minV, maxV Real) {
		ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
		Ax := scale.X * ax
		Ay := scale.Y * ay
		Az := scale.Z * az
		Aw := scale.W * aw
		off := majorRadius*math.Sqrt(Ax*Ax+Ay*Ay) +
			middleRadius*math.Sqrt(Ax*Ax+Ay*Ay+Az*Az) +
			minorRadius*math.Sqrt(Ax*Ax+Ay*Ay+Az*Az+Aw*Aw)
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
	d.AABBMin = Point4{minX, minY, minZ, minW}
	d.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	d.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	d.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	d.colorArr = [3]Real{color.R, color.G, color.B}
	d.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	d.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - d.refl[i] - d.refr[i] - d.diff[i]
		if p < 0 {
			p = 0
		}
		d.pAbs[i] = p
		d.iorInv[i] = 1 / d.iorArr[i]
		n := d.iorArr[i]
		r0 := (n - 1) / (n + 1)
		d.f0[i] = r0 * r0
	}

	DebugLog("Created ditorus: %+v", d)
	return d, nil
}

func intersectRayDitorus(O Point4, D Vector4, tor *Ditorus) (hit objectHit, ok bool) {
	Op := Vector4{O.X - tor.Center.X, O.Y - tor.Center.Y, O.Z - tor.Center.Z, O.W - tor.Center.W}
	Ol := tor.RT.MulVec(Op)
	Dl := tor.RT.MulVec(D)

	R := tor.MajorRadius
	M := tor.MiddleRadius
	r := tor.MinorRadius
	sx, sy, sz, sw := tor.Scale.X, tor.Scale.Y, tor.Scale.Z, tor.Scale.W

	box, ok := boxIntervalLocal(Ol, Dl, sx*(R+M+r), sy*(R+M+r), sz*(M+r), sw*r)
	if !ok {
		return objectHit{}, false
	}

	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	// In the scaled space the solid is the set of points within r of the
	// torus surface in XYZ: the distance to that surface minus r is an exact
	// signed distance. e is the signed distance to it inside XYZ.
	dist := func(tt Real) rayDist {
		q := A.Add(B.Mul(tt))
		rho := math.Hypot(q.X, q.Y)
		p := rho - R
		E := math.Hypot(p, q.Z)
		e := E - M
		delta := math.Hypot(e, q.W)
		rd := rayDist{d: delta - r}
		if rho > 0 && E > 0 && delta > 0 {
			rd.rate = (e/E*(p/rho*(q.X*B.X+q.Y*B.Y)+q.Z*B.Z) + q.W*B.W) / delta
			rd.smooth = math.Min(delta, math.Min(E, rho))
		}
		return rd
	}

	tHit, inv, ok := traceSolid(dist, B.Len(), box.t0, box.t1)
	if !ok {
		return objectHit{}, false
	}

	q := A.Add(B.Mul(tHit))
	rho := math.Hypot(q.X, q.Y)
	p := rho - R
	E := math.Hypot(p, q.Z)
	kz, kxy := Real(0), Real(0)
	if E > 0 {
		kz = (E - M) / E
		if rho > 0 {
			kxy = kz * p / rho
		}
	}
	Nl := Vector4{kxy * q.X / sx, kxy * q.Y / sy, kz * q.Z / sz, q.W / sw}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	Nw := tor.R.MulVec(Nl).Norm()
	return objectHit{t: tHit, Nw: Nw, mat: tor, inv: inv}, true
}

func (d *Ditorus) PAbsCh(i int) Real   { return d.pAbs[i] }
func (d *Ditorus) DiffCh(i int) Real   { return d.diff[i] }
func (d *Ditorus) ColorCh(i int) Real  { return d.colorArr[i] }
func (d *Ditorus) F0Ch(i int) Real     { return d.f0[i] }
func (d *Ditorus) ReflCh(i int) Real   { return d.refl[i] }
func (d *Ditorus) RefrCh(i int) Real   { return d.refr[i] }
func (d *Ditorus) IORCh(i int) Real    { return d.iorArr[i] }
func (d *Ditorus) IORInvCh(i int) Real { return d.iorInv[i] }
