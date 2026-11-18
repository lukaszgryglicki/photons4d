package photons4d

import (
	"fmt"
	"math"
)

const DuotorusConservativeAABB = false

// Duotorus / Clifford-torus tube (local space):
//
//	( sqrt(u² + v²) - Rxy )² + ( sqrt(z² + s²) - Rzw )² <= r²
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
//	Scale         = (sx, sy, sz, sw)
//	MajorRadiusXY = Rxy
//	MajorRadiusZW = Rzw
//	MinorRadius   = r
//
// This is a tube around the product of two orthogonal circles:
//
//	u² + v² = Rxy²
//	z² + s² = Rzw²
//
// Topologically this is a T² × D² solid, i.e. a genuinely 4D torus-like object.
// The local object is then rotated by R and translated to Center.
type Duotorus struct {
	Center        Point4
	Scale         Vector4
	MajorRadiusXY Real
	MajorRadiusZW Real
	MinorRadius   Real
	R             Mat4
	RT            Mat4

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

func NewDuotorus(center Point4, scale Vector4, majorRadiusXY, majorRadiusZW, minorRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Duotorus, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("duotorus scale must be >0 on all axes, got %+v", scale)
	}
	if !(majorRadiusXY > 0) {
		return nil, fmt.Errorf("duotorus majorRadiusXY must be > 0, got %.6g", majorRadiusXY)
	}
	if !(majorRadiusZW > 0) {
		return nil, fmt.Errorf("duotorus majorRadiusZW must be > 0, got %.6g", majorRadiusZW)
	}
	if !(minorRadius > 0) {
		return nil, fmt.Errorf("duotorus minorRadius must be > 0, got %.6g", minorRadius)
	}
	if !(minorRadius < majorRadiusXY && minorRadius < majorRadiusZW) {
		return nil, fmt.Errorf("duotorus requires minorRadius < min(majorRadiusXY, majorRadiusZW), got Rxy=%.6g Rzw=%.6g r=%.6g", majorRadiusXY, majorRadiusZW, minorRadius)
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
	d := &Duotorus{
		Center:        center,
		Scale:         scale,
		MajorRadiusXY: majorRadiusXY,
		MajorRadiusZW: majorRadiusZW,
		MinorRadius:   minorRadius,
		R:             R,
		RT:            R.Transpose(),
		Color:         color,
		Diffuse:       diffuse,
		Reflect:       reflectivity,
		Refract:       refractivity,
		IOR:           ior,
	}

	if DuotorusConservativeAABB {
		// Conservative analytic AABB from the rotated bounding box of the
		// axis-scaled canonical duotorus.
		hx := scale.X * (majorRadiusXY + minorRadius)
		hy := scale.Y * (majorRadiusXY + minorRadius)
		hz := scale.Z * (majorRadiusZW + minorRadius)
		hw := scale.W * (majorRadiusZW + minorRadius)

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
		d.AABBMin = Point4{minX, minY, minZ, minW}
		d.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	} else {
		// Exact support-based AABB.
		// In scaled local space this shape is:
		//   (circle in XY of radius Rxy) × (circle in ZW of radius Rzw),
		// thickened by a 4D ball of radius r.
		// So support along local direction a is:
		//   Rxy * sqrt((sx*ax)^2 + (sy*ay)^2)
		// + Rzw * sqrt((sz*az)^2 + (sw*aw)^2)
		// + r * sqrt((sx*ax)^2 + (sy*ay)^2 + (sz*az)^2 + (sw*aw)^2)
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			Ax := scale.X * ax
			Ay := scale.Y * ay
			Az := scale.Z * az
			Aw := scale.W * aw
			A := math.Sqrt(Ax*Ax + Ay*Ay)
			B := math.Sqrt(Az*Az + Aw*Aw)
			off := majorRadiusXY*A + majorRadiusZW*B + minorRadius*math.Sqrt(A*A+B*B)
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
	}

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

	DebugLog("Created duotorus: %+v", d)
	return d, nil
}

func intersectRayDuotorus(O Point4, D Vector4, tor *Duotorus) (hit objectHit, ok bool) {
	Op := Vector4{O.X - tor.Center.X, O.Y - tor.Center.Y, O.Z - tor.Center.Z, O.W - tor.Center.W}
	Ol := tor.RT.MulVec(Op)
	Dl := tor.RT.MulVec(D)

	Rxy := tor.MajorRadiusXY
	Rzw := tor.MajorRadiusZW
	r := tor.MinorRadius
	sx, sy, sz, sw := tor.Scale.X, tor.Scale.Y, tor.Scale.Z, tor.Scale.W

	box, ok := boxIntervalLocal(Ol, Dl, sx*(Rxy+r), sy*(Rxy+r), sz*(Rzw+r), sw*(Rzw+r))
	if !ok {
		return objectHit{}, false
	}

	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	// In the scaled space the solid is the set of points within r of the flat
	// torus (the product of the two circles): the distance to that torus minus
	// r is an exact signed distance.
	dist := func(tt Real) rayDist {
		q := A.Add(B.Mul(tt))
		a := math.Hypot(q.X, q.Y)
		b := math.Hypot(q.Z, q.W)
		pa := a - Rxy
		pb := b - Rzw
		delta := math.Hypot(pa, pb)
		rd := rayDist{d: delta - r}
		if a > 0 && b > 0 && delta > 0 {
			rd.rate = (pa/a*(q.X*B.X+q.Y*B.Y) + pb/b*(q.Z*B.Z+q.W*B.W)) / delta
			rd.smooth = math.Min(delta, math.Min(a, b))
		}
		return rd
	}

	tHit, inv, ok := traceSolid(dist, B.Len(), box.t0, box.t1)
	if !ok {
		return objectHit{}, false
	}

	q := A.Add(B.Mul(tHit))
	ka, kb := Real(0), Real(0)
	if a := math.Hypot(q.X, q.Y); a > 0 {
		ka = (a - Rxy) / a
	}
	if b := math.Hypot(q.Z, q.W); b > 0 {
		kb = (b - Rzw) / b
	}
	Nl := Vector4{ka * q.X / sx, ka * q.Y / sy, kb * q.Z / sz, kb * q.W / sw}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}
	Nw := tor.R.MulVec(Nl).Norm()
	return objectHit{t: tHit, Nw: Nw, mat: tor, inv: inv}, true
}

func (d *Duotorus) PAbsCh(i int) Real   { return d.pAbs[i] }
func (d *Duotorus) DiffCh(i int) Real   { return d.diff[i] }
func (d *Duotorus) ColorCh(i int) Real  { return d.colorArr[i] }
func (d *Duotorus) F0Ch(i int) Real     { return d.f0[i] }
func (d *Duotorus) ReflCh(i int) Real   { return d.refl[i] }
func (d *Duotorus) RefrCh(i int) Real   { return d.refr[i] }
func (d *Duotorus) IORCh(i int) Real    { return d.iorArr[i] }
func (d *Duotorus) IORInvCh(i int) Real { return d.iorInv[i] }
