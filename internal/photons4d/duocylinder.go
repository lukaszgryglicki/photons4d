package photons4d

import (
	"fmt"
	"math"
)

const DuocylinderConservativeAABB = false

// Duocylinder (local space):
//
//	x²/rx² + y²/ry² <= 1
//	z²/rz² + w²/rw² <= 1
//
// with Scale = (rx, ry, rz, rw), where:
//   - (rx, ry) are the semi-axes of the first 2D ellipse in the XY plane
//   - (rz, rw) are the semi-axes of the second 2D ellipse in the ZW plane
//
// So this is D² × D² in 4D: the Cartesian product of two disks/ellipses.
// The local object is then rotated by R and translated to Center.
type Duocylinder struct {
	Center Point4
	Scale  Vector4
	R      Mat4
	RT     Mat4

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

func NewDuocylinder(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Duocylinder, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("duocylinder scale must be >0 on all axes, got %+v", scale)
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
	d := &Duocylinder{
		Center:  center,
		Scale:   scale,
		R:       R,
		RT:      R.Transpose(),
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}

	if DuocylinderConservativeAABB {
		// Conservative rotated bounding box of the local extents.
		hx, hy, hz, hw := scale.X, scale.Y, scale.Z, scale.W
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
		// For a world-axis row a=(ax,ay,az,aw), support is:
		//   sqrt((rx*ax)^2 + (ry*ay)^2) + sqrt((rz*az)^2 + (rw*aw)^2)
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			off1 := math.Sqrt((scale.X*ax)*(scale.X*ax) + (scale.Y*ay)*(scale.Y*ay))
			off2 := math.Sqrt((scale.Z*az)*(scale.Z*az) + (scale.W*aw)*(scale.W*aw))
			off := off1 + off2
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

	DebugLog("Created duocylinder: %+v", d)
	return d, nil
}

func intervalEllipse2DLocal(ou, ov, du, dv, ru, rv Real) (tInterval, bool) {
	const eps = 1e-12
	invRu2 := 1 / (ru * ru)
	invRv2 := 1 / (rv * rv)

	A := du*du*invRu2 + dv*dv*invRv2
	B := 2 * (ou*du*invRu2 + ov*dv*invRv2)
	C := ou*ou*invRu2 + ov*ov*invRv2 - 1

	if math.Abs(A) < eps {
		if C <= 0 {
			return tInterval{Real(math.Inf(-1)), Real(math.Inf(1))}, true
		}
		return tInterval{}, false
	}

	disc := B*B - 4*A*C
	if disc < 0 {
		if disc > -1e-15 {
			disc = 0
		} else {
			return tInterval{}, false
		}
	}
	sqrtD := math.Sqrt(disc)
	t0 := (-B - sqrtD) / (2 * A)
	t1 := (-B + sqrtD) / (2 * A)
	if t0 > t1 {
		t0, t1 = t1, t0
	}
	return tInterval{t0, t1}, true
}

func intersectRayDuocylinder(O Point4, D Vector4, dc *Duocylinder) (hit objectHit, ok bool) {
	Op := Vector4{O.X - dc.Center.X, O.Y - dc.Center.Y, O.Z - dc.Center.Z, O.W - dc.Center.W}
	Ol := dc.RT.MulVec(Op)
	Dl := dc.RT.MulVec(D)

	rx, ry, rz, rw := dc.Scale.X, dc.Scale.Y, dc.Scale.Z, dc.Scale.W

	ivXY, okXY := intervalEllipse2DLocal(Ol.X, Ol.Y, Dl.X, Dl.Y, rx, ry)
	if !okXY {
		return objectHit{}, false
	}
	ivZW, okZW := intervalEllipse2DLocal(Ol.Z, Ol.W, Dl.Z, Dl.W, rz, rw)
	if !okZW {
		return objectHit{}, false
	}

	tEnter := ivXY.t0
	if ivZW.t0 > tEnter {
		tEnter = ivZW.t0
	}
	tExit := ivXY.t1
	if ivZW.t1 < tExit {
		tExit = ivZW.t1
	}

	const eps = 1e-12
	if tEnter > tExit || (tExit <= eps && tEnter <= eps) {
		return objectHit{}, false
	}

	inv := tEnter <= eps
	t := tEnter
	if inv {
		t = tExit
	}
	if t <= eps {
		return objectHit{}, false
	}

	P := Ol.Add(Dl.Mul(t))
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)
	invRw2 := 1 / (rw * rw)

	fXY := P.X*P.X*invRx2 + P.Y*P.Y*invRy2 - 1
	fZW := P.Z*P.Z*invRz2 + P.W*P.W*invRw2 - 1
	tol := 1e-8

	activeXY := math.Abs(fXY) <= tol
	activeZW := math.Abs(fZW) <= tol
	if !activeXY && !activeZW {
		if math.Abs(fXY) <= math.Abs(fZW) {
			activeXY = true
		} else {
			activeZW = true
		}
	}

	Nl := Vector4{}
	if activeXY {
		Nl = Nl.Add(Vector4{P.X * invRx2, P.Y * invRy2, 0, 0})
	}
	if activeZW {
		Nl = Nl.Add(Vector4{0, 0, P.Z * invRz2, P.W * invRw2})
	}
	if Nl.Len() == 0 {
		if math.Abs(P.X/rx)+math.Abs(P.Y/ry) >= math.Abs(P.Z/rz)+math.Abs(P.W/rw) {
			Nl = Vector4{1, 0, 0, 0}
		} else {
			Nl = Vector4{0, 0, 1, 0}
		}
	}

	Nw := dc.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: dc, inv: inv}, true
}

func (d *Duocylinder) PAbsCh(i int) Real   { return d.pAbs[i] }
func (d *Duocylinder) DiffCh(i int) Real   { return d.diff[i] }
func (d *Duocylinder) ColorCh(i int) Real  { return d.colorArr[i] }
func (d *Duocylinder) F0Ch(i int) Real     { return d.f0[i] }
func (d *Duocylinder) ReflCh(i int) Real   { return d.refl[i] }
func (d *Duocylinder) RefrCh(i int) Real   { return d.refr[i] }
func (d *Duocylinder) IORCh(i int) Real    { return d.iorArr[i] }
func (d *Duocylinder) IORInvCh(i int) Real { return d.iorInv[i] }
