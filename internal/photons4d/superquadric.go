package photons4d

import (
	"fmt"
	"math"
)

const SuperquadricConservativeAABB = false

// Superquadric / Lp hyperellipsoid (local space):
//
//	|x/sx|^p + |y/sy|^p + |z/sz|^p + |w/sw|^p <= 1
//
// with:
//
//	Scale = (sx, sy, sz, sw)
//	Power = p > 1
//
// For p = 2 this is the usual hyperellipsoid.
// As p grows, the shape becomes more box-like.
// The local object is then rotated by R and translated to Center.
type Superquadric struct {
	Center Point4
	Scale  Vector4
	Power  Real
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

func NewSuperquadric(center Point4, scale Vector4, power Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Superquadric, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("superquadric scale must be >0 on all axes, got %+v", scale)
	}
	if !(power > 1) {
		return nil, fmt.Errorf("superquadric power must be > 1, got %.6g", power)
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
	sq := &Superquadric{
		Center:  center,
		Scale:   scale,
		Power:   power,
		R:       R,
		RT:      R.Transpose(),
		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}

	if SuperquadricConservativeAABB {
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
		sq.AABBMin = Point4{minX, minY, minZ, minW}
		sq.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	} else {
		// Exact support-based AABB for the Lp ball under anisotropic scale.
		// Dual exponent q = p/(p-1).
		q := power / (power - 1)
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			tx := math.Pow(math.Abs(scale.X*ax), q)
			ty := math.Pow(math.Abs(scale.Y*ay), q)
			tz := math.Pow(math.Abs(scale.Z*az), q)
			tw := math.Pow(math.Abs(scale.W*aw), q)
			off := math.Pow(tx+ty+tz+tw, 1/q)
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
		sq.AABBMin = Point4{minX, minY, minZ, minW}
		sq.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	}

	sq.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	sq.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	sq.colorArr = [3]Real{color.R, color.G, color.B}
	sq.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	sq.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - sq.refl[i] - sq.refr[i] - sq.diff[i]
		if p < 0 {
			p = 0
		}
		sq.pAbs[i] = p
		sq.iorInv[i] = 1 / sq.iorArr[i]
		n := sq.iorArr[i]
		r0 := (n - 1) / (n + 1)
		sq.f0[i] = r0 * r0
	}

	DebugLog("Created superquadric: %+v", sq)
	return sq, nil
}

func superquadricImplicitLocal(p Vector4, scale Vector4, power Real) Real {
	tx := math.Pow(math.Abs(p.X/scale.X), power)
	ty := math.Pow(math.Abs(p.Y/scale.Y), power)
	tz := math.Pow(math.Abs(p.Z/scale.Z), power)
	tw := math.Pow(math.Abs(p.W/scale.W), power)
	return tx + ty + tz + tw - 1
}

func superquadricDerivAlong(Ol, Dl Vector4, scale Vector4, power, t Real) Real {
	term := func(o, d, s Real) Real {
		x := o + t*d
		ax := math.Abs(x / s)
		if ax == 0 {
			return 0
		}
		return power * math.Copysign(math.Pow(ax, power-1), x) * d / s
	}
	return term(Ol.X, Dl.X, scale.X) +
		term(Ol.Y, Dl.Y, scale.Y) +
		term(Ol.Z, Dl.Z, scale.Z) +
		term(Ol.W, Dl.W, scale.W)
}

func bisectSignRoot(f func(Real) Real, a, b Real) Real {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a
	}
	if fb == 0 {
		return b
	}
	lo, hi := a, b
	for i := 0; i < 100; i++ {
		m := 0.5 * (lo + hi)
		fm := f(m)
		if fm == 0 {
			return m
		}
		if (fa < 0 && fm > 0) || (fa > 0 && fm < 0) {
			hi = m
			fb = fm
		} else {
			lo = m
			fa = fm
		}
		_ = fb
	}
	return 0.5 * (lo + hi)
}

func superquadricMinimizer(Ol, Dl Vector4, scale Vector4, power, t0, t1 Real) Real {
	d0 := superquadricDerivAlong(Ol, Dl, scale, power, t0)
	if d0 >= 0 {
		return t0
	}
	d1 := superquadricDerivAlong(Ol, Dl, scale, power, t1)
	if d1 <= 0 {
		return t1
	}
	lo, hi := t0, t1
	for i := 0; i < 100; i++ {
		m := 0.5 * (lo + hi)
		dm := superquadricDerivAlong(Ol, Dl, scale, power, m)
		if dm > 0 {
			hi = m
		} else {
			lo = m
		}
	}
	return 0.5 * (lo + hi)
}

func intersectRaySuperquadric(O Point4, D Vector4, sq *Superquadric) (hit objectHit, ok bool) {
	Op := Vector4{O.X - sq.Center.X, O.Y - sq.Center.Y, O.Z - sq.Center.Z, O.W - sq.Center.W}
	Ol := sq.RT.MulVec(Op)
	Dl := sq.RT.MulVec(D)

	box, ok := boxIntervalLocal(Ol, Dl, sq.Scale.X, sq.Scale.Y, sq.Scale.Z, sq.Scale.W)
	if !ok {
		return objectHit{}, false
	}

	f := func(t Real) Real {
		p := Ol.Add(Dl.Mul(t))
		return superquadricImplicitLocal(p, sq.Scale, sq.Power)
	}

	const eps = 1e-12
	const tol = 1e-10

	t0, t1 := box.t0, box.t1
	minT := superquadricMinimizer(Ol, Dl, sq.Scale, sq.Power, t0, t1)
	fmin := f(minT)
	if fmin > tol {
		return objectHit{}, false
	}

	fAt0 := f(0)
	inv := fAt0 <= tol
	var tHit Real
	found := false

	if inv {
		if t1 <= eps {
			return objectHit{}, false
		}
		f1 := f(t1)
		if f1 <= tol {
			tHit = t1
			found = true
		} else {
			tHit = bisectSignRoot(f, 0, t1)
			found = tHit > eps
		}
	} else {
		f0 := f(t0)
		if f0 <= tol {
			if t0 > eps {
				tHit = t0
				found = true
			}
		} else if fmin <= tol {
			tHit = bisectSignRoot(f, t0, minT)
			found = tHit > eps
		}
	}

	if !found || tHit <= eps {
		return objectHit{}, false
	}

	Pq := Ol.Add(Dl.Mul(tHit))
	gradComp := func(x, s Real) Real {
		ax := math.Abs(x / s)
		if ax == 0 {
			return 0
		}
		return sq.Power * math.Copysign(math.Pow(ax, sq.Power-1), x) / s
	}
	Nl := Vector4{
		gradComp(Pq.X, sq.Scale.X),
		gradComp(Pq.Y, sq.Scale.Y),
		gradComp(Pq.Z, sq.Scale.Z),
		gradComp(Pq.W, sq.Scale.W),
	}
	if Nl.Len() == 0 {
		return objectHit{}, false
	}

	Nw := sq.R.MulVec(Nl).Norm()
	return objectHit{t: tHit, Nw: Nw, mat: sq, inv: inv}, true
}

func (s *Superquadric) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *Superquadric) DiffCh(i int) Real   { return s.diff[i] }
func (s *Superquadric) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *Superquadric) F0Ch(i int) Real     { return s.f0[i] }
func (s *Superquadric) ReflCh(i int) Real   { return s.refl[i] }
func (s *Superquadric) RefrCh(i int) Real   { return s.refr[i] }
func (s *Superquadric) IORCh(i int) Real    { return s.iorArr[i] }
func (s *Superquadric) IORInvCh(i int) Real { return s.iorInv[i] }
