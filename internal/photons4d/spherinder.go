package photons4d

import (
	"fmt"
	"math"
)

const SpherinderConservativeAABB = false

// Spherinder (local space, axis along local W):
//
//	x²/rx² + y²/ry² + z²/rz² <= 1
//	|w| <= h
//
// with Scale = (rx, ry, rz, h), where:
//   - rx, ry, rz are cross-section radii
//   - h is the half-length along local W
//
// The local object is then rotated by R and translated to Center.
// Spherinder: a 4D cylinder with a 3-ball/ellipsoidal cross-section extruded along local W.
// Local-space definition:
//
//	x^2/rx^2 + y^2/ry^2 + z^2/rz^2 <= 1
//	|w| <= h
//
// where Scale=(rx, ry, rz, h). The object is then rotated by R and translated to Center.
type Spherinder struct {
	Center Point4
	Scale  Vector4 // X/Y/Z = cross-section radii, W = half-length along local W
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

func NewSpherinder(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Spherinder, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("spherinder scale must be >0 on all axes, got %+v", scale)
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
	s := &Spherinder{
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

	if SpherinderConservativeAABB {
		abs := func(x Real) Real {
			if x < 0 {
				return -x
			}
			return x
		}
		extent := func(row int) (min, max Real) {
			off := abs(R.M[row][0])*scale.X + abs(R.M[row][1])*scale.Y + abs(R.M[row][2])*scale.Z + abs(R.M[row][3])*scale.W
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
		s.AABBMin = Point4{minX, minY, minZ, minW}
		s.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	} else {
		// Exact support-based AABB.
		// For a world-axis row a=(ax,ay,az,aw), the local support is:
		//   sqrt((ax*rx)^2 + (ay*ry)^2 + (az*rz)^2) + h*|aw|
		// because this is an ellipsoidal 3-ball extruded symmetrically along W.
		rx, ry, rz, h := scale.X, scale.Y, scale.Z, scale.W
		axisExtent := func(row int) (min, max Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			off := math.Sqrt((ax*rx)*(ax*rx)+(ay*ry)*(ay*ry)+(az*rz)*(az*rz)) + math.Abs(aw)*h
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
		s.AABBMin = Point4{minX, minY, minZ, minW}
		s.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	}

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

	DebugLog("Created spherinder: %+v", s)
	return s, nil
}

func intersectRaySpherinder(O Point4, D Vector4, s *Spherinder) (hit objectHit, ok bool) {
	Op := Vector4{O.X - s.Center.X, O.Y - s.Center.Y, O.Z - s.Center.Z, O.W - s.Center.W}
	Ol := s.RT.MulVec(Op)
	Dl := s.RT.MulVec(D)

	rx, ry, rz, h := s.Scale.X, s.Scale.Y, s.Scale.Z, s.Scale.W
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)

	const eps = 1e-12

	sideEnter, sideExit := Real(math.Inf(-1)), Real(math.Inf(1))
	A := Dl.X*Dl.X*invRx2 + Dl.Y*Dl.Y*invRy2 + Dl.Z*Dl.Z*invRz2
	B := 2 * (Ol.X*Dl.X*invRx2 + Ol.Y*Dl.Y*invRy2 + Ol.Z*Dl.Z*invRz2)
	C := Ol.X*Ol.X*invRx2 + Ol.Y*Ol.Y*invRy2 + Ol.Z*Ol.Z*invRz2 - 1
	if math.Abs(float64(A)) < eps {
		if C > 0 {
			return objectHit{}, false
		}
	} else {
		disc := B*B - 4*A*C
		if disc < 0 {
			return objectHit{}, false
		}
		sqrtD := math.Sqrt(disc)
		inv2A := 1 / (2 * A)
		t0 := (-B - sqrtD) * inv2A
		t1 := (-B + sqrtD) * inv2A
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		sideEnter, sideExit = t0, t1
	}

	capEnter, capExit := Real(math.Inf(-1)), Real(math.Inf(1))
	if math.Abs(float64(Dl.W)) < eps {
		if Ol.W < -h || Ol.W > h {
			return objectHit{}, false
		}
	} else {
		t0 := (-h - Ol.W) / Dl.W
		t1 := (h - Ol.W) / Dl.W
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		capEnter, capExit = t0, t1
	}

	tEnter := sideEnter
	if capEnter > tEnter {
		tEnter = capEnter
	}
	tExit := sideExit
	if capExit < tExit {
		tExit = capExit
	}
	if tEnter > tExit {
		return objectHit{}, false
	}
	if tExit <= eps && tEnter <= eps {
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
	capTol := 1e-9 * (1 + h)
	if math.Abs(float64(math.Abs(P.W)-h)) <= capTol {
		Nl := Vector4{0, 0, 0, 1}
		if P.W < 0 {
			Nl = Nl.Mul(-1)
		}
		Nw := s.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: s, inv: inv}, true
	}

	Nl := Vector4{P.X * invRx2, P.Y * invRy2, P.Z * invRz2, 0}
	Nw := s.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: s, inv: inv}, true
}

func (s *Spherinder) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *Spherinder) DiffCh(i int) Real   { return s.diff[i] }
func (s *Spherinder) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *Spherinder) F0Ch(i int) Real     { return s.f0[i] }
func (s *Spherinder) ReflCh(i int) Real   { return s.refl[i] }
func (s *Spherinder) RefrCh(i int) Real   { return s.refr[i] }
func (s *Spherinder) IORCh(i int) Real    { return s.iorArr[i] }
func (s *Spherinder) IORInvCh(i int) Real { return s.iorInv[i] }
