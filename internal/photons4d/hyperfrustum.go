package photons4d

import (
	"fmt"
	"math"
	"sort"
)

const HyperFrustumConservativeAABB = false

// HyperFrustum / truncated hypercone (local space, axis along local W):
//
//	x²/rx² + y²/ry² + z²/rz² <= m(w)²
//	-h <= w <= h
//
//	m(w) = a + b w
//	a = (lower + upper) / 2
//	b = (upper - lower) / (2 h)
//
// with Scale = (rx, ry, rz, h), where:
//   - lower is the cross-section multiplier at local w = -h
//   - upper is the cross-section multiplier at local w = +h
//
// So the cross-section is an ellipsoidal 3-ball whose overall size changes
// linearly along W, but keeps the same XYZ aspect ratio throughout.
//
// Special cases:
//   - lower == upper  -> spherinder-like elliptical cylinder
//   - lower < upper   -> widening frustum
//   - lower > upper   -> narrowing frustum
//
// The local object is then rotated by R and translated to Center.
type HyperFrustum struct {
	Center      Point4
	Scale       Vector4 // X/Y/Z = base shape radii, W = half-length along local W
	LowerRadius Real
	UpperRadius Real
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

func NewHyperFrustum(center Point4, scale Vector4, lowerRadius, upperRadius Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*HyperFrustum, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("hyperfrustum scale must be >0 on all axes, got %+v", scale)
	}
	if !(lowerRadius > 0) {
		return nil, fmt.Errorf("hyperfrustum lowerRadius must be > 0, got %.6g", lowerRadius)
	}
	if !(upperRadius > 0) {
		return nil, fmt.Errorf("hyperfrustum upperRadius must be > 0, got %.6g", upperRadius)
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
	hf := &HyperFrustum{
		Center:      center,
		Scale:       scale,
		LowerRadius: lowerRadius,
		UpperRadius: upperRadius,
		R:           R,
		RT:          R.Transpose(),
		Color:       color,
		Diffuse:     diffuse,
		Reflect:     reflectivity,
		Refract:     refractivity,
		IOR:         ior,
	}

	if HyperFrustumConservativeAABB {
		// Conservative rotated box of the larger end radius.
		m := math.Max(lowerRadius, upperRadius)
		hx := scale.X * m
		hy := scale.Y * m
		hz := scale.Z * m
		hw := scale.W

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
		hf.AABBMin = Point4{minX, minY, minZ, minW}
		hf.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	} else {
		// Exact support-based AABB.
		// For a world-axis row a=(ax,ay,az,aw), support is linear in w,
		// so the maximum is attained at one of the two end caps:
		//   lower * S - h*aw
		//   upper * S + h*aw
		// where S = sqrt((rx*ax)^2 + (ry*ay)^2 + (rz*az)^2)
		rx, ry, rz, h := scale.X, scale.Y, scale.Z, scale.W
		axisExtent := func(row int) (minV, maxV Real) {
			ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
			S := math.Sqrt((rx*ax)*(rx*ax) + (ry*ay)*(ry*ay) + (rz*az)*(rz*az))

			maxLocal := math.Max(lowerRadius*S-h*aw, upperRadius*S+h*aw)
			minLocal := math.Min(-lowerRadius*S-h*aw, -upperRadius*S+h*aw)

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
			return c + minLocal, c + maxLocal
		}
		minX, maxX := axisExtent(0)
		minY, maxY := axisExtent(1)
		minZ, maxZ := axisExtent(2)
		minW, maxW := axisExtent(3)
		hf.AABBMin = Point4{minX, minY, minZ, minW}
		hf.AABBMax = Point4{maxX, maxY, maxZ, maxW}
	}

	hf.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	hf.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	hf.colorArr = [3]Real{color.R, color.G, color.B}
	hf.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	hf.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - hf.refl[i] - hf.refr[i] - hf.diff[i]
		if p < 0 {
			p = 0
		}
		hf.pAbs[i] = p
		hf.iorInv[i] = 1 / hf.iorArr[i]
		n := hf.iorArr[i]
		r0 := (n - 1) / (n + 1)
		hf.f0[i] = r0 * r0
	}

	DebugLog("Created hyperfrustum: %+v", hf)
	return hf, nil
}

func intersectRayHyperFrustum(O Point4, D Vector4, f *HyperFrustum) (hit objectHit, ok bool) {
	Op := Vector4{O.X - f.Center.X, O.Y - f.Center.Y, O.Z - f.Center.Z, O.W - f.Center.W}
	Ol := f.RT.MulVec(Op)
	Dl := f.RT.MulVec(D)

	rx, ry, rz, h := f.Scale.X, f.Scale.Y, f.Scale.Z, f.Scale.W
	lower, upper := f.LowerRadius, f.UpperRadius

	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)

	// m(w) = a + b w
	a := 0.5 * (lower + upper)
	b := (upper - lower) / (2 * h)

	const eps = 1e-12

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

	crossO := Ol.X*Ol.X*invRx2 + Ol.Y*Ol.Y*invRy2 + Ol.Z*Ol.Z*invRz2
	crossOD := Ol.X*Dl.X*invRx2 + Ol.Y*Dl.Y*invRy2 + Ol.Z*Dl.Z*invRz2
	crossD := Dl.X*Dl.X*invRx2 + Dl.Y*Dl.Y*invRy2 + Dl.Z*Dl.Z*invRz2

	M0 := a + b*Ol.W
	M1 := b * Dl.W

	A := crossD - M1*M1
	B := 2 * (crossOD - M0*M1)
	C := crossO - M0*M0

	cands := []Real{capEnter, capExit}

	if math.Abs(float64(A)) < eps {
		if math.Abs(float64(B)) >= eps {
			t := -C / B
			if t >= capEnter-1e-9 && t <= capExit+1e-9 {
				cands = append(cands, t)
			}
		}
	} else {
		disc := B*B - 4*A*C
		if disc >= 0 {
			sqrtD := math.Sqrt(disc)
			t0 := (-B - sqrtD) / (2 * A)
			t1 := (-B + sqrtD) / (2 * A)
			if t0 >= capEnter-1e-9 && t0 <= capExit+1e-9 {
				cands = append(cands, t0)
			}
			if t1 >= capEnter-1e-9 && t1 <= capExit+1e-9 {
				cands = append(cands, t1)
			}
		}
	}

	sort.Slice(cands, func(i, j int) bool { return cands[i] < cands[j] })
	dedup := make([]Real, 0, len(cands))
	for _, x := range cands {
		if len(dedup) == 0 || math.Abs(float64(x-dedup[len(dedup)-1])) > 1e-10 {
			dedup = append(dedup, x)
		}
	}
	cands = dedup

	insideAt := func(t Real) bool {
		p := Ol.Add(Dl.Mul(t))
		if p.W < -h-1e-9 || p.W > h+1e-9 {
			return false
		}
		m := a + b*p.W
		if m <= 0 {
			return false
		}
		val := p.X*p.X*invRx2 + p.Y*p.Y*invRy2 + p.Z*p.Z*invRz2 - m*m
		return val <= 1e-9
	}

	delta := 1 + rx + ry + rz + h + lower + upper
	probePoint := func(a, b Real) Real {
		switch {
		case math.IsInf(a, -1) && math.IsInf(b, 1):
			return 0
		case math.IsInf(a, -1):
			return b - delta
		case math.IsInf(b, 1):
			return a + delta
		default:
			return 0.5 * (a + b)
		}
	}

	var tEnter, tExit Real
	found := false
	for i := 0; i+1 < len(cands); i++ {
		a, b := cands[i], cands[i+1]
		if b-a <= eps {
			continue
		}
		m := probePoint(a, b)
		if insideAt(m) {
			tEnter, tExit = a, b
			found = true
			break
		}
	}
	if !found {
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
	if math.Abs(float64(P.W-h)) <= capTol {
		Nl := Vector4{0, 0, 0, 1}
		Nw := f.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: f, inv: inv}, true
	}
	if math.Abs(float64(P.W+h)) <= capTol {
		Nl := Vector4{0, 0, 0, -1}
		Nw := f.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: f, inv: inv}, true
	}

	m := a + b*P.W
	Nl := Vector4{
		2 * P.X * invRx2,
		2 * P.Y * invRy2,
		2 * P.Z * invRz2,
		-2 * m * b,
	}
	if Nl.Len() == 0 {
		Nl = Vector4{0, 0, 0, 1}
	}
	Nw := f.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: f, inv: inv}, true
}

func (f *HyperFrustum) PAbsCh(i int) Real   { return f.pAbs[i] }
func (f *HyperFrustum) DiffCh(i int) Real   { return f.diff[i] }
func (f *HyperFrustum) ColorCh(i int) Real  { return f.colorArr[i] }
func (f *HyperFrustum) F0Ch(i int) Real     { return f.f0[i] }
func (f *HyperFrustum) ReflCh(i int) Real   { return f.refl[i] }
func (f *HyperFrustum) RefrCh(i int) Real   { return f.refr[i] }
func (f *HyperFrustum) IORCh(i int) Real    { return f.iorArr[i] }
func (f *HyperFrustum) IORInvCh(i int) Real { return f.iorInv[i] }
