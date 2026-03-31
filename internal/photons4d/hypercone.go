package photons4d

import (
	"fmt"
	"math"
	"sort"
)

// HyperCone (local space, axis along local W):
//
//	x²/rx² + y²/ry² + z²/rz² <= u²
//	u = (w + h) / (2h)
//	-h <= w <= h
//
// with Scale = (rx, ry, rz, h), where:
//   - the apex is at local w = -h
//   - the base is at local w = +h
//   - the base radii are (rx, ry, rz)
//
// So the cross-section shrinks linearly to a point at the apex and grows
// linearly to an ellipsoidal 3-ball at the base. The local object is then
// rotated by R and translated to Center.
type HyperCone struct {
	Center Point4
	Scale  Vector4 // X/Y/Z = base radii, W = half-length along local W
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

func NewHyperCone(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*HyperCone, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("hypercone scale must be >0 on all axes, got %+v", scale)
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
	hc := &HyperCone{
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

	// Exact support-based AABB.
	// For a world-axis row a=(ax,ay,az,aw), local support is attained either:
	//   - at the apex:           -h * aw
	//   - on the base ellipsoid: sqrt((ax*rx)^2 + (ay*ry)^2 + (az*rz)^2) + h * aw
	// because the support is linear in the interpolation parameter from apex to base.
	rx, ry, rz, h := scale.X, scale.Y, scale.Z, scale.W
	axisExtent := func(row int) (minV, maxV Real) {
		ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
		S := math.Sqrt((ax*rx)*(ax*rx) + (ay*ry)*(ay*ry) + (az*rz)*(az*rz))
		maxLocal := math.Max(-h*aw, S+h*aw)
		minLocal := math.Min(-h*aw, -S+h*aw)
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
	hc.AABBMin = Point4{minX, minY, minZ, minW}
	hc.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	hc.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	hc.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	hc.colorArr = [3]Real{color.R, color.G, color.B}
	hc.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	hc.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - hc.refl[i] - hc.refr[i] - hc.diff[i]
		if p < 0 {
			p = 0
		}
		hc.pAbs[i] = p
		hc.iorInv[i] = 1 / hc.iorArr[i]
		n := hc.iorArr[i]
		r0 := (n - 1) / (n + 1)
		hc.f0[i] = r0 * r0
	}

	DebugLog("Created hypercone: %+v", hc)
	return hc, nil
}

func intersectRayHyperCone(O Point4, D Vector4, c *HyperCone) (hit objectHit, ok bool) {
	Op := Vector4{O.X - c.Center.X, O.Y - c.Center.Y, O.Z - c.Center.Z, O.W - c.Center.W}
	Ol := c.RT.MulVec(Op)
	Dl := c.RT.MulVec(D)

	rx, ry, rz, h := c.Scale.X, c.Scale.Y, c.Scale.Z, c.Scale.W
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)
	const eps = 1e-12

	sideValue := func(p Vector4) Real {
		u := (p.W + h) / (2 * h)
		return p.X*p.X*invRx2 + p.Y*p.Y*invRy2 + p.Z*p.Z*invRz2 - u*u
	}
	insideAt := func(t Real) bool {
		p := Ol.Add(Dl.Mul(t))
		if p.W < -h-1e-9 || p.W > h+1e-9 {
			return false
		}
		return sideValue(p) <= 1e-9
	}

	// Cap interval from -h <= w <= h
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

	// Side surface roots of:
	//   x²/rx² + y²/ry² + z²/rz² - ((w+h)/(2h))² = 0
	u0 := (Ol.W + h) / (2 * h)
	du := Dl.W / (2 * h)
	crossO := Ol.X*Ol.X*invRx2 + Ol.Y*Ol.Y*invRy2 + Ol.Z*Ol.Z*invRz2
	crossOD := Ol.X*Dl.X*invRx2 + Ol.Y*Dl.Y*invRy2 + Ol.Z*Dl.Z*invRz2
	crossD := Dl.X*Dl.X*invRx2 + Dl.Y*Dl.Y*invRy2 + Dl.Z*Dl.Z*invRz2
	A := crossD - du*du
	B := 2 * (crossOD - u0*du)
	C := crossO - u0*u0

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

	delta := 1 + rx + ry + rz + h
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
		Nw := c.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
	}

	// Side normal from gradient of:
	//   F = x²/rx² + y²/ry² + z²/rz² - ((w+h)/(2h))²
	u := (P.W + h) / (2 * h)
	Nl := Vector4{
		2 * P.X * invRx2,
		2 * P.Y * invRy2,
		2 * P.Z * invRz2,
		-u / h,
	}
	if Nl.Len() == 0 {
		// Singular apex fallback
		Nl = Vector4{0, 0, 0, -1}
	}
	Nw := c.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
}

func (c *HyperCone) PAbsCh(i int) Real   { return c.pAbs[i] }
func (c *HyperCone) DiffCh(i int) Real   { return c.diff[i] }
func (c *HyperCone) ColorCh(i int) Real  { return c.colorArr[i] }
func (c *HyperCone) F0Ch(i int) Real     { return c.f0[i] }
func (c *HyperCone) ReflCh(i int) Real   { return c.refl[i] }
func (c *HyperCone) RefrCh(i int) Real   { return c.refr[i] }
func (c *HyperCone) IORCh(i int) Real    { return c.iorArr[i] }
func (c *HyperCone) IORInvCh(i int) Real { return c.iorInv[i] }
