package photons4d

import (
	"fmt"
	"math"
	"sort"
)

// HyperCapsule (local space, axis along local W):
//
//	min over t in [-h, h] of:
//	  x²/rx² + y²/ry² + z²/rz² + (w - t)²/rw² <= 1
//
// Equivalently:
//
//	x²/rx² + y²/ry² + z²/rz² + max(0, |w|-h)²/rw² <= 1
//
// with:
//
//	Scale      = (rx, ry, rz, rw)
//	HalfLength = h
//
// So this is a line segment along local W, swept by a 4D ellipsoid with
// radii (rx, ry, rz, rw). The local object is then rotated by R and
// translated to Center.
type HyperCapsule struct {
	Center     Point4
	Scale      Vector4 // ellipsoid radii of the swept profile: (rx, ry, rz, rw)
	HalfLength Real    // segment half-length along local W
	R          Mat4
	RT         Mat4

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

type tInterval struct {
	t0, t1 Real
}

func NewHyperCapsule(center Point4, scale Vector4, halfLength Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*HyperCapsule, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("hypercapsule scale must be >0 on all axes, got %+v", scale)
	}
	if !(halfLength > 0) {
		return nil, fmt.Errorf("hypercapsule halfLength must be > 0, got %.6g", halfLength)
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
	hc := &HyperCapsule{
		Center:     center,
		Scale:      scale,
		HalfLength: halfLength,
		R:          R,
		RT:         R.Transpose(),
		Color:      color,
		Diffuse:    diffuse,
		Reflect:    reflectivity,
		Refract:    refractivity,
		IOR:        ior,
	}

	// Exact support-based AABB for:
	//   segment along local W of half-length h
	//   Minkowski-summed with ellipsoid radii (rx, ry, rz, rw)
	//
	// Support along row-vector a is:
	//   h*|aw| + sqrt((ax*rx)^2 + (ay*ry)^2 + (az*rz)^2 + (aw*rw)^2)
	rx, ry, rz, rw := scale.X, scale.Y, scale.Z, scale.W
	h := halfLength
	axisExtent := func(row int) (minV, maxV Real) {
		ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
		off := h*math.Abs(aw) + math.Sqrt((ax*rx)*(ax*rx)+(ay*ry)*(ay*ry)+(az*rz)*(az*rz)+(aw*rw)*(aw*rw))
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

	DebugLog("Created hypercapsule: %+v", hc)
	return hc, nil
}

func intervalSpherinderLocal(Ol, Dl Vector4, rx, ry, rz, h Real) (tInterval, bool) {
	const eps = 1e-12
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)

	sideEnter, sideExit := Real(math.Inf(-1)), Real(math.Inf(1))
	A := Dl.X*Dl.X*invRx2 + Dl.Y*Dl.Y*invRy2 + Dl.Z*Dl.Z*invRz2
	B := 2 * (Ol.X*Dl.X*invRx2 + Ol.Y*Dl.Y*invRy2 + Ol.Z*Dl.Z*invRz2)
	C := Ol.X*Ol.X*invRx2 + Ol.Y*Ol.Y*invRy2 + Ol.Z*Ol.Z*invRz2 - 1
	if math.Abs(float64(A)) < eps {
		if C > 0 {
			return tInterval{}, false
		}
	} else {
		disc := B*B - 4*A*C
		if disc < 0 {
			return tInterval{}, false
		}
		sqrtD := math.Sqrt(disc)
		t0 := (-B - sqrtD) / (2 * A)
		t1 := (-B + sqrtD) / (2 * A)
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		sideEnter, sideExit = t0, t1
	}

	slabEnter, slabExit := Real(math.Inf(-1)), Real(math.Inf(1))
	if math.Abs(float64(Dl.W)) < eps {
		if Ol.W < -h || Ol.W > h {
			return tInterval{}, false
		}
	} else {
		t0 := (-h - Ol.W) / Dl.W
		t1 := (h - Ol.W) / Dl.W
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		slabEnter, slabExit = t0, t1
	}

	t0 := sideEnter
	if slabEnter > t0 {
		t0 = slabEnter
	}
	t1 := sideExit
	if slabExit < t1 {
		t1 = slabExit
	}
	if t0 > t1 {
		return tInterval{}, false
	}
	return tInterval{t0, t1}, true
}

func intervalEllipsoidLocal(Ol, Dl Vector4, centerW, rx, ry, rz, rw Real) (tInterval, bool) {
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)
	invRw2 := 1 / (rw * rw)

	ow := Ol.W - centerW
	A := Dl.X*Dl.X*invRx2 + Dl.Y*Dl.Y*invRy2 + Dl.Z*Dl.Z*invRz2 + Dl.W*Dl.W*invRw2
	B := 2 * (Ol.X*Dl.X*invRx2 + Ol.Y*Dl.Y*invRy2 + Ol.Z*Dl.Z*invRz2 + ow*Dl.W*invRw2)
	C := Ol.X*Ol.X*invRx2 + Ol.Y*Ol.Y*invRy2 + Ol.Z*Ol.Z*invRz2 + ow*ow*invRw2 - 1

	const eps = 1e-12
	if math.Abs(float64(A)) < eps {
		if C > 0 {
			return tInterval{}, false
		}
		return tInterval{Real(math.Inf(-1)), Real(math.Inf(1))}, true
	}
	disc := B*B - 4*A*C
	if disc < 0 {
		return tInterval{}, false
	}
	sqrtD := math.Sqrt(disc)
	t0 := (-B - sqrtD) / (2 * A)
	t1 := (-B + sqrtD) / (2 * A)
	if t0 > t1 {
		t0, t1 = t1, t0
	}
	return tInterval{t0, t1}, true
}

func mergeIntervals(in []tInterval) []tInterval {
	if len(in) == 0 {
		return nil
	}
	sort.Slice(in, func(i, j int) bool {
		if in[i].t0 == in[j].t0 {
			return in[i].t1 < in[j].t1
		}
		return in[i].t0 < in[j].t0
	})
	out := make([]tInterval, 0, len(in))
	cur := in[0]
	for i := 1; i < len(in); i++ {
		if in[i].t0 <= cur.t1+1e-10 {
			if in[i].t1 > cur.t1 {
				cur.t1 = in[i].t1
			}
			continue
		}
		out = append(out, cur)
		cur = in[i]
	}
	out = append(out, cur)
	return out
}

func intersectRayHyperCapsule(O Point4, D Vector4, c *HyperCapsule) (hit objectHit, ok bool) {
	Op := Vector4{O.X - c.Center.X, O.Y - c.Center.Y, O.Z - c.Center.Z, O.W - c.Center.W}
	Ol := c.RT.MulVec(Op)
	Dl := c.RT.MulVec(D)

	rx, ry, rz, rw := c.Scale.X, c.Scale.Y, c.Scale.Z, c.Scale.W
	h := c.HalfLength
	const eps = 1e-12

	parts := make([]tInterval, 0, 3)
	if iv, ok := intervalSpherinderLocal(Ol, Dl, rx, ry, rz, h); ok {
		parts = append(parts, iv)
	}
	if iv, ok := intervalEllipsoidLocal(Ol, Dl, +h, rx, ry, rz, rw); ok {
		parts = append(parts, iv)
	}
	if iv, ok := intervalEllipsoidLocal(Ol, Dl, -h, rx, ry, rz, rw); ok {
		parts = append(parts, iv)
	}
	if len(parts) == 0 {
		return objectHit{}, false
	}

	union := mergeIntervals(parts)
	if len(union) == 0 {
		return objectHit{}, false
	}

	inv := false
	var t Real
	found := false

	for _, iv := range union {
		if iv.t0 <= eps && iv.t1 > eps {
			inv = true
			t = iv.t1
			found = true
			break
		}
	}
	if !found {
		for _, iv := range union {
			if iv.t0 > eps {
				t = iv.t0
				found = true
				break
			}
		}
	}
	if !found || t <= eps {
		return objectHit{}, false
	}

	P := Ol.Add(Dl.Mul(t))
	invRx2 := 1 / (rx * rx)
	invRy2 := 1 / (ry * ry)
	invRz2 := 1 / (rz * rz)
	invRw2 := 1 / (rw * rw)
	tol := 1e-8

	// Classify the surface piece.
	if P.W > h+tol {
		Nl := Vector4{P.X * invRx2, P.Y * invRy2, P.Z * invRz2, (P.W - h) * invRw2}
		Nw := c.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
	}
	if P.W < -h-tol {
		Nl := Vector4{P.X * invRx2, P.Y * invRy2, P.Z * invRz2, (P.W + h) * invRw2}
		Nw := c.R.MulVec(Nl).Norm()
		return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
	}

	// Central spherinder side or seam.
	Nl := Vector4{P.X * invRx2, P.Y * invRy2, P.Z * invRz2, 0}
	if Nl.Len() == 0 {
		// Degenerate seam fallback.
		if P.W >= 0 {
			Nl = Vector4{0, 0, 0, 1}
		} else {
			Nl = Vector4{0, 0, 0, -1}
		}
	}
	Nw := c.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
}

func (c *HyperCapsule) PAbsCh(i int) Real   { return c.pAbs[i] }
func (c *HyperCapsule) DiffCh(i int) Real   { return c.diff[i] }
func (c *HyperCapsule) ColorCh(i int) Real  { return c.colorArr[i] }
func (c *HyperCapsule) F0Ch(i int) Real     { return c.f0[i] }
func (c *HyperCapsule) ReflCh(i int) Real   { return c.refl[i] }
func (c *HyperCapsule) RefrCh(i int) Real   { return c.refr[i] }
func (c *HyperCapsule) IORCh(i int) Real    { return c.iorArr[i] }
func (c *HyperCapsule) IORInvCh(i int) Real { return c.iorInv[i] }
