package photons4d

import (
	"fmt"
	"math"
	"sort"
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

func torusImplicitScaled(sx, sy, sz, sw, R, r Real) Real {
	sumsq := sx*sx + sy*sy + sz*sz + sw*sw
	a := sumsq + R*R - r*r
	return a*a - 4*R*R*(sx*sx+sw*sw)
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

func polyEvalQuartic(c4, c3, c2, c1, c0, x Real) Real {
	return ((((c4*x)+c3)*x+c2)*x+c1)*x + c0
}

func polyEvalCubic(c3, c2, c1, c0, x Real) Real {
	return (((c3*x)+c2)*x+c1)*x + c0
}

func uniqueSortedRoots(xs []Real, tol Real) []Real {
	if len(xs) == 0 {
		return nil
	}
	sort.Slice(xs, func(i, j int) bool { return xs[i] < xs[j] })
	out := make([]Real, 0, len(xs))
	for _, x := range xs {
		if len(out) == 0 || math.Abs(float64(x-out[len(out)-1])) > float64(tol) {
			out = append(out, x)
		}
	}
	return out
}

func solveQuadraticReal(a, b, c Real) []Real {
	const eps = 1e-12
	if math.Abs(float64(a)) < eps {
		if math.Abs(float64(b)) < eps {
			return nil
		}
		return []Real{-c / b}
	}
	disc := b*b - 4*a*c
	if disc < 0 {
		return nil
	}
	if disc == 0 {
		return []Real{-b / (2 * a)}
	}
	sqrtD := math.Sqrt(disc)
	r0 := (-b - sqrtD) / (2 * a)
	r1 := (-b + sqrtD) / (2 * a)
	return uniqueSortedRoots([]Real{r0, r1}, 1e-10)
}

func realRootsCubic(a, b, c, d Real) []Real {
	const eps = 1e-12

	if math.Abs(float64(a)) < eps {
		return solveQuadraticReal(b, c, d)
	}

	// Normalize to x^3 + A x^2 + B x + C = 0
	A := b / a
	B := c / a
	C := d / a

	// Depressed cubic y^3 + p y + q = 0, x = y - A/3
	p := B - A*A/3
	q := 2*A*A*A/27 - A*B/3 + C

	disc := q*q/4 + p*p*p/27
	shift := A / 3

	if disc > eps {
		sqrtD := math.Sqrt(disc)
		u := math.Cbrt(-q/2 + sqrtD)
		v := math.Cbrt(-q/2 - sqrtD)
		return []Real{u + v - shift}
	}

	if math.Abs(float64(disc)) <= eps {
		if math.Abs(float64(q)) <= eps {
			return []Real{-shift}
		}
		u := math.Cbrt(-q / 2)
		return uniqueSortedRoots([]Real{2*u - shift, -u - shift}, 1e-10)
	}

	// disc < 0: three real roots
	m := 2 * math.Sqrt(-p/3)
	arg := -q / (2 * math.Sqrt(-(p*p*p)/27))
	if arg < -1 {
		arg = -1
	}
	if arg > 1 {
		arg = 1
	}
	phi := math.Acos(arg)

	r0 := m*math.Cos(phi/3) - shift
	r1 := m*math.Cos((phi+2*math.Pi)/3) - shift
	r2 := m*math.Cos((phi+4*math.Pi)/3) - shift
	return uniqueSortedRoots([]Real{r0, r1, r2}, 1e-10)
}

func quarticTol(c4, c3, c2, c1, c0 Real) Real {
	m := math.Abs(c4)
	for _, v := range []Real{c3, c2, c1, c0} {
		av := math.Abs(v)
		if av > m {
			m = av
		}
	}
	if m < 1 {
		m = 1
	}
	return 1e-10 * m
}

func bisectQuarticRoot(c4, c3, c2, c1, c0, a, b Real) Real {
	fa := polyEvalQuartic(c4, c3, c2, c1, c0, a)
	fb := polyEvalQuartic(c4, c3, c2, c1, c0, b)
	if fa == 0 {
		return a
	}
	if fb == 0 {
		return b
	}
	lo, hi := a, b
	for i := 0; i < 80; i++ {
		m := 0.5 * (lo + hi)
		fm := polyEvalQuartic(c4, c3, c2, c1, c0, m)
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

func quarticRootsInInterval(c4, c3, c2, c1, c0, tmin, tmax Real) []Real {
	tol := quarticTol(c4, c3, c2, c1, c0)

	crit := realRootsCubic(4*c4, 3*c3, 2*c2, c1)
	points := make([]Real, 0, len(crit)+2)
	points = append(points, tmin)
	for _, x := range crit {
		if x > tmin+1e-12 && x < tmax-1e-12 {
			points = append(points, x)
		}
	}
	points = append(points, tmax)
	points = uniqueSortedRoots(points, 1e-12)

	roots := make([]Real, 0, 4)

	for _, p := range points {
		if math.Abs(float64(polyEvalQuartic(c4, c3, c2, c1, c0, p))) <= float64(tol) {
			roots = append(roots, p)
		}
	}

	for i := 0; i+1 < len(points); i++ {
		a, b := points[i], points[i+1]
		if b-a <= 1e-12 {
			continue
		}
		fa := polyEvalQuartic(c4, c3, c2, c1, c0, a)
		fb := polyEvalQuartic(c4, c3, c2, c1, c0, b)
		if (fa < 0 && fb > 0) || (fa > 0 && fb < 0) {
			roots = append(roots, bisectQuarticRoot(c4, c3, c2, c1, c0, a, b))
		}
	}

	roots = uniqueSortedRoots(roots, 1e-9)
	out := make([]Real, 0, len(roots))
	for _, x := range roots {
		if x >= tmin-1e-10 && x <= tmax+1e-10 {
			out = append(out, x)
		}
	}
	return out
}

func intersectRaySpheritorus(O Point4, D Vector4, t *Spheritorus) (hit objectHit, ok bool) {
	Op := Vector4{O.X - t.Center.X, O.Y - t.Center.Y, O.Z - t.Center.Z, O.W - t.Center.W}
	Ol := t.RT.MulVec(Op)
	Dl := t.RT.MulVec(D)

	R := t.MajorRadius
	r := t.MinorRadius
	sx, sy, sz, sw := t.Scale.X, t.Scale.Y, t.Scale.Z, t.Scale.W

	// Conservative local bounding box for the torus, used to restrict the quartic solve.
	hx := sx * (R + r)
	hy := sy * r
	hz := sz * r
	hw := sw * (R + r)

	box, ok := boxIntervalLocal(Ol, Dl, hx, hy, hz, hw)
	if !ok {
		return objectHit{}, false
	}

	A := Vector4{Ol.X / sx, Ol.Y / sy, Ol.Z / sz, Ol.W / sw}
	B := Vector4{Dl.X / sx, Dl.Y / sy, Dl.Z / sz, Dl.W / sw}

	M := A.Dot(A) + R*R - r*r
	N := 2 * A.Dot(B)
	P := B.Dot(B)

	Q := A.X*A.X + A.W*A.W
	Tw := 2 * (A.X*B.X + A.W*B.W)
	U := B.X*B.X + B.W*B.W

	c4 := P * P
	c3 := 2 * P * N
	c2 := N*N + 2*P*M - 4*R*R*U
	c1 := 2*N*M - 4*R*R*Tw
	c0 := M*M - 4*R*R*Q

	roots := quarticRootsInInterval(c4, c3, c2, c1, c0, box.t0, box.t1)
	if len(roots) == 0 {
		return objectHit{}, false
	}

	points := make([]Real, 0, len(roots)+2)
	points = append(points, box.t0)
	points = append(points, roots...)
	points = append(points, box.t1)
	points = uniqueSortedRoots(points, 1e-9)

	const eps = 1e-12
	inv := false
	var tHit Real
	found := false

	for i := 0; i+1 < len(points); i++ {
		a, b := points[i], points[i+1]
		if b-a <= 1e-10 {
			continue
		}
		m := 0.5 * (a + b)
		if polyEvalQuartic(c4, c3, c2, c1, c0, m) <= 0 {
			// Ray starts inside this interval -> exiting hit.
			if a <= eps && b > eps {
				inv = true
				tHit = b
				found = true
				break
			}
			// First interior interval ahead of the ray -> entering hit.
			if a > eps {
				inv = false
				tHit = a
				found = true
				break
			}
		}
	}

	if !found || tHit <= eps {
		return objectHit{}, false
	}

	Pq := Ol.Add(Dl.Mul(tHit))
	Sx := Pq.X / sx
	Sy := Pq.Y / sy
	Sz := Pq.Z / sz
	Sw := Pq.W / sw

	sum := Sx*Sx + Sy*Sy + Sz*Sz + Sw*Sw
	Aval := sum + R*R - r*r
	Bval := sum - R*R - r*r

	gradS := Vector4{
		4 * Sx * Bval,
		4 * Sy * Aval,
		4 * Sz * Aval,
		4 * Sw * Bval,
	}
	gradL := Vector4{
		gradS.X / sx,
		gradS.Y / sy,
		gradS.Z / sz,
		gradS.W / sw,
	}
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
