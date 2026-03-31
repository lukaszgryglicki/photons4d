package photons4d

import (
	"fmt"
	"math"
	"math/cmplx"
)

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

	// Conservative analytic AABB from the rotated bounding box of the
	// axis-scaled canonical duotorus.
	//
	// In local coordinates:
	//   x,y are bounded by (Rxy + r)
	//   z,w are bounded by (Rzw + r)
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

func polyTrimCF(p []Real) []Real {
	n := len(p)
	for n > 1 && math.Abs(float64(p[n-1])) < 1e-15 {
		n--
	}
	out := make([]Real, n)
	copy(out, p[:n])
	return out
}

func polyAddCF(a, b []Real) []Real {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]Real, n)
	for i := range a {
		out[i] += a[i]
	}
	for i := range b {
		out[i] += b[i]
	}
	return polyTrimCF(out)
}

func polySubCF(a, b []Real) []Real {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]Real, n)
	for i := range a {
		out[i] += a[i]
	}
	for i := range b {
		out[i] -= b[i]
	}
	return polyTrimCF(out)
}

func polyScaleCF(a []Real, s Real) []Real {
	out := make([]Real, len(a))
	for i := range a {
		out[i] = a[i] * s
	}
	return polyTrimCF(out)
}

func polyMulCF(a, b []Real) []Real {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	out := make([]Real, len(a)+len(b)-1)
	for i := range a {
		for j := range b {
			out[i+j] += a[i] * b[j]
		}
	}
	return polyTrimCF(out)
}

func polyEvalCF(p []Real, x Real) Real {
	acc := Real(0)
	for i := len(p) - 1; i >= 0; i-- {
		acc = acc*x + p[i]
	}
	return acc
}

func polyEvalCFCplx(p []Real, z complex128) complex128 {
	acc := complex(0.0, 0.0)
	for i := len(p) - 1; i >= 0; i-- {
		acc = acc*z + complex(p[i], 0)
	}
	return acc
}

func durandKernerRealRootsCF(coeffs []Real) []Real {
	p := polyTrimCF(coeffs)
	deg := len(p) - 1
	switch deg {
	case -1, 0:
		return nil
	case 1:
		return []Real{-p[0] / p[1]}
	case 2:
		return solveQuadraticReal(p[2], p[1], p[0])
	case 3:
		return realRootsCubic(p[3], p[2], p[1], p[0])
	}

	lead := p[deg]
	monic := make([]Real, len(p))
	for i := range p {
		monic[i] = p[i] / lead
	}

	radius := 1.0
	for i := 0; i < deg; i++ {
		v := 1 + math.Abs(float64(monic[i]))
		if v > radius {
			radius = v
		}
	}

	roots := make([]complex128, deg)
	for k := 0; k < deg; k++ {
		ang := 2 * math.Pi * float64(k) / float64(deg)
		roots[k] = complex(radius*math.Cos(ang), radius*math.Sin(ang))
	}

	const maxIter = 200
	for it := 0; it < maxIter; it++ {
		converged := true
		for i := range roots {
			zi := roots[i]
			denom := complex(1.0, 0.0)
			for j := range roots {
				if i == j {
					continue
				}
				diff := zi - roots[j]
				if cmplx.Abs(diff) < 1e-14 {
					diff += complex(1e-12*float64(i+1), 1e-12*float64(j+1))
				}
				denom *= diff
			}
			delta := polyEvalCFCplx(monic, zi) / denom
			roots[i] = zi - delta
			if cmplx.Abs(delta) > 1e-12 {
				converged = false
			}
		}
		if converged {
			break
		}
	}

	scale := 1.0
	for _, v := range p {
		av := math.Abs(float64(v))
		if av > scale {
			scale = av
		}
	}

	out := make([]Real, 0, deg)
	for _, z := range roots {
		if math.Abs(imag(z)) > 1e-7 {
			continue
		}
		x := Real(real(z))
		if math.Abs(float64(polyEvalCF(p, x))) > 1e-6*scale {
			continue
		}
		out = append(out, x)
	}
	return uniqueSortedRoots(out, 1e-8)
}

func duotorusImplicitLocal(p Vector4, scale Vector4, Rxy, Rzw, r Real) Real {
	u := p.X / scale.X
	v := p.Y / scale.Y
	z := p.Z / scale.Z
	s := p.W / scale.W

	a := math.Hypot(u, v)
	b := math.Hypot(z, s)

	da := a - Rxy
	db := b - Rzw
	return da*da + db*db - r*r
}

func intersectRayDuotorus(O Point4, D Vector4, tor *Duotorus) (hit objectHit, ok bool) {
	Op := Vector4{O.X - tor.Center.X, O.Y - tor.Center.Y, O.Z - tor.Center.Z, O.W - tor.Center.W}
	Ol := tor.RT.MulVec(Op)
	Dl := tor.RT.MulVec(D)

	Rxy := tor.MajorRadiusXY
	Rzw := tor.MajorRadiusZW
	r := tor.MinorRadius
	sx, sy, sz, sw := tor.Scale.X, tor.Scale.Y, tor.Scale.Z, tor.Scale.W

	hx := sx * (Rxy + r)
	hy := sy * (Rxy + r)
	hz := sz * (Rzw + r)
	hw := sw * (Rzw + r)

	box, ok := boxIntervalLocal(Ol, Dl, hx, hy, hz, hw)
	if !ok {
		return objectHit{}, false
	}

	// A(t) = (x/sx)^2 + (y/sy)^2
	a0 := Ol.X*Ol.X/(sx*sx) + Ol.Y*Ol.Y/(sy*sy)
	a1 := 2 * (Ol.X*Dl.X/(sx*sx) + Ol.Y*Dl.Y/(sy*sy))
	a2 := Dl.X*Dl.X/(sx*sx) + Dl.Y*Dl.Y/(sy*sy)
	A := []Real{a0, a1, a2}

	// B(t) = (z/sz)^2 + (w/sw)^2
	b0 := Ol.Z*Ol.Z/(sz*sz) + Ol.W*Ol.W/(sw*sw)
	b1 := 2 * (Ol.Z*Dl.Z/(sz*sz) + Ol.W*Dl.W/(sw*sw))
	b2 := Dl.Z*Dl.Z/(sz*sz) + Dl.W*Dl.W/(sw*sw)
	B := []Real{b0, b1, b2}

	// C = A + B + Rxy^2 + Rzw^2 - r^2
	C := polyAddCF(A, B)
	C[0] += Rxy*Rxy + Rzw*Rzw - r*r

	// Squaring twice:
	//   (sqrt(A)-Rxy)^2 + (sqrt(B)-Rzw)^2 = r^2
	// => C = 2 Rxy sqrt(A) + 2 Rzw sqrt(B)
	// => (C^2 - 4Rxy^2 A - 4Rzw^2 B)^2 = 64 Rxy^2 Rzw^2 A B
	C2 := polyMulCF(C, C)
	E := polySubCF(polySubCF(C2, polyScaleCF(A, 4*Rxy*Rxy)), polyScaleCF(B, 4*Rzw*Rzw))
	P := polySubCF(polyMulCF(E, E), polyScaleCF(polyMulCF(A, B), 64*Rxy*Rzw*Rxy*Rzw))

	roots := durandKernerRealRootsCF(P)
	if len(roots) == 0 {
		return objectHit{}, false
	}

	origTol := 1e-7 * (1 + Rxy + Rzw + r)
	validRoots := make([]Real, 0, len(roots))
	for _, x := range roots {
		if x < box.t0-1e-9 || x > box.t1+1e-9 {
			continue
		}
		p := Ol.Add(Dl.Mul(x))
		if math.Abs(float64(duotorusImplicitLocal(p, tor.Scale, Rxy, Rzw, r))) <= float64(origTol) {
			validRoots = append(validRoots, x)
		}
	}
	if len(validRoots) == 0 {
		return objectHit{}, false
	}
	validRoots = uniqueSortedRoots(validRoots, 1e-8)

	insideAt := func(t Real) bool {
		p := Ol.Add(Dl.Mul(t))
		return duotorusImplicitLocal(p, tor.Scale, Rxy, Rzw, r) <= 1e-9
	}

	points := make([]Real, 0, len(validRoots)+2)
	points = append(points, box.t0)
	points = append(points, validRoots...)
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
		if insideAt(m) {
			if a <= eps && b > eps {
				inv = true
				tHit = b
				found = true
				break
			}
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
	u := Pq.X / sx
	v := Pq.Y / sy
	z := Pq.Z / sz
	s := Pq.W / sw

	a := math.Hypot(u, v)
	b := math.Hypot(z, s)

	const gEps = 1e-12
	ga := Real(0)
	gb := Real(0)
	if a > gEps {
		ga = 2 * (a - Rxy) / a
	}
	if b > gEps {
		gb = 2 * (b - Rzw) / b
	}

	Nl := Vector4{
		ga * Pq.X / (sx * sx),
		ga * Pq.Y / (sy * sy),
		gb * Pq.Z / (sz * sz),
		gb * Pq.W / (sw * sw),
	}
	if Nl.Len() == 0 {
		if math.Abs(Pq.X)+math.Abs(Pq.Y) >= math.Abs(Pq.Z)+math.Abs(Pq.W) {
			Nl = Vector4{1, 0, 0, 0}
		} else {
			Nl = Vector4{0, 0, 1, 0}
		}
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
