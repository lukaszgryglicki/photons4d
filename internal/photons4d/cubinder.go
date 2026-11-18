package photons4d

import (
	"fmt"
	"math"
)

// Cubinder (local space):
//
//	x²/rx² + y²/ry² <= 1
//	|z| <= hz
//	|w| <= hw
//
// with Scale = (rx, ry, hz, hw), where:
//   - (rx, ry) are the semi-axes of the 2D ellipse in the XY plane
//   - (hz, hw) are the half-sizes of the rectangle in the ZW plane
//
// So this is D² × I² in 4D: the Cartesian product of a disk/ellipse and a
// square/rectangle. It lies between the 8-cell (I² × I²) and the duocylinder
// (D² × D²). The local object is then rotated by R and translated to Center.
type Cubinder struct {
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

func NewCubinder(center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*Cubinder, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("cubinder scale must be >0 on all axes, got %+v", scale)
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
	c := &Cubinder{
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
	// For a world-axis row a=(ax,ay,az,aw), support is:
	//   sqrt((rx*ax)^2 + (ry*ay)^2) + hz*|az| + hw*|aw|
	axisExtent := func(row int) (minV, maxV Real) {
		ax, ay, az, aw := R.M[row][0], R.M[row][1], R.M[row][2], R.M[row][3]
		off := math.Sqrt((scale.X*ax)*(scale.X*ax)+(scale.Y*ay)*(scale.Y*ay)) + scale.Z*math.Abs(az) + scale.W*math.Abs(aw)
		var ctr Real
		switch row {
		case 0:
			ctr = center.X
		case 1:
			ctr = center.Y
		case 2:
			ctr = center.Z
		default:
			ctr = center.W
		}
		return ctr - off, ctr + off
	}
	minX, maxX := axisExtent(0)
	minY, maxY := axisExtent(1)
	minZ, maxZ := axisExtent(2)
	minW, maxW := axisExtent(3)
	c.AABBMin = Point4{minX, minY, minZ, minW}
	c.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	c.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	c.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	c.colorArr = [3]Real{color.R, color.G, color.B}
	c.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	c.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - c.refl[i] - c.refr[i] - c.diff[i]
		if p < 0 {
			p = 0
		}
		c.pAbs[i] = p
		c.iorInv[i] = 1 / c.iorArr[i]
		n := c.iorArr[i]
		r0 := (n - 1) / (n + 1)
		c.f0[i] = r0 * r0
	}

	DebugLog("Created cubinder: %+v", c)
	return c, nil
}

func intersectRayCubinder(O Point4, D Vector4, c *Cubinder) (hit objectHit, ok bool) {
	Op := Vector4{O.X - c.Center.X, O.Y - c.Center.Y, O.Z - c.Center.Z, O.W - c.Center.W}
	Ol := c.RT.MulVec(Op)
	Dl := c.RT.MulVec(D)

	rx, ry, hz, hw := c.Scale.X, c.Scale.Y, c.Scale.Z, c.Scale.W

	ivXY, okXY := intervalEllipse2DLocal(Ol.X, Ol.Y, Dl.X, Dl.Y, rx, ry)
	if !okXY {
		return objectHit{}, false
	}

	const eps = 1e-12

	// Piece of the boundary met on entry and on exit: 0 = round side,
	// 1 = the Z faces, 2 = the W faces.
	tEnter, tExit := ivXY.t0, ivXY.t1
	enter, exit := 0, 0
	slab := func(o, d, h Real, piece int) bool {
		if math.Abs(d) < eps {
			return !(o < -h || o > h)
		}
		t0 := (-h - o) / d
		t1 := (h - o) / d
		if t0 > t1 {
			t0, t1 = t1, t0
		}
		if t0 > tEnter {
			tEnter, enter = t0, piece
		}
		if t1 < tExit {
			tExit, exit = t1, piece
		}
		return true
	}
	if !slab(Ol.Z, Dl.Z, hz, 1) || !slab(Ol.W, Dl.W, hw, 2) {
		return objectHit{}, false
	}
	if tEnter > tExit || (tExit <= eps && tEnter <= eps) {
		return objectHit{}, false
	}

	inv := tEnter <= eps
	t, piece := tEnter, enter
	if inv {
		t, piece = tExit, exit
	}
	if t <= eps {
		return objectHit{}, false
	}

	P := Ol.Add(Dl.Mul(t))
	var Nl Vector4
	switch piece {
	case 1:
		Nl = Vector4{0, 0, math.Copysign(1, P.Z), 0}
	case 2:
		Nl = Vector4{0, 0, 0, math.Copysign(1, P.W)}
	default:
		Nl = Vector4{P.X / (rx * rx), P.Y / (ry * ry), 0, 0}
		if Nl.Len() == 0 {
			return objectHit{}, false
		}
	}

	Nw := c.R.MulVec(Nl).Norm()
	return objectHit{t: t, Nw: Nw, mat: c, inv: inv}, true
}

func (c *Cubinder) PAbsCh(i int) Real   { return c.pAbs[i] }
func (c *Cubinder) DiffCh(i int) Real   { return c.diff[i] }
func (c *Cubinder) ColorCh(i int) Real  { return c.colorArr[i] }
func (c *Cubinder) F0Ch(i int) Real     { return c.f0[i] }
func (c *Cubinder) ReflCh(i int) Real   { return c.refl[i] }
func (c *Cubinder) RefrCh(i int) Real   { return c.refr[i] }
func (c *Cubinder) IORCh(i int) Real    { return c.iorArr[i] }
func (c *Cubinder) IORInvCh(i int) Real { return c.iorInv[i] }
