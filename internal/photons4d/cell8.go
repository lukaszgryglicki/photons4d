package photons4d

import (
	"fmt"
	"math"
)

// Cell8: axis-aligned box in local space, then rotated about origin and translated to Center.
// 'Half' are half-lengths along the local X,Y,Z,W axes.
// pAbs = 1 - (reflect + refract) (per channel). That’s your absorption knob.
// The non-absorption budget (avail = reflect + refract) gets split each hit by
// Fresnel–Schlick (angle-dependent), biased by your reflect/refract to steer it.
// color tints whatever survives (both reflected and refracted) per interaction.
// ior (per channel) sets Fresnel strength and dispersion (B > G > R ⇒ stronger “prism”).

type Cell8 struct {
	Center Point4
	Half   Vector4 // half-sizes: Lx/2, Ly/2, Lz/2, Lw/2
	R      Mat4    // local->world rotation
	RT     Mat4    // world->local rotation (R^T, since R is orthonormal)

	// Material (per-channel):
	Color   RGB // tint/filter applied to reflected & refracted energy
	Diffuse RGB // Lambertian probability mass per channel
	Reflect RGB // fraction reflected per channel (0..1)
	Refract RGB // fraction refracted per channel (0..1)
	IOR     RGB // index of refraction inside cell8 per channel

	// cached
	Normals  [4]Vector4 // world-space unit normals for +X,+Y,+Z,+W faces
	AABBMin  Point4     // conservative 4D AABB (min)
	AABBMax  Point4     // conservative 4D AABB (max)
	refl     [3]Real    // [R,G,B]
	refr     [3]Real    // [R,G,B]
	colorArr [3]Real    // [R,G,B]
	diff     [3]Real    // [R,G,B]
	iorArr   [3]Real    // [R,G,B]
	iorInv   [3]Real    // [1/R, 1/G, 1/B]
	pAbs     [3]Real    // 1 - refl - refr (clamped ≥ 0)
	f0       [3]Real    // Schlick F0 per channel = ((ior-1)/(ior+1))^2
}

// NewCell8 constructs a cell8 and precomputes caches.
func NewCell8(
	center Point4,
	scale Vector4, // full edge lengths
	angles Rot4, // radians
	color, diffuse, reflectivity, refractivity, ior RGB,
) (*Cell8, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("cell8 scale must be >0 on all axes, got %+v", scale)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }

	type ch struct {
		n       string
		r, t, d Real
	}
	chk := []ch{
		{"R", reflectivity.R, refractivity.R, diffuse.R},
		{"G", reflectivity.G, refractivity.G, diffuse.G},
		{"B", reflectivity.B, refractivity.B, diffuse.B},
	}
	for _, c := range chk {
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
	hc := Cell8{
		Center: center,
		Half:   Vector4{scale.X * 0.5, scale.Y * 0.5, scale.Z * 0.5, scale.W * 0.5},
		R:      R,
		RT:     R.Transpose(),

		Color:   color,
		Diffuse: diffuse,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}

	// world normals = columns of R (already unit)
	hc.Normals[0] = Vector4{R.M[0][0], R.M[1][0], R.M[2][0], R.M[3][0]}
	hc.Normals[1] = Vector4{R.M[0][1], R.M[1][1], R.M[2][1], R.M[3][1]}
	hc.Normals[2] = Vector4{R.M[0][2], R.M[1][2], R.M[2][2], R.M[3][2]}
	hc.Normals[3] = Vector4{R.M[0][3], R.M[1][3], R.M[2][3], R.M[3][3]}

	// precompute AABB by projecting extents
	abs := func(x Real) Real {
		if x < 0 {
			return -x
		}
		return x
	}
	extent := func(row int) (min, max Real) {
		off := abs(R.M[row][0])*hc.Half.X + abs(R.M[row][1])*hc.Half.Y + abs(R.M[row][2])*hc.Half.Z + abs(R.M[row][3])*hc.Half.W
		var c Real
		switch row {
		case 0:
			c = hc.Center.X
		case 1:
			c = hc.Center.Y
		case 2:
			c = hc.Center.Z
		default:
			c = hc.Center.W
		}
		return c - off, c + off
	}
	minX, maxX := extent(0)
	minY, maxY := extent(1)
	minZ, maxZ := extent(2)
	minW, maxW := extent(3)
	hc.AABBMin = Point4{minX, minY, minZ, minW}
	hc.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	// material arrays
	hc.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	hc.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	hc.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	hc.colorArr = [3]Real{color.R, color.G, color.B}
	hc.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - hc.refl[i] - hc.refr[i] - hc.diff[i]
		if p < 0 {
			p = 0
		}
		hc.pAbs[i] = p
		hc.iorInv[i] = 1 / hc.iorArr[i]
		// Fresnel-Schlick F0 for air(1.0) <-> material(ior)
		n := hc.iorArr[i]
		r0 := (n - 1) / (n + 1)
		hc.f0[i] = r0 * r0
	}

	// DebugLog("Created Cell8: Center=%+v, Half=%+v, R=%+v, Color=%+v, Reflectivity=%+v, Refractivity=%+v, IOR=%+v", hc.Center, hc.Half, hc.R, hc.Color, hc.Reflect, hc.Refract, hc.IOR)
	DebugLog("Created 8-cell: %+v", hc)
	return &hc, nil
}

func intersectRayCell8(O Point4, D Vector4, h *Cell8) (hit objectHit, ok bool) {
	// world -> local
	Op := Vector4{O.X - h.Center.X, O.Y - h.Center.Y, O.Z - h.Center.Z, O.W - h.Center.W}
	Ol := h.RT.MulVec(Op)
	Dl := h.RT.MulVec(D)

	// Use real infinities for clarity.
	ninf := math.Inf(-1)
	inf := math.Inf(1)
	tmin, tmax := ninf, inf
	enterAxis, enterSign := -1, 0
	exitAxis, exitSign := -1, 0

	axis := func(o, d, half Real, ax int) (bool, Real, Real, int, int, int, int) {
		const eps = 1e-12
		if math.Abs(d) < eps {
			if math.Abs(o) > half {
				return false, 0, 0, 0, 0, 0, 0
			}
			// Inside this slab and parallel: do not tighten; preserve axis id.
			// Use +1 as a neutral sign so later flip logic is consistent.
			return true, ninf, inf, ax, +1, ax, +1
		}
		t1 := (-half - o) / d
		t2 := (half - o) / d
		sgnEnter := -1
		sgnExit := +1
		if t1 > t2 {
			t1, t2 = t2, t1
			sgnEnter, sgnExit = +1, -1
		}
		return true, t1, t2, ax, sgnEnter, ax, sgnExit
	}

	okX, a1, a2, eAx1, eSg1, xAx1, xSg1 := axis(Ol.X, Dl.X, h.Half.X, 0)
	if !okX {
		return objectHit{}, false
	}
	okY, b1, b2, eAx2, eSg2, xAx2, xSg2 := axis(Ol.Y, Dl.Y, h.Half.Y, 1)
	if !okY {
		return objectHit{}, false
	}
	okZ, c1, c2, eAx3, eSg3, xAx3, xSg3 := axis(Ol.Z, Dl.Z, h.Half.Z, 2)
	if !okZ {
		return objectHit{}, false
	}
	okW, d1, d2, eAx4, eSg4, xAx4, xSg4 := axis(Ol.W, Dl.W, h.Half.W, 3)
	if !okW {
		return objectHit{}, false
	}

	type slab struct {
		t1, t2             Real
		eAx, eSg, xAx, xSg int
	}
	sl := []slab{
		{a1, a2, eAx1, eSg1, xAx1, xSg1},
		{b1, b2, eAx2, eSg2, xAx2, xSg2},
		{c1, c2, eAx3, eSg3, xAx3, xSg3},
		{d1, d2, eAx4, eSg4, xAx4, xSg4},
	}
	// Small tolerance to make face choice deterministic under ties.
	const tie = 1e-15
	for _, s := range sl {
		if s.t1 > tmin+tie {
			tmin, enterAxis, enterSign = s.t1, s.eAx, s.eSg
		}
		if s.t2 < tmax-tie {
			tmax, exitAxis, exitSign = s.t2, s.xAx, s.xSg
		}
	}
	if tmax < 0 || tmin > tmax {
		return objectHit{}, false
	}

	// More stable inside/outside test near faces.
	const epsInside = 1e-12
	inv := tmin < epsInside && tmax > epsInside

	// If nothing tightened (all slabs parallel), we cannot define a face.
	if enterAxis < 0 && exitAxis < 0 {
		return objectHit{}, false
	}

	var t Real
	ax, sg := 0, 0
	if !inv {
		t, ax, sg = tmin, enterAxis, enterSign // entering
	} else {
		t, ax, sg = tmax, exitAxis, exitSign // exiting
	}

	// pick cached world normal (unit), flip by sign
	Nw := h.Normals[ax]
	if sg < 0 {
		Nw = Nw.Mul(-1)
	}

	return objectHit{t: t, Nw: Nw, mat: h, inv: inv}, true
}

func (c *Cell8) PAbsCh(i int) Real   { return c.pAbs[i] }
func (c *Cell8) DiffCh(i int) Real   { return c.diff[i] }
func (c *Cell8) ColorCh(i int) Real  { return c.colorArr[i] }
func (c *Cell8) F0Ch(i int) Real     { return c.f0[i] }
func (c *Cell8) ReflCh(i int) Real   { return c.refl[i] }
func (c *Cell8) RefrCh(i int) Real   { return c.refr[i] }
func (c *Cell8) IORCh(i int) Real    { return c.iorArr[i] }
func (c *Cell8) IORInvCh(i int) Real { return c.iorInv[i] }
