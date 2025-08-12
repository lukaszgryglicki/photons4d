package photons4d

import (
	"fmt"
	"math"
)

// HyperSphere: start from a unit 4D sphere in local space,
// scale along local axes (Radii.X/Y/Z/W), rotate by R about origin,
// then translate to Center.
type HyperSphere struct {
	Center Point4
	R      Mat4 // local->world rotation
	RT     Mat4 // world->local rotation (R^T)
	Radii  Vector4

	// Material (per-channel):
	Color   RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	// cached
	AABBMin  Point4
	AABBMax  Point4
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func NewHyperSphere(
	center Point4,
	radii Vector4, // semi-axes after scaling (≥0)
	angles Rot4,
	color, reflectivity, refractivity, ior RGB,
) (*HyperSphere, error) {
	if !(radii.X > 0 && radii.Y > 0 && radii.Z > 0 && radii.W > 0) {
		return nil, fmt.Errorf("hypersphere radii must be >0 on all axes, got %+v", radii)
	}
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	type ch struct {
		n    string
		r, t Real
	}
	chk := []ch{
		{"R", reflectivity.R, refractivity.R},
		{"G", reflectivity.G, refractivity.G},
		{"B", reflectivity.B, refractivity.B},
	}
	for _, c := range chk {
		if !in01(c.r) || !in01(c.t) {
			return nil, fmt.Errorf("reflect/refract must be in [0,1]; channel %s got reflect=%.6g refract=%.6g", c.n, c.r, c.t)
		}
		if c.r+c.t > 1+1e-12 {
			return nil, fmt.Errorf("per-channel reflect+refract must be ≤1; channel %s got %.6g", c.n, c.r+c.t)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	R := rotFromAngles(angles)
	hs := &HyperSphere{
		Center: center,
		R:      R,
		RT:     R.Transpose(),
		Radii:  radii,

		Color:   color,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}

	// AABB for rotated ellipsoid: project axis extents
	abs := func(x Real) Real {
		if x < 0 {
			return -x
		}
		return x
	}
	axis := radii
	extent := func(row int) (min, max Real) {
		off := abs(R.M[row][0])*axis.X + abs(R.M[row][1])*axis.Y + abs(R.M[row][2])*axis.Z + abs(R.M[row][3])*axis.W
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
	hs.AABBMin = Point4{minX, minY, minZ, minW}
	hs.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	// material caches
	hs.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	hs.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	hs.colorArr = [3]Real{color.R, color.G, color.B}
	hs.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - hs.refl[i] - hs.refr[i]
		if p < 0 {
			p = 0
		}
		hs.pAbs[i] = p
		hs.iorInv[i] = 1 / hs.iorArr[i]
		n := hs.iorArr[i]
		r0 := (n - 1) / (n + 1)
		hs.f0[i] = r0 * r0
	}

	DebugLog("Created hypersphere: %+v", hs)
	return hs, nil
}

// Ray/ellipsoid intersection via unit-sphere transform.
// Local: y = RT*(x - C). Unit-sphere coords: s = y / Radii.
// Solve ||s_o + t s_d||^2 = 1.
func intersectRayHyperSphere(O Point4, D Vector4, h *HyperSphere) (hit objectHit, ok bool) {
	// world -> local
	Op := Vector4{O.X - h.Center.X, O.Y - h.Center.Y, O.Z - h.Center.Z, O.W - h.Center.W}
	Ol := h.RT.MulVec(Op)
	Dl := h.RT.MulVec(D)

	// scale to unit-sphere space
	ax := h.Radii
	Os := Vector4{Ol.X / ax.X, Ol.Y / ax.Y, Ol.Z / ax.Z, Ol.W / ax.W}
	Ds := Vector4{Dl.X / ax.X, Dl.Y / ax.Y, Dl.Z / ax.Z, Dl.W / ax.W}

	a := Ds.Dot(Ds)
	b := 2 * Os.Dot(Ds)
	c := Os.Dot(Os) - 1
	disc := b*b - 4*a*c
	if disc < 0 {
		return objectHit{}, false
	}
	sqrtD := math.Sqrt(disc)
	inv2a := 1 / (2 * a)
	t0 := (-b - sqrtD) * inv2a
	t1 := (-b + sqrtD) * inv2a
	if t0 > t1 {
		t0, t1 = t1, t0
	}

	const eps = 1e-12
	inv := false
	t := t0
	if t <= eps {
		t = t1
		inv = true // origin was inside ⇒ first positive is the exit
	}
	if t <= eps {
		return objectHit{}, false
	}

	// Hit point in unit-sphere coords
	S := Os.Add(Ds.Mul(t)) // on unit sphere
	// World-space normal: N ∝ R * (S / Radii)  (then normalize)
	nl := Vector4{S.X / ax.X, S.Y / ax.Y, S.Z / ax.Z, S.W / ax.W}
	Nw := h.R.MulVec(nl).Norm()

	return objectHit{t: t, Nw: Nw, hs: h, inv: inv}, true
}
