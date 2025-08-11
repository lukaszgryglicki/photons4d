package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"image"
	"image/color/palette"
	"image/draw"
	"image/gif"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sync"
	"sync/atomic"
	"time"
)

type Real = float64

// Channel indices for readability.
const (
	ChR        = 0
	ChG        = 1
	ChB        = 2
	SceneResX  = 256
	SceneResY  = 256
	SceneResZ  = 256
	ProbeRays  = 100_000
	Spp        = 128 // samples per voxel target
	GIFOut     = "volume.gif"
	GIFDelay   = 5 // 100ths of a second per frame
	Gamma      = 0.75
	MaxBounces = 32
	NumShards  = 1024
)

type shardLocks struct{ mu [NumShards]sync.Mutex }

// ---------------------------------------------------------
// Geometry & color
// ---------------------------------------------------------

// Point4 represents a point in 4-dimensional space.
type Point4 struct {
	X, Y, Z, W Real
}

// Vector4 represents a direction (not a position) in 4D space.
type Vector4 struct {
	X, Y, Z, W Real
}

// 4×4 matrix (row-major)
type Mat4 struct {
	M [4][4]Real
}

// Angles in radians for rotations in coordinate planes.
type Rot4 struct {
	XY, XZ, XW, YZ, YW, ZW Real
}

// RGB stores color components; each should be in [0,1].
type RGB struct {
	R, G, B Real
}

// ---------------------------------------------------------
// Light
// ---------------------------------------------------------

// ConeLight4 is a 4D "spot" light: rays originate at Origin and
// are restricted to a cone around Direction with half-angle Angle (radians).
type ConeLight4 struct {
	Origin    Point4
	Direction Vector4 // should be unit-length; constructor normalizes
	Color     RGB     // each in [0,1]
	Angle     Real    // half-angle in radians, (0, π]
	// cached
	cosAngle    Real
	oneMinusCos Real
	colorSum    Real
	thrR, thrG  Real    // cumulative thresholds in [0,colorSum]
	U, V, W     Vector4 // orthonormal basis s.t. {U,V,W} ⟂ Direction and unit
}

// Hypercube4: axis-aligned box in local space, then rotated about origin and translated to Center.
// 'Half' are half-lengths along the local X,Y,Z,W axes.
type Hypercube4 struct {
	Center Point4
	Half   Vector4 // half-sizes: Lx/2, Ly/2, Lz/2, Lw/2
	R      Mat4    // local->world rotation
	RT     Mat4    // world->local rotation (R^T, since R is orthonormal)
	// Material (per-channel):
	Color   RGB // tint/filter applied to reflected & refracted energy
	Reflect RGB // fraction reflected per channel (0..1)
	Refract RGB // fraction refracted per channel (0..1)
	IOR     RGB // index of refraction inside cube per channel (e.g., {1.5,1.5,1.6})
	// cached
	Normals          [4]Vector4 // world-space unit face normals for axes X,Y,Z,W
	AABBMin, AABBMax Point4     // conservative world AABB
	refl             [3]Real
	refr             [3]Real
	pAbs             [3]Real
	ior              [3]Real
	iorInv           [3]Real
}

// ---------------------------------------------------------
// Scene: 3D cube at fixed W with flat voxel buffer
// ---------------------------------------------------------

// Scene3D stores a 3D volume (axis-aligned in X,Y,Z) embedded at W=Center.W.
// Voxels are in a flat buffer: len = Nx*Ny*Nz*3 (RGB).
type Scene3D struct {
	Center               Point4
	Width, Height, Depth Real
	Nx, Ny, Nz           int
	MaxBounces           int
	Buf                  []Real // flat: (((i*Ny)+j)*Nz + k)*3 + c
	Hypercubes           []*Hypercube4
	// cached
	MinX, MaxX, MinY, MaxY, MinZ, MaxZ Real
	InvSpanX, InvSpanY, InvSpanZ       Real
	StrideX, StrideY                   int // for flat index: base = i*StrideX + j*StrideY + k*3
}

// ---------- JSON config ----------

type SceneCfg struct {
	Center     Point4 `json:"center"`               // { "X":0,"Y":0,"Z":0,"W":0 }
	Width      Real   `json:"width"`                // e.g. 2.0
	Height     Real   `json:"height"`               // e.g. 2.0
	Depth      Real   `json:"depth"`                // e.g. 2.0
	MaxBounces int    `json:"maxBounces,omitempty"` // e.g. 10
}

type LightCfg struct {
	Origin    Point4  `json:"origin"`    // { "X":..., "Y":..., "Z":..., "W":... }
	Direction Vector4 `json:"direction"` // { "X":..., "Y":..., "Z":..., "W":... }
	Color     RGB     `json:"color"`     // { "R":..., "G":..., "B":... } in [0,1]
	AngleDeg  Real    `json:"angleDeg"`  // degrees half angle
}

type Config struct {
	SceneResX  int            `json:"sceneResX"`
	SceneResY  int            `json:"sceneResY"`
	SceneResZ  int            `json:"sceneResZ"`
	ProbeRays  int            `json:"probeRays"`          // trials per light for p_hit
	Spp        int            `json:"spp"`                // samples per voxel target
	GIFOut     string         `json:"gifOut"`             // output filename
	GIFDelay   int            `json:"gifDelay,omitempty"` // 100ths of a sec per frame; default 5
	Gamma      Real           `json:"gamma,omitempty"`    // tone map gamma; default 0.75
	Scene      SceneCfg       `json:"scene"`
	Lights     []LightCfg     `json:"lights"`
	Hypercubes []HypercubeCfg `json:"hypercubes,omitempty"`
}

// Rotation in degrees for JSON (friendlier than radians).
type Rot4Deg struct {
	XY Real `json:"xy"`
	XZ Real `json:"xz"`
	XW Real `json:"xw"`
	YZ Real `json:"yz"`
	YW Real `json:"yw"`
	ZW Real `json:"zw"`
}

type HypercubeCfg struct {
	Center Point4  `json:"center"`
	Size   Vector4 `json:"size"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`   // required
	Reflect RGB `json:"reflect"` // required
	Refract RGB `json:"refract"` // required
	IOR     RGB `json:"ior"`     // required (per-channel)
}

type cubeHit struct {
	t   Real    // param distance along ray
	Nw  Vector4 // world-space unit normal at hit
	hc  *Hypercube4
	inv bool // true if we were inside and are exiting
}

// Build validates and constructs the runtime object (no defaults).
func (hc HypercubeCfg) Build() (*Hypercube4, error) {
	rad := hc.RotDeg.Radians()

	// Optional sanity check: warn if all angles are ~zero.
	if math.Abs(rad.XY)+math.Abs(rad.XZ)+math.Abs(rad.XW)+
		math.Abs(rad.YZ)+math.Abs(rad.YW)+math.Abs(rad.ZW) < 1e-12 {
		debugLog("Hypercube rotation is ~zero; check JSON 'rotDeg' keys (xy,xz,xw,yz,yw,zw).")
	}

	return NewHypercube4(
		hc.Center,
		hc.Size,
		rad,
		hc.Color,
		hc.Reflect,
		hc.Refract,
		hc.IOR,
	)
}

func (r Rot4Deg) Radians() Rot4 {
	const k = math.Pi / 180
	return Rot4{
		XY: r.XY * k, XZ: r.XZ * k, XW: r.XW * k,
		YZ: r.YZ * k, YW: r.YW * k, ZW: r.ZW * k,
	}
}

func (sl *shardLocks) lock(idx int)   { sl.mu[idx&(NumShards-1)].Lock() }
func (sl *shardLocks) unlock(idx int) { sl.mu[idx&(NumShards-1)].Unlock() }

// Add lets you translate a Point4 by a Vector4.
func (p Point4) Add(v Vector4) Point4 {
	point := Point4{p.X + v.X, p.Y + v.Y, p.Z + v.Z, p.W + v.W}
	// debugLog("Adding Vector4 to Point4: (%f, %f, %f, %f) + (%f, %f, %f, %f) = (%f, %f, %f, %f)", p.X, p.Y, p.Z, p.W, v.X, v.Y, v.Z, v.W, point.X, point.Y, point.Z, point.W)
	return point
}

// Vector functions
func (a Vector4) Add(b Vector4) Vector4 { return Vector4{a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W} }
func (a Vector4) Sub(b Vector4) Vector4 { return Vector4{a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W} }
func (v Vector4) Mul(s Real) Vector4    { return Vector4{v.X * s, v.Y * s, v.Z * s, v.W * s} }

// Dot returns the dot product between two 4D vectors.
func (a Vector4) Dot(b Vector4) Real {
	dot := a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W
	// debugLog("Dot product: (%f, %f, %f, %f) . (%f, %f, %f, %f) = %f", a.X, a.Y, a.Z, a.W, b.X, b.Y, b.Z, b.W, dot)
	return dot
}

// Len returns the Euclidean length of the vector.
func (v Vector4) Len() Real {
	return math.Sqrt(v.Dot(v))
}

// Norm returns a unit-length version of the vector.
// If the vector is (near) zero, it returns the input unchanged.
func (v Vector4) Norm() Vector4 {
	l := v.Len()
	if l == 0 {
		return v
	}
	return Vector4{v.X / l, v.Y / l, v.Z / l, v.W / l}
}

func I4() Mat4 {
	return Mat4{M: [4][4]Real{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}}
}

func (A Mat4) Mul(B Mat4) Mat4 {
	var R Mat4
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				sum += A.M[r][k] * B.M[k][c]
			}
			R.M[r][c] = sum
		}
	}
	return R
}

func (A Mat4) Transpose() Mat4 {
	var R Mat4
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			R.M[r][c] = A.M[c][r]
		}
	}
	return R
}

func (A Mat4) MulVec(v Vector4) Vector4 {
	return Vector4{
		A.M[0][0]*v.X + A.M[0][1]*v.Y + A.M[0][2]*v.Z + A.M[0][3]*v.W,
		A.M[1][0]*v.X + A.M[1][1]*v.Y + A.M[1][2]*v.Z + A.M[1][3]*v.W,
		A.M[2][0]*v.X + A.M[2][1]*v.Y + A.M[2][2]*v.Z + A.M[2][3]*v.W,
		A.M[3][0]*v.X + A.M[3][1]*v.Y + A.M[3][2]*v.Z + A.M[3][3]*v.W,
	}
}

func rotXY(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][1] = c, -s
	M.M[1][0], M.M[1][1] = s, c
	return M
}
func rotXZ(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][2] = c, -s
	M.M[2][0], M.M[2][2] = s, c
	return M
}
func rotXW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][3] = c, -s
	M.M[3][0], M.M[3][3] = s, c
	return M
}
func rotYZ(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[1][1], M.M[1][2] = c, -s
	M.M[2][1], M.M[2][2] = s, c
	return M
}
func rotYW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[1][1], M.M[1][3] = c, -s
	M.M[3][1], M.M[3][3] = s, c
	return M
}
func rotZW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[2][2], M.M[2][3] = c, -s
	M.M[3][2], M.M[3][3] = s, c
	return M
}

// Compose rotation from angles.
// Order matters; here we apply right-to-left: first XY, then XZ, XW, YZ, YW, ZW.
func rotFromAngles(r Rot4) Mat4 {
	R := I4()
	R = rotZW(r.ZW).Mul(R)
	R = rotYW(r.YW).Mul(R)
	R = rotYZ(r.YZ).Mul(R)
	R = rotXW(r.XW).Mul(R)
	R = rotXZ(r.XZ).Mul(R)
	R = rotXY(r.XY).Mul(R)
	return R
}

func NewHypercube4(
	center Point4,
	size Vector4, // full edge lengths in world units
	angles Rot4, // rotation in radians
	color, reflectivity, refractivity, ior RGB,
) (*Hypercube4, error) {

	// Size must be positive on all axes
	if !(size.X > 0 && size.Y > 0 && size.Z > 0 && size.W > 0) {
		return nil, fmt.Errorf("hypercube size must be >0 on all axes, got %+v", size)
	}

	// Helper: in [0,1]
	in01 := func(x Real) bool { return x >= 0 && x <= 1 }

	// Validate [0,1] ranges for reflect/refract; and per-channel sum ≤ 1
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

	// IOR must be strictly positive
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	// Build rotation
	R := rotFromAngles(angles)

	hc := Hypercube4{
		Center: center,
		Half:   Vector4{size.X * 0.5, size.Y * 0.5, size.Z * 0.5, size.W * 0.5},
		R:      R,
		RT:     R.Transpose(),

		Color:   color,
		Reflect: reflectivity,
		Refract: refractivity,
		IOR:     ior,
	}

	// world normals = columns of R (already unit)
	hc.Normals[0] = Vector4{R.M[0][0], R.M[1][0], R.M[2][0], R.M[3][0]}
	hc.Normals[1] = Vector4{R.M[0][1], R.M[1][1], R.M[2][1], R.M[3][1]}
	hc.Normals[2] = Vector4{R.M[0][2], R.M[1][2], R.M[2][2], R.M[3][2]}
	hc.Normals[3] = Vector4{R.M[0][3], R.M[1][3], R.M[2][3], R.M[3][3]}

	// precompute material arrays
	hc.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	hc.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	hc.ior = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - hc.refl[i] - hc.refr[i]
		if p < 0 {
			p = 0
		}
		hc.pAbs[i] = p
		hc.iorInv[i] = 1 / hc.ior[i]
	}

	// static world AABB (conservative) of rotated 4D box projected into XYZ and W:
	// For each world axis α ∈ {X,Y,Z,W}, extents are sum_j |R_{αj}| * Half_j
	abs := func(x Real) Real {
		if x < 0 {
			return -x
		}
		return x
	}
	ext := func(row int) (min, max Real) {
		off := abs(R.M[row][0])*hc.Half.X + abs(R.M[row][1])*hc.Half.Y + abs(R.M[row][2])*hc.Half.Z + abs(R.M[row][3])*hc.Half.W
		c := []Real{hc.Center.X, hc.Center.Y, hc.Center.Z, hc.Center.W}[row]
		return c - off, c + off
	}
	minX, maxX := ext(0)
	minY, maxY := ext(1)
	minZ, maxZ := ext(2)
	minW, maxW := ext(3)
	hc.AABBMin = Point4{minX, minY, minZ, minW}
	hc.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	// debugLog("Created Hypercube4: Center=%+v, Half=%+v, R=%+v, Color=%+v, Reflectivity=%+v, Refractivity=%+v, IOR=%+v", hc.Center, hc.Half, hc.R, hc.Color.clamp01(), hc.Reflect.clamp01(), hc.Refract.clamp01(), hc.IOR)
	debugLog("Created Hypercube4: Center=%+v, Half=%+v, R=%+v, Color=%+v, Reflectivity=%+v, Refractivity=%+v, IOR=%+v", hc.Center, hc.Half, hc.R, hc.Color, hc.Reflect, hc.Refract, hc.IOR)
	return &hc, nil
}

// clamp01 clamps each channel to [0,1] (useful for validation).
func (c RGB) clamp01() RGB {
	cl := func(x Real) Real {
		if x < 0 {
			return 0
		}
		if x > 1 {
			return 1
		}
		return x
	}
	return RGB{cl(c.R), cl(c.G), cl(c.B)}
}

// Build an orthonormal basis {U,V,W} in the 3D subspace orthogonal to A (unit).
func orthonormal3(a Vector4) (u, v, w Vector4) {
	// pick a helper not (near-)parallel to a
	h := Vector4{1, 0, 0, 0}
	if math.Abs(a.X) > 0.9 {
		h = Vector4{0, 1, 0, 0}
	}
	u = (h.Sub(a.Mul(h.Dot(a)))).Norm()
	// any vector not collinear with u or a
	h2 := Vector4{0, 0, 1, 0}
	v = (h2.Sub(a.Mul(h2.Dot(a))).Sub(u.Mul(h2.Dot(u)))).Norm()
	// third via Gram–Schmidt
	// pick a different seed to reduce degeneracy risk
	h3 := Vector4{0, 0, 0, 1}
	w = (h3.Sub(a.Mul(h3.Dot(a))).Sub(u.Mul(h3.Dot(u))).Sub(v.Mul(h3.Dot(v)))).Norm()
	return
}

// NewConeLight4 constructs a cone light, normalizing direction and clamping color.
// Angle must be in (0, π]. Returns an error if the direction is zero or angle is invalid.
func NewConeLight4(origin Point4, dir Vector4, color RGB, angle Real) (*ConeLight4, error) {
	if angle <= 0 || angle > math.Pi {
		return nil, errors.New("angle must be in (0, π]")
	}
	n := dir.Norm()
	if n.Len() == 0 {
		return nil, errors.New("direction must be non-zero")
	}
	cosA := math.Cos(angle)
	c := color.clamp01()
	L := &ConeLight4{
		Origin: origin, Direction: n, Color: c, Angle: angle,
		cosAngle:    cosA,
		oneMinusCos: 1 - cosA,
		colorSum:    c.R + c.G + c.B,
	}
	L.thrR = c.R                    // [0, colorSum)
	L.thrG = c.R + c.G              // [0, colorSum)
	L.U, L.V, L.W = orthonormal3(n) // build once
	return L, nil
}

func randNormal(rng *rand.Rand) Real {
	// Box–Muller
	u1 := rng.Float64()
	u2 := rng.Float64()
	r := math.Sqrt(-2 * math.Log(math.Max(u1, 1e-12)))
	return r * math.Cos(2*math.Pi*u2)
}

func (l *ConeLight4) SampleDir(rng *rand.Rand) Vector4 {
	axis := l.Direction // already unit
	// Sample orthonormal direction in the 3D subspace orthogonal to axis
	var u Vector4
	for {
		r := Vector4{randNormal(rng), randNormal(rng), randNormal(rng), randNormal(rng)}
		ortho := r.Sub(axis.Mul(r.Dot(axis)))
		if ortho.Len() > 1e-12 {
			u = ortho.Norm()
			break
		}
	}
	// Uniform over spherical cap: cosφ is linear in U
	// cosφ ∈ [cosθ, 1], with θ = l.Angle
	u1 := rng.Float64()
	cosPhi := 1 - u1*(l.oneMinusCos)
	sinPhi := math.Sqrt(1 - cosPhi*cosPhi)
	// Rotate in the axis/u 2D plane
	return axis.Mul(cosPhi).Add(u.Mul(sinPhi)).Norm()
}

// NewScene3D allocates a zero-initialized flat voxel grid for the given size and resolution.
func NewScene3D(center Point4, width, height, depth Real, nx, ny, nz, maxBounces int) *Scene3D {
	if nx <= 0 || ny <= 0 || nz <= 0 {
		panic("voxel resolution must be positive")
	}
	total := nx * ny * nz * 3
	halfX, halfY, halfZ := width*0.5, height*0.5, depth*0.5
	minX, maxX := center.X-halfX, center.X+halfX
	minY, maxY := center.Y-halfY, center.Y+halfY
	minZ, maxZ := center.Z-halfZ, center.Z+halfZ
	invSpanX := 1 / (maxX - minX)
	invSpanY := 1 / (maxY - minY)
	invSpanZ := 1 / (maxZ - minZ)

	strideY := nz * 3
	strideX := ny * strideY
	s := &Scene3D{
		Center: center,
		Width:  width, Height: height, Depth: depth,
		Nx: nx, Ny: ny, Nz: nz,
		Buf:        make([]Real, total),
		MaxBounces: maxBounces,
		MinX:       minX, MaxX: maxX, MinY: minY, MaxY: maxY, MinZ: minZ, MaxZ: maxZ,
		InvSpanX: invSpanX, InvSpanY: invSpanY, InvSpanZ: invSpanZ,
		StrideX: strideX, StrideY: strideY,
	}
	debugLog("Created scene center=%+v, size=(%.2f, %.2f, %.2f), resolution=(%d, %d, %d), maxBounces=%d", center, width, height, depth, nx, ny, nz, maxBounces)
	return s
}

// VoxelSize returns the physical size of each voxel along X,Y,Z.
func (s *Scene3D) VoxelSize() (dx, dy, dz Real) {
	dx, dy, dz = s.Width/Real(s.Nx), s.Height/Real(s.Ny), s.Depth/Real(s.Nz)
	debugLogOnce("Voxel size: (%.5f, %.5f, %.5f)", dx, dy, dz)
	return
}

// Bounds returns the min/max corners (in X,Y,Z; W is fixed) of the cube.
func (s *Scene3D) Bounds() (minX, maxX, minY, maxY, minZ, maxZ, W Real) {
	halfX := s.Width * 0.5
	halfY := s.Height * 0.5
	halfZ := s.Depth * 0.5
	minX, maxX, minY, maxY, minZ, maxZ, W = s.Center.X-halfX, s.Center.X+halfX,
		s.Center.Y-halfY, s.Center.Y+halfY,
		s.Center.Z-halfZ, s.Center.Z+halfZ,
		s.Center.W
	debugLogOnce("Scene bounds: X=(%.5f, %.5f), Y=(%.5f, %.5f), Z=(%.5f, %.5f), W=%.5f", minX, maxX, minY, maxY, minZ, maxZ, W)
	return
}

func (s *Scene3D) VoxelIndexOf(p Point4) (ok bool, i, j, k int, ux, uy, uz Real) {
	if p.X < s.MinX || p.X >= s.MaxX || p.Y < s.MinY || p.Y >= s.MaxY || p.Z < s.MinZ || p.Z >= s.MaxZ {
		return false, 0, 0, 0, 0, 0, 0
	}
	ux = (p.X - s.MinX) * s.InvSpanX
	uy = (p.Y - s.MinY) * s.InvSpanY
	uz = (p.Z - s.MinZ) * s.InvSpanZ
	i = int(ux * Real(s.Nx))
	j = int(uy * Real(s.Ny))
	k = int(uz * Real(s.Nz))
	if i == s.Nx {
		i = s.Nx - 1
	}
	if j == s.Ny {
		j = s.Ny - 1
	}
	if k == s.Nz {
		k = s.Nz - 1
	}
	return true, i, j, k, ux, uy, uz
}

func (s *Scene3D) idx(i, j, k, c int) int {
	return i*s.StrideX + j*s.StrideY + k*3 + c
}

func (s *Scene3D) AddHypercube(h *Hypercube4) {
	s.Hypercubes = append(s.Hypercubes, h)
}

// castSingleRay traces one photon from `light` through the scene.
//   - Picks ONE RGB channel (R/G/B) with probability ∝ light.Color.
//   - Marches the ray, always taking the nearest event: a hypercube hit or the W-plane hit.
//   - On hypercube, Russian-roulette per channel: absorb / reflect / refract.
//   - Energy is tinted by cube.Color[ch] on each surface interaction.
//   - Travel distance is accumulated across all segments; deposit uses 1/(dist^2 + eps).
//   - If deposit==true, it writes into scene.Buf with shard locks; otherwise it only tests
//     and returns whether the photon would have deposited.
//
// Returns true iff it deposits into any voxel on the plane.
func castSingleRay(light *ConeLight4, scene *Scene3D, rng *rand.Rand, locks *shardLocks, deposit bool) bool {
	// Choose the photon channel; if all zero, no photon.
	ch := pickChannelFast(light.colorSum, light.thrR, light.thrG, rng)
	// ch := pickChannelByColor(light.Color, rng)
	if ch < 0 {
		return false
	}

	// Initial throughput/energy: sum of channels so that E[deposit] = Color[ch]*w
	throughput := light.Color.R + light.Color.G + light.Color.B

	// Start at the light with a sampled unit direction.
	O := light.Origin
	D := light.SampleDir(rng) // unit

	const eps = 1e-6
	const bump = 1e-6

	totalDist := 0.0

	for bounce := 0; bounce < scene.MaxBounces; bounce++ {
		// Next plane hit (W == scene.Center.W)
		tPlane := planeHitT(scene, O, D)

		// Next hypercube hit
		hit, okCube := nearestCube(scene, O, D)

		// If no cube ahead or plane is closer → try to deposit on the plane.
		if !okCube || tPlane < hit.t {
			if !isFinite(tPlane) {
				return false // parallel to plane and no cube to stop us
			}

			P := O.Add(D.Mul(tPlane))
			if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
				totalDist += tPlane
				w := 1.0 / (totalDist*totalDist + eps)

				if deposit {
					base := scene.idx(i, j, k, ChR)
					if locks != nil {
						locks.lock(base)
					}
					scene.Buf[base+ch] += Real(throughput * w)
					if locks != nil {
						locks.unlock(base)
					}
				}
				return true
			}
			// Missed the volume bounds on the plane.
			return false
		}

		// Otherwise: hypercube is the first event.
		P := O.Add(D.Mul(hit.t))
		totalDist += hit.t

		hc := hit.hc
		// Per-channel probabilities.
		pRefl := rgbAt(hc.Reflect, ch)
		pRefr := rgbAt(hc.Refract, ch)
		pAbs := 1.0 - pRefl - pRefr
		if pAbs < 0 {
			pAbs = 0 // numeric clamp
		}

		// Russian roulette: absorb / reflect / refract.
		u := rng.Float64()
		if u < pAbs {
			return false // absorbed
		}

		// Tint throughput by the cube color (per channel).
		throughput *= rgbAt(hc.Color, ch)

		N := hit.Nw // outward normal (unit)

		if u < pAbs+pRefl {
			// Reflect
			D = reflect4(D, N)
			O = Point4{
				P.X + D.X*bump,
				P.Y + D.Y*bump,
				P.Z + D.Z*bump,
				P.W + D.W*bump,
			}
			continue
		}

		// Refract
		nOutside := 1.0
		nInside := rgbAt(hc.IOR, ch)
		var n1, n2 Real
		if hit.inv {
			// We were inside, exiting to outside.
			n1, n2 = nInside, nOutside
		} else {
			// Entering the cube.
			n1, n2 = nOutside, nInside
		}
		eta := n1 / n2

		if T, ok := refract4(D, N, eta); ok {
			D = T
			O = Point4{
				P.X + D.X*bump,
				P.Y + D.Y*bump,
				P.Z + D.Z*bump,
				P.W + D.W*bump,
			}
			continue
		}

		// Total internal reflection fallback.
		D = reflect4(D, N)
		O = Point4{
			P.X + D.X*bump,
			P.Y + D.Y*bump,
			P.Z + D.Z*bump,
			P.W + D.W*bump,
		}
	}

	// Safety bail-out (path didn't terminate sanely).
	return false
}

// raysPerLight must match len(lights). This gives you full control (e.g., after measuring p_hit per light).
func fireRaysFromLights(lights []*ConeLight4, scene *Scene3D, raysPerLight []int) {
	if len(lights) == 0 || len(raysPerLight) != len(lights) {
		return
	}

	// Total rays (for progress).
	totalRays := 0
	for _, n := range raysPerLight {
		if n > 0 {
			totalRays += n
		}
	}
	if totalRays == 0 {
		return
	}

	workers := runtime.NumCPU()
	if workers < 1 {
		workers = 1
	}

	// Distribute each light's rays across workers (evenly, with remainder spread).
	per := make([][]int, len(lights)) // [light][worker] -> count
	for li, n := range raysPerLight {
		per[li] = make([]int, workers)
		base, rem := n/workers, n%workers
		for w := 0; w < workers; w++ {
			per[li][w] = base
			if w < rem {
				per[li][w]++
			}
		}
	}

	var counter int64
	nextPrint := int64(1)
	if totalRays >= 100 {
		nextPrint = int64(totalRays / 100) // ~1%
	}

	locks := &shardLocks{}
	var wg sync.WaitGroup
	wg.Add(workers)

	for w := 0; w < workers; w++ {
		wid := w
		go func() {
			defer wg.Done()
			seed := time.Now().UnixNano() ^ int64(uint64(wid)*0x9e3779b97f4a7c15)
			rng := rand.New(rand.NewSource(seed))
			for li, L := range lights {
				n := per[li][wid]
				for s := 0; s < n; s++ {
					_ = castSingleRay(L, scene, rng, locks, true)
					fired := atomic.AddInt64(&counter, 1)
					if fired%nextPrint == 0 {
						fmt.Printf("[PROGRESS] %.2f%%\n", Real(fired)*100/Real(totalRays))
					}
				}
			}
		}()
	}

	wg.Wait()
}

func estimateHitProb(light *ConeLight4, scene *Scene3D, trials int) Real {
	if trials <= 0 {
		return 0
	}
	workers := runtime.NumCPU()
	if workers < 1 {
		workers = 1
	}
	if workers > trials {
		workers = trials
	}

	per, rem := trials/workers, trials%workers
	var wg sync.WaitGroup
	hitsCh := make(chan int, workers)

	for w := 0; w < workers; w++ {
		n := per
		if w < rem {
			n++
		}
		if n == 0 {
			continue
		}
		wg.Add(1)
		go func(wid, n int) {
			defer wg.Done()
			// independent RNG per worker
			seed := time.Now().UnixNano() ^ int64(uint64(wid)*0x9e3779b97f4a7c15)
			rng := rand.New(rand.NewSource(seed))

			localHits := 0
			for i := 0; i < n; i++ {
				if castSingleRay(light, scene, rng, nil, false) {
					localHits++
				}
			}
			hitsCh <- localHits
		}(w, n)
	}

	wg.Wait()
	close(hitsCh)

	totalHits := 0
	for h := range hitsCh {
		totalHits += h
	}
	return Real(totalHits) / Real(trials)
}

func intersectRayHypercube(O Point4, D Vector4, h *Hypercube4) (hit cubeHit, ok bool) {
	// world -> local
	Op := Vector4{O.X - h.Center.X, O.Y - h.Center.Y, O.Z - h.Center.Z, O.W - h.Center.W}
	Ol := h.RT.MulVec(Op)
	Dl := h.RT.MulVec(D)

	const inf = 1e300
	tmin, tmax := -inf, inf
	enterAxis, enterSign := -1, 0
	exitAxis, exitSign := -1, 0

	axis := func(o, d, half Real, ax int) (bool, Real, Real, int, int, int, int) {
		const eps = 1e-12
		if math.Abs(d) < eps {
			if math.Abs(o) > half {
				return false, 0, 0, 0, 0, 0, 0
			}
			return true, -inf, inf, 0, 0, 0, 0
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
		return cubeHit{}, false
	}
	okY, b1, b2, eAx2, eSg2, xAx2, xSg2 := axis(Ol.Y, Dl.Y, h.Half.Y, 1)
	if !okY {
		return cubeHit{}, false
	}
	okZ, c1, c2, eAx3, eSg3, xAx3, xSg3 := axis(Ol.Z, Dl.Z, h.Half.Z, 2)
	if !okZ {
		return cubeHit{}, false
	}
	okW, d1, d2, eAx4, eSg4, xAx4, xSg4 := axis(Ol.W, Dl.W, h.Half.W, 3)
	if !okW {
		return cubeHit{}, false
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
	for _, s := range sl {
		if s.t1 > tmin {
			tmin, enterAxis, enterSign = s.t1, s.eAx, s.eSg
		}
		if s.t2 < tmax {
			tmax, exitAxis, exitSign = s.t2, s.xAx, s.xSg
		}
	}
	if tmax < 0 || tmin > tmax {
		return cubeHit{}, false
	}

	inv := tmin < 0 && tmax > 0
	var t Real
	ax, sg := 0, 0
	if !inv {
		t, ax, sg = tmin, enterAxis, enterSign // entering
	} else {
		t, ax, sg = tmax, exitAxis, exitSign // exiting
	}

	// normal in local space
	var Nl Vector4
	switch ax {
	case 0:
		Nl = Vector4{Real(sg), 0, 0, 0}
	case 1:
		Nl = Vector4{0, Real(sg), 0, 0}
	case 2:
		Nl = Vector4{0, 0, Real(sg), 0}
	case 3:
		Nl = Vector4{0, 0, 0, Real(sg)}
	}
	Nw := h.R.MulVec(Nl).Norm()

	return cubeHit{t: t, Nw: Nw, hc: h, inv: inv}, true
}

func rayAABBParametric(O Point4, D Vector4, minP, maxP Point4) (bool, Real) {
	// 4D slab; we only need tNear for sorting
	tmin, tmax := -1e300, 1e300
	// inline per-axis to avoid loops and divisions by zero:
	upd := func(o, d, mn, mx Real) {
		const eps = 1e-12
		if math.Abs(d) < eps {
			if o < mn || o > mx {
				tmin = 1
				tmax = 0
				return
			} // force miss
			return
		}
		t1 := (mn - o) / d
		t2 := (mx - o) / d
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
	}
	upd(O.X, D.X, minP.X, maxP.X)
	upd(O.Y, D.Y, minP.Y, maxP.Y)
	upd(O.Z, D.Z, minP.Z, maxP.Z)
	upd(O.W, D.W, minP.W, maxP.W)
	if tmax < 0 || tmin > tmax {
		return false, 0
	}
	return true, tmin
}

func nearestCube(scene *Scene3D, O Point4, D Vector4) (cubeHit, bool) {
	best := cubeHit{}
	bestT := 1e300
	okAny := false
	for _, h := range scene.Hypercubes {
		if ok, tNear := rayAABBParametric(O, D, h.AABBMin, h.AABBMax); !ok || tNear > bestT {
			continue
		}
		hit, ok := intersectRayHypercube(O, D, h)
		if ok && hit.t > 1e-12 && hit.t < bestT {
			bestT, best, okAny = hit.t, hit, true
		}
	}
	return best, okAny
}

// reflection & refraction in 4D
func reflect4(I, N Vector4) Vector4 {
	return I.Sub(N.Mul(2 * I.Dot(N))).Norm()
}

func refract4(I, N Vector4, eta Real) (Vector4, bool) {
	cosi := -I.Dot(N)
	if cosi < -1 {
		cosi = -1
	}
	if cosi > 1 {
		cosi = 1
	}
	k := 1 - eta*eta*(1-cosi*cosi)
	if k < 0 {
		return Vector4{}, false // TIR
	}
	T := I.Mul(eta).Add(N.Mul(eta*cosi - math.Sqrt(k)))
	return T.Norm(), true
}

// plane distance to W = scene.Center.W (from O along D), +Inf if none
func planeHitT(scene *Scene3D, O Point4, D Vector4) Real {
	if math.Abs(D.W) < 1e-12 {
		return math.Inf(1)
	}
	t := (scene.Center.W - O.W) / D.W
	if t <= 1e-9 {
		return math.Inf(1)
	}
	return t
}

func isFinite(x Real) bool { return !math.IsInf(x, 0) && !math.IsNaN(x) }

func pickChannelFast(csum, thrR, thrG Real, rng *rand.Rand) int {
	if csum <= 0 {
		return -1
	}
	u := rng.Float64() * csum
	if u < thrR {
		return ChR
	}
	if u < thrG {
		return ChG
	}
	return ChB
}

func rgbAt(c RGB, ch int) Real {
	switch ch {
	case ChR:
		return c.R
	case ChG:
		return c.G
	default:
		return c.B
	}
}

// SaveAnimatedGIF writes a GIF with one frame per Z slice (k = 0..Nz-1).
// delay is in 100ths of a second (e.g., 5 => 20 fps).
// per-slice normalization + optional gamma (e.g., 0.7 brightens).
func SaveAnimatedGIF(scene *Scene3D, path string, delay int, gamma Real) error {
	Nx, Ny, Nz := scene.Nx, scene.Ny, scene.Nz

	out := &gif.GIF{
		Image:     make([]*image.Paletted, 0, Nz),
		Delay:     make([]int, 0, Nz),
		LoopCount: 0,
	}
	rgba := image.NewNRGBA(image.Rect(0, 0, Nx, Ny))

	// helper: scalar → 0..255 with gamma
	toByte := func(v, scale Real) uint8 {
		if v <= 0 {
			return 0
		}
		n := v * scale // 0..1 ideally
		if n > 1 {
			n = 1
		}
		if gamma != 1 {
			n = math.Pow(n, 1.0/gamma)
		}
		return uint8(math.Round(n * 255))
	}

	for k := 0; k < Nz; k++ {
		// Progress logging
		if k%max(1, Nz/100) == 0 { // ~1% steps
			percent := Real(k+1) * 100 / Real(Nz)
			fmt.Printf("[GIF] %.2f%%\n", percent)
		}
		// 1) find max over this slice
		sliceMax := 0.0
		for j := 0; j < Ny; j++ {
			for i := 0; i < Nx; i++ {
				idx := scene.idx(i, j, k, ChR)
				// use luminance-ish max across channels
				r := Real(scene.Buf[idx])
				g := Real(scene.Buf[idx+1])
				b := Real(scene.Buf[idx+2])
				// peak across channels
				if r > sliceMax {
					sliceMax = r
				}
				if g > sliceMax {
					sliceMax = g
				}
				if b > sliceMax {
					sliceMax = b
				}
			}
		}
		if sliceMax == 0 {
			sliceMax = 1 // avoid div-by-zero, will be black anyway
		}
		scale := 1.0 / sliceMax

		// 2) fill RGBA (flip Y so up is up)
		for j := 0; j < Ny; j++ {
			y := Ny - 1 - j
			for i := 0; i < Nx; i++ {
				idx := scene.idx(i, j, k, ChR)
				r := toByte(Real(scene.Buf[idx]), scale)
				g := toByte(Real(scene.Buf[idx+1]), scale)
				b := toByte(Real(scene.Buf[idx+2]), scale)
				p := rgba.PixOffset(i, y)
				rgba.Pix[p] = r
				rgba.Pix[p+1] = g
				rgba.Pix[p+2] = b
				rgba.Pix[p+3] = 255
			}
		}

		// 3) Quantize to paletted for GIF
		pimg := image.NewPaletted(rgba.Bounds(), palette.Plan9)
		draw.FloydSteinberg.Draw(pimg, pimg.Bounds(), rgba, image.Point{})

		out.Image = append(out.Image, pimg)
		out.Delay = append(out.Delay, delay)
	}

	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return gif.EncodeAll(f, out)
}

func loadConfig(path string) (*Config, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var cfg Config
	if err := json.Unmarshal(data, &cfg); err != nil {
		return nil, err
	}
	// Defaults / validation
	if cfg.SceneResX <= 0 {
		cfg.SceneResX = SceneResX
	}
	if cfg.SceneResY <= 0 {
		cfg.SceneResY = SceneResY
	}
	if cfg.SceneResZ <= 0 {
		cfg.SceneResZ = SceneResZ
	}
	if cfg.ProbeRays <= 0 {
		cfg.ProbeRays = ProbeRays
	}
	if cfg.Spp <= 0 {
		cfg.Spp = Spp
	}
	if cfg.GIFOut == "" {
		cfg.GIFOut = GIFOut
	}
	if cfg.GIFDelay <= 0 {
		cfg.GIFDelay = GIFDelay
	}
	if cfg.Gamma <= 0 {
		cfg.Gamma = Gamma
	}
	if len(cfg.Lights) == 0 {
		return nil, fmt.Errorf("config has no lights")
	}
	if cfg.Scene.MaxBounces <= 0 {
		cfg.Scene.MaxBounces = MaxBounces
	}
	return &cfg, nil
}

// ---------------------------------------------------------
// main
// ---------------------------------------------------------
func main() {
	rand.Seed(time.Now().UnixNano())
	var (
		cfg *Config
		err error
	)
	if len(os.Args) < 2 {
		cfg, err = loadConfig("config.json")
	} else {
		cfg, err = loadConfig(os.Args[1])
	}
	if err != nil {
		panic(err)
	}

	// Resolution
	Nx, Ny, Nz := cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ

	// Scene from JSON
	scene := NewScene3D(cfg.Scene.Center, cfg.Scene.Width, cfg.Scene.Height, cfg.Scene.Depth, Nx, Ny, Nz, cfg.Scene.MaxBounces)

	// Lights from JSON
	lights := make([]*ConeLight4, 0, len(cfg.Lights))
	for i, Lc := range cfg.Lights {
		angle := Lc.AngleDeg * math.Pi / 180.0
		L, err := NewConeLight4(Lc.Origin, Lc.Direction, Lc.Color, angle)
		if err != nil {
			panic(fmt.Errorf("light[%d] invalid: %w", i, err))
		}
		lights = append(lights, L)
	}

	// Build hypercubes from JSON (if any)
	for i, hc := range cfg.Hypercubes {
		h, err := hc.Build()
		if err != nil {
			fmt.Printf("cannot add hypercube[%d]: %+v, skipping", i, err)
			continue
		}
		scene.AddHypercube(h)
	}

	// Rays per light to hit target SPP
	Nvox := Nx * Ny * Nz
	needRays := make([]int, len(lights))
	totalRays := 0
	for i, L := range lights {
		p := estimateHitProb(L, scene, cfg.ProbeRays)
		if p < 1e-7 {
			p = 1e-7 // avoid div-by-zero; treat as extremely low hit prob
		}
		need := int(3 * Real(cfg.Spp) * Real(Nvox) / p)
		if need < cfg.ProbeRays {
			need = cfg.ProbeRays
		} // at least probe rays
		debugLog("light[%d] p_hit≈%.3f → rays=%d for %d spp @ %dx%dx%d",
			i, p, need, cfg.Spp, Nx, Ny, Nz)
		needRays[i] = need
		totalRays += need
	}

	// Fire rays (multi-light, sharded writes)
	start := time.Now()
	fireRaysFromLights(lights, scene, needRays)
	runtime.GC()
	elapsed := time.Since(start)

	// Diagnostics
	ci, cj, ck := Nx/2, Ny/2, Nz/2
	debugLog("Rays: %d, time: %s", totalRays, elapsed)
	debugLog("Center voxel RGB: (%.6g, %.6g, %.6g)",
		Real(scene.Buf[scene.idx(ci, cj, ck, ChR)]),
		Real(scene.Buf[scene.idx(ci, cj, ck, ChG)]),
		Real(scene.Buf[scene.idx(ci, cj, ck, ChB)]),
	)

	// Save GIF
	if err := SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma); err != nil {
		panic(err)
	}
	fmt.Println("Saved animated GIF:", cfg.GIFOut)
}
