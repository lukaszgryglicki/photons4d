package main

import (
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
	ChR      = 0
	ChG      = 1
	ChB      = 2
	SceneRes = 400
	Spp      = 64
)

const numShards = 4096 // power of two is nice
type shardLocks struct{ mu [numShards]sync.Mutex }

func (sl *shardLocks) lock(idx int)   { sl.mu[idx&(numShards-1)].Lock() }
func (sl *shardLocks) unlock(idx int) { sl.mu[idx&(numShards-1)].Unlock() }

// ---------------------------------------------------------
// Geometry & color
// ---------------------------------------------------------

// Point4 represents a point in 4-dimensional space.
type Point4 struct {
	X, Y, Z, W float64
}

// Add lets you translate a Point4 by a Vector4.
func (p Point4) Add(v Vector4) Point4 {
	point := Point4{p.X + v.X, p.Y + v.Y, p.Z + v.Z, p.W + v.W}
	// debugLog("Adding Vector4 to Point4: (%f, %f, %f, %f) + (%f, %f, %f, %f) = (%f, %f, %f, %f)", p.X, p.Y, p.Z, p.W, v.X, v.Y, v.Z, v.W, point.X, point.Y, point.Z, point.W)
	return point
}

// Vector4 represents a direction (not a position) in 4D space.
type Vector4 struct {
	X, Y, Z, W float64
}

// Vector functions
func (a Vector4) Add(b Vector4) Vector4 { return Vector4{a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W} }
func (a Vector4) Sub(b Vector4) Vector4 { return Vector4{a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W} }
func (v Vector4) Mul(s float64) Vector4 { return Vector4{v.X * s, v.Y * s, v.Z * s, v.W * s} }

// Dot returns the dot product between two 4D vectors.
func (a Vector4) Dot(b Vector4) float64 {
	dot := a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W
	// debugLog("Dot product: (%f, %f, %f, %f) . (%f, %f, %f, %f) = %f", a.X, a.Y, a.Z, a.W, b.X, b.Y, b.Z, b.W, dot)
	return dot
}

// Len returns the Euclidean length of the vector.
func (v Vector4) Len() float64 {
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

// RGB stores color components; each should be in [0,1].
type RGB struct {
	R, G, B float64
}

// clamp01 clamps each channel to [0,1] (useful for validation).
func (c RGB) clamp01() RGB {
	cl := func(x float64) float64 {
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

// ---------------------------------------------------------
// Light
// ---------------------------------------------------------

// ConeLight4 is a 4D "spot" light: rays originate at Origin and
// are restricted to a cone around Direction with half-angle Angle (radians).
type ConeLight4 struct {
	Origin    Point4
	Direction Vector4 // should be unit-length; constructor normalizes
	Color     RGB     // each in [0,1]
	Angle     float64 // half-angle in radians, (0, π]

	// cached
	cosAngle float64
}

// NewConeLight4 constructs a cone light, normalizing direction and clamping color.
// Angle must be in (0, π]. Returns an error if the direction is zero or angle is invalid.
func NewConeLight4(origin Point4, dir Vector4, color RGB, angle float64) (*ConeLight4, error) {
	if angle <= 0 || angle > math.Pi {
		return nil, errors.New("angle must be in (0, π]")
	}
	n := dir.Norm()
	if n.Len() == 0 {
		return nil, errors.New("direction must be non-zero")
	}
	light := &ConeLight4{
		Origin:    origin,
		Direction: n,
		Color:     color.clamp01(),
		Angle:     angle,
		cosAngle:  math.Cos(angle),
	}
	debugLog("Created light %+v", light)
	return light, nil
}

// ContainsDir reports whether a given direction lies inside the cone.
// The input does not have to be unit-length.
func (l *ConeLight4) ContainsDir(v Vector4) bool {
	u := v.Norm()
	// angle between unit vectors via dot product threshold
	contains := l.Direction.Dot(u) >= l.cosAngle
	debugLog("ContainsDir: light direction (%f, %f, %f, %f) vs input (%f, %f, %f, %f) => %v", l.Direction.X, l.Direction.Y, l.Direction.Z, l.Direction.W, u.X, u.Y, u.Z, u.W, contains)
	return contains
}

func (l *ConeLight4) ContainsPoint(p Point4) bool {
	v := Vector4{
		p.X - l.Origin.X,
		p.Y - l.Origin.Y,
		p.Z - l.Origin.Z,
		p.W - l.Origin.W,
	}
	dist := v.Len()
	if dist == 0 {
		debugLog("ContainsPoint: point is at the origin of the cone: %+v, returning true", l.Origin)
		return true // point at origin of cone
	}
	u := v.Mul(1.0 / dist) // normalize
	insideAngle := l.Direction.Dot(u) >= l.cosAngle
	debugLog("ContainsPoint: light origin (%f, %f, %f, %f) vs point (%f, %f, %f, %f) => dist: %f, u: %+v, inside angle: %v", l.Origin.X, l.Origin.Y, l.Origin.Z, l.Origin.W, p.X, p.Y, p.Z, p.W, dist, u, insideAngle)
	return insideAngle // && dist <= maxLength (if you want range limit)
}

func randNormalGlobal() float64 {
	// Box–Muller
	u1 := rand.Float64()
	u2 := rand.Float64()
	r := math.Sqrt(-2 * math.Log(math.Max(u1, 1e-12)))
	r2 := r * math.Cos(2*math.Pi*u2)
	// LG:skip
	// debugLog("Generated random normal: u1=%.5f, u2=%.5f, r=%.5f, r2=%.5f", u1, u2, r, r2)
	return r2
}

func randNormal(rng *rand.Rand) float64 {
	// Box–Muller
	u1 := rng.Float64()
	u2 := rng.Float64()
	r := math.Sqrt(-2 * math.Log(math.Max(u1, 1e-12)))
	return r * math.Cos(2*math.Pi*u2)
}

// SampleDir returns a unit direction inside the cone around l.Direction with half-angle l.Angle.
// Notes:
//   - This picks phi uniformly in [0, Angle]. That’s simple, but not uniform in solid angle.
//     If you want solid-angle-uniform sampling on S^3’s spherical cap later, we can tweak the phi distribution.
func (l *ConeLight4) SampleDirGlobal() Vector4 {
	axis := l.Direction.Norm() // should already be unit, but be safe

	// 1) pick an angle within the cone (simple uniform in [0, Angle])
	phi := rand.Float64() * l.Angle
	c, s := math.Cos(phi), math.Sin(phi)
	// LG:skip
	// debugLog("SampleDir: angle phi=%.5f (cos=%.5f, sin=%.5f)", phi, c, s)

	// 2) sample a random unit vector orthogonal to axis (lives in the 3D subspace ⟂ to axis)
	//    Do: draw a random 4D normal vector r, remove its projection onto axis, normalize.
	for {
		r := Vector4{randNormalGlobal(), randNormalGlobal(), randNormalGlobal(), randNormalGlobal()}
		ortho := r.Sub(axis.Mul(r.Dot(axis)))
		if ortho.Len() > 1e-12 {
			u := ortho.Norm()
			// 3) rotate off the axis by phi inside the cone
			rd := axis.Mul(c).Add(u.Mul(s)).Norm()
			// LG:skip
			// debugLog("SampleDir: sampled direction (%.5f, %.5f, %.5f, %.5f)", rd.X, rd.Y, rd.Z, rd.W)
			return rd
		}
		// if degenerate, resample
	}
}

func (l *ConeLight4) SampleDir(rng *rand.Rand) Vector4 {
	axis := l.Direction.Norm()
	phi := rng.Float64() * l.Angle
	c, s := math.Cos(phi), math.Sin(phi)

	for {
		r := Vector4{
			randNormal(rng),
			randNormal(rng),
			randNormal(rng),
			randNormal(rng),
		}
		ortho := r.Sub(axis.Mul(r.Dot(axis)))
		if ortho.Len() > 1e-12 {
			u := ortho.Norm()
			return axis.Mul(c).Add(u.Mul(s)).Norm()
		}
	}
}

// ---------------------------------------------------------
// Scene: 3D cube at fixed W with flat voxel buffer
// ---------------------------------------------------------

// Scene3D stores a 3D volume (axis-aligned in X,Y,Z) embedded at W=Center.W.
// Voxels are in a flat buffer: len = Nx*Ny*Nz*3 (RGB).
type Scene3D struct {
	Center               Point4
	Width, Height, Depth float64
	Nx, Ny, Nz           int
	Buf                  []Real // flat: (((i*Ny)+j)*Nz + k)*3 + c
}

// NewScene3D allocates a zero-initialized flat voxel grid for the given size and resolution.
func NewScene3D(center Point4, width, height, depth float64, nx, ny, nz int) *Scene3D {
	if nx <= 0 || ny <= 0 || nz <= 0 {
		panic("voxel resolution must be positive")
	}
	total := nx * ny * nz * 3
	s := &Scene3D{
		Center: center,
		Width:  width, Height: height, Depth: depth,
		Nx: nx, Ny: ny, Nz: nz,
		Buf: make([]Real, total),
	}
	debugLog("Created scene center=%+v, size=(%.2f, %.2f, %.2f), resolution=(%d, %d, %d)", center, width, height, depth, nx, ny, nz)
	return s
}

// VoxelSize returns the physical size of each voxel along X,Y,Z.
func (s *Scene3D) VoxelSize() (dx, dy, dz float64) {
	dx, dy, dz = s.Width/float64(s.Nx), s.Height/float64(s.Ny), s.Depth/float64(s.Nz)
	debugLogOnce("Voxel size: (%.5f, %.5f, %.5f)", dx, dy, dz)
	return
}

// Bounds returns the min/max corners (in X,Y,Z; W is fixed) of the cube.
func (s *Scene3D) Bounds() (minX, maxX, minY, maxY, minZ, maxZ, W float64) {
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

// VoxelIndexOf maps a 4D point to voxel indices and also returns normalized coords (u,v,w) in [0,1].
func (s *Scene3D) VoxelIndexOf(p Point4) (ok bool, i, j, k int, ux, uy, uz float64) {
	minX, maxX, minY, maxY, minZ, maxZ, _ := s.Bounds()
	if p.X < minX || p.X >= maxX || p.Y < minY || p.Y >= maxY || p.Z < minZ || p.Z >= maxZ {
		// LG:skip
		// debugLog("VoxelIndexOf: point (%f, %f, %f, %f) is outside bounds", p.X, p.Y, p.Z, p.W)
		return false, 0, 0, 0, 0, 0, 0
	}
	ux = (p.X - minX) / (maxX - minX)
	uy = (p.Y - minY) / (maxY - minY)
	uz = (p.Z - minZ) / (maxZ - minZ)
	i = int(ux * float64(s.Nx))
	j = int(uy * float64(s.Ny))
	k = int(uz * float64(s.Nz))
	if i == s.Nx {
		i = s.Nx - 1
	}
	if j == s.Ny {
		j = s.Ny - 1
	}
	if k == s.Nz {
		k = s.Nz - 1
	}
	// LG:skip
	// debugLog("VoxelIndexOf: point (%f, %f, %f, %f) maps to voxel i=%d j=%d k=%d (u=%.3f v=%.3f w=%.3f)", p.X, p.Y, p.Z, p.W, i, j, k, ux, uy, uz)
	return true, i, j, k, ux, uy, uz
}

// Flat buffer index helper (c ∈ {ChR,ChG,ChB}).
func (s *Scene3D) idx(i, j, k, c int) int {
	return (((i*s.Ny)+j)*s.Nz+k)*3 + c
}

// ---------------------------------------------------------
// Ray cast: one ray → deposit
// ---------------------------------------------------------

// castOneRay fires a single ray. If it hits the scene (W-plane + inside XYZ bounds),
// it accumulates inverse-square intensity into the hit voxel and prints debug info.
func castOneRay(light *ConeLight4, scene *Scene3D) {
	// sample direction
	D := light.SampleDirGlobal() // unit
	// LG:skip
	// debugLog("Ray dir: (%.5f, %.5f, %.5f, %.5f)", D.X, D.Y, D.Z, D.W)

	// intersect with hyperplane W = scene.Center.W
	den := D.W
	if math.Abs(den) < 1e-12 {
		// LG:skip
		// debugLog("MISS: ray parallel to scene hyperplane (W constant).")
		return
	}
	t := (scene.Center.W - light.Origin.W) / den
	if t <= 0 {
		// LG:skip
		// debugLog("MISS: intersection behind origin (t=%.6f)", t)
		return
	}
	P := light.Origin.Add(D.Mul(t))
	// LG:skip
	// debugLog("Intersection t=%.6f at P=(%.5f, %.5f, %.5f, %.5f)", t, P.X, P.Y, P.Z, P.W)

	// check inside bounds and map to voxel
	// _, _, _ => ux, uy, uz
	if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
		// inverse-square falloff by traveled distance (since D is unit, distance = t)
		const eps = 1e-6
		w := 1.0 / (t*t + eps)

		// deposit (per-channel scale by light color)
		addR := light.Color.R * w
		addG := light.Color.G * w
		addB := light.Color.B * w

		scene.Buf[scene.idx(i, j, k, ChR)] += Real(addR)
		scene.Buf[scene.idx(i, j, k, ChG)] += Real(addG)
		scene.Buf[scene.idx(i, j, k, ChB)] += Real(addB)

		// LG:skip
		// debugLog("HIT: voxel i=%d j=%d k=%d  (u=%.3f v=%.3f w=%.3f)", i, j, k, ux, uy, uz)
		// debugLog("     added RGB = (%.6g, %.6g, %.6g)", addR, addG, addB)
	} else {
		// LG:skip
		// debugLog("MISS: intersection outside scene XYZ bounds.")
	}
}

// Same as castOneRay, but deposits into 'buf' instead of scene.Buf.
// Safe to call from multiple goroutines with different 'buf's.
func castOneRayInto(light *ConeLight4, scene *Scene3D, buf []Real, rng *rand.Rand) {
	D := light.SampleDir(rng) // unit

	den := D.W
	if math.Abs(den) < 1e-12 {
		return
	}
	t := (scene.Center.W - light.Origin.W) / den
	if t <= 0 {
		return
	}
	P := light.Origin.Add(D.Mul(t))

	if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
		const eps = 1e-6
		w := 1.0 / (t*t + eps)
		addR := light.Color.R * w
		addG := light.Color.G * w
		addB := light.Color.B * w

		idxR := scene.idx(i, j, k, ChR)
		buf[idxR+0] += Real(addR)
		buf[idxR+1] += Real(addG)
		buf[idxR+2] += Real(addB)
	}
}

// Write directly into scene.Buf
func castOneRayShard(light *ConeLight4, scene *Scene3D, rng *rand.Rand, locks *shardLocks) {
	D := light.SampleDir(rng)
	den := D.W
	if math.Abs(den) < 1e-12 {
		return
	}
	t := (scene.Center.W - light.Origin.W) / den
	if t <= 0 {
		return
	}
	P := light.Origin.Add(D.Mul(t))
	if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
		const eps = 1e-6
		w := 1.0 / (t*t + eps)
		base := scene.idx(i, j, k, ChR)
		locks.lock(base)
		scene.Buf[base+0] += Real(light.Color.R * w)
		scene.Buf[base+1] += Real(light.Color.G * w)
		scene.Buf[base+2] += Real(light.Color.B * w)
		locks.unlock(base)
	}
}

// shoot many rays
func fireRays(light *ConeLight4, scene *Scene3D, n int) {
	for s := 0; s < n; s++ {
		castOneRay(light, scene)
	}
}

// Uses all CPU cores. Spawns NumCPU workers; each gets its own local buffer.
// After all workers finish, reduces local buffers into scene.Buf.
func fireRaysParallel(light *ConeLight4, scene *Scene3D, totalRays int) {
	workers := runtime.NumCPU()
	if workers < 1 {
		workers = 1
	}
	raysPer := totalRays / workers
	rem := totalRays % workers

	debugLogOnce("Launching %d workers (rays: %d each, +1 for first %d workers)", workers, raysPer, rem)

	// local buffers per worker
	locals := make([][]Real, workers)
	for w := 0; w < workers; w++ {
		locals[w] = make([]Real, len(scene.Buf))
	}

	var counter int64
	nextPrint := int64(totalRays / 100) // every ~1%
	if nextPrint < 1 {
		nextPrint = 1 // at least once
	}
	var wg sync.WaitGroup
	wg.Add(workers)
	for w := 0; w < workers; w++ {
		// shards: first 'rem' workers do one extra ray
		count := raysPer
		if w < rem {
			count++
		}

		// Each worker runs independently and writes only to its local buffer
		go func(wid, n int) {
			defer wg.Done()
			// (Optional) Per-worker RNG for speed:
			// r := rand.New(rand.NewSource(time.Now().UnixNano() + int64(wid)))
			// If you switch SampleDir/randNormal to use 'r', contention drops further.
			local := locals[wid]
			rng := rand.New(rand.NewSource(time.Now().UnixNano() + int64(wid)))
			for s := 0; s < n; s++ {
				castOneRayInto(light, scene, local, rng)
				// update counter atomically
				fired := atomic.AddInt64(&counter, 1)
				if fired%nextPrint == 0 {
					percent := float64(fired) / float64(totalRays) * 100
					fmt.Printf("[PROGRESS] %.2f%%\n", percent)
				}
			}
		}(w, count)
	}
	wg.Wait()

	// Reduce locals → global
	// Single-threaded reduction is often fine; if len is huge, you can parallelize this too.
	buf := scene.Buf
	for w := 0; w < workers; w++ {
		src := locals[w]
		for i := 0; i < len(buf); i++ {
			buf[i] += src[i]
		}
	}
}

func fireRaysParallelShard(light *ConeLight4, scene *Scene3D, totalRays int) {
	workers := runtime.NumCPU()
	if workers < 1 {
		workers = 1
	}
	raysPer, rem := totalRays/workers, totalRays%workers

	var counter int64
	nextPrint := int64(1)
	if totalRays >= 100 {
		nextPrint = int64(totalRays / 100)
	}

	locks := &shardLocks{}
	var wg sync.WaitGroup
	wg.Add(workers)
	for w := 0; w < workers; w++ {
		count := raysPer
		if w < rem {
			count++
		}
		go func(wid, n int) {
			defer wg.Done()
			rng := rand.New(rand.NewSource(time.Now().UnixNano() ^ int64(wid)))
			for s := 0; s < n; s++ {
				castOneRayShard(light, scene, rng, locks)
				fired := atomic.AddInt64(&counter, 1)
				if fired%nextPrint == 0 {
					fmt.Printf("[PROGRESS] %.2f%%\n", float64(fired)*100/float64(totalRays))
				}
			}
		}(w, count)
	}
	wg.Wait()
}

func estimateHitProb(light *ConeLight4, scene *Scene3D, trials int) float64 {
	hits := 0
	for i := 0; i < trials; i++ {
		D := light.SampleDirGlobal()
		den := D.W
		if math.Abs(den) < 1e-12 {
			continue
		}
		t := (scene.Center.W - light.Origin.W) / den
		if t <= 0 {
			continue
		}
		P := light.Origin.Add(D.Mul(t))
		if ok, _, _, _, _, _, _ := scene.VoxelIndexOf(P); ok {
			hits++
		}
	}
	return float64(hits) / float64(trials)
}

// find the global max channel value across the whole buffer (for consistent brightness)
func (s *Scene3D) maxChannel() float64 {
	maxv := 0.0
	for i := 0; i < len(s.Buf); i++ {
		v := float64(s.Buf[i])
		if v > maxv {
			maxv = v
		}
	}
	if maxv == 0 {
		maxv = 1 // avoid div-by-zero
	}
	return maxv
}

// SaveAnimatedGIF writes a GIF with one frame per Z slice (k = 0..Nz-1).
// delay is in 100ths of a second (e.g., 5 => 20 fps).
// per-slice normalization + optional gamma (e.g., 0.7 brightens).
func SaveAnimatedGIF(scene *Scene3D, path string, delay int, gamma float64) error {
	Nx, Ny, Nz := scene.Nx, scene.Ny, scene.Nz

	out := &gif.GIF{
		Image:     make([]*image.Paletted, 0, Nz),
		Delay:     make([]int, 0, Nz),
		LoopCount: 0,
	}
	rgba := image.NewNRGBA(image.Rect(0, 0, Nx, Ny))

	// helper: scalar → 0..255 with gamma
	toByte := func(v, scale float64) uint8 {
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
			percent := float64(k) * 100 / float64(Nz)
			fmt.Printf("[GIF] %.2f%%\n", percent)
		}
		// 1) find max over this slice
		sliceMax := 0.0
		for j := 0; j < Ny; j++ {
			for i := 0; i < Nx; i++ {
				idx := scene.idx(i, j, k, ChR)
				// use luminance-ish max across channels
				r := float64(scene.Buf[idx+0])
				g := float64(scene.Buf[idx+1])
				b := float64(scene.Buf[idx+2])
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
				r := toByte(float64(scene.Buf[idx+0]), scale)
				g := toByte(float64(scene.Buf[idx+1]), scale)
				b := toByte(float64(scene.Buf[idx+2]), scale)
				p := rgba.PixOffset(i, y)
				rgba.Pix[p+0] = r
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

// ---------------------------------------------------------
// main
// ---------------------------------------------------------

func main() {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("4d")

	// Choose resolution. Start with 256^3 to check speed/memory first.
	// 512^3 with float32 ≈ 1.6 GB; ensure you have enough RAM.
	Nx, Ny, Nz := SceneRes, SceneRes, SceneRes
	// Nx, Ny, Nz := 512, 512, 512 // <-- enable if you have >2 GB free

	// Example scene: 2×2×2 cube at W=0
	scene := NewScene3D(Point4{0, 0, 0, 0}, 2.0, 2.0, 2.0, Nx, Ny, Nz)

	// Example light: at W=+1, pointing toward -W with a warm color, ~25° cone
	light, err := NewConeLight4(
		Point4{0, 0, 0, 1.0},
		Vector4{0, 0, 0, -1}, // aim toward decreasing W
		RGB{1.0, 0.8, 0.2},   // color in 0..1
		25*math.Pi/180.0,     // half-angle
	)
	if err != nil {
		panic(err)
	}

	p := estimateHitProb(light, scene, 100000) // 1e5 is plenty quick
	need := int(float64(Spp) * float64(SceneRes*SceneRes*SceneRes) / p)
	debugLog("p_hit≈%.3f, rays needed for %d spp @ N=%d: %d\n", p, Spp, SceneRes, need)

	// Fire a LOT of rays to saturate the scene
	totalRays := int(need)
	start := time.Now()
	// fireRays(light, scene, totalRays)
	// fireRaysParallel(light, scene, totalRays)
	fireRaysParallelShard(light, scene, totalRays)
	runtime.GC()
	elapsed := time.Since(start)

	// Diagnostics
	ci, cj, ck := Nx/2, Ny/2, Nz/2
	debugLog("Rays: %d, time: %s", totalRays, elapsed)
	debugLog("Center voxel RGB: (%.6g, %.6g, %.6g)",
		float64(scene.Buf[scene.idx(ci, cj, ck, ChR)]),
		float64(scene.Buf[scene.idx(ci, cj, ck, ChG)]),
		float64(scene.Buf[scene.idx(ci, cj, ck, ChB)]),
	)

	// Sum total energy (slow for huge grids; ok as a quick check)
	var sumR, sumG, sumB float64
	voxels := len(scene.Buf)
	for i := 0; i < voxels; i += 3 {
		sumR += float64(scene.Buf[i+ChR])
		sumG += float64(scene.Buf[i+ChG])
		sumB += float64(scene.Buf[i+ChB])
	}
	voxels /= 3
	sumR /= float64(voxels)
	sumG /= float64(voxels)
	sumB /= float64(voxels)
	debugLog("Average energy  R: %.6g  G: %.6g  B: %.6g", sumR, sumG, sumB)

	// Save animation: one frame per Z slice
	outPath := "volume.gif"
	delay := 5    // 5 × 1/100s = 50 ms per frame ≈ 20 FPS
	gamma := 0.75 // tweak 0.6..1.0 (lower = brighter)
	if err := SaveAnimatedGIF(scene, outPath, delay, gamma); err != nil {
		panic(err)
	}
	fmt.Println("Saved animated GIF:", outPath)
}
