package photons4d

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sync"
	"sync/atomic"
	"time"
)

func castSingleRay(light *Light, scene *Scene, rng *rand.Rand, locks *shardLocks, deposit bool) bool {
	// Choose the photon channel.
	ch := pickChannel(light.colorSum, light.thrR, light.thrG, rng)

	// Initial throughput/energy: sum of channels so that E[deposit] = Color[ch]*w
	throughput := light.Color.R + light.Color.G + light.Color.B

	// Start at the light with a sampled unit direction.
	O := light.Origin
	D := light.SampleDir(rng) // approx unit
	// normalize for safety: starts
	l2 := D.Dot(D)
	if l2 != 0 {
		D = D.Mul(1 / math.Sqrt(l2))
	}
	// normalize for safety: ends

	totalDist, w := 0.0, 0.0

	for bounce := 0; bounce < scene.MaxBounces; bounce++ {
		// optional periodic re-normalization for numerical stability
		if (bounce & 3) == 3 {
			l2 := D.Dot(D)
			if l2 > 0 {
				D = D.Mul(1 / math.Sqrt(l2))
			}
		}

		// Next plane hit (W == scene.Center.W)
		tPlane := planeHit(scene, O, D)

		// Next object hit (cube/ellipsoid) with AABB pre-cull
		hit, okObj := nearestHit(scene, O, D, tPlane)

		// If no cube ahead or plane is closer → try to deposit on the plane.
		if !okObj || tPlane < hit.t {
			if !isFinite(tPlane) {
				if Debug {
					logRay("parallel_to_scene", Miss, O, D, Point4{}, bounce, totalDist)
				}
				// parallel to plane and no cube to stop us
				return false
			}

			P := O.Add(D.Mul(tPlane))
			if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
				totalDist += tPlane
				if AttenuateD2 {
					w = 1.0 / (totalDist*totalDist + epsDist)
				} else {
					w = 1.0
				}

				if deposit && UseLocks {
					base := scene.idx(i, j, k, ChR)
					if locks != nil {
						locks.lock(base)
					}
					scene.Buf[base+ch] += Real(throughput * w * light.voidCoeff)
					if locks != nil {
						locks.unlock(base)
					}
				} else if deposit {
					base := scene.idx(i, j, k, ChR)
					scene.Buf[base+ch] += Real(throughput * w * light.voidCoeff)
				}
				if Debug {
					logRay("hit_scene", Hit, O, D, P, bounce, totalDist)
				}
				return true
			}
			if Debug {
				logRay("miss_scene", Miss, O, D, P, bounce, totalDist)
			}
			return false
		}

		// Otherwise: hypercube is the first event.
		P := O.Add(D.Mul(hit.t))
		totalDist += hit.t

		// Outward normal (unit)
		N := hit.Nw

		// --- Absorption & Fresnel split (per channel) ---
		pAbs := hit.pAbsCh(ch) // object's absorption knob
		avail := 1 - pAbs      // budget to split between reflection/refraction
		if avail <= 0 {
			if Debug {
				logRay("absorbed", Absorb, O, D, P, bounce, totalDist)
			}
			return false
		}

		// cosθ for Schlick (ensure in [0,1])
		var cosTheta Real
		if hit.inv {
			cosTheta = D.Dot(N) // exiting: incident vs outward normal
		} else {
			cosTheta = -D.Dot(N) // entering: flip sign
		}
		if cosTheta < 0 {
			cosTheta = 0
		} else if cosTheta > 1 {
			cosTheta = 1
		}

		// Schlick Fresnel using cached F0 (air ↔ material)
		F0 := hit.f0Ch(ch)
		x := 1 - cosTheta
		x2 := x * x
		x5 := x2 * x2 * x
		F := F0 + (1-F0)*x5

		// Bias Fresnel by your reflect/refract knobs to control the split.
		rW := hit.reflCh(ch) * F
		tW := hit.refrCh(ch) * (1 - F)
		f := F
		if sum := rW + tW; sum > 0 {
			f = rW / sum
		}

		pReflDyn := avail * f // final reflect probability this hit

		// --- Russian roulette ---
		u := rng.Float64()
		if u < pAbs {
			if Debug {
				logRay("absorbed", Absorb, O, D, P, bounce, totalDist)
			}
			return false
		}

		// Survived → tint throughput by cube color.
		throughput *= hit.colorCh(ch)

		if u < pAbs+pReflDyn {
			// Reflect
			D = reflect4(D, N)
			// keep direction unit-length
			if l2 := D.Dot(D); l2 > 0 {
				D = D.Mul(1 / math.Sqrt(l2))
			}
			O = Point4{
				P.X + D.X*bumpShift,
				P.Y + D.Y*bumpShift,
				P.Z + D.Z*bumpShift,
				P.W + D.W*bumpShift,
			}
			if Debug {
				logRay("reflected", Reflect, O, D, P, bounce, totalDist)
			}
			continue
		}

		// Refract
		var eta Real
		if hit.inv {
			// exiting: inside -> outside
			eta = hit.iorCh(ch)
		} else {
			// entering: outside -> inside
			eta = hit.iorInvCh(ch)
		}
		if T, ok := refract4(D, N, eta); ok {
			D = T
			// keep direction unit-length
			if l2 := D.Dot(D); l2 > 0 {
				D = D.Mul(1 / math.Sqrt(l2))
			}
			O = Point4{
				P.X + D.X*bumpShift,
				P.Y + D.Y*bumpShift,
				P.Z + D.Z*bumpShift,
				P.W + D.W*bumpShift,
			}
			if Debug {
				logRay("refracted", Refract, O, D, P, bounce, totalDist)
			}
			continue
		}

		// Total internal reflection fallback.
		D = reflect4(D, N)
		if l2 := D.Dot(D); l2 > 0 {
			D = D.Mul(1 / math.Sqrt(l2))
		}
		O = Point4{
			P.X + D.X*bumpShift,
			P.Y + D.Y*bumpShift,
			P.Z + D.Z*bumpShift,
			P.W + D.W*bumpShift,
		}
		if Debug {
			logRay("total_internal_refraction", TIR, O, D, P, bounce, totalDist)
		}
	}

	if Debug {
		logRay("too_complex", RecurenceLimit, O, D, Point4{}, scene.MaxBounces, totalDist)
	}
	return false
}

// raysPerLight must match len(lights). This gives you full control (e.g., after measuring p_hit per light).
func castRays(lights []*Light, scene *Scene, raysPerLight []int) {
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
