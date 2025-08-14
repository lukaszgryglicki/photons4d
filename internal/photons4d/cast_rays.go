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
	throughput := light.Intensity * (light.Color.R + light.Color.G + light.Color.B)

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

		// Next object hit (cell/ellipsoid) with AABB pre-cull
		hit, okObj := NearestHitFunc(scene, O, D, tPlane)

		// If no cell/ellipsoid ahead or plane is closer → try to deposit on the plane.
		if !okObj || tPlane < hit.t {
			if !isFinite(tPlane) {
				if Debug && deposit {
					logRay("parallel_to_scene", Miss, O, D, Point4{}, bounce, totalDist)
				}
				// parallel to plane and no object to stop us
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
					scene.Buf[base+ch] += Real(throughput * w)
					if locks != nil {
						locks.unlock(base)
					}
				} else if deposit {
					base := scene.idx(i, j, k, ChR)
					scene.Buf[base+ch] += Real(throughput * w)
				}
				if Debug && deposit {
					logRay("hit_scene", Hit, O, D, P, bounce, totalDist)
				}
				return true
			}
			if Debug && deposit {
				logRay("miss_scene", Miss, O, D, P, bounce, totalDist)
			}
			return false
		}

		// Otherwise: cell8 is the first event.
		P := O.Add(D.Mul(hit.t))
		totalDist += hit.t

		// Outward normal (unit)
		N := hit.Nw

		// --- Absorption & Fresnel split (per channel) ---
		pAbs := hit.pAbsCh(ch) // object's absorption knob
		avail := 1 - pAbs      // budget to split among reflection/refraction/diffuse
		if avail <= 0 {
			if Debug && deposit {
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

		// Bias Fresnel by reflect/refract knobs; diffuse is angle-independent.
		pReflDyn, pRefrDyn, pDiffDyn := splitEventWeights(hit.reflCh(ch), hit.refrCh(ch), hit.diffCh(ch), F, pAbs)

		// --- Russian roulette ---
		u := rng.Float64()
		// Note: thresholds are explicit and use all three dynamic probabilities.
		if u < pAbs {
			if Debug && deposit {
				logRay("absorbed", Absorb, O, D, P, bounce, totalDist)
			}
			return false
		}

		// Survived → tint throughput by object color.
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
			if Debug && deposit {
				logRay("reflected", Reflect, O, D, P, bounce, totalDist)
			}
			continue
		}
		if u < pAbs+pReflDyn+pRefrDyn {
			// Refract
			Nn := N
			var eta Real
			if D.Dot(N) < 0 {
				eta = hit.iorInvCh(ch) // air -> material
			} else {
				Nn = N.Mul(-1)
				eta = hit.iorCh(ch) // material -> air
			}
			if T, ok := refract4(D, Nn, eta); ok {
				D = T
				if l2 := D.Dot(D); l2 > 0 {
					D = D.Mul(1 / math.Sqrt(l2))
				}
				O = Point4{P.X + D.X*bumpShift, P.Y + D.Y*bumpShift, P.Z + D.Z*bumpShift, P.W + D.W*bumpShift}
				if Debug && deposit {
					logRay("refracted", Refract, O, D, P, bounce, totalDist)
				}
				continue
			}
			// TIR → fall through to specular reflection
			D = reflect4(D, N)
			if l2 := D.Dot(D); l2 > 0 {
				D = D.Mul(1 / math.Sqrt(l2))
			}
			O = Point4{P.X + D.X*bumpShift, P.Y + D.Y*bumpShift, P.Z + D.Z*bumpShift, P.W + D.W*bumpShift}
			if Debug && deposit {
				logRay("total_internal_refraction", TIR, O, D, P, bounce, totalDist)
			}
			continue
		}
		if u < pAbs+pReflDyn+pRefrDyn+pDiffDyn {
			// Diffuse scatter (Lambertian-style on S^3 hemisphere around N)
			D = sampleDiffuseDir(N, rng)
			O = Point4{P.X + D.X*bumpShift, P.Y + D.Y*bumpShift, P.Z + D.Z*bumpShift, P.W + D.W*bumpShift}
			if Debug && deposit {
				logRay("diffuse", Diffuse, O, D, P, bounce, totalDist)
			}
			continue
		}
		// Should not happen (u in [0,1)), but keep safe fallback.
		return false
	}

	if Debug && deposit {
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

	// before launching workers
	start := time.Now()

	var counter int64 // incremented by workers
	// Zero-padding width based on number of slices.
	width := 1
	if totalRays > 1 {
		width = int(math.Log10(Real(totalRays-1))) + 1
	}
	fmt.Printf("[PROGRESS] %5.2f%% | rays=%*d/%*d | ETA=—\n", 0.0, width, 0, width, totalRays)

	// progress monitor (tick every 1s)
	done := make(chan struct{})
	monitorDone := make(chan struct{})
	go func() {
		defer close(monitorDone)
		if totalRays == 0 {
			return
		}
		ticker := time.NewTicker(2 * time.Second)
		defer ticker.Stop()

		lastPrintedAt := time.Now()
		lastPctBucketPrinted := -1 // 0..10 (for 0..100%), -1 means none yet

		for {
			select {
			case <-ticker.C:
				fired := atomic.LoadInt64(&counter)
				if fired >= int64(totalRays) {
					// Final line
					elapsed := time.Since(start)
					rate := float64(fired) / elapsed.Seconds()
					if rate < 1e-9 {
						rate = 0
					}
					fmt.Printf("[PROGRESS] 100.00%% | rays=%*d/%*d | rate=%.0f/s | elapsed=%s\n",
						width, fired, width, totalRays, rate, elapsed.Truncate(time.Millisecond))
					return
				}

				// Compute stats
				elapsed := time.Since(start)
				pct := float64(fired) * 100.0 / float64(totalRays)
				rate := float64(fired) / math.Max(elapsed.Seconds(), 1e-9)
				eta := time.Duration(0)
				if rate > 0 {
					rem := float64(totalRays-int(fired)) / rate
					eta = time.Duration(rem * float64(time.Second))
				}

				// Print if 10s since last print OR we crossed another 10% bucket
				tenSecPassed := time.Since(lastPrintedAt) >= 10*time.Second
				currentBucket := int(pct / 10.0) // 0..10
				crossed10pct := currentBucket != lastPctBucketPrinted

				if tenSecPassed || crossed10pct {
					// When printing due to bucket crossing, update the bucket.
					// When printing due to time, do NOT update the bucket so we still
					// get a message exactly at 10%, 20%, ... later.
					line := fmt.Sprintf("[PROGRESS] %5.2f%% | rays=%*d/%*d | rate=%.0f/s | ETA=%02d:%02d:%02d\n",
						pct, width, fired, width, totalRays, rate,
						int(eta.Hours()), int(eta.Minutes())%60, int(eta.Seconds())%60)
					fmt.Print(line)
					lastPrintedAt = time.Now()
					if crossed10pct {
						lastPctBucketPrinted = currentBucket
					}
				}

			case <-done:
				// One last summary before exiting (in case we didn’t hit 100% in the ticker loop)
				fired := atomic.LoadInt64(&counter)
				elapsed := time.Since(start)
				pct := float64(fired) * 100.0 / float64(totalRays)
				rate := float64(fired) / math.Max(elapsed.Seconds(), 1e-9)
				fmt.Printf("[PROGRESS] %5.2f%% | rays=%*d/%*d | rate=%.0f/s | elapsed=%s\n",
					pct, width, fired, width, totalRays, rate, elapsed.Truncate(time.Millisecond))
				return
			}
		}
	}()

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
					atomic.AddInt64(&counter, 1)
				}
			}
		}()
	}

	wg.Wait()
	close(done)   // stop the monitor
	<-monitorDone // wait for final progress line

}
