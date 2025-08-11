
func castSingleRay(light *ConeLight4, scene *Scene3D, rng *rand.Rand, locks *shardLocks, deposit bool) bool {
	// Choose the photon channel; if all zero, no photon.
	ch := pickChannelFast(light.colorSum, light.thrR, light.thrG, rng)

	// Initial throughput/energy: sum of channels so that E[deposit] = Color[ch]*w
	throughput := light.Color.R + light.Color.G + light.Color.B

	// Start at the light with a sampled unit direction.
	O := light.Origin
	D := light.SampleDir(rng) // approx unit

	totalDist := 0.0

	for bounce := 0; bounce < scene.MaxBounces; bounce++ {
		// optional periodic re-normalization for numerical stability
		if (bounce & 3) == 3 {
			l2 := D.Dot(D)
			if l2 > 0 {
				D = D.Mul(1 / math.Sqrt(l2))
			}
		}

		// Next plane hit (W == scene.Center.W)
		tPlane := planeHitT(scene, O, D)

		// Next hypercube hit (with AABB pre-cull)
		hit, okCube := nearestCube(scene, O, D)

		// If no cube ahead or plane is closer → try to deposit on the plane.
		if !okCube || tPlane < hit.t {
			if !isFinite(tPlane) {
				if Debug {
					logRay("parallel_to_scene", Miss, O, D, Point4{}, bounce, totalDist)
				}
				return false // parallel to plane and no cube to stop us
			}

			P := O.Add(D.Mul(tPlane))
			if ok, i, j, k, _, _, _ := scene.VoxelIndexOf(P); ok {
				totalDist += tPlane
				w := 1.0 / (totalDist*totalDist + epsDist)

				if deposit {
					base := scene.idx(i, j, k, ChR)
					if locks != nil {
						locks.lock(base)
					}
					scene.Buf[base+ch] += Real(throughput * w * light.voidCoeff)
					if locks != nil {
						locks.unlock(base)
					}
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

		hc := hit.hc

		// Precompute outward normal (unit)
		N := hit.Nw

		// Per-channel base absorption (cached).
		pAbs := hc.pAbs[ch]

		// Fresnel-Schlick per channel, scaled by available non-absorption energy.
		// cosθ depends on whether we're entering (outside->inside) or exiting (inside->outside).
		var cosTheta Real
		if hit.inv {
			// exiting: angle between incident ray and outward normal
			cosTheta = D.Dot(N)
		} else {
			// entering: use -D·N
			cosTheta = -D.Dot(N)
		}
		if cosTheta < 0 {
			cosTheta = 0
		} else if cosTheta > 1 {
			cosTheta = 1
		}
		F := hc.f0[ch] + (1-hc.f0[ch])*math.Pow(1-cosTheta, 5)

		avail := 1 - pAbs     // total energy available to split between refl/refr
		pReflDyn := avail * F // dynamic reflect prob via Fresnel
		// pRefrDyn := avail - pReflDyn // not needed explicitly for roulette

		// Russian roulette: absorb / reflect / refract.
		u := rng.Float64()
		if u < pAbs {
			if Debug {
				logRay("absorbed", Absorb, O, D, P, bounce, totalDist)
			}
			return false // absorbed
		}

		// Tint throughput by the cube color (per channel) after survival.
		throughput *= hc.colorArr[ch]

		if u < pAbs+pReflDyn {
			// Reflect (no renorm needed if inputs are unit)
			D = reflect4(D, N)
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
			eta = hc.iorArr[ch]
		} else {
			// entering: outside -> inside
			eta = hc.iorInv[ch]
		}

		if T, ok := refract4(D, N, eta); ok {
			D = T
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
