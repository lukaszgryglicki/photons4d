package photons4d

import (
	"math"
	"strings"
	"time"
)

func Run(cfgPath string) error {
	cfg, err := loadConfig(cfgPath)
	if err != nil {
		return err
	}

	if Debug {
		cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ, cfg.Spp, cfg.Scene.MaxBounces =
			imax(cfg.SceneResX>>3, 4),
			imax(cfg.SceneResY>>3, 4),
			imax(cfg.SceneResZ>>2, 1),
			imax(cfg.Spp>>2, 4),
			imax(cfg.Scene.MaxBounces>>1, 8)
		cfg.GIFOut = strings.Replace(cfg.GIFOut, ".gif", "_debug.gif", 1)
		DebugLog("Debug mode: reduced resolution to %d x %d x %d, spp to %d and max bounces to %d", cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ, cfg.Spp, cfg.Scene.MaxBounces)
	}

	Nx, Ny, Nz := cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ
	scene := NewScene(cfg.Scene.Center, cfg.Scene.Width, cfg.Scene.Height, cfg.Scene.Depth, Nx, Ny, Nz, cfg.Scene.MaxBounces, cfg.Scene.EnvHypersphere)

	lights := make([]*Light, 0, len(cfg.Lights))
	for _, Lc := range cfg.Lights {
		angle := Lc.AngleDeg * math.Pi / 180.0
		L, err := NewLight(Lc.Origin, Lc.Direction, Lc.Color, angle, Lc.Intensity)
		if err != nil {
			return err
		}
		lights = append(lights, L)
	}

	for _, hc := range cfg.Cells8 {
		h, err := hc.Build()
		if err != nil {
			continue
		}
		scene.AddCell8(h)
	}

	for _, scfg := range cfg.Hyperspheres {
		h, err := scfg.Build()
		if err != nil {
			continue
		}
		scene.AddHyperSphere(h)
	}

	for _, scfg := range cfg.Cells5 {
		sx, err := scfg.Build()
		if err != nil {
			continue
		}
		scene.AddCell5(sx)
	}

	for _, scfg := range cfg.Cells16 {
		obj, err := scfg.Build()
		if err != nil {
			continue
		}
		scene.AddCell16(obj)
	}

	for _, scfg := range cfg.Cells24 {
		obj, err := scfg.Build()
		if err != nil {
			continue
		}
		scene.AddCell24(obj)
	}

	for _, cc := range cfg.Cells120 {
		h, err := cc.Build()
		if err != nil {
			continue
		}
		scene.AddCell120(h)
	}

	for _, cc := range cfg.Cells600 {
		h, err := cc.Build()
		if err != nil {
			continue
		}
		scene.AddCell600(h)
	}
	nObjects := scene.NObjects()
	DebugLog("Scene created with %d objects", nObjects)
	if AlwaysBVH && NeverBVH {
		// Prefer explicit “always on” if both are set
		DebugLog("Both ALWAYS_BVH and NEVER_BVH are set; ALWAYS_BVH wins")
		NeverBVH = false
	}
	if AlwaysBVH {
		DebugLog("AlwaysBVH is set, using BVH of AABBs")
		NearestHitFunc = nearestHitBVH
	} else if NeverBVH {
		DebugLog("NeverBVH is set, using nearestHit function")
		NearestHitFunc = nearestHit
	} else {
		if nObjects < AABBBVHFromNObjects {
			NearestHitFunc = nearestHit
			DebugLog("Using nearestHit function (instead of BVH of AABB) for %d objects", nObjects)
		} else {
			NearestHitFunc = nearestHitBVH
			DebugLog("Using BVH of AABBs for %d objects", nObjects)
		}
	}

	Nvox := Nx * Ny * Nz
	needRays := make([]int, len(lights))
	totalRays := 0
	for i, L := range lights {
		p := estimateHitProb(L, scene, cfg.ProbeRays)
		if p < 1e-7 {
			DebugLog("Light #%d, hit probability too low: %.12f, setting to 1e-7", i, p)
			p = 1e-7
		}
		need := int(3 * Real(cfg.Spp) * Real(Nvox) / p)
		if need < cfg.ProbeRays {
			need = cfg.ProbeRays
		}
		needRays[i] = need
		totalRays += need
		DebugLog("Light #%d, needs: %d rays, scene hit probability %.12f", i, need, p)
	}
	DebugLog("Total rays needed: %d", totalRays)
	if AlwaysBVH || nObjects >= AABBBVHFromNObjects {
		DumpAABBBVH(scene, false)
	}

	start := time.Now()
	castRays(lights, scene, needRays)
	elapsed := time.Since(start)
	DebugLog("Rays: %d, time: %s", totalRays, elapsed)

	if Debug {
		raysStats()
		if !DumpAABBBVH(scene, false) {
			DebugLog("Current scene does not use AABB BVH, here is what it will look like if used:")
			DumpAABBBVH(scene, true)
		}
	}

	if err := SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma); err != nil {
		panic(err)
	}
	DebugLog("Saved animated GIF: %s", cfg.GIFOut)
	if PNG {
		prefix := strings.Replace(cfg.GIFOut, ".gif", "", 1)
		prefix = strings.Replace(prefix, "gifs/", "pngs/", 1)
		if err := SavePNGSequence16(scene, prefix, cfg.Gamma); err != nil {
			panic(err)
		}
		DebugLog("Saved PNG sequence with prefix: %s", prefix)
	}

	if RAW {
		fn := strings.Replace(cfg.GIFOut, ".gif", ".raw", 1)
		fn = strings.Replace(fn, "gifs/", "raws/", 1)
		if err := scene.SaveRawRGB64(fn); err != nil {
			panic(err)
		}
		DebugLog("Saved RAW scene: %s", fn)
	}
	return nil
}
