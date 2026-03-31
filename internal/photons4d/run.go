package photons4d

import (
	"fmt"
	"math"
	"path/filepath"
	"strings"
	"time"
)

func recoveryRawPathFor(out string) string {
	base := filepath.Base(out)
	ext := filepath.Ext(base)
	name := strings.TrimSuffix(base, ext)
	if name == "" {
		name = "scene"
	}
	return name + ".recovery.raw"
}

func saveWithRetry(label string, attempts int, fn func() error) error {
	var lastErr error
	for attempt := 1; attempt <= attempts; attempt++ {
		if err := fn(); err == nil {
			return nil
		} else {
			lastErr = err
			DebugLog("Save attempt %d/%d for %s failed: %v", attempt, attempts, label, err)
		}
		if attempt < attempts {
			time.Sleep(time.Duration(attempt) * 250 * time.Millisecond)
		}
	}
	return lastErr
}

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

	for i, hc := range cfg.Cells8 {
		h, err := hc.Build()
		if err != nil {
			panic(fmt.Errorf("cells8[%d]: %w", i, err))
		}
		scene.AddCell8(h)
	}

	for i, scfg := range cfg.Hyperspheres {
		h, err := scfg.Build()
		if err != nil {
			panic(fmt.Errorf("hyperspheres[%d]: %w", i, err))
		}
		scene.AddHyperSphere(h)
	}

	for i, scfg := range cfg.Cells5 {
		sx, err := scfg.Build()
		if err != nil {
			panic(fmt.Errorf("cells5[%d]: %w", i, err))
		}
		scene.AddCell5(sx)
	}

	for i, scfg := range cfg.Cells16 {
		obj, err := scfg.Build()
		if err != nil {
			panic(fmt.Errorf("cells16[%d]: %w", i, err))
		}
		scene.AddCell16(obj)
	}

	for i, scfg := range cfg.Cells24 {
		obj, err := scfg.Build()
		if err != nil {
			panic(fmt.Errorf("cells24[%d]: %w", i, err))
		}
		scene.AddCell24(obj)
	}

	for i, cc := range cfg.Cells120 {
		h, err := cc.Build()
		if err != nil {
			panic(fmt.Errorf("cells120[%d]: %w", i, err))
		}
		scene.AddCell120(h)
	}

	for i, cc := range cfg.Cells600 {
		h, err := cc.Build()
		if err != nil {
			panic(fmt.Errorf("cells600[%d]: %w", i, err))
		}
		scene.AddCell600(h)
	}
	for i, scfg := range cfg.Spherinders {
		h, err := scfg.Build()
		if err != nil {
			panic(fmt.Errorf("spherinders[%d]: %w", i, err))
		}
		scene.AddSpherinder(h)
	}
	for i, hcfg := range cfg.HyperCones {
		h, err := hcfg.Build()
		if err != nil {
			panic(fmt.Errorf("hypercones[%d]: %w", i, err))
		}
		scene.AddHyperCone(h)
	}
	for i, hcfg := range cfg.HyperCapsules {
		h, err := hcfg.Build()
		if err != nil {
			panic(fmt.Errorf("hypercapsules[%d]: %w", i, err))
		}
		scene.AddHyperCapsule(h)
	}
	for i, tcfg := range cfg.Spheritori {
		h, err := tcfg.Build()
		if err != nil {
			panic(fmt.Errorf("spheritori[%d]: %w", i, err))
		}
		scene.AddSpheritorus(h)
	}
	nObjects := scene.NObjects()
	DebugLog("Scene created with %d objects", nObjects)
	if AlwaysBVH && NeverBVH {
		// Prefer explicit “always on” if both are set
		DebugLog("Both ALWAYS_BVH and NEVER_BVH are set; ALWAYS_BVH wins")
		NeverBVH = false
	}
	if AlwaysBVH {
		root := getOrBuildBVH(scene)
		DebugLog("AlwaysBVH is set, using BVH of AABBs")
		// NearestHitFunc = nearestHitBVH
		NearestHitFunc = func(_ *Scene, O Point4, D Vector4, tMax Real) (objectHit, bool) {
			return traverseNearest(root, O, D, tMax)
		}
	} else if NeverBVH {
		DebugLog("NeverBVH is set, using nearestHit function")
		NearestHitFunc = nearestHit
	} else {
		if nObjects < AABBBVHFromNObjects {
			NearestHitFunc = nearestHit
			DebugLog("Using nearestHit function (instead of BVH of AABB) for %d objects", nObjects)
		} else {
			//NearestHitFunc = nearestHitBVH
			root := getOrBuildBVH(scene)
			NearestHitFunc = func(_ *Scene, O Point4, D Vector4, tMax Real) (objectHit, bool) {
				return traverseNearest(root, O, D, tMax)
			}
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

	recoveryRaw := recoveryRawPathFor(cfg.GIFOut)
	if err := saveWithRetry("animated GIF", 3, func() error {
		return SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma)
	}); err != nil {
		if rawErr := scene.SaveRawRGB64(recoveryRaw); rawErr != nil {
			return fmt.Errorf("save animated GIF %q: %w; recovery raw save %q also failed: %v", cfg.GIFOut, err, recoveryRaw, rawErr)
		}
		return fmt.Errorf("save animated GIF %q: %w; recovery raw saved to %q", cfg.GIFOut, err, recoveryRaw)
	}
	DebugLog("Saved animated GIF: %s", cfg.GIFOut)
	if PNG {
		prefix := strings.Replace(cfg.GIFOut, ".gif", "", 1)
		prefix = strings.Replace(prefix, "gifs/", "pngs/", 1)
		if err := saveWithRetry("PNG sequence", 3, func() error {
			return SavePNGSequence16(scene, prefix, cfg.Gamma)
		}); err != nil {
			if rawErr := scene.SaveRawRGB64(recoveryRaw); rawErr != nil {
				return fmt.Errorf("save PNG sequence %q: %w; recovery raw save %q also failed: %v", prefix, err, recoveryRaw, rawErr)
			}
			return fmt.Errorf("save PNG sequence %q: %w; recovery raw saved to %q", prefix, err, recoveryRaw)
		}
		DebugLog("Saved PNG sequence with prefix: %s", prefix)
	}

	if RAW {
		fn := strings.Replace(cfg.GIFOut, ".gif", ".raw", 1)
		fn = strings.Replace(fn, "gifs/", "raws/", 1)
		if err := saveWithRetry("raw scene", 3, func() error {
			return scene.SaveRawRGB64(fn)
		}); err != nil {
			fallback := recoveryRawPathFor(fn)
			if rawErr := scene.SaveRawRGB64(fallback); rawErr != nil {
				return fmt.Errorf("save raw scene %q: %w; fallback raw save %q also failed: %v", fn, err, fallback, rawErr)
			}
			return fmt.Errorf("save raw scene %q: %w; fallback raw saved to %q", fn, err, fallback)
		}
		DebugLog("Saved RAW scene: %s", fn)
	}
	return nil
}
