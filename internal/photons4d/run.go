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
		DebugLog("Debug mode: reduced resolution to %d x %d x %d, spp to %d and max bounces to %d", cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ, cfg.Spp, cfg.Scene.MaxBounces)
	}

	Nx, Ny, Nz := cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ
	scene := NewScene(cfg.Scene.Center, cfg.Scene.Width, cfg.Scene.Height, cfg.Scene.Depth, Nx, Ny, Nz, cfg.Scene.MaxBounces)

	lights := make([]*Light, 0, len(cfg.Lights))
	for _, Lc := range cfg.Lights {
		angle := Lc.AngleDeg * math.Pi / 180.0
		L, err := NewLight(Lc.Origin, Lc.Direction, Lc.Color, angle, Lc.VoidLight)
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
		scene.AddTwentyFourCell(obj)
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

	start := time.Now()
	castRays(lights, scene, needRays)
	elapsed := time.Since(start)
	DebugLog("Rays: %d, time: %s", totalRays, elapsed)

	if Debug {
		raysStats()
	}

	if PNG {
		prefix := strings.Replace(cfg.GIFOut, ".gif", "", 1)
		prefix = strings.Replace(prefix, "gifs/", "pngs/", 1)
		if err := SavePNGSequence16(scene, prefix, cfg.Gamma); err != nil {
			panic(err)
		}
		DebugLog("Saved PNG sequence with prefix: %s", prefix)
	} else {
		if err := SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma); err != nil {
			panic(err)
		}
		DebugLog("Saved animated GIF: %s", cfg.GIFOut)
	}
	return nil
}
