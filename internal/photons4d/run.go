package photons4d

import (
	"math"
	"time"
)

func Run(cfgPath string) error {
	cfg, err := loadConfig(cfgPath)
	if err != nil {
		return err
	}

	Nx, Ny, Nz := cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ
	scene := NewScene(cfg.Scene.Center, cfg.Scene.Width, cfg.Scene.Height, cfg.Scene.Depth, Nx, Ny, Nz, cfg.Scene.MaxBounces)

	lights := make([]*Light, 0, len(cfg.Lights))
	for i, Lc := range cfg.Lights {
		angle := Lc.AngleDeg * math.Pi / 180.0
		L, err := NewLight(Lc.Origin, Lc.Direction, Lc.Color, angle, Lc.VoidLight)
		if err != nil {
			return err
		}
		_ = i
		lights = append(lights, L)
	}

	for i, hc := range cfg.Hypercubes {
		h, err := hc.Build()
		if err != nil {
			_ = i // keep if you print debug
			continue
		}
		scene.AddHypercube(h)
	}

	Nvox := Nx * Ny * Nz
	needRays := make([]int, len(lights))
	totalRays := 0
	for i, L := range lights {
		_ = i
		p := estimateHitProb(L, scene, cfg.ProbeRays)
		if p < 1e-7 {
			p = 1e-7
		}
		need := int(3 * Real(cfg.Spp) * Real(Nvox) / p)
		if need < cfg.ProbeRays {
			need = cfg.ProbeRays
		}
		needRays[i] = need
		totalRays += need
	}

	start := time.Now()
	castRays(lights, scene, needRays)
	elapsed := time.Since(start)
	DebugLog("Rays: %d, time: %s", totalRays, elapsed)

	if Debug {
		raysStats()
	}

	// Save GIF
	if err := SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma); err != nil {
		panic(err)
	}
	DebugLog("Saved animated GIF: %s", cfg.GIFOut)
	return nil
}
