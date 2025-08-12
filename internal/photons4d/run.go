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

	for _, hc := range cfg.Hypercubes {
		h, err := hc.Build()
		if err != nil {
			continue
		}
		scene.AddHypercube(h)
	}

	for i, scfg := range cfg.Hyperspheres {
		h, err := scfg.Build()
		if err != nil {
			_ = i // keep if you log
			continue
		}
		scene.AddHyperSphere(h)
	}

	Nvox := Nx * Ny * Nz
	needRays := make([]int, len(lights))
	totalRays := 0
	for i, L := range lights {
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
		DebugLog("Light #%d, needs: %d rays", i, need)
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
		// Save GIF
		if err := SaveAnimatedGIF(scene, cfg.GIFOut, cfg.GIFDelay, cfg.Gamma); err != nil {
			panic(err)
		}
		DebugLog("Saved animated GIF: %s", cfg.GIFOut)
	}
	return nil
}
