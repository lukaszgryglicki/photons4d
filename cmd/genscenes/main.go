// Command genscenes generates the photons4d "showcase" scene files.
//
// Each per-object scene frames a single 4D primitive in its best light: glassy /
// dispersive material, several wide cone lights above the W=0 deposition slab,
// and a genuine 4D rotation (with xw/yw/zw components) so that the per-Z-slice
// "video" sweeps through structurally different cross-sections of the object.
//
// A handful of "mixed" scenes pack every supported primitive into one frame for
// a grand tour, and one scene uses the environment-hypersphere (escape) mode.
//
// Resolution and samples-per-voxel are scaled up for a large machine: every
// scene is emitted at ResMultXY*ResMultXY*ResMultZ the base voxel resolution and
// SppMult the base samples-per-voxel (see the *Mult constants below).
//
// Usage:
//
//	go run ./cmd/genscenes              # write scenes/showcase/*.json
//	go run ./cmd/genscenes -outdir DIR  # write somewhere else
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
)

// Scale factors applied to every emitted scene (relative to the base values
// baked into each scene below). Tuned for a 512-core / 1-2 TB machine.
const (
	ResMultXY = 2 // multiply sceneResX and sceneResY
	ResMultZ  = 2 // multiply sceneResZ (number of video frames)
	SppMult   = 6 // multiply spp
)

// ---------- JSON value helpers ----------------------------------------------

type obj = map[string]any

func r5(x float64) float64 { return math.Round(x*1e5) / 1e5 }

func p4(x, y, z, w float64) obj { return obj{"X": r5(x), "Y": r5(y), "Z": r5(z), "W": r5(w)} }
func v4(x, y, z, w float64) obj { return obj{"X": r5(x), "Y": r5(y), "Z": r5(z), "W": r5(w)} }
func rgb(r, g, b float64) obj   { return obj{"R": r5(r), "G": r5(g), "B": r5(b)} }

func rot(xy, xz, xw, yz, yw, zw float64) obj {
	return obj{"xy": xy, "xz": xz, "xw": xw, "yz": yz, "yw": yw, "zw": zw}
}

// merge returns a fresh object combining geometry fields and material fields.
func merge(maps ...obj) obj {
	out := obj{}
	for _, m := range maps {
		for k, v := range m {
			out[k] = v
		}
	}
	return out
}

func material(color, diffuse, reflect, refract, ior obj) obj {
	return obj{"color": color, "diffuse": diffuse, "reflect": reflect, "refract": refract, "ior": ior}
}

// ---------- material presets -------------------------------------------------
// Per channel, reflect+refract+diffuse must be <= 1 (a little absorption left).

// glass: dispersive glass; blue bends more than red -> prismatic caustics.
func glass(tint obj, n, spread, refr, refl float64) obj {
	return material(tint, rgb(0, 0, 0), rgb(refl, refl, refl), rgb(refr, refr, refr),
		rgb(n-spread, n, n+spread))
}
func glassDefault(tint obj, n, spread float64) obj { return glass(tint, n, spread, 0.9, 0.07) }

func coloredGlass(tint obj, n, spread float64) obj { return glass(tint, n, spread, 0.88, 0.08) }

func mirror(tint obj) obj {
	return material(tint, rgb(0.04, 0.04, 0.04), rgb(0.82, 0.85, 0.9), rgb(0, 0, 0),
		rgb(1.6, 1.6, 1.6))
}

// ---------- light rigs -------------------------------------------------------

// studioLights: several wide cone lights from above (W+), only modestly tilted.
// All cones are wide enough that the bulk of each light's photons reach the W=0
// deposition slab directly (high hit-probability => bounded ray budget for a
// given samples-per-voxel target), regardless of object shape. The object in
// the beam carves a coloured caustic out of that flood of light. Tints
// (white / warm / cool / green) blend on the slab for rich colour.
func studioLights(mainIntensity, wide float64) []obj {
	return []obj{
		{"origin": p4(0, 0, 0, 1.25), "direction": v4(0, 0, 0, -1),
			"color": rgb(1, 1, 1), "angleDeg": wide, "intensity": mainIntensity},
		{"origin": p4(0, 0, 0, 1.18), "direction": v4(0.16, -0.10, 0.05, -1),
			"color": rgb(1.0, 0.58, 0.30), "angleDeg": wide - 6, "intensity": 0.5},
		{"origin": p4(0, 0, 0, 1.18), "direction": v4(-0.16, 0.10, -0.05, -1),
			"color": rgb(0.36, 0.64, 1.0), "angleDeg": wide - 6, "intensity": 0.5},
		{"origin": p4(0, 0, 0, 1.15), "direction": v4(0.02, 0.16, 0.12, -1),
			"color": rgb(0.5, 1.0, 0.55), "angleDeg": wide - 8, "intensity": 0.34},
	}
}

// ---------- scene writer -----------------------------------------------------

type sceneOpts struct {
	name    string
	objects obj // map of objectKey -> []obj
	lights  []obj
	gamma   float64
	bounces int
	baseXY  int
	baseZ   int
	baseSpp int
	width   float64
	height  float64
	depth   float64
	escape  bool
	probe   int
	delay   int
}

func (o sceneOpts) withDefaults() sceneOpts {
	if o.gamma == 0 {
		o.gamma = 0.7
	}
	if o.bounces == 0 {
		o.bounces = 42
	}
	if o.baseXY == 0 {
		o.baseXY = 300
	}
	if o.baseZ == 0 {
		o.baseZ = 64
	}
	if o.baseSpp == 0 {
		o.baseSpp = 128
	}
	if o.width == 0 {
		o.width = 2.2
	}
	if o.height == 0 {
		o.height = 2.2
	}
	if o.depth == 0 {
		o.depth = 2.2
	}
	if o.probe == 0 {
		o.probe = 200000
	}
	if o.delay == 0 {
		o.delay = 6
	}
	return o
}

func writeScene(outdir string, in sceneOpts) error {
	o := in.withDefaults()
	scn := obj{
		"center": p4(0, 0, 0, 0),
		"width":  o.width, "height": o.height, "depth": o.depth,
		"maxBounces": o.bounces,
	}
	if o.escape {
		scn["escape"] = true
	}
	cfg := obj{
		"sceneResX": o.baseXY * ResMultXY,
		"sceneResY": o.baseXY * ResMultXY,
		"sceneResZ": o.baseZ * ResMultZ,
		"probeRays": o.probe,
		"spp":       o.baseSpp * SppMult,
		"gifOut":    "./gifs/" + o.name + ".gif",
		"gifDelay":  o.delay,
		"gamma":     o.gamma,
		"scene":     scn,
		"lights":    o.lights,
	}
	for k, v := range o.objects {
		cfg[k] = v
	}

	data, err := json.MarshalIndent(cfg, "", "  ")
	if err != nil {
		return err
	}
	path := filepath.Join(outdir, o.name+".json")
	if err := os.WriteFile(path, append(data, '\n'), 0o644); err != nil {
		return err
	}
	keys := make([]string, 0)
	for k := range o.objects {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	fmt.Printf("wrote %-44s %v\n", path, keys)
	return nil
}

// ---------- hero placement / rotation ---------------------------------------

const heroW = 0.5

func heroRot(i int) obj {
	table := []obj{
		rot(18, -12, 22, 14, -16, 28),
		rot(0, 24, 18, 10, 26, -14),
		rot(28, 8, -20, -12, 18, 24),
		rot(12, -22, 30, 8, -10, 16),
	}
	return table[((i%len(table))+len(table))%len(table)]
}

// ring returns n positions on a circle in the XZ plane at height y and depth W=w.
func ring(n int, radius, y, w, phase float64) []obj {
	out := make([]obj, 0, n)
	for i := 0; i < n; i++ {
		a := phase + 2*math.Pi*float64(i)/float64(n)
		out = append(out, p4(radius*math.Cos(a), y, radius*math.Sin(a), w))
	}
	return out
}

func single(key string, o obj) obj { return obj{key: []obj{o}} }

// ---------- main -------------------------------------------------------------

func main() {
	outdir := flag.String("outdir", "scenes/showcase", "directory to write showcase scene JSON files")
	flag.Parse()

	if err := os.MkdirAll(*outdir, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, "mkdir:", err)
		os.Exit(1)
	}

	var scenes []sceneOpts

	add := func(s sceneOpts) { scenes = append(scenes, s) }

	// === Per-object showcase scenes ==========================================

	add(sceneOpts{name: "showcase_hypersphere", lights: studioLights(1.0, 34),
		objects: single("hyperspheres", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.58, 0.44, 0.62, 0.34), "rotDeg": heroRot(0),
		}, glassDefault(rgb(0.96, 0.99, 1.0), 1.52, 0.06)))})

	add(sceneOpts{name: "showcase_cell5", lights: studioLights(1.0, 34),
		objects: single("cells5", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.66, 0.66, 0.66, 0.66), "rotDeg": heroRot(1),
		}, glassDefault(rgb(1.0, 0.93, 0.85), 1.5, 0.05)))})

	add(sceneOpts{name: "showcase_cell8", lights: studioLights(1.0, 34),
		objects: single("cells8", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.74, 0.74, 0.74, 0.56), "rotDeg": heroRot(2),
		}, glassDefault(rgb(0.9, 0.97, 1.0), 1.5, 0.06)))})

	add(sceneOpts{name: "showcase_cell16", lights: studioLights(1.0, 34),
		objects: single("cells16", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.7, 0.7, 0.7, 0.7), "rotDeg": heroRot(3),
		}, glassDefault(rgb(1.0, 0.9, 0.95), 1.51, 0.055)))})

	add(sceneOpts{name: "showcase_cell24", lights: studioLights(1.0, 34),
		objects: single("cells24", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.64, 0.64, 0.64, 0.64), "rotDeg": heroRot(0),
		}, glassDefault(rgb(0.88, 1.0, 0.92), 1.52, 0.06)))})

	add(sceneOpts{name: "showcase_cell120", lights: studioLights(1.0, 34), baseSpp: 160,
		objects: single("cells120", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.62, 0.62, 0.62, 0.62), "rotDeg": heroRot(1),
		}, glassDefault(rgb(0.95, 0.92, 1.0), 1.54, 0.07)))})

	add(sceneOpts{name: "showcase_cell600", lights: studioLights(1.0, 34), baseSpp: 160,
		objects: single("cells600", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.62, 0.62, 0.62, 0.62), "rotDeg": heroRot(2),
		}, glassDefault(rgb(1.0, 0.95, 0.9), 1.54, 0.07)))})

	add(sceneOpts{name: "showcase_starcell", lights: studioLights(1.0, 34), baseSpp: 150, width: 2.6,
		objects: obj{"starCells": []obj{
			merge(obj{"kind": "star-cell-16", "center": p4(-0.62, 0.30, 0, heroW),
				"scale": v4(0.2, 0.2, 0.2, 0.2), "coreRadius": 0.4, "spikeLength": 0.62, "sharpness": 6.0,
				"rotDeg": rot(12, -14, 18, 16, -10, 22)}, coloredGlass(rgb(1.0, 0.55, 0.3), 1.52, 0.06)),
			merge(obj{"kind": "star-cell-24", "center": p4(0, -0.34, 0, heroW),
				"scale": v4(0.22, 0.22, 0.22, 0.22), "coreRadius": 0.42, "spikeLength": 0.6, "sharpness": 6.4,
				"rotDeg": rot(20, 10, -16, -12, 18, 26)}, coloredGlass(rgb(0.4, 0.85, 1.0), 1.52, 0.06)),
			merge(obj{"kind": "star-cell-600", "center": p4(0.62, 0.30, 0, heroW),
				"scale": v4(0.2, 0.2, 0.2, 0.2), "coreRadius": 0.46, "spikeLength": 0.56, "sharpness": 7.0,
				"rotDeg": rot(8, 22, 12, 18, 26, 14)}, coloredGlass(rgb(0.7, 1.0, 0.55), 1.52, 0.06)),
		}}})

	add(sceneOpts{name: "showcase_spherinder", lights: studioLights(1.0, 34),
		objects: single("spherinders", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.46, 0.46, 0.46, 0.42), "rotDeg": heroRot(3),
		}, glassDefault(rgb(0.92, 0.98, 1.0), 1.5, 0.05)))})

	// Cone shrunk a touch vs the others: a centred cone scatters more of the
	// beam, so a smaller cone keeps its hit-probability in line with the rest.
	add(sceneOpts{name: "showcase_hypercone", lights: studioLights(1.0, 36),
		objects: single("hypercones", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.44, 0.44, 0.44, 0.46), "rotDeg": heroRot(0),
		}, glassDefault(rgb(1.0, 0.9, 0.82), 1.52, 0.06)))})

	add(sceneOpts{name: "showcase_hypercapsule", lights: studioLights(1.0, 34),
		objects: single("hypercapsules", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.34, 0.34, 0.34, 0.3), "halfLength": 0.46,
			"rotDeg": heroRot(1),
		}, glassDefault(rgb(0.9, 0.95, 1.0), 1.5, 0.05)))})

	add(sceneOpts{name: "showcase_spheritorus", lights: studioLights(1.0, 34),
		objects: single("spheritori", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(1, 1, 1, 1),
			"majorRadius": 0.5, "minorRadius": 0.17, "rotDeg": heroRot(2),
		}, glassDefault(rgb(1.0, 0.92, 0.98), 1.52, 0.06)))})

	add(sceneOpts{name: "showcase_duotorus", lights: studioLights(1.0, 34), baseSpp: 150,
		objects: single("duotori", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(1, 1, 1, 1),
			"majorRadiusXY": 0.42, "majorRadiusZW": 0.32, "minorRadius": 0.12, "rotDeg": heroRot(3),
		}, glassDefault(rgb(0.9, 1.0, 0.95), 1.53, 0.065)))})

	add(sceneOpts{name: "showcase_duocylinder", lights: studioLights(1.0, 34),
		objects: single("duocylinders", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.46, 0.46, 0.46, 0.46), "rotDeg": heroRot(0),
		}, glassDefault(rgb(0.95, 0.95, 1.0), 1.5, 0.05)))})

	add(sceneOpts{name: "showcase_torisphere", lights: studioLights(1.0, 34), baseSpp: 150,
		objects: single("torispheres", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(1, 1, 1, 1),
			"majorRadius": 0.46, "minorRadius": 0.16, "rotDeg": heroRot(1),
		}, glassDefault(rgb(1.0, 0.96, 0.9), 1.52, 0.06)))})

	add(sceneOpts{name: "showcase_superquadric", lights: studioLights(1.0, 34),
		objects: single("superquadrics", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.56, 0.56, 0.56, 0.46), "power": 4.0,
			"rotDeg": heroRot(2),
		}, glassDefault(rgb(0.93, 0.97, 1.0), 1.51, 0.055)))})

	add(sceneOpts{name: "showcase_hyperfrustum", lights: studioLights(1.0, 34),
		objects: single("hyperfrustums", merge(obj{
			"center": p4(0, 0, 0, heroW), "scale": v4(0.5, 0.5, 0.5, 0.5),
			"lowerRadius": 0.34, "upperRadius": 0.78, "rotDeg": heroRot(3),
		}, glassDefault(rgb(1.0, 0.93, 0.88), 1.52, 0.06)))})

	// === Mixed scenes ========================================================

	add(sceneOpts{name: "showcase_mixed_glass_garden", lights: studioLights(1.0, 28),
		gamma: 0.72, bounces: 48, baseSpp: 128, width: 2.4, height: 2.4, depth: 2.4,
		objects: glassGarden()})

	add(sceneOpts{name: "showcase_mixed_prism_chamber", lights: prismLights(),
		gamma: 0.68, bounces: 56, baseXY: 320, baseZ: 72, baseSpp: 160,
		width: 2.4, height: 2.4, depth: 2.4, objects: prismChamber()})

	add(sceneOpts{name: "showcase_mixed_grand",
		lights: append(studioLights(1.05, 30), obj{
			"origin": p4(0, -0.48, -0.40, 1.10), "direction": v4(0, 0.30, 0.26, -1),
			"color": rgb(0.9, 0.5, 1.0), "angleDeg": 12.0, "intensity": 0.28}),
		gamma: 0.74, bounces: 56, baseZ: 72, baseSpp: 140,
		width: 2.5, height: 2.5, depth: 2.5, objects: grandTour()})

	add(sceneOpts{name: "showcase_escape_environment", lights: escapeLights(), escape: true,
		gamma: 0.72, bounces: 48, baseSpp: 160, width: 2.0, height: 2.0, depth: 2.0,
		objects: escapeObjs()})

	for _, s := range scenes {
		if err := writeScene(*outdir, s); err != nil {
			fmt.Fprintln(os.Stderr, "write", s.name, ":", err)
			os.Exit(1)
		}
	}
	fmt.Printf("\nDone: %d scenes -> %s  (res x%d/%d/%d, spp x%d)\n",
		len(scenes), *outdir, ResMultXY, ResMultXY, ResMultZ, SppMult)
}

// ---------- mixed-scene object sets -----------------------------------------

func glassGarden() obj {
	outer := ring(9, 0.78, 0.22, heroW, 0.1)
	inner := ring(8, 0.40, -0.24, heroW, 0.5)
	tints := []obj{
		rgb(1.0, 0.55, 0.4), rgb(0.4, 0.8, 1.0), rgb(0.6, 1.0, 0.55), rgb(1.0, 0.85, 0.4),
		rgb(0.85, 0.5, 1.0), rgb(0.4, 1.0, 0.85), rgb(1.0, 0.45, 0.65), rgb(0.5, 0.7, 1.0),
		rgb(0.95, 0.97, 1.0),
	}
	gm := func(i int) obj { return coloredGlass(tints[i%len(tints)], 1.5+0.02*float64(i%4), 0.06) }

	objs := obj{}
	push := func(key string, o obj) {
		if arr, ok := objs[key].([]obj); ok {
			objs[key] = append(arr, o)
		} else {
			objs[key] = []obj{o}
		}
	}

	push("cells8", merge(obj{"center": outer[0], "scale": v4(0.3, 0.3, 0.3, 0.26), "rotDeg": heroRot(0)}, gm(0)))
	push("cells16", merge(obj{"center": outer[1], "scale": v4(0.32, 0.32, 0.32, 0.32), "rotDeg": heroRot(1)}, gm(1)))
	push("cells24", merge(obj{"center": outer[2], "scale": v4(0.3, 0.3, 0.3, 0.3), "rotDeg": heroRot(2)}, gm(2)))
	push("cells120", merge(obj{"center": outer[3], "scale": v4(0.3, 0.3, 0.3, 0.3), "rotDeg": heroRot(3)}, gm(3)))
	push("cells600", merge(obj{"center": outer[4], "scale": v4(0.3, 0.3, 0.3, 0.3), "rotDeg": heroRot(0)}, gm(4)))
	push("spheritori", merge(obj{"center": outer[5], "scale": v4(1, 1, 1, 1), "majorRadius": 0.24, "minorRadius": 0.09, "rotDeg": heroRot(1)}, gm(5)))
	push("torispheres", merge(obj{"center": outer[6], "scale": v4(1, 1, 1, 1), "majorRadius": 0.24, "minorRadius": 0.09, "rotDeg": heroRot(2)}, gm(6)))
	push("duotori", merge(obj{"center": outer[7], "scale": v4(1, 1, 1, 1), "majorRadiusXY": 0.22, "majorRadiusZW": 0.16, "minorRadius": 0.07, "rotDeg": heroRot(3)}, gm(7)))
	push("cells5", merge(obj{"center": outer[8], "scale": v4(0.34, 0.34, 0.34, 0.34), "rotDeg": heroRot(0)}, gm(8)))

	push("hyperspheres", merge(obj{"center": inner[0], "scale": v4(0.22, 0.18, 0.24, 0.16), "rotDeg": heroRot(1)}, gm(1)))
	push("spherinders", merge(obj{"center": inner[1], "scale": v4(0.2, 0.2, 0.2, 0.18), "rotDeg": heroRot(2)}, gm(2)))
	push("hypercones", merge(obj{"center": inner[2], "scale": v4(0.22, 0.22, 0.22, 0.22), "rotDeg": heroRot(3)}, gm(3)))
	push("hypercapsules", merge(obj{"center": inner[3], "scale": v4(0.16, 0.16, 0.16, 0.14), "halfLength": 0.2, "rotDeg": heroRot(0)}, gm(4)))
	push("duocylinders", merge(obj{"center": inner[4], "scale": v4(0.2, 0.2, 0.2, 0.2), "rotDeg": heroRot(1)}, gm(5)))
	push("superquadrics", merge(obj{"center": inner[5], "scale": v4(0.22, 0.22, 0.22, 0.2), "power": 4.0, "rotDeg": heroRot(2)}, gm(6)))
	push("hyperfrustums", merge(obj{"center": inner[6], "scale": v4(0.22, 0.22, 0.22, 0.22), "lowerRadius": 0.3, "upperRadius": 0.7, "rotDeg": heroRot(3)}, gm(7)))
	push("starCells", merge(obj{"kind": "star-cell-24", "center": inner[7], "scale": v4(0.12, 0.12, 0.12, 0.12),
		"coreRadius": 0.42, "spikeLength": 0.58, "sharpness": 6.4, "rotDeg": heroRot(0)}, gm(0)))
	return objs
}

func prismChamber() obj {
	objs := obj{}
	objs["hyperspheres"] = []obj{merge(obj{"center": p4(0, 0, 0, 0.52), "scale": v4(0.5, 0.5, 0.5, 0.3),
		"rotDeg": heroRot(0)}, glass(rgb(0.98, 0.99, 1.0), 1.6, 0.1, 0.92, 0.05))}
	objs["cells8"] = []obj{merge(obj{"center": p4(-0.7, 0.25, 0.1, 0.55), "scale": v4(0.34, 0.34, 0.34, 0.3),
		"rotDeg": heroRot(1)}, mirror(rgb(1.0, 0.85, 0.7)))}
	objs["cells24"] = []obj{merge(obj{"center": p4(0.7, 0.25, -0.1, 0.55), "scale": v4(0.32, 0.32, 0.32, 0.32),
		"rotDeg": heroRot(2)}, mirror(rgb(0.7, 0.85, 1.0)))}
	objs["superquadrics"] = []obj{merge(obj{"center": p4(-0.55, -0.4, 0.2, 0.5), "scale": v4(0.24, 0.24, 0.24, 0.2),
		"power": 6.0, "rotDeg": heroRot(3)}, coloredGlass(rgb(1.0, 0.4, 0.45), 1.55, 0.09))}
	objs["spheritori"] = []obj{merge(obj{"center": p4(0.55, -0.4, -0.2, 0.5), "scale": v4(1, 1, 1, 1),
		"majorRadius": 0.26, "minorRadius": 0.09, "rotDeg": heroRot(0)}, coloredGlass(rgb(0.45, 1.0, 0.65), 1.55, 0.09))}
	objs["duotori"] = []obj{merge(obj{"center": p4(0, 0.5, 0, 0.5), "scale": v4(1, 1, 1, 1),
		"majorRadiusXY": 0.22, "majorRadiusZW": 0.18, "minorRadius": 0.07, "rotDeg": heroRot(1)},
		coloredGlass(rgb(0.5, 0.6, 1.0), 1.55, 0.09))}
	objs["hypercones"] = []obj{merge(obj{"center": p4(0, -0.05, 0.55, 0.5), "scale": v4(0.22, 0.22, 0.22, 0.24),
		"rotDeg": heroRot(2)}, coloredGlass(rgb(1.0, 0.8, 0.35), 1.55, 0.09))}
	return objs
}

// prismLights: tight beams aimed through the central lens for sharp, dramatic
// dispersion fans (lower hit-probability => heavier, but that is the point).
func prismLights() []obj {
	return []obj{
		{"origin": p4(0, 0, 0, 1.25), "direction": v4(0, 0, 0, -1), "color": rgb(1, 1, 1), "angleDeg": 12.0, "intensity": 1.0},
		{"origin": p4(0.9, 0, 0, 1.0), "direction": v4(-0.5, 0, 0, -1), "color": rgb(1.0, 0.3, 0.3), "angleDeg": 9.0, "intensity": 0.55},
		{"origin": p4(-0.9, 0, 0, 1.0), "direction": v4(0.5, 0, 0, -1), "color": rgb(0.3, 0.5, 1.0), "angleDeg": 9.0, "intensity": 0.55},
		{"origin": p4(0, 0.9, 0, 1.0), "direction": v4(0, -0.5, 0, -1), "color": rgb(0.35, 1.0, 0.4), "angleDeg": 9.0, "intensity": 0.5},
		{"origin": p4(0, -0.9, 0, 1.0), "direction": v4(0, 0.5, 0, -1), "color": rgb(1.0, 0.85, 0.3), "angleDeg": 9.0, "intensity": 0.5},
	}
}

func grandTour() obj {
	objs := glassGarden()
	push := func(key string, o obj) {
		if arr, ok := objs[key].([]obj); ok {
			objs[key] = append(arr, o)
		} else {
			objs[key] = []obj{o}
		}
	}
	push("cells120", merge(obj{"center": p4(0, 0, 0, 0.52), "scale": v4(0.34, 0.34, 0.34, 0.34),
		"rotDeg": heroRot(2)}, glassDefault(rgb(0.97, 0.95, 1.0), 1.56, 0.08)))
	push("cells600", merge(obj{"center": p4(0, 0, 0, 0.86), "scale": v4(0.42, 0.42, 0.42, 0.22),
		"rotDeg": heroRot(3)}, glassDefault(rgb(0.9, 0.95, 1.0), 1.5, 0.05)))
	push("starCells", merge(obj{"kind": "star-cell-120", "center": p4(0, 0, 0, 0.2),
		"scale": v4(0.16, 0.16, 0.16, 0.16), "coreRadius": 0.44, "spikeLength": 0.56, "sharpness": 7.0,
		"rotDeg": heroRot(1)}, coloredGlass(rgb(1.0, 0.6, 0.85), 1.55, 0.09)))
	return objs
}

func escapeObjs() obj {
	return obj{
		"cells600": []obj{merge(obj{"center": p4(0, 0, 0, 0), "scale": v4(0.5, 0.5, 0.5, 0.5),
			"rotDeg": heroRot(0)}, glassDefault(rgb(0.95, 0.97, 1.0), 1.55, 0.08))},
		"hyperspheres": []obj{merge(obj{"center": p4(0.55, 0.2, -0.15, 0.1), "scale": v4(0.22, 0.22, 0.22, 0.22),
			"rotDeg": heroRot(1)}, coloredGlass(rgb(1.0, 0.5, 0.4), 1.5, 0.08))},
		"starCells": []obj{merge(obj{"kind": "star-cell-24", "center": p4(-0.5, -0.2, 0.2, -0.1),
			"scale": v4(0.16, 0.16, 0.16, 0.16), "coreRadius": 0.42, "spikeLength": 0.58, "sharpness": 6.4,
			"rotDeg": heroRot(2)}, coloredGlass(rgb(0.45, 0.85, 1.0), 1.52, 0.06))},
		"duotori": []obj{merge(obj{"center": p4(0.1, 0.5, 0.3, 0), "scale": v4(1, 1, 1, 1),
			"majorRadiusXY": 0.2, "majorRadiusZW": 0.16, "minorRadius": 0.07, "rotDeg": heroRot(3)},
			coloredGlass(rgb(0.5, 1.0, 0.6), 1.52, 0.06))},
	}
}

func escapeLights() []obj {
	return []obj{
		{"origin": p4(0, 0, 0, 1.4), "direction": v4(0, 0, 0, -1), "color": rgb(1, 1, 1), "angleDeg": 80.0, "intensity": 1.0},
		{"origin": p4(1.1, 0.6, 0, 0.4), "direction": v4(-1.0, -0.5, 0, -0.3), "color": rgb(1.0, 0.5, 0.35), "angleDeg": 40.0, "intensity": 0.5},
		{"origin": p4(-1.1, -0.6, 0.4, 0.4), "direction": v4(1.0, 0.5, -0.3, -0.3), "color": rgb(0.4, 0.6, 1.0), "angleDeg": 40.0, "intensity": 0.5},
	}
}
