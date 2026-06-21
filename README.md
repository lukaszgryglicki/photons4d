# photons4d

A 4D forward photon-tracing renderer. Lights emit photons into a 4D scene; the
energy that lands on the 3D hyperplane slab at `W = scene.center.W` is recorded
as a voxel volume and exported as an animated GIF / 16-bit PNG sequence / MP4,
one frame per Z slice.

## Build

```
make build       # debug build (verbose logging, used by tests)
make release     # optimized binary
make test        # go test ./internal/photons4d
```

## Showcase scenes

`scenes/showcase/*.json` are generated, ready-to-render scenes — one per
supported 4D primitive (hypersphere, 5/8/16/24/120/600-cell, star cells,
spherinder, hypercone, hypercapsule, spheritorus, duotorus, duocylinder,
torisphere, superquadric, hyperfrustum) plus four "mixed" scenes
(`glass_garden`, `prism_chamber`, `grand`, and an `escape_environment` scene
that uses the environment-hypersphere mode). Each frames its object(s) in
dispersive glass under several coloured cone lights, so the per-Z-slice video
sweeps through the object's changing 4D cross sections as prismatic caustics.

Regenerate them with the Go generator (edit `cmd/genscenes/main.go` to tweak):

```
go run ./cmd/genscenes               # writes scenes/showcase/*.json
```

Resolution and samples-per-voxel are scaled for a large machine via the
`ResMultXY` / `ResMultZ` / `SppMult` constants in `cmd/genscenes/main.go`
(currently 2× resolution on each axis and 6× samples-per-voxel).

## Rendering the showcase to MP4

```
./run_showcase_mp4.sh                 # render every showcase scene -> mp4/<name>.mp4
FPS=12 ./run_showcase_mp4.sh          # set MP4 frame rate
SCENES="scenes/showcase/showcase_cell8.json" ./run_showcase_mp4.sh   # one scene
DEBUG=1 ./run_showcase_mp4.sh         # quick low-res preview of all scenes
KEEP_PNGS=1 ./run_showcase_mp4.sh     # keep the intermediate PNG frames
```

These are heavy scenes by design (4D caustics need many photons); run a subset
or `DEBUG=1` first to preview. Requires `ffmpeg` (libx265).
