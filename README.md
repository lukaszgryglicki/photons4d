# photons4d

A 4D forward photon-tracing renderer. Lights emit photons into a 4D scene; the
energy that lands on the 3D hyperplane slab at `W = scene.center.W` is recorded
as a voxel volume and exported as an animated GIF / 16-bit PNG sequence / MP4,
one frame per Z slice.

## Build

```
make build       # debug build (verbose logging)
make release     # optimized binary ./photons4d
make static      # optimized binary without cgo ./photons4d.static
make install     # release + static, photons4d.static -> $(INSTALL_DIR)/photons4d (default /data/scripts)
make test        # go test ./internal/photons4d
make clean
```

Builds on FreeBSD, Linux and macOS. `make static` gives a static executable on
Linux and FreeBSD; macOS has none, there the binary needs only the libraries of
the system. On x86-64 an optimized build uses the highest level (`GOAMD64`) that
the building machine runs; pass `GOAMD64=v3` (or `v2`, `v1`) for a binary that
has to run on older processors.

## Objects

A scene lists its objects under the keys below. Every object has `center`,
`scale` (per axis, 1 when left out), `rotDeg` (`xy`, `xz`, `xw`, `yz`, `yw`,
`zw`) and a material (`color`, `diffuse`, `reflect`, `refract`, `ior`); the
last column names what else it takes.

| key | solid | fields |
| --- | --- | --- |
| `hyperspheres` | ball, ellipsoid | |
| `cells5` `cells8` `cells16` `cells24` `cells120` `cells600` | the six convex regular polychora | |
| `polychora` | convex uniform polychoron, or its dual | `kind`, `dual` |
| `starPolychora` | regular star polychoron | `kind` |
| `polyhedralPrisms` `polyhedralPyramids` `polyhedralBipyramids` | prism, pyramid, bipyramid over a uniform polyhedron | `base` |
| `duoprisms` `duopyramids` | product and hull of two regular polygons | `p`, `q` |
| `starCells` | ball with the spikes of a regular polychoron | `kind`, `coreRadius`, `spikeLength`, `sharpness` |
| `spherinders` `cubinders` `duocylinders` | ball × segment, disk × square, disk × disk | |
| `prismicCylinders` | disk × regular polygon | `sides` |
| `hypercones` `hyperfrustums` | cone and frustum over a ball | `lowerRadius`, `upperRadius` (frustum) |
| `coninders` `cylindrones` `dicones` | cone × segment, cone over a cylinder, cone over a cone | |
| `hypercapsules` | ball swept along a segment | `halfLength` |
| `spheritori` `torispheres` | tube around a circle, around a sphere | `majorRadius`, `minorRadius` |
| `duotori` | tube around a flat torus (tiger) | `majorRadiusXY`, `majorRadiusZW`, `minorRadius` |
| `ditori` | tube around a torus | `majorRadius`, `middleRadius`, `minorRadius` |
| `torinders` | torus × segment | `majorRadius`, `minorRadius` |
| `kleinBottles` | wall around a Klein bottle | `majorRadius`, `minorRadius`, `thickness`, `twist` |
| `projectivePlanes` | wall around the real projective plane | `radius`, `thickness` |
| `mobiusBands` | wall around a Möbius band | `majorRadius`, `halfWidth`, `thickness`, `twist` |
| `paraboloids` `hyperboloids` | paraboloid, hyperboloid of one sheet over a ball | `waist` (hyperboloid) |
| `sphericalSegments` | ball between two parallel hyperplanes | `lower`, `upper` |
| `lenses` | part common to two balls | `thickness` |
| `reuleauxSolids` | body of constant width, Reuleaux simplex | `kind` |
| `superquadrics` | `\|x\|^p + \|y\|^p + \|z\|^p + \|w\|^p <= 1` | `power` |
| `bracketopes` | body of a bracket notation | `shape` |

`polychora`: `kind` is an operation on a regular polychoron — `rectified`,
`truncated`, `bitruncated`, `cantellated`, `runcinated`, `cantitruncated`,
`runcitruncated` or `omnitruncated` before `5-cell`, `tesseract`, `16-cell`,
`24-cell`, `120-cell` or `600-cell`, as in `"truncated 120-cell"` — or
`snub 24-cell` or `grand antiprism`: the 47 convex uniform polychora that are
not prisms. It can also be a Coxeter diagram: `x5o3o3x` is the runcinated
120-cell (`x` is a ringed node, `o` a plain one, the numbers are the branches;
`2` stands for no branch, so `x3o2x7o` is a duoprism). `"dual": true` gives the
dual. Before `scale` the vertices lie at the distance 1 from the centre (the
farthest ones in the case of a dual).

`starPolychora`: `kind` is one of the ten regular star polychora, by name
(`icosahedral 120-cell`, `small stellated 120-cell`, `great 120-cell`,
`grand 120-cell`, `great stellated 120-cell`, `grand stellated 120-cell`,
`great grand 120-cell`, `great icosahedral 120-cell`, `grand 600-cell`,
`great grand stellated 120-cell`) or by Schläfli symbol (`{5/2,5,3}`). The solid
is all that the polychoron winds around; its vertices lie at the distance 1.

`base` is one of the 5 Platonic and 13 Archimedean solids by name
(`truncated icosahedron`, `snub cube`), a prism or an antiprism
(`7-gonal prism`, `pentagonal antiprism`) or a diagram (`x5x3o`); its
circumradius is 1 and it lies in XYZ. The prism is as long as an edge of the
base, which makes it uniform; the pyramid has its apex at `w = -1` and its base
at `w = +1`, the bipyramid its base at `w = 0` and its apexes at `w = ±1`.

The cones follow `hypercones`: the apex at -1 and the base at +1 of the axis,
which is Z for the coninder and W for the cylindrone; the dicone joins the disk
at `w = +1` to the segment along Z at `w = -1`. The polygons of `duoprisms`,
`duopyramids` and `prismicCylinders` have the circumradius 1.

`bracketopes`: `shape` joins the four axes (`I`) with round brackets (the root
of the sum of the squares), square brackets (the largest) and angle brackets
(the sum): `(IIII)` is the ball, `[IIII]` the tesseract, `<IIII>` the 16-cell,
`[(II)(II)]` the duocylinder, `<(II)(II)>` its dual and `<(III)I>` the
bipyramid over a ball. The axes may be named instead (`x`, `y`, `z`, `w`), and
then groups may share them: `[(xyz)(xyw)]` is the part common to two
spherinders, `[(xy)(xz)(xw)(yz)(yw)(zw)]` the part common to six cubinders.

`kleinBottles`, `projectivePlanes`, `mobiusBands`: surfaces without an inside.
A ray would miss a surface, so each is the wall of all points within
`thickness`/2 of it. A circle (of `minorRadius`, around a point at
`majorRadius` from the centre) or, for the band, a segment (`halfWidth` to both
sides of such a point) sweeps the surface: its plane goes once around XY and
makes `twist` half-turns in ZW on the way,
`(a cos u, a sin u, b cos(twist·u/2), b sin(twist·u/2))`. One half-turn (when
`twist` is left out) gives the Klein bottle and the Möbius band, which do not
cross themselves in 4D; an even `twist` gives a torus or a ring twisted in 4D.
The projective plane is the Klein bottle whose circles (of `radius`) all pass
through the centre.

The paraboloid has its tip at `w = -1` and the radius 1 at `w = +1`. The
hyperboloid has the radius 1 at `w = ±1` and `waist` in the middle (0.5 when
left out; 0 gives two cones tip to tip). `sphericalSegments` keep the part of
the unit ball between `w = lower` and `w = upper` (0 and 1 when left out: a
half of the ball). The rim of a lens is the sphere of the radius 1 in `w = 0`
and `thickness` is its size along W (0.6 when left out; 2 gives the ball).
`reuleauxSolids`: the `rotor` (when `kind` is left out) is the Reuleaux
triangle of the width 2 turned around W, a body of the same width in every
direction that rolls between two hyperplanes as a ball does; the `simplex` is
the part common to five balls around the corners of a 5-cell.

## Showcase scenes

`scenes/showcase/*.json` are generated, ready-to-render scenes — one per
supported 4D primitive (hypersphere, 5/8/16/24/120/600-cell, star cells,
spherinder, hypercone, hypercapsule, spheritorus, duotorus, duocylinder,
torisphere, superquadric, hyperfrustum, cubinder, ditorus, duoprism), scenes
with uniform and star polychora, polyhedral solids, a duopyramid, bracketopes,
cone solids, a prismic cylinder and a torinder, a Klein bottle, a projective
plane and a Möbius band (alone, together and with one to three half-turns),
quadrics, a lens and rolling bodies, plus five "mixed" scenes
(`glass_garden`, `prism_chamber`, `grand`, `menagerie`, and an
`escape_environment` scene
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
