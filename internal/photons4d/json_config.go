package photons4d

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

type SceneCfg struct {
	Center     Point4 `json:"center"`
	Width      Real   `json:"width"`
	Height     Real   `json:"height"`
	Depth      Real   `json:"depth"`
	MaxBounces int    `json:"maxBounces,omitempty"`
	// When true, the "scene" is an infinite hypersphere (environment).
	// Rays that do not hit any object "escape" and deposit into angle bins
	// derived from their canonical 4D direction (α,β,γ).
	EnvHypersphere bool `json:"escape,omitempty"`
}

type LightCfg struct {
	Origin    Point4  `json:"origin"`
	Direction Vector4 `json:"direction"`
	Intensity Real    `json:"intensity"`
	Color     RGB     `json:"color"`
	AngleDeg  Real    `json:"angleDeg"`
}

type Config struct {
	SceneResX     int               `json:"sceneResX"`
	SceneResY     int               `json:"sceneResY"`
	SceneResZ     int               `json:"sceneResZ"`
	ProbeRays     int               `json:"probeRays"`
	Spp           int               `json:"spp"`
	GIFOut        string            `json:"gifOut"`
	GIFDelay      int               `json:"gifDelay,omitempty"`
	Gamma         Real              `json:"gamma,omitempty"`
	Scene         SceneCfg          `json:"scene"`
	Lights        []LightCfg        `json:"lights"`
	Hyperspheres  []HyperSphereCfg  `json:"hyperspheres,omitempty"`
	Cells5        []Cell5Cfg        `json:"cells5,omitempty"`
	Cells8        []Cell8Cfg        `json:"cells8,omitempty"`
	Cells16       []Cell16Cfg       `json:"cells16,omitempty"`
	Cells24       []Cell24Cfg       `json:"cells24,omitempty"`
	Cells120      []Cell120Cfg      `json:"cells120,omitempty"`
	Cells600      []Cell600Cfg      `json:"cells600,omitempty"`
	StarCells     []StarCellCfg     `json:"starCells,omitempty"`
	Spherinders   []SpherinderCfg   `json:"spherinders,omitempty"`
	HyperCones    []HyperConeCfg    `json:"hypercones,omitempty"`
	HyperCapsules []HyperCapsuleCfg `json:"hypercapsules,omitempty"`
	Spheritori    []SpheritorusCfg  `json:"spheritori,omitempty"`
	Duotori       []DuotorusCfg     `json:"duotori,omitempty"`
	Duocylinders  []DuocylinderCfg  `json:"duocylinders,omitempty"`
	Torispheres   []TorisphereCfg   `json:"torispheres,omitempty"`
	Superquadrics []SuperquadricCfg `json:"superquadrics,omitempty"`
	HyperFrustums []HyperFrustumCfg `json:"hyperfrustums,omitempty"`
	Cubinders     []CubinderCfg     `json:"cubinders,omitempty"`
	Ditori        []DitorusCfg      `json:"ditori,omitempty"`
	Duoprisms     []DuoprismCfg     `json:"duoprisms,omitempty"`

	Polychora            []PolychoronCfg      `json:"polychora,omitempty"`
	StarPolychora        []StarPolychoronCfg  `json:"starPolychora,omitempty"`
	PolyhedralPrisms     []PolyhedralCfg      `json:"polyhedralPrisms,omitempty"`
	PolyhedralPyramids   []PolyhedralCfg      `json:"polyhedralPyramids,omitempty"`
	PolyhedralBipyramids []PolyhedralCfg      `json:"polyhedralBipyramids,omitempty"`
	Duopyramids          []DuopyramidCfg      `json:"duopyramids,omitempty"`
	Bracketopes          []BracketopeCfg      `json:"bracketopes,omitempty"`
	Coninders            []ConeSolidCfg       `json:"coninders,omitempty"`
	Cylindrones          []ConeSolidCfg       `json:"cylindrones,omitempty"`
	Dicones              []ConeSolidCfg       `json:"dicones,omitempty"`
	PrismicCylinders     []PrismicCylinderCfg `json:"prismicCylinders,omitempty"`
	Torinders            []TorinderCfg        `json:"torinders,omitempty"`

	KleinBottles      []KleinBottleCfg      `json:"kleinBottles,omitempty"`
	ProjectivePlanes  []ProjectivePlaneCfg  `json:"projectivePlanes,omitempty"`
	MobiusBands       []MobiusBandCfg       `json:"mobiusBands,omitempty"`
	Paraboloids       []ParaboloidCfg       `json:"paraboloids,omitempty"`
	Hyperboloids      []HyperboloidCfg      `json:"hyperboloids,omitempty"`
	SphericalSegments []SphericalSegmentCfg `json:"sphericalSegments,omitempty"`
	Lenses            []LensCfg             `json:"lenses,omitempty"`
	ReuleauxSolids    []ReuleauxSolidCfg    `json:"reuleauxSolids,omitempty"`
}

// Rotation in degrees for JSON (friendlier than radians).
type Rot4Deg struct {
	XY Real `json:"xy"`
	XZ Real `json:"xz"`
	XW Real `json:"xw"`
	YZ Real `json:"yz"`
	YW Real `json:"yw"`
	ZW Real `json:"zw"`
}

type HyperSphereCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"` // optional per-axis scale; defaults 1
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type SpherinderCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell5Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"` // defaults 1 on each axis
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell8Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell16Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell24Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell120Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell600Cfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type StarCellCfg struct {
	Kind        string  `json:"kind"`
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	CoreRadius  Real    `json:"coreRadius,omitempty"`
	SpikeLength Real    `json:"spikeLength,omitempty"`
	Sharpness   Real    `json:"sharpness,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type HyperConeCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type HyperCapsuleCfg struct {
	Center     Point4  `json:"center"`
	Scale      Vector4 `json:"scale,omitempty"`
	HalfLength Real    `json:"halfLength,omitempty"`
	RotDeg     Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type SpheritorusCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	MajorRadius Real    `json:"majorRadius,omitempty"`
	MinorRadius Real    `json:"minorRadius,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type DuotorusCfg struct {
	Center        Point4  `json:"center"`
	Scale         Vector4 `json:"scale,omitempty"`
	MajorRadiusXY Real    `json:"majorRadiusXY,omitempty"`
	MajorRadiusZW Real    `json:"majorRadiusZW,omitempty"`
	MinorRadius   Real    `json:"minorRadius,omitempty"`
	RotDeg        Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type DuocylinderCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type TorisphereCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	MajorRadius Real    `json:"majorRadius,omitempty"`
	MinorRadius Real    `json:"minorRadius,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type SuperquadricCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	Power  Real    `json:"power,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type HyperFrustumCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	LowerRadius Real    `json:"lowerRadius,omitempty"`
	UpperRadius Real    `json:"upperRadius,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type CubinderCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type DitorusCfg struct {
	Center       Point4  `json:"center"`
	Scale        Vector4 `json:"scale,omitempty"`
	MajorRadius  Real    `json:"majorRadius,omitempty"`
	MiddleRadius Real    `json:"middleRadius,omitempty"`
	MinorRadius  Real    `json:"minorRadius,omitempty"`
	RotDeg       Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type DuoprismCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	P      int     `json:"p,omitempty"`
	Q      int     `json:"q,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// PolychoronCfg: a convex uniform polychoron by name ("truncated 120-cell",
// "snub 24-cell", "grand antiprism", ...) or by Coxeter diagram ("x5o3o3x"),
// or its dual.
type PolychoronCfg struct {
	Kind   string  `json:"kind"`
	Dual   bool    `json:"dual,omitempty"`
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// StarPolychoronCfg: a regular star polychoron by name ("grand 600-cell")
// or by Schläfli symbol ("{3,3,5/2}").
type StarPolychoronCfg struct {
	Kind   string  `json:"kind"`
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// PolyhedralCfg: the prism, pyramid or bipyramid over a uniform polyhedron
// ("truncated icosahedron", "snub cube", "7-gonal antiprism", "x5x3o", ...).
type PolyhedralCfg struct {
	Base   string  `json:"base"`
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type DuopyramidCfg DuoprismCfg

// BracketopeCfg: a shape in the bracket notation, such as "<(II)(II)>".
type BracketopeCfg struct {
	Shape  string  `json:"shape"`
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// ConeSolidCfg serves the coninders, cylindrones and dicones.
type ConeSolidCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type PrismicCylinderCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	Sides  int     `json:"sides,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type TorinderCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	MajorRadius Real    `json:"majorRadius,omitempty"`
	MinorRadius Real    `json:"minorRadius,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// KleinBottleCfg: the wall of the given thickness around a Klein bottle. An
// even twist (the number of half-turns, 1 if left out) gives a twisted torus.
type KleinBottleCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	MajorRadius Real    `json:"majorRadius,omitempty"`
	MinorRadius Real    `json:"minorRadius,omitempty"`
	Thickness   Real    `json:"thickness,omitempty"`
	Twist       int     `json:"twist,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type ProjectivePlaneCfg struct {
	Center    Point4  `json:"center"`
	Scale     Vector4 `json:"scale,omitempty"`
	Radius    Real    `json:"radius,omitempty"`
	Thickness Real    `json:"thickness,omitempty"`
	RotDeg    Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type MobiusBandCfg struct {
	Center      Point4  `json:"center"`
	Scale       Vector4 `json:"scale,omitempty"`
	MajorRadius Real    `json:"majorRadius,omitempty"`
	HalfWidth   Real    `json:"halfWidth,omitempty"`
	Thickness   Real    `json:"thickness,omitempty"`
	Twist       int     `json:"twist,omitempty"`
	RotDeg      Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type ParaboloidCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// HyperboloidCfg: waist is the radius in the middle (0.5 if left out, 0 for
// two cones tip to tip); the radius at both ends is 1 before the scale.
type HyperboloidCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	Waist  *Real   `json:"waist,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// SphericalSegmentCfg: the part of the ball between w = lower and w = upper
// (0 and 1 if left out: a half of the ball).
type SphericalSegmentCfg struct {
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	Lower  *Real   `json:"lower,omitempty"`
	Upper  *Real   `json:"upper,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type LensCfg struct {
	Center    Point4  `json:"center"`
	Scale     Vector4 `json:"scale,omitempty"`
	Thickness Real    `json:"thickness,omitempty"`
	RotDeg    Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

// ReuleauxSolidCfg: kind is "rotor" (if left out) or "simplex".
type ReuleauxSolidCfg struct {
	Kind   string  `json:"kind,omitempty"`
	Center Point4  `json:"center"`
	Scale  Vector4 `json:"scale,omitempty"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Diffuse RGB `json:"diffuse"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

func (r Rot4Deg) Radians() Rot4 {
	const k = math.Pi / 180
	return Rot4{
		XY: r.XY * k, XZ: r.XZ * k, XW: r.XW * k,
		YZ: r.YZ * k, YW: r.YW * k, ZW: r.ZW * k,
	}
}

// Build validates and constructs the runtime object (no defaults).
func (hc Cell8Cfg) Build() (*Cell8, error) {
	rad := hc.RotDeg.Radians()
	if math.Abs(rad.XY)+math.Abs(rad.XZ)+math.Abs(rad.XW)+
		math.Abs(rad.YZ)+math.Abs(rad.YW)+math.Abs(rad.ZW) < 1e-12 {
		// DebugLog("Cell8 rotation is ~zero; check JSON 'rotDeg' keys (xy,xz,xw,yz,yw,zw).")
	}
	sc := hc.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell8(hc.Center, sc, rad, hc.Color, hc.Diffuse, hc.Reflect, hc.Refract, hc.IOR)
}

func (hs HyperSphereCfg) Build() (*HyperSphere, error) {
	rad := hs.RotDeg.Radians()
	sc := hs.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	if sc.X <= 0 || sc.Y <= 0 || sc.Z <= 0 || sc.W <= 0 {
		return nil, fmt.Errorf("scale must be > 0 on all axes, got %+v", sc)
	}
	radii := Vector4{sc.X, sc.Y, sc.Z, sc.W}
	return NewHyperSphere(hs.Center, radii, rad, hs.Color, hs.Diffuse, hs.Reflect, hs.Refract, hs.IOR)
}

func (s Cell5Cfg) Build() (*Cell5, error) {
	rad := s.RotDeg.Radians()
	sc := s.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell5(s.Center, sc, rad, s.Color, s.Diffuse, s.Reflect, s.Refract, s.IOR)
}

func (s Cell16Cfg) Build() (*Cell16, error) {
	rad := s.RotDeg.Radians()
	sc := s.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell16(s.Center, sc, rad, s.Color, s.Diffuse, s.Reflect, s.Refract, s.IOR)
}

func (s Cell24Cfg) Build() (*Cell24, error) {
	rad := s.RotDeg.Radians()
	sc := s.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell24(s.Center, sc, rad, s.Color, s.Diffuse, s.Reflect, s.Refract, s.IOR)
}

func (c Cell120Cfg) Build() (*Cell120, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell120(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}
func (c Cell600Cfg) Build() (*Cell600, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCell600(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c StarCellCfg) Build() (*StarCell, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	core := c.CoreRadius
	if core == 0 {
		core = defaultStarCellCoreRadius
	}
	spike := c.SpikeLength
	if spike == 0 {
		spike = defaultStarCellSpikeLength
	}
	sharpness := c.Sharpness
	if sharpness == 0 {
		sharpness = defaultStarCellSharpness
	}
	return NewStarCell(c.Kind, c.Center, sc, core, spike, sharpness, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c SpherinderCfg) Build() (*Spherinder, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewSpherinder(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c HyperConeCfg) Build() (*HyperCone, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewHyperCone(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c HyperCapsuleCfg) Build() (*HyperCapsule, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	hl := c.HalfLength
	if hl == 0 {
		hl = 1
	}
	return NewHyperCapsule(c.Center, sc, hl, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c SpheritorusCfg) Build() (*Spheritorus, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	R := c.MajorRadius
	if R == 0 {
		R = 1
	}
	r := c.MinorRadius
	if r == 0 {
		r = 0.25
	}
	return NewSpheritorus(c.Center, sc, R, r, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c DuotorusCfg) Build() (*Duotorus, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	rxy := c.MajorRadiusXY
	if rxy == 0 {
		rxy = 1
	}
	rzw := c.MajorRadiusZW
	if rzw == 0 {
		rzw = 1
	}
	r := c.MinorRadius
	if r == 0 {
		r = 0.25
	}
	return NewDuotorus(c.Center, sc, rxy, rzw, r, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c DuocylinderCfg) Build() (*Duocylinder, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewDuocylinder(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c TorisphereCfg) Build() (*Torisphere, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	R := c.MajorRadius
	if R == 0 {
		R = 1
	}
	r := c.MinorRadius
	if r == 0 {
		r = 0.25
	}
	return NewTorisphere(c.Center, sc, R, r, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c SuperquadricCfg) Build() (*Superquadric, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	p := c.Power
	if p == 0 {
		p = 2
	}
	return NewSuperquadric(c.Center, sc, p, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c HyperFrustumCfg) Build() (*HyperFrustum, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	l := c.LowerRadius
	if l == 0 {
		l = 0.5
	}
	u := c.UpperRadius
	if u == 0 {
		u = 1.0
	}
	return NewHyperFrustum(c.Center, sc, l, u, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c CubinderCfg) Build() (*Cubinder, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return NewCubinder(c.Center, sc, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c DitorusCfg) Build() (*Ditorus, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	R := c.MajorRadius
	if R == 0 {
		R = 1
	}
	m := c.MiddleRadius
	if m == 0 {
		m = 0.4
	}
	r := c.MinorRadius
	if r == 0 {
		r = 0.15
	}
	return NewDitorus(c.Center, sc, R, m, r, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c DuoprismCfg) Build() (*Duoprism, error) {
	sc := c.Scale
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	p := c.P
	if p == 0 {
		p = 3
	}
	q := c.Q
	if q == 0 {
		q = 3
	}
	return NewDuoprism(c.Center, sc, p, q, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

// scaleOrOne replaces the missing (zero) parts of a scale by 1.
func scaleOrOne(sc Vector4) Vector4 {
	if sc.X == 0 {
		sc.X = 1
	}
	if sc.Y == 0 {
		sc.Y = 1
	}
	if sc.Z == 0 {
		sc.Z = 1
	}
	if sc.W == 0 {
		sc.W = 1
	}
	return sc
}

func (c PolychoronCfg) Build() (*Polytope, error) {
	return NewPolychoron(c.Kind, c.Dual, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c StarPolychoronCfg) Build() (*StarPolychoron, error) {
	return NewStarPolychoron(c.Kind, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

// Build makes the prism, the pyramid or the bipyramid (shape) over the base.
func (c PolyhedralCfg) Build(shape string) (*Polytope, error) {
	return NewPolyhedralSolid(shape, c.Base, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c DuopyramidCfg) Build() (*Polytope, error) {
	p := c.P
	if p == 0 {
		p = 3
	}
	q := c.Q
	if q == 0 {
		q = 3
	}
	return NewDuopyramid(c.Center, scaleOrOne(c.Scale), p, q, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c BracketopeCfg) Build() (*ConvexSolid, error) {
	return NewBracketope(c.Shape, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

// Build makes the coninder, the cylindrone or the dicone (kind).
func (c ConeSolidCfg) Build(kind string) (*ConvexSolid, error) {
	build := NewConinder
	switch kind {
	case "cylindrone":
		build = NewCylindrone
	case "dicone":
		build = NewDicone
	}
	return build(c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c PrismicCylinderCfg) Build() (*ConvexSolid, error) {
	n := c.Sides
	if n == 0 {
		n = 3
	}
	return NewPrismicCylinder(n, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c TorinderCfg) Build() (*Torinder, error) {
	R := c.MajorRadius
	if R == 0 {
		R = 1
	}
	r := c.MinorRadius
	if r == 0 {
		r = 0.25
	}
	return NewTorinder(c.Center, scaleOrOne(c.Scale), R, r, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

// orDefault replaces a missing (zero) value by its default.
func orDefault(v, def Real) Real {
	if v == 0 {
		return def
	}
	return v
}

func (c KleinBottleCfg) Build() (*KleinSurface, error) {
	twist := c.Twist
	if twist == 0 {
		twist = 1
	}
	return NewKleinBottle(c.Center, scaleOrOne(c.Scale), orDefault(c.MajorRadius, 1), orDefault(c.MinorRadius, 0.4), orDefault(c.Thickness, 0.1), twist,
		c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c ProjectivePlaneCfg) Build() (*KleinSurface, error) {
	return NewProjectivePlane(c.Center, scaleOrOne(c.Scale), orDefault(c.Radius, 0.5), orDefault(c.Thickness, 0.1),
		c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c MobiusBandCfg) Build() (*KleinSurface, error) {
	twist := c.Twist
	if twist == 0 {
		twist = 1
	}
	return NewMobiusBand(c.Center, scaleOrOne(c.Scale), orDefault(c.MajorRadius, 1), orDefault(c.HalfWidth, 0.4), orDefault(c.Thickness, 0.1), twist,
		c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c ParaboloidCfg) Build() (*RevolvedQuadric, error) {
	return NewParaboloid(c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c HyperboloidCfg) Build() (*RevolvedQuadric, error) {
	waist := Real(0.5)
	if c.Waist != nil {
		waist = *c.Waist
	}
	return NewHyperboloid(c.Center, scaleOrOne(c.Scale), waist, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c SphericalSegmentCfg) Build() (*RevolvedQuadric, error) {
	lower, upper := Real(0), Real(1)
	if c.Lower != nil {
		lower = *c.Lower
	}
	if c.Upper != nil {
		upper = *c.Upper
	}
	return NewSphericalSegment(c.Center, scaleOrOne(c.Scale), lower, upper, c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c LensCfg) Build() (*ConvexSolid, error) {
	return NewLens(c.Center, scaleOrOne(c.Scale), orDefault(c.Thickness, 0.6), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func (c ReuleauxSolidCfg) Build() (*ConvexSolid, error) {
	return NewReuleauxSolid(c.Kind, c.Center, scaleOrOne(c.Scale), c.RotDeg.Radians(), c.Color, c.Diffuse, c.Reflect, c.Refract, c.IOR)
}

func loadConfig(path string) (*Config, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var cfg Config
	if err := json.Unmarshal(data, &cfg); err != nil {
		return nil, err
	}
	// Defaults / validation
	if cfg.SceneResX <= 0 {
		cfg.SceneResX = SceneResX
	}
	if cfg.SceneResY <= 0 {
		cfg.SceneResY = SceneResY
	}
	if cfg.SceneResZ <= 0 {
		cfg.SceneResZ = SceneResZ
	}
	if cfg.ProbeRays <= 0 {
		cfg.ProbeRays = ProbeRays
	}
	if cfg.Spp <= 0 {
		cfg.Spp = Spp
	}
	if cfg.GIFOut == "" {
		cfg.GIFOut = GIFOut
	}
	if cfg.GIFDelay <= 0 {
		cfg.GIFDelay = GIFDelay
	}
	if cfg.Gamma <= 0 {
		cfg.Gamma = Gamma
	}
	if len(cfg.Lights) == 0 {
		return nil, fmt.Errorf("config has no lights")
	}
	if cfg.Scene.MaxBounces <= 0 {
		cfg.Scene.MaxBounces = MaxBounces
	}
	if ForceEscape {
		cfg.Scene.EnvHypersphere = true
	}
	if EscapeSPPAdjust && cfg.Scene.EnvHypersphere {
		cfg.Spp *= EscapeSppFactor
	}
	DebugLog("Loaded config from %s: size=(%d, %d, %d), probe=%d, SPP=%d, gamma=%f", path, cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ, cfg.ProbeRays, cfg.Spp, cfg.Gamma)
	return &cfg, nil
}
