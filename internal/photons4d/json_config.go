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
}

type LightCfg struct {
	Origin    Point4  `json:"origin"`
	Direction Vector4 `json:"direction"`
	VoidLight bool    `json:"voidLight,omitempty"`
	Color     RGB     `json:"color"`
	AngleDeg  Real    `json:"angleDeg"`
}

type Config struct {
	SceneResX    int                 `json:"sceneResX"`
	SceneResY    int                 `json:"sceneResY"`
	SceneResZ    int                 `json:"sceneResZ"`
	ProbeRays    int                 `json:"probeRays"`
	Spp          int                 `json:"spp"`
	GIFOut       string              `json:"gifOut"`
	GIFDelay     int                 `json:"gifDelay,omitempty"`
	Gamma        Real                `json:"gamma,omitempty"`
	Scene        SceneCfg            `json:"scene"`
	Lights       []LightCfg          `json:"lights"`
	Hyperspheres []HyperSphereCfg    `json:"hyperspheres,omitempty"`
	Cells5       []Cell5Cfg          `json:"cells5,omitempty"`
	Cells8       []Cell8Cfg          `json:"cells8,omitempty"`
	Cells16      []SixteenCellCfg    `json:"cells16,omitempty"`
	Cells24      []TwentyFourCellCfg `json:"cells24,omitempty"`
	Cells120     []Cell120Cfg        `json:"cells120,omitempty"`
	Cells600     []Cell600Cfg        `json:"cells600,omitempty"`
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

type Cell8Cfg struct {
	Center Point4  `json:"center"`
	Size   Vector4 `json:"size"`
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type HyperSphereCfg struct {
	Center Point4  `json:"center"`
	Radius Real    `json:"radius"`          // base radius
	Scale  Vector4 `json:"scale,omitempty"` // optional per-axis scale; defaults 1
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type Cell5Cfg struct {
	Center Point4  `json:"center"`
	Side   Real    `json:"side"`            // edge length before per-axis Scale
	Scale  Vector4 `json:"scale,omitempty"` // defaults 1 on each axis
	RotDeg Rot4Deg `json:"rotDeg"`

	Color   RGB `json:"color"`
	Reflect RGB `json:"reflect"`
	Refract RGB `json:"refract"`
	IOR     RGB `json:"ior"`
}

type SixteenCellCfg struct {
	Center  Point4  `json:"center"`
	Side    Real    `json:"side"`
	Scale   Vector4 `json:"scale,omitempty"`
	RotDeg  Rot4Deg `json:"rotDeg"`
	Color   RGB     `json:"color"`
	Reflect RGB     `json:"reflect"`
	Refract RGB     `json:"refract"`
	IOR     RGB     `json:"ior"`
}

type TwentyFourCellCfg struct {
	Center  Point4  `json:"center"`
	Side    Real    `json:"side"`
	Scale   Vector4 `json:"scale,omitempty"`
	RotDeg  Rot4Deg `json:"rotDeg"`
	Color   RGB     `json:"color"`
	Reflect RGB     `json:"reflect"`
	Refract RGB     `json:"refract"`
	IOR     RGB     `json:"ior"`
}

type Cell120Cfg struct {
	Center  Point4  `json:"center"`
	Radius  Real    `json:"radius"`
	Scale   Vector4 `json:"scale,omitempty"`
	RotDeg  Rot4Deg `json:"rotDeg"`
	Color   RGB     `json:"color"`
	Reflect RGB     `json:"reflect"`
	Refract RGB     `json:"refract"`
	IOR     RGB     `json:"ior"`
}

type Cell600Cfg struct {
	Center  Point4  `json:"center"`
	Radius  Real    `json:"radius"`
	Scale   Vector4 `json:"scale,omitempty"`
	RotDeg  Rot4Deg `json:"rotDeg"`
	Color   RGB     `json:"color"`
	Reflect RGB     `json:"reflect"`
	Refract RGB     `json:"refract"`
	IOR     RGB     `json:"ior"`
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

	return NewCell8(
		hc.Center,
		hc.Size,
		rad,
		hc.Color,
		hc.Reflect,
		hc.Refract,
		hc.IOR,
	)
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
	radii := Vector4{hs.Radius * sc.X, hs.Radius * sc.Y, hs.Radius * sc.Z, hs.Radius * sc.W}
	return NewHyperSphere(hs.Center, radii, rad, hs.Color, hs.Reflect, hs.Refract, hs.IOR)
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
	return NewCell5(s.Center, s.Side, sc, rad, s.Color, s.Reflect, s.Refract, s.IOR)
}

func (s SixteenCellCfg) Build() (*SixteenCell, error) {
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
	return NewSixteenCell(s.Center, s.Side, sc, rad, s.Color, s.Reflect, s.Refract, s.IOR)
}

func (s TwentyFourCellCfg) Build() (*TwentyFourCell, error) {
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
	return NewTwentyFourCell(s.Center, s.Side, sc, rad, s.Color, s.Reflect, s.Refract, s.IOR)
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
	return NewCell120(c.Center, c.Radius, sc, c.RotDeg.Radians(), c.Color, c.Reflect, c.Refract, c.IOR)
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
	return NewCell600(c.Center, c.Radius, sc, c.RotDeg.Radians(), c.Color, c.Reflect, c.Refract, c.IOR)
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
	return &cfg, nil
}
