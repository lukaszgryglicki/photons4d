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
	SceneResX    int              `json:"sceneResX"`
	SceneResY    int              `json:"sceneResY"`
	SceneResZ    int              `json:"sceneResZ"`
	ProbeRays    int              `json:"probeRays"`
	Spp          int              `json:"spp"`
	GIFOut       string           `json:"gifOut"`
	GIFDelay     int              `json:"gifDelay,omitempty"`
	Gamma        Real             `json:"gamma,omitempty"`
	Scene        SceneCfg         `json:"scene"`
	Lights       []LightCfg       `json:"lights"`
	Hyperspheres []HyperSphereCfg `json:"hyperspheres,omitempty"`
	Cells5       []Cell5Cfg       `json:"cells5,omitempty"`
	Cells8       []Cell8Cfg       `json:"cells8,omitempty"`
	Cells16      []Cell16Cfg      `json:"cells16,omitempty"`
	Cells24      []Cell24Cfg      `json:"cells24,omitempty"`
	Cells120     []Cell120Cfg     `json:"cells120,omitempty"`
	Cells600     []Cell600Cfg     `json:"cells600,omitempty"`
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
	DebugLog("Loaded config from %s: size=(%d, %d, %d), probe=%d, SPP=%d, gamma=%f", path, cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ, cfg.ProbeRays, cfg.Spp, cfg.Gamma)
	return &cfg, nil
}
