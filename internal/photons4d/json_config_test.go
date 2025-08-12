package photons4d

import (
	"math"
	"testing"
)

func TestRot4DegRadians(t *testing.T) {
	r := Rot4Deg{XY: 90, XZ: 180, XW: 0, YZ: 30, YW: -45, ZW: 10}.Radians()
	if math.Abs(float64(r.XY-math.Pi/2)) > 1e-12 || math.Abs(float64(r.XZ-math.Pi)) > 1e-12 {
		t.Fatal("degree->radian conversion wrong")
	}
}

func TestCell8CfgBuildAndValidation(t *testing.T) {
	// Valid
	h, err := (Cell8Cfg{
		Center: Point4{0, 0, 0, 0.5},
		Scale:  Vector4{1, 1, 1, 1},
		RotDeg: Rot4Deg{},
		Color:  RGB{1, 1, 1}, Reflect: RGB{0, 0, 0}, Refract: RGB{1, 1, 1}, IOR: RGB{1.2, 1.2, 1.2},
	}).Build()
	if err != nil || h == nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Invalid size
	h2, err := (Cell8Cfg{
		Center: Point4{}, Scale: Vector4{0, 1, 1, 1},
		RotDeg: Rot4Deg{}, Color: RGB{1, 1, 1}, Reflect: RGB{0, 0, 0}, Refract: RGB{1, 1, 1}, IOR: RGB{1.1, 1.1, 1.1},
	}).Build()
	if h2.Half != (Vector4{.5, .5, .5, .5}) {
		t.Fatalf("Scale defaults failed: %+v", h2.Half)
	}
}

func TestHyperSphereCfgScaleDefaults(t *testing.T) {
	hs, err := (HyperSphereCfg{
		Center: Point4{0, 0, 0, 1},
		Scale:  Vector4{.5, .5, .5, .5},
		RotDeg: Rot4Deg{},
		Color:  RGB{1, 1, 1}, Reflect: RGB{0, 0, 0}, Refract: RGB{1, 1, 1}, IOR: RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if hs.Radii != (Vector4{0.5, 0.5, 0.5, 0.5}) {
		t.Fatalf("Scale defaults failed: %+v", hs.Radii)
	}
}
