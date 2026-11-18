package photons4d

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
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
		Center:  Point4{0, 0, 0, 0.5},
		Scale:   Vector4{1, 1, 1, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.2, 1.2, 1.2},
	}).Build()
	if err != nil || h == nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Invalid size
	h2, err := (Cell8Cfg{
		Center: Point4{}, Scale: Vector4{0, 1, 1, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.1, 1.1, 1.1},
	}).Build()
	if h2.Half != (Vector4{.5, .5, .5, .5}) {
		t.Fatalf("Scale defaults failed: %+v", h2.Half)
	}
}

func TestHyperSphereCfgScaleDefaults(t *testing.T) {
	hs, err := (HyperSphereCfg{
		Center:  Point4{0, 0, 0, 1},
		Scale:   Vector4{.5, .5, .5, .5},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if hs.Radii != (Vector4{0.5, 0.5, 0.5, 0.5}) {
		t.Fatalf("Scale defaults failed: %+v", hs.Radii)
	}
}

func TestHyperConeCfgScaleDefaults(t *testing.T) {
	hc, err := (HyperConeCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if hc.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", hc.Scale)
	}
}

func TestHyperCapsuleCfgDefaults(t *testing.T) {
	hc, err := (HyperCapsuleCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if hc.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", hc.Scale)
	}
	if hc.HalfLength != 1 {
		t.Fatalf("HalfLength default failed: %v", hc.HalfLength)
	}
}

func TestSpheritorusCfgDefaults(t *testing.T) {
	st, err := (SpheritorusCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if st.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", st.Scale)
	}
	if st.MajorRadius != 1 || st.MinorRadius != 0.25 {
		t.Fatalf("Radius defaults failed: major=%v minor=%v", st.MajorRadius, st.MinorRadius)
	}
}

func TestDuotorusCfgDefaults(t *testing.T) {
	dt, err := (DuotorusCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if dt.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", dt.Scale)
	}
	if dt.MajorRadiusXY != 1 || dt.MajorRadiusZW != 1 || dt.MinorRadius != 0.25 {
		t.Fatalf("Radius defaults failed: Rxy=%v Rzw=%v r=%v", dt.MajorRadiusXY, dt.MajorRadiusZW, dt.MinorRadius)
	}
}

func TestSpherinderCfgDefaults(t *testing.T) {
	sp, err := (SpherinderCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if sp.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", sp.Scale)
	}
}

func TestDuocylinderCfgDefaults(t *testing.T) {
	dc, err := (DuocylinderCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if dc.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", dc.Scale)
	}
}

func TestTorisphereCfgDefaults(t *testing.T) {
	ts, err := (TorisphereCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if ts.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", ts.Scale)
	}
	if ts.MajorRadius != 1 || ts.MinorRadius != 0.25 {
		t.Fatalf("Radius defaults failed: major=%v minor=%v", ts.MajorRadius, ts.MinorRadius)
	}
}

func TestSuperquadricCfgDefaults(t *testing.T) {
	sq, err := (SuperquadricCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if sq.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", sq.Scale)
	}
	if sq.Power != 2 {
		t.Fatalf("Power default failed: %v", sq.Power)
	}
}

func TestHyperFrustumCfgDefaults(t *testing.T) {
	hf, err := (HyperFrustumCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if hf.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", hf.Scale)
	}
	if hf.LowerRadius != 0.5 || hf.UpperRadius != 1.0 {
		t.Fatalf("Radius defaults failed: lower=%v upper=%v", hf.LowerRadius, hf.UpperRadius)
	}
}

func TestStarCellCfgDefaults(t *testing.T) {
	sc, err := (StarCellCfg{
		Kind:    "star-cell-24",
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if sc.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", sc.Scale)
	}
	if sc.CoreRadius != defaultStarCellCoreRadius {
		t.Fatalf("CoreRadius default failed: %v", sc.CoreRadius)
	}
	if sc.SpikeLength != defaultStarCellSpikeLength {
		t.Fatalf("SpikeLength default failed: %v", sc.SpikeLength)
	}
	if sc.Sharpness != defaultStarCellSharpness {
		t.Fatalf("Sharpness default failed: %v", sc.Sharpness)
	}
	if sc.Kind != "star-cell-24" {
		t.Fatalf("kind canonicalization failed: %q", sc.Kind)
	}
}

func TestCubinderCfgDefaults(t *testing.T) {
	c, err := (CubinderCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if c.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", c.Scale)
	}
}

func TestDitorusCfgDefaults(t *testing.T) {
	d, err := (DitorusCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if d.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", d.Scale)
	}
	if d.MajorRadius != 1 || d.MiddleRadius != 0.4 || d.MinorRadius != 0.15 {
		t.Fatalf("Radius defaults failed: major=%v middle=%v minor=%v", d.MajorRadius, d.MiddleRadius, d.MinorRadius)
	}
}

func TestDuoprismCfgDefaults(t *testing.T) {
	d, err := (DuoprismCfg{
		Center:  Point4{0, 0, 0, 1},
		RotDeg:  Rot4Deg{},
		Color:   RGB{1, 1, 1},
		Diffuse: RGB{0, 0, 0},
		Reflect: RGB{0, 0, 0},
		Refract: RGB{1, 1, 1},
		IOR:     RGB{1.4, 1.4, 1.4},
	}).Build()
	if err != nil {
		t.Fatal(err)
	}
	if d.Scale != (Vector4{1, 1, 1, 1}) {
		t.Fatalf("Scale defaults failed: %+v", d.Scale)
	}
	if d.P != 3 || d.Q != 3 {
		t.Fatalf("Polygon defaults failed: p=%d q=%d", d.P, d.Q)
	}
}

// The three new kinds are read from JSON under their keys and built into the scene.
func TestNewObjectKindsLoadFromJSON(t *testing.T) {
	const mat = `"color": {"R":1,"G":1,"B":1}, "diffuse": {"R":0,"G":0,"B":0}, "reflect": {"R":0.1,"G":0.1,"B":0.1},
		"refract": {"R":0.8,"G":0.8,"B":0.8}, "ior": {"R":1.4,"G":1.5,"B":1.6}`
	path := filepath.Join(t.TempDir(), "scene.json")
	data := `{
		"sceneResX": 4, "sceneResY": 4, "sceneResZ": 2, "probeRays": 100, "spp": 1, "gifOut": "x.gif",
		"scene": {"center": {"X":0,"Y":0,"Z":0,"W":0}, "width": 2, "height": 2, "depth": 2},
		"lights": [{"origin": {"X":0,"Y":0,"Z":0,"W":1}, "direction": {"X":0,"Y":0,"Z":0,"W":-1},
			"color": {"R":1,"G":1,"B":1}, "angleDeg": 30, "intensity": 1}],
		"cubinders": [{"center": {"X":0.4,"Y":0,"Z":0,"W":0.5}, "scale": {"X":0.2,"Y":0.2,"Z":0.1,"W":0.15},
			"rotDeg": {"xy":10,"xw":20}, ` + mat + `}],
		"ditori": [{"center": {"X":-0.4,"Y":0,"Z":0,"W":0.5}, "majorRadius": 0.3, "middleRadius": 0.1, "minorRadius": 0.04,
			"rotDeg": {"yw":30}, ` + mat + `}],
		"duoprisms": [{"center": {"X":0,"Y":0.4,"Z":0,"W":0.5}, "scale": {"X":0.2,"Y":0.2,"Z":0.2,"W":0.2}, "p": 5, "q": 6,
			"rotDeg": {"zw":15}, ` + mat + `}]
	}`
	if err := os.WriteFile(path, []byte(data), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg, err := loadConfig(path)
	if err != nil {
		t.Fatal(err)
	}
	scene, _, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(scene.Cubinders) != 1 || len(scene.Ditori) != 1 || len(scene.Duoprisms) != 1 || scene.NObjects() != 3 {
		t.Fatalf("objects: cubinders=%d ditori=%d duoprisms=%d total=%d", len(scene.Cubinders), len(scene.Ditori), len(scene.Duoprisms), scene.NObjects())
	}
	if scene.Cubinders[0].Scale != (Vector4{0.2, 0.2, 0.1, 0.15}) {
		t.Fatalf("cubinder scale: %+v", scene.Cubinders[0].Scale)
	}
	if d := scene.Ditori[0]; d.MajorRadius != 0.3 || d.MiddleRadius != 0.1 || d.MinorRadius != 0.04 {
		t.Fatalf("ditorus radii: %v %v %v", d.MajorRadius, d.MiddleRadius, d.MinorRadius)
	}
	if d := scene.Duoprisms[0]; d.P != 5 || d.Q != 6 {
		t.Fatalf("duoprism polygons: p=%d q=%d", d.P, d.Q)
	}

	bad := strings.Replace(data, `"p": 5`, `"p": 2`, 1)
	if err := os.WriteFile(path, []byte(bad), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg, err = loadConfig(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, _, err = buildScene(cfg); err == nil || !strings.Contains(err.Error(), "duoprisms[0]") {
		t.Fatalf("a 2-gon must be rejected with the place of the error, got %v", err)
	}
}

// The polytopes, star polychora, convex solids and torinders are read from
// JSON under their keys, with their defaults, and built into the scene.
func TestPolytopeAndSolidKindsLoadFromJSON(t *testing.T) {
	const mat = `"color": {"R":1,"G":1,"B":1}, "diffuse": {"R":0,"G":0,"B":0}, "reflect": {"R":0.1,"G":0.1,"B":0.1},
		"refract": {"R":0.8,"G":0.8,"B":0.8}, "ior": {"R":1.4,"G":1.5,"B":1.6}`
	const at = `"center": {"X":0.1,"Y":0.2,"Z":0.3,"W":0.4}, `
	path := filepath.Join(t.TempDir(), "scene.json")
	data := `{
		"sceneResX": 4, "sceneResY": 4, "sceneResZ": 2, "probeRays": 100, "spp": 1, "gifOut": "x.gif",
		"scene": {"center": {"X":0,"Y":0,"Z":0,"W":0}, "width": 2, "height": 2, "depth": 2},
		"lights": [{"origin": {"X":0,"Y":0,"Z":0,"W":1}, "direction": {"X":0,"Y":0,"Z":0,"W":-1},
			"color": {"R":1,"G":1,"B":1}, "angleDeg": 30, "intensity": 1}],
		"polychora": [
			{"kind": "truncated 120-cell", ` + at + `"scale": {"X":0.2,"Y":0.3,"Z":0.2,"W":0.2}, "rotDeg": {"xy":10}, ` + mat + `},
			{"kind": "x3o4o3x", "dual": true, ` + at + mat + `},
			{"kind": "grand antiprism", ` + at + mat + `}],
		"starPolychora": [
			{"kind": "Great Grand Stellated 120-cell", ` + at + mat + `},
			{"kind": "{5,5/2,5}", ` + at + `"scale": {"X":0.5}, ` + mat + `}],
		"polyhedralPrisms": [{"base": "snub cube", ` + at + mat + `}],
		"polyhedralPyramids": [{"base": "icosahedron", ` + at + mat + `}, {"base": "5-gonal antiprism", ` + at + mat + `}],
		"polyhedralBipyramids": [{"base": "x5x3x", ` + at + mat + `}],
		"duopyramids": [{` + at + mat + `}, {"p": 4, "q": 6, ` + at + mat + `}],
		"bracketopes": [{"shape": "<(II)(II)>", ` + at + mat + `}, {"shape": "duospindle", ` + at + mat + `}],
		"coninders": [{` + at + mat + `}],
		"cylindrones": [{` + at + mat + `}],
		"dicones": [{` + at + mat + `}],
		"prismicCylinders": [{` + at + mat + `}, {"sides": 7, ` + at + mat + `}],
		"torinders": [{` + at + mat + `}, {"majorRadius": 0.5, "minorRadius": 0.2, "scale": {"W":3}, ` + at + mat + `}]
	}`
	if err := os.WriteFile(path, []byte(data), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg, err := loadConfig(path)
	if err != nil {
		t.Fatal(err)
	}
	scene, _, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(scene.Polytopes) != 9 || len(scene.StarPolychora) != 2 || len(scene.ConvexSolids) != 7 || len(scene.Torinders) != 2 || scene.NObjects() != 20 {
		t.Fatalf("objects: polytopes=%d stars=%d convex=%d torinders=%d total=%d",
			len(scene.Polytopes), len(scene.StarPolychora), len(scene.ConvexSolids), len(scene.Torinders), scene.NObjects())
	}
	var names []string
	for _, p := range scene.Polytopes {
		names = append(names, fmt.Sprintf("%s/%d", p.Name, len(p.N)))
		if p.Center != (Point4{0.1, 0.2, 0.3, 0.4}) {
			t.Errorf("%s: center %v", p.Name, p.Center)
		}
	}
	// name and number of cells, in the order of the keys
	if got, want := strings.Join(names, ", "), "truncated 120-cell/720, dual x3o4o3x/144, grand antiprism/320, "+
		"snub cube prism/40, icosahedron pyramid/21, 5-gonal antiprism pyramid/13, x5x3x bipyramid/124, 3-3 duopyramid/9, 4-6 duopyramid/24"; got != want {
		t.Errorf("polytopes:\n got %s\nwant %s", got, want)
	}
	if p := scene.Polytopes[0]; p.Scale != (Vector4{0.2, 0.3, 0.2, 0.2}) || p.R == I4() {
		t.Errorf("truncated 120-cell: scale %v, rotation %v", p.Scale, p.R)
	}
	if s := scene.StarPolychora; s[0].Name != "{5/2,3,3}" || s[1].Name != "{5,5/2,5}" || s[0].Scale != (Vector4{1, 1, 1, 1}) || s[1].Scale != (Vector4{0.5, 1, 1, 1}) {
		t.Errorf("star polychora: %s %v, %s %v", s[0].Name, s[0].Scale, s[1].Name, s[1].Scale)
	}
	names = names[:0]
	for _, c := range scene.ConvexSolids {
		names = append(names, c.Name)
	}
	if got, want := strings.Join(names, ", "), "bracketope <(II)(II)>, bracketope <(II)(II)>, coninder, cylindrone, dicone, "+
		"3-gonal prismic cylinder, 7-gonal prismic cylinder"; got != want {
		t.Errorf("convex solids:\n got %s\nwant %s", got, want)
	}
	if a, b := scene.Torinders[0], scene.Torinders[1]; a.MajorRadius != 1 || a.MinorRadius != 0.25 || a.Scale != (Vector4{1, 1, 1, 1}) ||
		b.MajorRadius != 0.5 || b.MinorRadius != 0.2 || b.Scale != (Vector4{1, 1, 1, 3}) {
		t.Errorf("torinders: %v/%v %v, %v/%v %v", a.MajorRadius, a.MinorRadius, a.Scale, b.MajorRadius, b.MinorRadius, b.Scale)
	}

	// a wrong entry is reported with its key and place
	for _, c := range []struct{ old, new, where string }{
		{`"kind": "grand antiprism"`, `"kind": "grand prism"`, "polychora[2]"},
		{`"kind": "{5,5/2,5}"`, `"kind": "{5,5/2,4}"`, "starPolychora[1]"},
		{`"base": "snub cube"`, `"base": "snub sphere"`, "polyhedralPrisms[0]"},
		{`"base": "5-gonal antiprism"`, `"base": "2-gonal antiprism"`, "polyhedralPyramids[1]"},
		{`"base": "x5x3x"`, `"base": "x5x3x3x"`, "polyhedralBipyramids[0]"},
		{`"p": 4, "q": 6`, `"p": 4, "q": 2`, "duopyramids[1]"},
		{`"shape": "duospindle"`, `"shape": "(III)"`, "bracketopes[1]"},
		{`"sides": 7`, `"sides": 2`, "prismicCylinders[1]"},
		{`"majorRadius": 0.5, "minorRadius": 0.2`, `"majorRadius": 0.2, "minorRadius": 0.5`, "torinders[1]"},
		{`"coninders": [{`, `"coninders": [{"scale": {"X":-1}, `, "coninders[0]"},
		{`"cylindrones": [{`, `"cylindrones": [{"scale": {"X":-1}, `, "cylindrones[0]"},
		{`"dicones": [{`, `"dicones": [{"scale": {"X":-1}, `, "dicones[0]"},
	} {
		if !strings.Contains(data, c.old) {
			t.Fatalf("%q is not in the scene", c.old)
		}
		if err := os.WriteFile(path, []byte(strings.Replace(data, c.old, c.new, 1)), 0o644); err != nil {
			t.Fatal(err)
		}
		cfg, err := loadConfig(path)
		if err != nil {
			t.Fatal(err)
		}
		if _, _, err = buildScene(cfg); err == nil || !strings.Contains(err.Error(), c.where) {
			t.Errorf("%s must be rejected with its place, got %v", c.new, err)
		}
	}
}

// The walls around surfaces and the round solids are read from JSON under
// their keys, with their defaults, and built into the scene.
func TestKleinAndRoundKindsLoadFromJSON(t *testing.T) {
	const mat = `"color": {"R":1,"G":1,"B":1}, "diffuse": {"R":0,"G":0,"B":0}, "reflect": {"R":0.1,"G":0.1,"B":0.1},
		"refract": {"R":0.8,"G":0.8,"B":0.8}, "ior": {"R":1.4,"G":1.5,"B":1.6}`
	const at = `"center": {"X":0.1,"Y":0.2,"Z":0.3,"W":0.4}, `
	path := filepath.Join(t.TempDir(), "scene.json")
	data := `{
		"sceneResX": 4, "sceneResY": 4, "sceneResZ": 2, "probeRays": 100, "spp": 1, "gifOut": "x.gif",
		"scene": {"center": {"X":0,"Y":0,"Z":0,"W":0}, "width": 2, "height": 2, "depth": 2},
		"lights": [{"origin": {"X":0,"Y":0,"Z":0,"W":1}, "direction": {"X":0,"Y":0,"Z":0,"W":-1},
			"color": {"R":1,"G":1,"B":1}, "angleDeg": 30, "intensity": 1}],
		"kleinBottles": [{` + at + mat + `},
			{"majorRadius": 0.6, "minorRadius": 0.25, "thickness": 0.06, "twist": 3, "scale": {"W":2}, "rotDeg": {"xw":30}, ` + at + mat + `},
			{"twist": 2, ` + at + mat + `}],
		"projectivePlanes": [{` + at + mat + `}, {"radius": 0.3, "thickness": 0.04, ` + at + mat + `}],
		"mobiusBands": [{` + at + mat + `}, {"majorRadius": 0.7, "halfWidth": 0.2, "thickness": 0.05, "twist": 4, ` + at + mat + `}],
		"paraboloids": [{"scale": {"X":0.5,"Y":0.5,"Z":0.5,"W":0.8}, ` + at + mat + `}],
		"hyperboloids": [{` + at + mat + `}, {"waist": 0, ` + at + mat + `}, {"waist": 0.8, ` + at + mat + `}],
		"sphericalSegments": [{` + at + mat + `}, {"lower": -0.4, "upper": 0, ` + at + mat + `}, {"lower": -1, ` + at + mat + `}],
		"lenses": [{` + at + mat + `}, {"thickness": 0.2, ` + at + mat + `}],
		"reuleauxSolids": [{` + at + mat + `}, {"kind": "simplex", ` + at + mat + `}],
		"bracketopes": [{"shape": "[(xyz)(xyw)]", ` + at + mat + `}]
	}`
	if err := os.WriteFile(path, []byte(data), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg, err := loadConfig(path)
	if err != nil {
		t.Fatal(err)
	}
	scene, _, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(scene.KleinSurfaces) != 7 || len(scene.RevolvedQuadrics) != 7 || len(scene.ConvexSolids) != 5 || scene.NObjects() != 19 {
		t.Fatalf("objects: klein=%d quadrics=%d convex=%d total=%d", len(scene.KleinSurfaces), len(scene.RevolvedQuadrics), len(scene.ConvexSolids), scene.NObjects())
	}
	var got []string
	for _, k := range scene.KleinSurfaces {
		got = append(got, fmt.Sprintf("%s %g/%g/%g/%d", k.Name, k.MajorRadius, k.MinorRadius, k.Thickness, k.Twist))
		if k.Center != (Point4{0.1, 0.2, 0.3, 0.4}) {
			t.Errorf("%s: center %v", k.Name, k.Center)
		}
	}
	// name, radii, thickness and twist in the order of the keys
	if got, want := strings.Join(got, ", "), "klein bottle 1/0.4/0.1/1, klein bottle 0.6/0.25/0.06/3, twisted torus 1/0.4/0.1/2, "+
		"projective plane 0.5/0.5/0.1/1, projective plane 0.3/0.3/0.04/1, mobius band 1/0.4/0.1/1, twisted ring 0.7/0.2/0.05/4"; got != want {
		t.Errorf("walls around surfaces:\n got %s\nwant %s", got, want)
	}
	if k := scene.KleinSurfaces[1]; k.Scale != (Vector4{1, 1, 1, 2}) || k.R == I4() || !scene.KleinSurfaces[5].Band || k.Band {
		t.Errorf("klein bottle: scale %v, rotation %v", k.Scale, k.R)
	}
	got = got[:0]
	for _, q := range scene.RevolvedQuadrics {
		got = append(got, fmt.Sprintf("%s %.6g/%.6g/%.6g %g..%g", q.Name, q.K0, q.K1, q.K2, q.Lower, q.Upper))
	}
	if got, want := strings.Join(got, ", "), "paraboloid 0.5/0.5/0 -1..1, hyperboloid 0.25/0/0.75 -1..1, hyperboloid 0/0/1 -1..1, "+
		"hyperboloid 0.64/0/0.36 -1..1, spherical segment 1/0/-1 0..1, spherical segment 1/0/-1 -0.4..0, spherical segment 1/0/-1 -1..1"; got != want {
		t.Errorf("quadrics:\n got %s\nwant %s", got, want)
	}
	if q := scene.RevolvedQuadrics[0]; q.Scale != (Vector4{0.5, 0.5, 0.5, 0.8}) {
		t.Errorf("paraboloid scale %v", q.Scale)
	}
	got = got[:0]
	for _, c := range scene.ConvexSolids {
		got = append(got, fmt.Sprintf("%s %.3f", c.Name, c.AABBMax.W-c.AABBMin.W))
	}
	// name and size along W: the lenses are as thick as told, the rotor is 2 wide
	if got, want := strings.Join(got, ", "), "bracketope [(xyz)(xyw)] 2.000, lens 0.600, lens 0.200, reuleaux rotor 2.000, reuleaux simplex 1.581"; got != want {
		t.Errorf("convex solids:\n got %s\nwant %s", got, want)
	}

	// a wrong entry is reported with its key and place
	for _, c := range []struct{ old, new, where string }{
		{`"majorRadius": 0.6, "minorRadius": 0.25`, `"majorRadius": 0.2, "minorRadius": 0.25`, "kleinBottles[1]"},
		{`"twist": 2`, `"twist": 99`, "kleinBottles[2]"},
		{`"radius": 0.3`, `"radius": -0.3`, "projectivePlanes[1]"},
		{`"halfWidth": 0.2`, `"halfWidth": -0.2`, "mobiusBands[1]"},
		{`"scale": {"X":0.5,"Y":0.5,"Z":0.5,"W":0.8}`, `"scale": {"X":0.5,"Y":-0.5,"Z":0.5,"W":0.8}`, "paraboloids[0]"},
		{`"waist": 0.8`, `"waist": 1.8`, "hyperboloids[2]"},
		{`"lower": -0.4, "upper": 0`, `"lower": 0.4, "upper": 0`, "sphericalSegments[1]"},
		{`"thickness": 0.2`, `"thickness": 2.5`, "lenses[1]"},
		{`"kind": "simplex"`, `"kind": "cube"`, "reuleauxSolids[1]"},
		{`"shape": "[(xyz)(xyw)]"`, `"shape": "[(xyz)(xy)]"`, "bracketopes[0]"},
	} {
		if !strings.Contains(data, c.old) {
			t.Fatalf("%q is not in the scene", c.old)
		}
		if err := os.WriteFile(path, []byte(strings.Replace(data, c.old, c.new, 1)), 0o644); err != nil {
			t.Fatal(err)
		}
		cfg, err := loadConfig(path)
		if err != nil {
			t.Fatal(err)
		}
		if _, _, err = buildScene(cfg); err == nil || !strings.Contains(err.Error(), c.where) {
			t.Errorf("%s must be rejected with its place, got %v", c.new, err)
		}
	}
}
