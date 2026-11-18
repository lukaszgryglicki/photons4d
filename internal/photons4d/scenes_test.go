package photons4d

import (
	"path/filepath"
	"testing"
)

// Every scene file of the repository must load and build.
func TestAllSceneFilesLoadAndBuild(t *testing.T) {
	files, err := filepath.Glob("../../scenes/*.json")
	if err != nil {
		t.Fatal(err)
	}
	more, err := filepath.Glob("../../scenes/showcase/*.json")
	if err != nil {
		t.Fatal(err)
	}
	files = append(files, more...)
	if len(files) < 100 {
		t.Fatalf("only %d scene files found", len(files))
	}
	for _, f := range files {
		cfg, err := loadConfig(f)
		if err != nil {
			t.Errorf("%s: load: %v", f, err)
			continue
		}
		cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ = 2, 2, 2
		scene, lights, err := buildScene(cfg)
		if err != nil {
			t.Errorf("%s: build: %v", f, err)
			continue
		}
		if len(lights) == 0 {
			t.Errorf("%s: no lights", f)
		}
		if scene.NObjects() != len(collectSceneObjects(scene)) {
			t.Errorf("%s: %d objects, %d of them in the BVH", f, scene.NObjects(), len(collectSceneObjects(scene)))
		}
	}
}
