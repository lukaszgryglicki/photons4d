package photons4d

import (
	"os"
	"path/filepath"
	"testing"
	"time"
)

func TestMeasurePackRatio(t *testing.T) {
	if os.Getenv("MEASURE") == "" {
		t.Skip("set MEASURE=1")
	}
	dir := t.TempDir()
	scenePath := writeDistTestScene(t, dir, filepath.Join(dir, "x.gif"))
	cfg, err := prepareConfig(scenePath)
	if err != nil {
		t.Fatal(err)
	}
	cfg.SceneResX, cfg.SceneResY, cfg.SceneResZ = 128, 128, 48
	scene, lights, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	castRays(lights, scene, []int{8000000})
	idx, vals := extractSparseAndZero(scene.Buf)
	n := len(vals)
	t.Logf("entries: %d (buffer %d, %.1f%% dense)", n, len(scene.Buf), 100*float64(n)/float64(len(scene.Buf)))
	t0 := time.Now()
	packed, err := packSparse(idx, vals)
	if err != nil {
		t.Fatal(err)
	}
	tPack := time.Since(t0)
	rawLen := 0
	{
		_, _, rl, err := unpackSparse(packed)
		if err != nil {
			t.Fatal(err)
		}
		rawLen = rl
	}
	t1 := time.Now()
	for i := 0; i < 3; i++ {
		_, _, _, _ = unpackSparse(packed)
	}
	tUnpack := time.Since(t1) / 3
	gobEst := 0
	for _, d := range idx {
		z := uint64(d)
		b := 1
		for z >= 0x80 {
			z >>= 7
			b++
		}
		gobEst += b
	}
	gobEst += 9 * n
	t.Logf("raw payload: %d B (%.1f B/entry) | gob-estimate: %d B (%.1f B/entry)", rawLen, float64(rawLen)/float64(n), gobEst, float64(gobEst)/float64(n))
	t.Logf("packed: %d B (%.2f B/entry) | ratio vs raw %.2fx, vs gob %.2fx", len(packed), float64(len(packed))/float64(n), float64(rawLen)/float64(len(packed)), float64(gobEst)/float64(len(packed)))
	t.Logf("pack: %s (%.0f MB/s) | unpack: %s (%.0f MB/s)", tPack, float64(rawLen)/tPack.Seconds()/1e6, tUnpack, float64(rawLen)/tUnpack.Seconds()/1e6)
}
