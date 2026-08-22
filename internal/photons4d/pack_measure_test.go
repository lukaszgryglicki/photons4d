package photons4d

import (
	"math"
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
	rawLen := 0
	type combo struct {
		name  string
		codec uint8
		level int
	}
	for _, cc := range []combo{
		{"flate-1", CodecFlate, 1}, {"flate-6", CodecFlate, 6}, {"flate-9", CodecFlate, 9},
		{"zstd-1", CodecZstd, 1}, {"zstd-3", CodecZstd, 3}, {"zstd-6", CodecZstd, 6}, {"zstd-9", CodecZstd, 9},
	} {
		t0 := time.Now()
		packed, err := packSparse(idx, vals, cc.codec, cc.level)
		if err != nil {
			t.Fatal(err)
		}
		tPack := time.Since(t0)
		gi, gv, rl, err := unpackSparse(packed, cc.codec)
		if err != nil {
			t.Fatal(err)
		}
		for i := range idx {
			if gi[i] != idx[i] || math.Float64bits(gv[i]) != math.Float64bits(vals[i]) {
				t.Fatalf("%s: lossless round-trip violated at %d", cc.name, i)
			}
		}
		rawLen = rl
		t1 := time.Now()
		for i := 0; i < 3; i++ {
			_, _, _, _ = unpackSparse(packed, cc.codec)
		}
		tUnpack := time.Since(t1) / 3
		t.Logf("%-8s packed: %9d B (%.2f B/entry) ratio %.2fx | pack %8s (%4.0f MB/s) unpack %8s (%4.0f MB/s)",
			cc.name, len(packed), float64(len(packed))/float64(n), float64(rawLen)/float64(len(packed)),
			tPack.Round(time.Millisecond), float64(rawLen)/tPack.Seconds()/1e6,
			tUnpack.Round(time.Millisecond), float64(rawLen)/tUnpack.Seconds()/1e6)
	}
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
}
