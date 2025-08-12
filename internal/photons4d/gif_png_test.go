package photons4d

import (
	"os"
	"path/filepath"
	"testing"
)

func tinyScene() *Scene {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 2, 2, 1, 4)
	// Put a single bright voxel so files are not empty/black
	i, j, k := 1, 1, 0
	base := s.idx(i, j, k, ChR)
	s.Buf[base+0] = 1.0
	s.Buf[base+1] = 0.5
	s.Buf[base+2] = 0.25
	return s
}

func TestSaveAnimatedGIF(t *testing.T) {
	s := tinyScene()
	tmp := filepath.Join(t.TempDir(), "out.gif")
	if err := SaveAnimatedGIF(s, tmp, 5, 0.8); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(tmp); err != nil {
		t.Fatalf("gif not written: %v", err)
	}
}

func TestSavePNGSequence16(t *testing.T) {
	s := tinyScene()
	prefix := filepath.Join(t.TempDir(), "frame")
	if err := SavePNGSequence16(s, prefix, 0.8); err != nil {
		t.Fatal(err)
	}
	// Only one slice => "_0.png"
	f := prefix + "_0.png"
	if _, err := os.Stat(f); err != nil {
		t.Fatalf("png not written: %v", err)
	}
}
