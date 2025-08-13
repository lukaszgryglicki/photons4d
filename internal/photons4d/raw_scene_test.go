package photons4d

import (
	"bufio"
	"encoding/binary"
	"os"
	"path/filepath"
	"testing"
)

// build a Buf laid out as (((i*Ny)+j)*Nz + k)*3 + c
func makeBuf(Nx, Ny, Nz int) []Real {
	n := Nx * Ny * Nz * 3
	buf := make([]Real, n)
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			for k := 0; k < Nz; k++ {
				base := float64(i*100 + j*10 + k)
				idx := (((i*Ny)+j)*Nz + k) * 3
				buf[idx+0] = Real(base + 0.1) // R
				buf[idx+1] = Real(base + 0.2) // G
				buf[idx+2] = Real(base + 0.3) // B
			}
		}
	}
	return buf
}

func TestSceneSaveRawRGB64_Basic(t *testing.T) {
	Nx, Ny, Nz := 2, 2, 3
	s := &Scene{
		Nx:  Nx,
		Ny:  Ny,
		Nz:  Nz,
		Buf: makeBuf(Nx, Ny, Nz),
	}

	path := filepath.Join(t.TempDir(), "scene.raw")
	if err := s.SaveRawRGB64(path); err != nil {
		t.Fatalf("SaveRawRGB64 error: %v", err)
	}

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open result file: %v", err)
	}
	defer f.Close()

	r := bufio.NewReader(f)

	// Header
	var hx, hy, hz int32
	if err := binary.Read(r, binary.LittleEndian, &hx); err != nil {
		t.Fatalf("read Nx: %v", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &hy); err != nil {
		t.Fatalf("read Ny: %v", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &hz); err != nil {
		t.Fatalf("read Nz: %v", err)
	}
	if int(hx) != Nx || int(hy) != Ny || int(hz) != Nz {
		t.Fatalf("header mismatch got (%d,%d,%d) want (%d,%d,%d)", hx, hy, hz, Nx, Ny, Nz)
	}

	// Body
	expCount := Nx * Ny * Nz * 3
	got := make([]float64, expCount)
	for i := 0; i < expCount; i++ {
		if err := binary.Read(r, binary.LittleEndian, &got[i]); err != nil {
			t.Fatalf("read body[%d]: %v", i, err)
		}
	}

	expected := make([]float64, expCount)
	for i := range s.Buf {
		expected[i] = float64(s.Buf[i])
	}

	if len(got) != len(expected) {
		t.Fatalf("value length mismatch got %d want %d", len(got), len(expected))
	}
	for i := range got {
		if got[i] != expected[i] {
			t.Fatalf("value %d mismatch got %v want %v", i, got[i], expected[i])
		}
	}

	st, err := f.Stat()
	if err != nil {
		t.Fatalf("stat file: %v", err)
	}
	wantSize := int64(12 + 8*expCount) // 3*int32 header + N*float64
	if st.Size() != wantSize {
		t.Fatalf("file size mismatch got %d want %d", st.Size(), wantSize)
	}
}

func TestSceneSaveRawRGB64_Errors(t *testing.T) {
	t.Run("negative dims", func(t *testing.T) {
		s := &Scene{Nx: -1, Ny: 1, Nz: 1, Buf: make([]Real, 3)}
		if err := s.SaveRawRGB64(filepath.Join(t.TempDir(), "neg.raw")); err == nil {
			t.Fatalf("expected error for negative dims, got nil")
		}
	})
	t.Run("buf length mismatch", func(t *testing.T) {
		s := &Scene{Nx: 1, Ny: 1, Nz: 1, Buf: make([]Real, 2)}
		if err := s.SaveRawRGB64(filepath.Join(t.TempDir(), "mismatch.raw")); err == nil {
			t.Fatalf("expected error for buf length mismatch, got nil")
		}
	})
}

func TestSceneSaveRawRGB64_Zero(t *testing.T) {
	s := &Scene{Nx: 0, Ny: 0, Nz: 0, Buf: nil}
	path := filepath.Join(t.TempDir(), "zero.raw")
	if err := s.SaveRawRGB64(path); err != nil {
		t.Fatalf("unexpected err: %v", err)
	}
	st, _ := os.Stat(path)
	if st.Size() != 12 { // header only
		t.Fatalf("want 12 bytes, got %d", st.Size())
	}
}
