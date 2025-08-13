package photons4d

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"os"
)

func (s *Scene) SaveRawRGB64(path string) error {
	// Sanity checks
	if s.Nx < 0 || s.Ny < 0 || s.Nz < 0 {
		return fmt.Errorf("negative dimensions: Nx=%d Ny=%d Nz=%d", s.Nx, s.Ny, s.Nz)
	}
	// Expect exactly Nx*Ny*Nz*3 values in Buf.
	exp := s.Nx * s.Ny * s.Nz * 3
	if len(s.Buf) != exp {
		return fmt.Errorf("Buf length mismatch: got %d, expected %d (Nx*Ny*Nz*3)", len(s.Buf), exp)
	}

	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	// Header: Nx, Ny, Nz as int32 (little-endian)
	if err := binary.Write(w, binary.LittleEndian, int32(s.Nx)); err != nil {
		return err
	}
	if err := binary.Write(w, binary.LittleEndian, int32(s.Ny)); err != nil {
		return err
	}
	if err := binary.Write(w, binary.LittleEndian, int32(s.Nz)); err != nil {
		return err
	}

	// Body: full Buf as float64, in current memory order
	for i := range s.Buf {
		if err := binary.Write(w, binary.LittleEndian, float64(s.Buf[i])); err != nil {
			return err
		}
	}

	if err := w.Flush(); err != nil {
		return err
	}
	_ = f.Sync() // optional

	return nil
}
