package photons4d

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"os"
	"path/filepath"
)

func (s *Scene) SaveRawRGB64(path string) error {
	// Sanity checks
	if s.Nx < 0 || s.Ny < 0 || s.Nz < 0 {
		return fmt.Errorf("negative dimensions: Nx=%d Ny=%d Nz=%d", s.Nx, s.Ny, s.Nz)
	}
	// Expect exactly Nx*Ny*Nz*3 values in Buf. Use 64-bit multiply to avoid overflow.
	exp64 := int64(s.Nx) * int64(s.Ny) * int64(s.Nz) * 3
	if int64(len(s.Buf)) != exp64 {
		return fmt.Errorf("Buf length mismatch: got %d, expected %d (Nx*Ny*Nz*3)", len(s.Buf), exp64)
	}

	// Make sure parent directory exists.
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
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

	// Body: write the full buffer in one shot as float64 (faster, fewer syscalls)
	if exp64 > 0 {
		// Convert to float64 if necessary
		// buf64 := make([]float64, len(s.Buf))
		// for i := range s.Buf {
		// 	buf64[i] = float64(s.Buf[i])
		// }
		if err := binary.Write(w, binary.LittleEndian, s.Buf); err != nil {
			return err
		}
	}

	if err := w.Flush(); err != nil {
		return err
	}
	_ = f.Sync() // optional

	return nil
}
