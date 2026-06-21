package photons4d

import (
	"fmt"
	"image"
	"image/color/palette"
	"image/draw"
	"image/gif"
	"math"
	"os"
	"path/filepath"
)

// SaveAnimatedGIF writes a GIF with one frame per Z slice (k = 0..Nz-1).
// delay is in 100ths of a second (e.g., 5 => 20 fps).
// per-slice normalization + optional gamma (e.g., 0.7 brightens).
func SaveAnimatedGIF(scene *Scene, path string, delay int, gamma Real) error {
	Nx, Ny, Nz := scene.Nx, scene.Ny, scene.Nz

	out := &gif.GIF{
		Image:     make([]*image.Paletted, 0, Nz),
		Delay:     make([]int, 0, Nz),
		LoopCount: 0,
	}
	rgba := image.NewNRGBA(image.Rect(0, 0, Nx, Ny))

	// helper: scalar → 0..255 with gamma
	toByte := func(v, scale Real) uint8 {
		if v <= 0 {
			return 0
		}
		n := v * scale // 0..1 ideally
		if n > 1 {
			n = 1
		}
		if gamma != 1 {
			// Documented contract (see func comment): gamma<1 brightens.
			// For n in (0,1), n^gamma with gamma<1 is > n (brighter); the
			// previous n^(1/gamma) inverted this and darkened instead.
			n = math.Pow(n, gamma)
		}
		return uint8(math.Round(n * 255))
	}

	for k := 0; k < Nz; k++ {
		// Progress logging
		// (relies on max(int,int) from your debug helpers)
		if k%imax(1, Nz/100) == 0 { // ~1% steps
			percent := Real(k+1) * 100 / Real(Nz)
			fmt.Printf("[GIF] %.2f%%\n", percent)
		}
		// 1) find max over this slice
		sliceMax := 0.0
		for j := 0; j < Ny; j++ {
			for i := 0; i < Nx; i++ {
				idx := scene.idx(i, j, k, ChR)
				// peak across channels
				if r := Real(scene.Buf[idx+0]); r > sliceMax {
					sliceMax = r
				}
				if g := Real(scene.Buf[idx+1]); g > sliceMax {
					sliceMax = g
				}
				if b := Real(scene.Buf[idx+2]); b > sliceMax {
					sliceMax = b
				}
			}
		}
		if sliceMax == 0 {
			sliceMax = 1 // avoid div-by-zero, will be black anyway
		}
		scale := 1.0 / sliceMax

		// 2) fill RGBA (flip Y so up is up)
		for j := 0; j < Ny; j++ {
			y := Ny - 1 - j
			rowOff := y * rgba.Stride
			for i := 0; i < Nx; i++ {
				idx := scene.idx(i, j, k, ChR)
				r := toByte(Real(scene.Buf[idx+0]), scale)
				g := toByte(Real(scene.Buf[idx+1]), scale)
				b := toByte(Real(scene.Buf[idx+2]), scale)
				p := rowOff + i*4
				rgba.Pix[p+0] = r
				rgba.Pix[p+1] = g
				rgba.Pix[p+2] = b
				rgba.Pix[p+3] = 255
			}
		}

		// 3) Quantize to paletted for GIF
		pimg := image.NewPaletted(rgba.Bounds(), palette.Plan9)
		draw.FloydSteinberg.Draw(pimg, pimg.Bounds(), rgba, image.Point{})

		out.Image = append(out.Image, pimg)
		out.Delay = append(out.Delay, delay)
	}

	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(path), filepath.Base(path)+".*.tmp")
	if err != nil {
		return err
	}
	tmpPath := tmp.Name()
	defer func() {
		_ = tmp.Close()
		_ = os.Remove(tmpPath)
	}()
	if err := gif.EncodeAll(tmp, out); err != nil {
		return err
	}
	if err := tmp.Sync(); err != nil {
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	if err := os.Rename(tmpPath, path); err != nil {
		return err
	}
	return nil
}
