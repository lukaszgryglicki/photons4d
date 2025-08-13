package photons4d

import (
	"fmt"
	"image"
	"image/png"
	"math"
	"os"
)

// SavePNGSequence16 writes one 16-bit PNG per Z slice (k = 0..Nz-1).
// Each frame is lossless PNG (DEFLATE) and uses 16 bits per channel.
//
// Note: PNG is lossless; the only quantization here is your own mapping
// from float radiance -> 16-bit via per-slice normalization and gamma.
func SavePNGSequence16(scene *Scene, prefix string, gamma Real) error {
	Nx, Ny, Nz := scene.Nx, scene.Ny, scene.Nz

	// Helper: map scalar -> [0..65535] with gamma.
	toU16 := func(v, scale Real) uint16 {
		if v <= 0 {
			return 0
		}
		n := v * scale // ideally in [0,1]
		if n > 1 {
			n = 1
		}
		if gamma != 1 {
			n = math.Pow(n, 1.0/gamma)
		}
		x := math.Round(n * 65535.0)
		if x < 0 {
			return 0
		}
		if x > 65535 {
			return 65535
		}
		return uint16(x)
	}

	// Zero-padding width based on number of slices.
	width := 1
	if Nz > 1 {
		width = int(math.Log10(Real(Nz-1))) + 1
	}

	// Progress print step (~1%).
	step := 1
	if Nz >= 100 {
		step = Nz / 100
	}

	for k := 0; k < Nz; k++ {
		if k%step == 0 {
			percent := Real(k+1) * 100 / Real(Nz)
			fmt.Printf("[PNG]  %.2f%%\n", percent)
		}

		// 1) Find max over this slice (peak across RGB channels).
		sliceMax := 0.0
		for j := 0; j < Ny; j++ {
			for i := 0; i < Nx; i++ {
				base := scene.idx(i, j, k, ChR)
				if r := Real(scene.Buf[base+0]); r > sliceMax {
					sliceMax = r
				}
				if g := Real(scene.Buf[base+1]); g > sliceMax {
					sliceMax = g
				}
				if b := Real(scene.Buf[base+2]); b > sliceMax {
					sliceMax = b
				}
			}
		}
		if sliceMax == 0 {
			sliceMax = 1 // avoid div-by-zero; the slice will be black
		}
		scale := 1.0 / sliceMax

		// 2) Fill a 16-bit image (flip Y so up is up).
		img := image.NewNRGBA64(image.Rect(0, 0, Nx, Ny))
		const pxBytes = 8 // 4 channels * 2 bytes/channel
		for j := 0; j < Ny; j++ {
			y := Ny - 1 - j
			rowOff := y * img.Stride
			for i := 0; i < Nx; i++ {
				base := scene.idx(i, j, k, ChR)
				r := toU16(Real(scene.Buf[base+0]), scale)
				g := toU16(Real(scene.Buf[base+1]), scale)
				b := toU16(Real(scene.Buf[base+2]), scale)
				a := uint16(0xFFFF)

				p := rowOff + i*pxBytes
				// NRGBA64 stores big-endian uint16 per channel: R,G, B, A.
				img.Pix[p+0] = uint8(r >> 8)
				img.Pix[p+1] = uint8(r)
				img.Pix[p+2] = uint8(g >> 8)
				img.Pix[p+3] = uint8(g)
				img.Pix[p+4] = uint8(b >> 8)
				img.Pix[p+5] = uint8(b)
				img.Pix[p+6] = uint8(a >> 8)
				img.Pix[p+7] = uint8(a)
			}
		}

		// 3) Write PNG (lossless). 16-bit is preserved because we used NRGBA64.
		full := fmt.Sprintf("%s_%0*d.png", prefix, width, k)
		f, err := os.Create(full)
		if err != nil {
			return err
		}

		enc := png.Encoder{CompressionLevel: png.BestCompression} // still lossless
		if err := enc.Encode(f, img); err != nil {
			f.Close()
			return err
		}
		if err := f.Close(); err != nil {
			return err
		}
	}

	return nil
}
