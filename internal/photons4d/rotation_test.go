package photons4d

import (
	"math"
	"testing"
)

func TestRotFromAngles_IsOrthonormal(t *testing.T) {
	R := rotFromAngles(Rot4{
		XY: math.Pi / 6,
		XZ: math.Pi / 7,
		XW: math.Pi / 5,
		YZ: math.Pi / 8,
		YW: math.Pi / 9,
		ZW: math.Pi / 10,
	})

	RT := R.Transpose()
	// Check R^T R ~ I
	P := RT.Mul(R)
	I := I4()
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			diff := math.Abs(float64(P.M[r][c] - I.M[r][c]))
			if diff > 1e-12 {
				t.Fatalf("R^T R != I at (%d,%d): %.3g", r, c, diff)
			}
		}
	}
}

func TestAxisRotations(t *testing.T) {
	// Single-plane rotations keep length and rotate only the intended coordinates.
	v := Vector4{1, 0, 0, 0}
	R := rotXY(math.Pi / 2)
	o := R.MulVec(v)
	// 90° in XY: (1,0,0,0) -> (0,1,0,0)
	if math.Abs(float64(o.X-0)) > 1e-12 || math.Abs(float64(o.Y-1)) > 1e-12 {
		t.Fatalf("rotXY failed: %+v", o)
	}
	if math.Abs(float64(o.Len()-1)) > 1e-12 {
		t.Fatalf("rotXY broke length: %.12g", o.Len())
	}
}
