package photons4d

import (
	"math"
	"testing"
)

func TestI4MulVec(t *testing.T) {
	I := I4()
	v := Vector4{1, 2, 3, 4}
	out := I.MulVec(v)
	if out != v {
		t.Fatalf("I*v != v: %+v", out)
	}
}

func TestTransposeAndMul(t *testing.T) {
	// simple nontrivial matrix
	M := Mat4{M: [4][4]Real{
		{1, 2, 3, 4},
		{0, 1, 0, 0.5},
		{2, 0, 1, -1},
		{0, 0, 0.25, 1},
	}}
	T := M.Transpose()
	// check transpose symmetry for a couple elements
	if T.M[0][1] != M.M[1][0] || T.M[3][2] != M.M[2][3] {
		t.Fatal("Transpose mismatch")
	}

	// (M^T M) should be symmetric positive-definite (just basic check)
	S := T.Mul(M)
	if math.Abs(float64(S.M[0][1]-S.M[1][0])) > 1e-12 {
		t.Fatal("M^T M not symmetric")
	}
}
