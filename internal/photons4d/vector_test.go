package photons4d

import (
	"math"
	"testing"
)

func TestVectorOps(t *testing.T) {
	v := Vector4{1, 2, 3, 4}
	w := Vector4{-1, 0.5, 2, -2}
	s := Real(3)

	add := v.Add(w)
	if add != (Vector4{0, 2.5, 5, 2}) {
		t.Fatalf("Add mismatch: %+v", add)
	}
	sub := v.Sub(w)
	if sub != (Vector4{2, 1.5, 1, 6}) {
		t.Fatalf("Sub mismatch: %+v", sub)
	}
	mul := v.Mul(s)
	if mul != (Vector4{3, 6, 9, 12}) {
		t.Fatalf("Mul mismatch: %+v", mul)
	}
	dot := v.Dot(w)
	wantDot := Real(1*(-1) + 2*0.5 + 3*2 + 4*(-2))
	if dot != wantDot {
		t.Fatalf("Dot mismatch: got %.12g want %.12g", dot, wantDot)
	}
	l := v.Len()
	if math.Abs(float64(l-math.Sqrt(30))) > 1e-12 {
		t.Fatalf("Len mismatch: %.12g", l)
	}
	n := v.Norm()
	if math.Abs(float64(n.Len()-1)) > 1e-12 {
		t.Fatalf("Norm not unit: %.12g", n.Len())
	}
}
