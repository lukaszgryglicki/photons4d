package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

func TestNewLightValidation(t *testing.T) {
	_, err := NewLight(Point4{}, Vector4{}, RGB{1, 1, 1}, 0.1, 1.0)
	if err == nil {
		t.Fatal("expected error for zero dir")
	}
	_, err = NewLight(Point4{}, Vector4{1, 0, 0, 0}, RGB{0, 0, 0}, math.Pi/4, 0.5)
	if err == nil {
		t.Fatal("expected error for zero color sum")
	}
	_, err = NewLight(Point4{}, Vector4{1, 0, 0, 0}, RGB{1, 1, 1}, 0, 0.0)
	if err == nil {
		t.Fatal("expected error for angle out of range")
	}
}

func TestSampleDirInsideConeUnit(t *testing.T) {
	L, err := NewLight(Point4{}, Vector4{0, 0, 0, 1}, RGB{1, 0, 0}, math.Pi/6, -1.0)
	if err != nil {
		t.Fatal(err)
	}
	rng := rand.New(rand.NewSource(1))
	for i := 0; i < 1000; i++ {
		d := L.SampleDir(rng)
		if math.Abs(float64(d.Len()-1)) > 1e-9 {
			t.Fatalf("SampleDir not unit: %.12g", d.Len())
		}
		// Inside cap: dot >= cosAngle
		if d.Dot(L.Direction) < L.cosAngle-1e-12 {
			t.Fatalf("Sample outside cone: dot=%.12g cos=%.12g", d.Dot(L.Direction), L.cosAngle)
		}
		// Check ortho basis is indeed orthogonal to Direction
		if math.Abs(float64(L.U.Dot(L.Direction))) > 1e-9 ||
			math.Abs(float64(L.V.Dot(L.Direction))) > 1e-9 ||
			math.Abs(float64(L.W.Dot(L.Direction))) > 1e-9 {
			t.Fatal("U/V/W not orthogonal to Direction")
		}
	}
}
