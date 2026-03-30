package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

type seqSource struct {
	vals []int64
	i    int
}

func (s *seqSource) Int63() int64 {
	v := s.vals[s.i%len(s.vals)]
	s.i++
	return v
}

func (s *seqSource) Seed(seed int64) {}

func TestUnitS3_SkipsZeroSecondDisc(t *testing.T) {
	const half = int64(1) << 62
	const threeQuarter = int64(3) << 61
	rng := rand.New(&seqSource{vals: []int64{
		half, half, half, half,
		half, half, threeQuarter, half,
	}})
	v := unitS3(rng)
	if !isFinite(v.X) || !isFinite(v.Y) || !isFinite(v.Z) || !isFinite(v.W) {
		t.Fatalf("unitS3 returned non-finite vector: %+v", v)
	}
	if math.Abs(float64(v.Len()-1)) > 1e-12 {
		t.Fatalf("unitS3 returned non-unit vector: %+v len=%.12g", v, v.Len())
	}
}

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
