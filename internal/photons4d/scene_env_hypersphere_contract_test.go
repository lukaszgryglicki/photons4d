package photons4d

import "testing"

func TestAngleIndexOf_OppositeDirectionsDiffer(t *testing.T) {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 10, 10, 12, 2, true)

	// Choose a generic direction not aligned to axes.
	D := Vector4{X: 0.3, Y: -0.7, Z: 0.2, W: 0.6}
	i1, j1, k1, _, _, _ := s.AngleIndexOf(D)
	i2, j2, k2, _, _, _ := s.AngleIndexOf(Vector4{X: -D.X, Y: -D.Y, Z: -D.Z, W: -D.W})

	if i1 == i2 && j1 == j2 && k1 == k2 {
		t.Fatalf("opposite directions fell into the same angular voxel: (%d,%d,%d)", i1, j1, k1)
	}
}

func TestAngleIndexOf_IndexRanges(t *testing.T) {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 7, 9, 11, 1, true)

	vecs := []Vector4{
		{1, 2, 3, 4},
		{-1, 0.1, 0.1, 0.1},
		{0, -1, 0.5, -0.25},
		{0, 0, 1, 0},
		{0, 0, 0, -1},
	}
	for _, v := range vecs {
		i, j, k, _, _, _ := s.AngleIndexOf(v)
		if i < 0 || i >= s.Nx || j < 0 || j >= s.Ny || k < 0 || k >= s.Nz {
			t.Fatalf("out-of-range angle index for v=%+v: got (%d,%d,%d) within (%d,%d,%d)",
				v, i, j, k, s.Nx, s.Ny, s.Nz)
		}
	}
}
