package photons4d

import "testing"

func TestSceneMapping(t *testing.T) {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 4, 6, 10, 20, 30, 8, false)
	if s.StrideY != s.Nz*3 {
		t.Fatalf("StrideY wrong: %d", s.StrideY)
	}
	if s.StrideX != s.Ny*s.StrideY {
		t.Fatalf("StrideX wrong: %d", s.StrideX)
	}
	// bounds: X∈[-1,1], Y∈[-2,2], Z∈[-3,3]
	ok, i, j, k, ux, uy, uz := s.VoxelIndexOf(Point4{0, 0, 0, 0})
	if !ok || i != 5 || j != 10 || k != 15 {
		t.Fatalf("IndexOf center wrong: ok=%v i=%d j=%d k=%d (ux=%.3f uy=%.3f uz=%.3f)", ok, i, j, k, ux, uy, uz)
	}
}

func TestIdxAndVoxelSize(t *testing.T) {
	s := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 4, 4, 4, 4, false)
	dx, dy, dz := s.VoxelSize()
	if dx != 0.5 || dy != 0.5 || dz != 0.5 {
		t.Fatalf("voxel size wrong: %.3f %.3f %.3f", dx, dy, dz)
	}
	idx := s.idx(1, 2, 3, ChG)
	// idx = 1*StrideX + 2*StrideY + 3*3 + 1
	if idx != 1*s.StrideX+2*s.StrideY+3*3+1 {
		t.Fatalf("idx wrong: %d", idx)
	}
}
