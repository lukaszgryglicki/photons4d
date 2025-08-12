package photons4d

import "testing"

func TestPointAdd(t *testing.T) {
	p := Point4{1, 2, 3, 4}
	v := Vector4{-1, 1, 0, 2}
	q := p.Add(v)
	if q != (Point4{0, 3, 3, 6}) {
		t.Fatalf("Add mismatch: %+v", q)
	}
}
