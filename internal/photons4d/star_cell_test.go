package photons4d

import (
	"math"
	"testing"
)

func TestStarCellDirectionCounts(t *testing.T) {
	cases := []struct {
		kind  string
		alias string
		want  int
	}{
		{"star-cell-5", "5", 5},
		{"star-cell-8", "8-cell", 8},
		{"star-cell-16", "star16", 16},
		{"star-cell-24", "24", 24},
		{"star-cell-120", "cell-120", 120},
		{"star-cell-600", "600", 600},
	}
	for _, tc := range cases {
		canon, dirs, err := starCellDirections(tc.alias)
		if err != nil {
			t.Fatalf("%s: unexpected error: %v", tc.alias, err)
		}
		if canon != tc.kind {
			t.Fatalf("%s: canonical kind mismatch: got %q want %q", tc.alias, canon, tc.kind)
		}
		if len(dirs) != tc.want {
			t.Fatalf("%s: direction count mismatch: got %d want %d", tc.kind, len(dirs), tc.want)
		}
		for i, d := range dirs {
			if math.Abs(d.Len()-1) > 1e-12 {
				t.Fatalf("%s: dir[%d] not unit: len=%.12g", tc.kind, i, d.Len())
			}
		}
	}
}

func TestNewStarCellAndIntersect(t *testing.T) {
	sc, err := NewStarCell(
		"8",
		Point4{0, 0, 0, 0.6},
		Vector4{0.2, 0.2, 0.2, 0.2},
		0.45,
		0.55,
		6,
		Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2},
	)
	if err != nil {
		t.Fatalf("NewStarCell: %v", err)
	}
	if sc.Kind != "star-cell-8" {
		t.Fatalf("kind canonicalization failed: %q", sc.Kind)
	}

	O := Point4{0, 0, 0, 0}
	D := Vector4{0, 0, 0, 1}
	h, ok := intersectRayStarCell(O, D, sc)
	if !ok {
		t.Fatalf("expected star cell hit")
	}
	if h.inv {
		t.Fatalf("expected entering hit")
	}
	if math.Abs(float64(h.t-0.4)) > 1e-9 {
		t.Fatalf("t wrong: %.12g", h.t)
	}
	if math.Abs(float64(h.Nw.W+1)) > 1e-9 {
		t.Fatalf("expected -W normal, got %+v", h.Nw)
	}

	O2 := Point4{0, 0, 0, 0.6}
	D2 := Vector4{0, 0, 0, 1}
	h2, ok2 := intersectRayStarCell(O2, D2, sc)
	if !ok2 || !h2.inv {
		t.Fatalf("expected inside->exit hit, ok=%v inv=%v", ok2, h2.inv)
	}
	if math.Abs(float64(h2.t-0.2)) > 1e-9 {
		t.Fatalf("exit t wrong: %.12g", h2.t)
	}
	if math.Abs(float64(h2.Nw.W-1)) > 1e-9 {
		t.Fatalf("expected +W normal, got %+v", h2.Nw)
	}
}

func TestSceneAABBWithStarCell(t *testing.T) {
	sc, err := NewStarCell(
		"star-cell-120",
		Point4{0.1, -0.2, 0.05, 0.3},
		Vector4{0.24, 0.18, 0.21, 0.16},
		0.44,
		0.56,
		6.5,
		Rot4{XY: 0.2, XZ: -0.15, XW: 0.1, YW: -0.2, ZW: 0.3},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0.1, 0.1, 0.1}, RGB{0.8, 0.8, 0.8}, RGB{1.3, 1.3, 1.3},
	)
	if err != nil {
		t.Fatal(err)
	}
	if !(sc.AABBMin.X <= sc.AABBMax.X && sc.AABBMin.Y <= sc.AABBMax.Y &&
		sc.AABBMin.Z <= sc.AABBMax.Z && sc.AABBMin.W <= sc.AABBMax.W) {
		t.Fatalf("AABB invalid: [%+v .. %+v]", sc.AABBMin, sc.AABBMax)
	}

	localSamples := make([]Vector4, 0, len(sc.Dirs)+8)
	for _, d := range sc.Dirs {
		localSamples = append(localSamples, Vector4{
			d.X * sc.OuterRadius * sc.Scale.X,
			d.Y * sc.OuterRadius * sc.Scale.Y,
			d.Z * sc.OuterRadius * sc.Scale.Z,
			d.W * sc.OuterRadius * sc.Scale.W,
		})
	}
	localSamples = append(localSamples,
		Vector4{+sc.CoreRadius * sc.Scale.X, 0, 0, 0},
		Vector4{-sc.CoreRadius * sc.Scale.X, 0, 0, 0},
		Vector4{0, +sc.CoreRadius * sc.Scale.Y, 0, 0},
		Vector4{0, -sc.CoreRadius * sc.Scale.Y, 0, 0},
		Vector4{0, 0, +sc.CoreRadius * sc.Scale.Z, 0},
		Vector4{0, 0, -sc.CoreRadius * sc.Scale.Z, 0},
		Vector4{0, 0, 0, +sc.CoreRadius * sc.Scale.W},
		Vector4{0, 0, 0, -sc.CoreRadius * sc.Scale.W},
	)

	for i, v := range localSamples {
		p := sc.R.MulVec(v)
		w := Point4{sc.Center.X + p.X, sc.Center.Y + p.Y, sc.Center.Z + p.Z, sc.Center.W + p.W}
		if w.X < sc.AABBMin.X-1e-12 || w.X > sc.AABBMax.X+1e-12 ||
			w.Y < sc.AABBMin.Y-1e-12 || w.Y > sc.AABBMax.Y+1e-12 ||
			w.Z < sc.AABBMin.Z-1e-12 || w.Z > sc.AABBMax.Z+1e-12 ||
			w.W < sc.AABBMin.W-1e-12 || w.W > sc.AABBMax.W+1e-12 {
			t.Fatalf("sample %d outside AABB: %+v not in [%+v .. %+v]", i, w, sc.AABBMin, sc.AABBMax)
		}
	}
}

// Rays that cut the +X spike of an axis-aligned star just below its tip: the
// chord gets as thin as 1e-6 of the size and must still be met at its entry.
func TestStarCellThinSpikeIsNotSkipped(t *testing.T) {
	const core, spike, sharp = 0.45, 0.55, 6.0
	sc, err := NewStarCell("8", Point4{0, 0, 0, 0}, Vector4{1, 1, 1, 1}, core, spike, sharp, Rot4{},
		RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1, 1, 1}, RGB{1.2, 1.2, 1.2})
	if err != nil {
		t.Fatal(err)
	}
	across := Vector4{0, 0.37, 0.22, -0.13}
	t0 := across.Len()
	D := across.Norm()
	for i := 0; i < 400; i++ {
		x0 := 1 - math.Pow(10, -2-4*Real(i)/400)
		// half-width of the spike at the height x0, from its definition
		out := func(tau Real) Real {
			r := math.Hypot(x0, tau)
			return r - (core + spike*math.Pow(x0/r, sharp))
		}
		lo, hi := Real(0), Real(0.5)
		for k := 0; k < 100; k++ {
			if m := 0.5 * (lo + hi); out(m) > 0 {
				hi = m
			} else {
				lo = m
			}
		}
		O := Point4{x0, -across.Y, -across.Z, -across.W}
		h, ok := intersectRayStarCell(O, D, sc)
		if !ok || h.inv {
			t.Fatalf("x0=%.9f: the spike (half-width %.3g) was missed: ok=%v inv=%v", x0, lo, ok, h.inv)
		}
		if math.Abs(h.t-(t0-lo)) > 1e-9 {
			t.Fatalf("x0=%.9f: t=%.12g, want %.12g", x0, h.t, t0-lo)
		}
	}
}
