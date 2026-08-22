package photons4d

import (
	"math"
	"sync"
)

// ---- helpers ----
func evenPerms4() [][]int {
	// All even permutations of (0,1,2,3).
	// There are 12 of them.
	base := [][]int{
		{0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
		{1, 0, 3, 2}, {1, 2, 0, 3}, {1, 3, 2, 0},
		{2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1},
		{3, 0, 2, 1}, {3, 1, 0, 2}, {3, 2, 1, 0},
	}
	return base
}

func allPerms4Distinct(a [4]Real) [][4]Real {
	// 24 permutations of distinct entries (we will dedup outside if needed)
	idx := [][]int{
		{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
		{1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
		{2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
		{3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0},
	}
	out := make([][4]Real, 0, 24)
	for _, p := range idx {
		out = append(out, [4]Real{a[p[0]], a[p[1]], a[p[2]], a[p[3]]})
	}
	return out
}

func signVariants(vals [4]Real, maskNonZero [4]bool, requireEvenMinus bool) [][4]Real {
	out := make([][4]Real, 0, 16)
	for s := 0; s < 16; s++ {
		v := vals

		// start with how many masked entries are already negative
		minus := 0
		for i := 0; i < 4; i++ {
			if maskNonZero[i] && v[i] < 0 {
				minus++
			}
		}

		// apply flips; each flip toggles sign ⇒ toggles minus parity
		for i := 0; i < 4; i++ {
			if !maskNonZero[i] {
				continue
			}
			if (s>>i)&1 == 1 {
				v[i] = -v[i]
				if v[i] < 0 {
					minus++
				} else {
					minus--
				}
			}
		}

		if requireEvenMinus && (minus&1) != 0 {
			continue
		}
		out = append(out, v)
	}
	return out
}

func pushUnique(set map[[4]int64]struct{}, out *[]Vector4, v [4]Real) {
	// key with 1e-12 quantization
	const q = 1e12
	k := [4]int64{
		int64(math.Round(Real(v[0] * q))),
		int64(math.Round(Real(v[1] * q))),
		int64(math.Round(Real(v[2] * q))),
		int64(math.Round(Real(v[3] * q))),
	}
	if _, ok := set[k]; ok {
		return
	}
	set[k] = struct{}{}
	*out = append(*out, Vector4{v[0], v[1], v[2], v[3]}.Norm())
}

func verts600Unit() []Vector4 {
	phi := (1 + math.Sqrt(5)) / 2
	inv := 1 / phi

	set := make(map[[4]int64]struct{}, 128)
	out := make([]Vector4, 0, 120)

	// 8 axes
	for a := 0; a < 4; a++ {
		for s := -1; s <= 1; s += 2 {
			v := [4]Real{0, 0, 0, 0}
			v[a] = Real(s)
			pushUnique(set, &out, v)
		}
	}

	// 16 of (±1/2, ±1/2, ±1/2, ±1/2)
	b := Real(0.5)
	for sx := -1; sx <= 1; sx += 2 {
		for sy := -1; sy <= 1; sy += 2 {
			for sz := -1; sz <= 1; sz += 2 {
				for sw := -1; sw <= 1; sw += 2 {
					pushUnique(set, &out, [4]Real{b * Real(sx), b * Real(sy), b * Real(sz), b * Real(sw)})
				}
			}
		}
	}

	// 96: EVEN permutations of (φ/2, 1/2, 1/(2φ), 0) with ALL sign
	// combinations of the nonzero coordinates (12 × 8 = 96).
	// The previous code used all 24 permutations restricted to an even number
	// of minus signs — a different set that is NOT the 600-cell snub family
	// and produced a deformed polytope.
	base := [4]Real{Real(0.5 * phi), Real(0.5), Real(0.5 * inv), 0}
	for _, p := range evenPerms4() {
		v := [4]Real{base[p[0]], base[p[1]], base[p[2]], base[p[3]]}
		mask := [4]bool{v[0] != 0, v[1] != 0, v[2] != 0, v[3] != 0}
		for _, vv := range allSignVariants(v, mask) {
			pushUnique(set, &out, vv)
		}
	}

	if len(out) != 120 {
		DebugLog("verts600Unit: expected 120, got %d", len(out))
	}
	return out
}

// allSignVariants returns every sign combination of the masked (nonzero) coordinates.
func allSignVariants(vals [4]Real, maskNonZero [4]bool) [][4]Real {
	return signVariants(vals, maskNonZero, false)
}

var (
	v120Once  sync.Once
	v120Cache []Vector4
)

// verts120Unit returns the 600 unit vertex directions of the regular
// 120-cell that is the polar dual of conv(verts600Unit()).
//
// They are computed as the (normalized) centers of the 600 tetrahedral
// facets of the 600-cell, i.e. of all 4-cliques of its edge graph. This
// guarantees EXACT dual alignment with verts600Unit(): the polytope
// {x : v·x <= r for all v in verts600Unit()} has its vertices precisely
// along these directions (all with the same support scale), which is what
// cellPoly relies on for its AABB.
//
// The previous hand-written coordinate table was not a regular 120-cell
// (wrong vertex families and sign/permutation rules), so the 600-cell
// built from it was irregular and both objects' AABBs missed geometry.
func verts120Unit() []Vector4 {
	v120Once.Do(func() {
		vs := verts600Unit()
		n := len(vs)

		// Edge length of the unit-circumradius 600-cell is 1/phi (~0.618).
		minD := Real(1e300)
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				if d := vs[i].Sub(vs[j]).Len(); d > 1e-9 && d < minD {
					minD = d
				}
			}
		}
		tol := minD * 1.0001

		adj := make([][]bool, n)
		nbr := make([][]int, n)
		for i := range adj {
			adj[i] = make([]bool, n)
		}
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				if vs[i].Sub(vs[j]).Len() <= tol {
					adj[i][j], adj[j][i] = true, true
					nbr[i] = append(nbr[i], j)
					nbr[j] = append(nbr[j], i)
				}
			}
		}

		out := make([]Vector4, 0, 600)
		for i := 0; i < n; i++ {
			for _, j := range nbr[i] {
				if j <= i {
					continue
				}
				for _, k := range nbr[i] {
					if k <= j || !adj[j][k] {
						continue
					}
					for _, l := range nbr[i] {
						if l <= k || !adj[j][l] || !adj[k][l] {
							continue
						}
						c := vs[i].Add(vs[j]).Add(vs[k]).Add(vs[l])
						out = append(out, c.Norm())
					}
				}
			}
		}
		if len(out) != 600 {
			DebugLog("verts120Unit: expected 600 tetrahedral cells, got %d", len(out))
		}
		v120Cache = out
	})
	return v120Cache
}
