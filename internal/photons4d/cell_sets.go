package photons4d

import "math"

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
	// Flip signs of the entries marked true in maskNonZero.
	// If requireEvenMinus, keep only even number of minus signs (among masked).
	out := make([][4]Real, 0, 16)
	for s := 0; s < 16; s++ {
		v := vals
		minus := 0
		ok := true
		for i := 0; i < 4; i++ {
			if !maskNonZero[i] {
				continue
			}
			if (s>>i)&1 == 1 {
				v[i] = -v[i]
				minus++
			}
		}
		if requireEvenMinus && (minus&1) == 1 {
			ok = false
		}
		if ok {
			out = append(out, v)
		}
	}
	return out
}

func pushUnique(set map[[4]int64]struct{}, out *[]Vector4, v [4]Real) {
	// key with 1e-12 quantization
	const q = 1e12
	k := [4]int64{
		int64(math.Round(float64(v[0] * q))),
		int64(math.Round(float64(v[1] * q))),
		int64(math.Round(float64(v[2] * q))),
		int64(math.Round(float64(v[3] * q))),
	}
	if _, ok := set[k]; ok {
		return
	}
	set[k] = struct{}{}
	*out = append(*out, Vector4{v[0], v[1], v[2], v[3]}.Norm())
}

// ---- 600-cell vertices (unit-radius) ----
// From Wikipedia 600-cell "Unit radius Cartesian coordinates":
// 8 of (0,0,0,±1) permuted; 16 of (±1/2,±1/2,±1/2,±1/2); and
// 96 even permutations of (0, ±1/2, ±φ/2, ±1/(2φ)) with an even number of minus signs.
func verts600Unit() []Vector4 {
	phi := (1 + math.Sqrt(5)) / 2
	inv := 1 / phi

	set := make(map[[4]int64]struct{}, 128)
	out := make([]Vector4, 0, 120)

	// A: 8 axis points
	for a := 0; a < 4; a++ {
		for s := -1; s <= 1; s += 2 {
			v := [4]Real{0, 0, 0, 0}
			v[a] = Real(s)
			pushUnique(set, &out, v)
		}
	}

	// B: 16 hypercube-like
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

	// C: 96 even perms of (0, 1/2, φ/2, 1/(2φ)) with even minus signs on the three nonzeros.
	a0 := Real(0)
	a1 := Real(0.5)
	a2 := Real(0.5 * phi)
	a3 := Real(0.5 * inv)
	base := [4]Real{a0, a1, a2, a3}
	eps := evenPerms4()
	for _, p := range eps {
		v := [4]Real{base[p[0]], base[p[1]], base[p[2]], base[p[3]]}
		mask := [4]bool{v[0] != 0, v[1] != 0, v[2] != 0, v[3] != 0}
		for _, vv := range signVariants(v, mask, true) {
			pushUnique(set, &out, vv)
		}
	}

	// sanity
	if len(out) != 120 {
		DebugLog("verts600Unit: expected 120, got %d", len(out))
	}
	return out
}

// ---- 120-cell vertices (unit-radius) ----
// From Wikipedia 120-cell "Unit radius coordinates" left column:
// 8 ({±1,0,0,0})
// 16 ({±1,±1,±1,±1}) / 2
// 96 ([0, ±(φ−1), ±1, ±φ]) / 2
// 32 ([±1, ±1, ±1, ±√5]) / √8
// 32 ([±(φ−1), ±(φ−1), ±(φ−1), ±φ^2]) / √8
// 96 ([0, ±(φ−1), ±φ, ±√5]) / √8
// 96 ([0, ±(φ−2), ±1, ±φ^2]) / √8
// 192 ([±(φ−1), ±1, ±φ, ±2]) / √8
func verts120Unit() []Vector4 {
	phi := (1 + math.Sqrt(5)) / 2
	inv := 1 / phi
	rt5 := math.Sqrt(5)
	rt8 := math.Sqrt(8)

	add := func(vals [][4]Real, out *[]Vector4, set map[[4]int64]struct{}) {
		for _, v := range vals {
			pushUnique(set, out, v)
		}
	}

	set := make(map[[4]int64]struct{}, 640)
	out := make([]Vector4, 0, 600)

	// 8 axes
	for a := 0; a < 4; a++ {
		for s := -1; s <= 1; s += 2 {
			add([][4]Real{{Real(s), 0, 0, 0}}, &out, set)
			out[len(out)-1] = Vector4{0, 0, 0, 0} // fixed below by pushUnique normalization trick, leave as is
		}
	}
	// overwrite properly (axis set above used to reserve; do directly):
	out = out[:0]
	for a := 0; a < 4; a++ {
		for s := -1; s <= 1; s += 2 {
			v := [4]Real{0, 0, 0, 0}
			v[a] = Real(s)
			pushUnique(set, &out, v)
		}
	}

	// 16: (±1/2, ±1/2, ±1/2, ±1/2)
	for sx := -1; sx <= 1; sx += 2 {
		for sy := -1; sy <= 1; sy += 2 {
			for sz := -1; sz <= 1; sz += 2 {
				for sw := -1; sw <= 1; sw += 2 {
					pushUnique(set, &out, [4]Real{
						0.5 * Real(sx), 0.5 * Real(sy), 0.5 * Real(sz), 0.5 * Real(sw),
					})
				}
			}
		}
	}

	// 96: (0, ±(φ−1)/2, ±1/2, ±φ/2), even minus on nonzeros, all perms of positions (not parity-restricted)
	c3 := [3]Real{Real(inv * 0.5), 0.5, Real(phi * 0.5)}
	for zero := 0; zero < 4; zero++ {
		// permutations of the other 3 distinct values = 6
		p3 := [][]int{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}
		for _, p := range p3 {
			val := [4]Real{}
			idx := 0
			for k := 0; k < 4; k++ {
				if k == zero {
					val[k] = 0
				} else {
					val[k] = c3[p[idx]]
					idx++
				}
			}
			mask := [4]bool{val[0] != 0, val[1] != 0, val[2] != 0, val[3] != 0}
			for _, vv := range signVariants(val, mask, true) {
				pushUnique(set, &out, vv)
			}
		}
	}

	// 32: ([±1, ±1, ±1, ±√5]) / √8 with even minus parity
	for pos := 0; pos < 4; pos++ {
		val := [4]Real{Real(1 / rt8), Real(1 / rt8), Real(1 / rt8), Real(1 / rt8)}
		val[pos] = Real(rt5 / rt8)
		mask := [4]bool{true, true, true, true}
		for _, vv := range signVariants(val, mask, true) {
			pushUnique(set, &out, vv)
		}
	}

	// 32: ([±(φ−1), ±(φ−1), ±(φ−1), ±φ^2]) / √8, even minus
	for pos := 0; pos < 4; pos++ {
		val := [4]Real{Real(inv / rt8), Real(inv / rt8), Real(inv / rt8), Real(inv / rt8)}
		val[pos] = Real((phi * phi) / rt8)
		mask := [4]bool{true, true, true, true}
		for _, vv := range signVariants(val, mask, true) {
			pushUnique(set, &out, vv)
		}
	}

	// 96: ([0, ±(φ−1), ±φ, ±√5]) / √8, even minus on the 3 nonzeros
	c3b := [3]Real{Real(inv / rt8), Real(phi / rt8), Real(rt5 / rt8)}
	for zero := 0; zero < 4; zero++ {
		p3 := [][]int{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}
		for _, p := range p3 {
			val := [4]Real{}
			idx := 0
			for k := 0; k < 4; k++ {
				if k == zero {
					val[k] = 0
				} else {
					val[k] = c3b[p[idx]]
					idx++
				}
			}
			mask := [4]bool{val[0] != 0, val[1] != 0, val[2] != 0, val[3] != 0}
			for _, vv := range signVariants(val, mask, true) {
				pushUnique(set, &out, vv)
			}
		}
	}

	// 96: ([0, ±(φ−2), ±1, ±φ^2]) / √8, even minus on the 3 nonzeros
	phim2 := phi - 2 // negative number
	c3c := [3]Real{Real(phim2 / rt8), Real(1 / rt8), Real((phi * phi) / rt8)}
	for zero := 0; zero < 4; zero++ {
		p3 := [][]int{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}
		for _, p := range p3 {
			val := [4]Real{}
			idx := 0
			for k := 0; k < 4; k++ {
				if k == zero {
					val[k] = 0
				} else {
					val[k] = c3c[p[idx]]
					idx++
				}
			}
			mask := [4]bool{val[0] != 0, val[1] != 0, val[2] != 0, val[3] != 0}
			for _, vv := range signVariants(val, mask, true) {
				pushUnique(set, &out, vv)
			}
		}
	}

	// 192: ([±(φ−1), ±1, ±φ, ±2]) / √8, even minus on all 4
	base := [4]Real{Real(inv / rt8), Real(1 / rt8), Real(phi / rt8), Real(2 / rt8)}
	for _, p := range allPerms4Distinct(base) {
		mask := [4]bool{true, true, true, true}
		for _, vv := range signVariants(p, mask, true) {
			pushUnique(set, &out, vv)
		}
	}

	if len(out) != 600 {
		DebugLog("verts120Unit: expected 600, got %d", len(out))
	}
	return out
}
