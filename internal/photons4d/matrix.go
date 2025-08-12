package photons4d

// 4×4 matrix (row-major)
type Mat4 struct {
	M [4][4]Real
}

func I4() Mat4 {
	return Mat4{M: [4][4]Real{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}}
}

func (A Mat4) Mul(B Mat4) Mat4 {
	var R Mat4
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				sum += A.M[r][k] * B.M[k][c]
			}
			R.M[r][c] = sum
		}
	}
	return R
}

func (A Mat4) Transpose() Mat4 {
	var R Mat4
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			R.M[r][c] = A.M[c][r]
		}
	}
	return R
}
