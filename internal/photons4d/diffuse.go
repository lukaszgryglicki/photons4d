package photons4d

import (
	"math"
	"math/rand"
)

// Fast deterministic orthonormal basis of the 3D tangent space to unit a in R^4.
// This avoids time.Now(), RNG construction and retries in the diffuse hot path.
func orthonormal3Fast(a Vector4) (u, v, w Vector4) {
	const eps = 1e-12

	proj := func(x Vector4) Vector4 {
		return x.Sub(a.Mul(x.Dot(a)))
	}

	helpers := [4]Vector4{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}

	absComp := [4]Real{math.Abs(a.X), math.Abs(a.Y), math.Abs(a.Z), math.Abs(a.W)}
	order := [4]int{0, 1, 2, 3}
	// Least-aligned helper first => best-conditioned first tangent vector.
	for i := 1; i < 4; i++ {
		x := order[i]
		j := i - 1
		for ; j >= 0 && absComp[order[j]] > absComp[x]; j-- {
			order[j+1] = order[j]
		}
		order[j+1] = x
	}

	var basis [3]Vector4
	n := 0
	for _, idx := range order {
		q := proj(helpers[idx])
		for i := 0; i < n; i++ {
			q = q.Sub(basis[i].Mul(q.Dot(basis[i])))
		}
		// Reproject once more to kill accumulated numerical drift back onto a.
		q = proj(q)
		l := q.Len()
		if l <= eps {
			continue
		}
		basis[n] = q.Mul(1 / l)
		n++
		if n == 3 {
			return basis[0], basis[1], basis[2]
		}
	}

	// Extremely unlikely fallback.
	return orthonormal3(a)
}

// sampleDiffuseDir returns a cosine-weighted unit direction on the S^3 hemisphere around unit N.
// Construction: pick a point uniformly in the unit 3-ball for the tangent part (U,V,W),
// then set the normal component to sqrt(1 - r^2). This generalizes the familiar 3D disk mapping.
func sampleDiffuseDir(N Vector4, rng *rand.Rand) Vector4 {
	// Build orthonormal basis of the tangent 3-space.
	// U, V, W := orthonormal3(N)
	U, V, W := orthonormal3Fast(N)

	// Uniform in unit 3-ball: radius r ~ U^(1/3), orientation from S^2.
	ux, uy, uz := sampleS2(rng) // unit on S^2
	r := math.Cbrt(rng.Float64())
	tx, ty, tz := r*ux, r*uy, r*uz

	// Normal component (hemisphere side aligned with +N).
	nn2 := 1 - (tx*tx + ty*ty + tz*tz)
	if nn2 < 0 {
		nn2 = 0
	}
	nn := math.Sqrt(nn2)

	// Assemble world-space direction.
	dir := U.Mul(tx).Add(V.Mul(ty)).Add(W.Mul(tz)).Add(N.Mul(nn))
	// Numerical safety: normalize.
	l2 := dir.Dot(dir)
	if l2 > 0 {
		dir = dir.Mul(1 / math.Sqrt(l2))
	}
	return dir
}
