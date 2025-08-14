package photons4d

import (
	"math"
	"math/rand"
)

// sampleDiffuseDir returns a cosine-weighted unit direction on the S^3 hemisphere around unit N.
// Construction: pick a point uniformly in the unit 3-ball for the tangent part (U,V,W),
// then set the normal component to sqrt(1 - r^2). This generalizes the familiar 3D disk mapping.
func sampleDiffuseDir(N Vector4, rng *rand.Rand) Vector4 {
	// Build orthonormal basis of the tangent 3-space.
	U, V, W := orthonormal3(N)

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
