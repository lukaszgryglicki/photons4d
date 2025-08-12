package photons4d

import "math"

// Vector4 represents a direction (not a position) in 4D space.
type Vector4 struct {
	X, Y, Z, W Real
}

// Vector functions
func (a Vector4) Add(b Vector4) Vector4 { return Vector4{a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W} }
func (a Vector4) Sub(b Vector4) Vector4 { return Vector4{a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W} }
func (v Vector4) Mul(s Real) Vector4    { return Vector4{v.X * s, v.Y * s, v.Z * s, v.W * s} }

// Dot returns the dot product between two 4D vectors.
func (a Vector4) Dot(b Vector4) Real {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W
}

// Len returns the Euclidean length of the vector.
func (v Vector4) Len() Real { return math.Sqrt(v.Dot(v)) }

// Norm returns a unit-length version of the vector.
func (v Vector4) Norm() Vector4 {
	l := v.Len()
	if l == 0 {
		return v
	}
	return Vector4{v.X / l, v.Y / l, v.Z / l, v.W / l}
}

func (A Mat4) MulVec(v Vector4) Vector4 {
	return Vector4{
		A.M[0][0]*v.X + A.M[0][1]*v.Y + A.M[0][2]*v.Z + A.M[0][3]*v.W,
		A.M[1][0]*v.X + A.M[1][1]*v.Y + A.M[1][2]*v.Z + A.M[1][3]*v.W,
		A.M[2][0]*v.X + A.M[2][1]*v.Y + A.M[2][2]*v.Z + A.M[2][3]*v.W,
		A.M[3][0]*v.X + A.M[3][1]*v.Y + A.M[3][2]*v.Z + A.M[3][3]*v.W,
	}
}
