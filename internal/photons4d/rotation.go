package photons4d

import "math"

// Angles in radians for rotations in coordinate planes.
type Rot4 struct {
	XY, XZ, XW, YZ, YW, ZW Real
}

func rotXY(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][1] = c, -s
	M.M[1][0], M.M[1][1] = s, c
	return M
}
func rotXZ(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][2] = c, -s
	M.M[2][0], M.M[2][2] = s, c
	return M
}
func rotXW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[0][0], M.M[0][3] = c, -s
	M.M[3][0], M.M[3][3] = s, c
	return M
}
func rotYZ(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[1][1], M.M[1][2] = c, -s
	M.M[2][1], M.M[2][2] = s, c
	return M
}
func rotYW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[1][1], M.M[1][3] = c, -s
	M.M[3][1], M.M[3][3] = s, c
	return M
}
func rotZW(a Real) Mat4 {
	c, s := math.Cos(a), math.Sin(a)
	M := I4()
	M.M[2][2], M.M[2][3] = c, -s
	M.M[3][2], M.M[3][3] = s, c
	return M
}

// Compose rotation from angles.
func rotFromAngles(r Rot4) Mat4 {
	R := I4()
	R = rotZW(r.ZW).Mul(R)
	R = rotYW(r.YW).Mul(R)
	R = rotYZ(r.YZ).Mul(R)
	R = rotXW(r.XW).Mul(R)
	R = rotXZ(r.XZ).Mul(R)
	R = rotXY(r.XY).Mul(R)
	return R
}
