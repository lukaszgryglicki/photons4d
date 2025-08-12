package photons4d

// Point4 represents a point in 4-dimensional space.
type Point4 struct {
	X, Y, Z, W Real
}

// Add lets you translate a Point4 by a Vector4.
func (p Point4) Add(v Vector4) Point4 {
	return Point4{p.X + v.X, p.Y + v.Y, p.Z + v.Z, p.W + v.W}
}
