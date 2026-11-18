package photons4d

import (
	"fmt"
	"math"
)

// ballsSupport returns the largest a·q over the part common to the balls of
// the radius rho around the centres. The point that reaches it lies on the
// spheres of some of the balls, which meet in a round sphere around the
// point of the flat through their centres that is equally far from them; it
// is the point of that sphere furthest along a.
func ballsSupport(centers []Vector4, rho Real, a Vector4) Real {
	best := math.Inf(-1)
	for mask := 1; mask < 1<<len(centers); mask++ {
		var first Vector4
		var e [4]Vector4
		m := -1
		for i, c := range centers {
			if mask>>i&1 == 0 {
				continue
			}
			if m < 0 {
				first = c
			} else if m < len(e) {
				e[m] = c.Sub(first)
			}
			m++
		}
		if m > len(e) {
			continue
		}
		// The middle is first + Σ l_i e_i with 2·G·l = diag(G), and a leaves
		// the flat by a - Σ k_i e_i with G·k = (a·e_i), G being the Gram matrix.
		var g [4][6]Real
		for i := 0; i < m; i++ {
			for j := 0; j < m; j++ {
				g[i][j] = e[i].Dot(e[j])
			}
			g[i][m] = g[i][i] / 2
			g[i][m+1] = a.Dot(e[i])
		}
		flat := true
		for c := 0; c < m && flat; c++ {
			p := c
			for r := c + 1; r < m; r++ {
				if math.Abs(g[r][c]) > math.Abs(g[p][c]) {
					p = r
				}
			}
			g[c], g[p] = g[p], g[c]
			if math.Abs(g[c][c]) < 1e-12 {
				flat = false
				break
			}
			for r := 0; r < m; r++ {
				if r == c {
					continue
				}
				f := g[r][c] / g[c][c]
				for j := c; j < m+2; j++ {
					g[r][j] -= f * g[c][j]
				}
			}
		}
		if !flat {
			continue
		}
		mid, out := first, a
		for i := 0; i < m; i++ {
			mid = mid.Add(e[i].Mul(g[i][m] / g[i][i]))
			out = out.Sub(e[i].Mul(g[i][m+1] / g[i][i]))
		}
		r2 := rho*rho - mid.Sub(first).Dot(mid.Sub(first))
		if r2 < 0 {
			continue
		}
		p := mid
		if l := out.Len(); l > 1e-12*a.Len() {
			p = mid.Add(out.Mul(math.Sqrt(r2) / l))
		}
		inside := true
		for _, c := range centers {
			if p.Sub(c).Len() > rho*(1+1e-9) {
				inside = false
				break
			}
		}
		if inside {
			best = math.Max(best, a.Dot(p))
		}
	}
	return best
}

// NewLens builds a lens with two spherical faces: the part common to two
// balls whose centres lie on the local W axis. Before Scale its rim is the
// sphere of the radius 1 in the hyperplane w = 0 and its faces reach
// thickness/2 along W on both sides; the thickness 2 gives the unit ball.
func NewLens(center Point4, scale Vector4, thickness Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	t := thickness / 2
	if !(t > 0 && t <= 1) {
		return nil, fmt.Errorf("lens thickness must be in (0,2], got %.6g", thickness)
	}
	// The balls have the radius t + k around (0,0,0,±k): with the rim at 1,
	// |q - c|² - radius² is |q|² - 1 ± 2·k·w.
	k := (1 - t*t) / (2 * t)
	f := func(q Vector4) (Real, Vector4) {
		side := Real(1)
		if q.W < 0 {
			side = -1
		}
		return q.Dot(q) - 1 + 2*k*side*q.W, Vector4{2 * q.X, 2 * q.Y, 2 * q.Z, 2 * (q.W + side*k)}
	}
	centers := []Vector4{{0, 0, 0, -k}, {0, 0, 0, k}}
	support := func(a Vector4) Real { return ballsSupport(centers, t+k, a) }
	return newConvexSolid("lens", f, support, 1, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
}

// NewReuleauxSolid builds a body of the Reuleaux kind, bounded by pieces of
// spheres around its own corners:
//
//   - "rotor": the Reuleaux triangle of the width 2 turned around its axis
//     of symmetry (the local W axis) in XYZ. It has the same width 2 in every
//     direction, as the unit ball has, so it rolls between two parallel
//     hyperplanes without moving them apart;
//   - "simplex": the part common to five balls around the corners of a
//     regular 5-cell (on the unit sphere) that reach the other four corners.
func NewReuleauxSolid(kind string, center Point4, scale Vector4, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*ConvexSolid, error) {
	switch shapeKey(kind) {
	case "", "rotor":
		top, base := 2/math.Sqrt(3), -1/math.Sqrt(3)
		f := func(q Vector4) (Real, Vector4) {
			// the ball around the corner on the axis
			u := Vector4{q.X, q.Y, q.Z, q.W - top}
			v, g := u.Len(), Vector4{}
			if v > 0 {
				g = u.Mul(1 / v)
			}
			// the balls around the corners of the base, of which the far one counts
			lat := math.Sqrt(q.X*q.X + q.Y*q.Y + q.Z*q.Z)
			if side := math.Hypot(lat+1, q.W-base); side > v {
				v, g = side, Vector4{0, 0, 0, (q.W - base) / side}
				if lat > 0 {
					k := (lat + 1) / (side * lat)
					g.X, g.Y, g.Z = k*q.X, k*q.Y, k*q.Z
				}
			}
			return v - 2, g
		}
		corners := []Vector4{{0, 0, 0, top}, {-1, 0, 0, base}, {1, 0, 0, base}}
		support := func(a Vector4) Real {
			return ballsSupport(corners, 2, Vector4{math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z), 0, 0, a.W})
		}
		return newConvexSolid("reuleaux rotor", f, support, top, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
	case "simplex":
		canon := canonicalCell5()
		corners := make([]Vector4, len(canon))
		for i, v := range canon {
			corners[i] = v.Norm()
		}
		edge := corners[0].Sub(corners[1]).Len()
		f := func(q Vector4) (Real, Vector4) {
			v, g := math.Inf(-1), Vector4{}
			for _, c := range corners {
				if d := q.Sub(c); d.Len() > v {
					v, g = d.Len(), d
				}
			}
			if v > 0 {
				g = g.Mul(1 / v)
			}
			return v - edge, g
		}
		support := func(a Vector4) Real { return ballsSupport(corners, edge, a) }
		return newConvexSolid("reuleaux simplex", f, support, 1, center, scale, angles, color, diffuse, reflectivity, refractivity, ior)
	}
	return nil, fmt.Errorf("unknown reuleaux solid %q (rotor, simplex)", kind)
}
