package photons4d

import (
	"fmt"
	"math"
	"strings"
)

const (
	defaultStarCellCoreRadius  Real = 0.45
	defaultStarCellSpikeLength Real = 0.55
	defaultStarCellSharpness   Real = 6.0
)

// StarCell is a non-convex, star-shaped 4D solid generated from the dual-cell
// directions of a regular polytope family. It is not an exact Schläfli-Hess
// star polychoron; instead it is a robust solid “stellation” that preserves a
// single inside/outside region, so reflection/refraction semantics remain
// compatible with the rest of the renderer.
//
// In normalized local space q, the boundary radius is:
//
//	rho(u) = CoreRadius + SpikeLength * max_i(max(0, u·dir_i))^Sharpness
//
// where u = q/|q| and dir_i are the family-specific dual directions.
//
// The normalized local shape is then scaled per axis by Scale, rotated by R,
// and translated by Center.
type StarCell struct {
	Kind        string
	Center      Point4
	Scale       Vector4
	CoreRadius  Real
	SpikeLength Real
	Sharpness   Real
	OuterRadius Real
	R           Mat4
	RT          Mat4
	Dirs        []Vector4
	slope       Real // largest change of the boundary radius per radian of direction

	Color   RGB
	Diffuse RGB
	Reflect RGB
	Refract RGB
	IOR     RGB

	AABBMin Point4
	AABBMax Point4

	// material caches
	refl     [3]Real
	refr     [3]Real
	colorArr [3]Real
	diff     [3]Real
	iorArr   [3]Real
	iorInv   [3]Real
	pAbs     [3]Real
	f0       [3]Real
}

func normalizeStarCellKind(kind string) string {
	k := strings.ToLower(strings.TrimSpace(kind))
	k = strings.NewReplacer("_", "-", " ", "-").Replace(k)
	for strings.Contains(k, "--") {
		k = strings.ReplaceAll(k, "--", "-")
	}
	return strings.Trim(k, "-")
}

func canonicalStarCellKind(kind string) (string, error) {
	switch normalizeStarCellKind(kind) {
	case "5", "5-cell", "cell5", "cell-5", "star5", "star-5", "star-cell-5":
		return "star-cell-5", nil
	case "8", "8-cell", "cell8", "cell-8", "star8", "star-8", "star-cell-8":
		return "star-cell-8", nil
	case "16", "16-cell", "cell16", "cell-16", "star16", "star-16", "star-cell-16":
		return "star-cell-16", nil
	case "24", "24-cell", "cell24", "cell-24", "star24", "star-24", "star-cell-24":
		return "star-cell-24", nil
	case "120", "120-cell", "cell120", "cell-120", "star120", "star-120", "star-cell-120":
		return "star-cell-120", nil
	case "600", "600-cell", "cell600", "cell-600", "star600", "star-600", "star-cell-600":
		return "star-cell-600", nil
	default:
		return "", fmt.Errorf("unsupported star cell kind %q (expected 5, 8, 16, 24, 120 or 600 family)", kind)
	}
}

func canonical5StarDirs() []Vector4 {
	base := canonicalCell5()
	out := make([]Vector4, 0, len(base))
	for _, v := range base {
		out = append(out, v.Norm())
	}
	return out
}

func canonical8StarDirs() []Vector4 {
	return []Vector4{
		{+1, 0, 0, 0},
		{-1, 0, 0, 0},
		{0, +1, 0, 0},
		{0, -1, 0, 0},
		{0, 0, +1, 0},
		{0, 0, -1, 0},
		{0, 0, 0, +1},
		{0, 0, 0, -1},
	}
}

func canonical16StarDirs() []Vector4 {
	out := make([]Vector4, 0, 16)
	k := 0.5
	for sx := -1; sx <= 1; sx += 2 {
		for sy := -1; sy <= 1; sy += 2 {
			for sz := -1; sz <= 1; sz += 2 {
				for sw := -1; sw <= 1; sw += 2 {
					out = append(out, Vector4{k * Real(sx), k * Real(sy), k * Real(sz), k * Real(sw)}.Norm())
				}
			}
		}
	}
	return out
}

func canonical24StarDirs() []Vector4 {
	base := canonical24Verts()
	out := make([]Vector4, 0, len(base))
	for _, v := range base {
		out = append(out, v.Norm())
	}
	return out
}

func starCellDirections(kind string) (string, []Vector4, error) {
	canon, err := canonicalStarCellKind(kind)
	if err != nil {
		return "", nil, err
	}
	switch canon {
	case "star-cell-5":
		return canon, canonical5StarDirs(), nil
	case "star-cell-8":
		return canon, canonical8StarDirs(), nil
	case "star-cell-16":
		return canon, canonical16StarDirs(), nil
	case "star-cell-24":
		return canon, canonical24StarDirs(), nil
	case "star-cell-120":
		return canon, verts600Unit(), nil // dual of the 120-cell
	case "star-cell-600":
		return canon, verts120Unit(), nil // dual of the 600-cell
	default:
		return "", nil, fmt.Errorf("internal: unhandled canonical star cell kind %q", canon)
	}
}

// starCellHoles returns the directions that lie farthest from the family
// directions: the vertices of the dual polytope.
func starCellHoles(canon string, dirs []Vector4) []Vector4 {
	switch canon {
	case "star-cell-8":
		return canonical16StarDirs()
	case "star-cell-16":
		return canonical8StarDirs()
	case "star-cell-24":
		seeds := dual24Dirs()
		out := make([]Vector4, 0, len(seeds))
		for _, v := range seeds {
			out = append(out, v.Norm())
		}
		return out
	case "star-cell-120":
		return verts120Unit()
	case "star-cell-600":
		return verts600Unit()
	default:
		out := make([]Vector4, 0, len(dirs))
		for _, v := range dirs {
			out = append(out, v.Mul(-1))
		}
		return out
	}
}

// starCellSlope bounds |d rho / d angle| of rho = core + spike*c^sharpness
// over the cosines c in [minCos, 1] that the nearest family direction can have.
func starCellSlope(spike, sharpness, minCos Real) Real {
	c := minCos
	if sharpness > 1 {
		if top := math.Sqrt((sharpness - 1) / sharpness); top > c {
			c = top
		}
	}
	return spike * sharpness * math.Pow(c, sharpness-1) * math.Sqrt(1-c*c)
}

func NewStarCell(kind string, center Point4, scale Vector4, coreRadius, spikeLength, sharpness Real, angles Rot4, color, diffuse, reflectivity, refractivity, ior RGB) (*StarCell, error) {
	if !(scale.X > 0 && scale.Y > 0 && scale.Z > 0 && scale.W > 0) {
		return nil, fmt.Errorf("star cell scale must be > 0 on all axes, got %+v", scale)
	}
	if coreRadius < 0 {
		return nil, fmt.Errorf("star cell coreRadius must be >= 0, got %.6g", coreRadius)
	}
	if spikeLength < 0 {
		return nil, fmt.Errorf("star cell spikeLength must be >= 0, got %.6g", spikeLength)
	}
	if coreRadius+spikeLength <= 0 {
		return nil, fmt.Errorf("star cell coreRadius + spikeLength must be > 0, got %.6g", coreRadius+spikeLength)
	}
	if !(sharpness > 0) {
		return nil, fmt.Errorf("star cell sharpness must be > 0, got %.6g", sharpness)
	}

	in01 := func(x Real) bool { return x >= 0 && x <= 1 }
	for _, c := range []struct {
		n       string
		r, t, d Real
	}{
		{"R", reflectivity.R, refractivity.R, diffuse.R},
		{"G", reflectivity.G, refractivity.G, diffuse.G},
		{"B", reflectivity.B, refractivity.B, diffuse.B},
	} {
		if !in01(c.r) || !in01(c.t) || !in01(c.d) {
			return nil, fmt.Errorf("reflect/refract/diffuse must be in [0,1]; channel %s got reflect=%.6g refract=%.6g diffuse=%.6g", c.n, c.r, c.t, c.d)
		}
		if c.r+c.t+c.d > 1+1e-12 {
			return nil, fmt.Errorf("per-channel reflect+refract+diffuse must be ≤1; channel %s got %.6g", c.n, c.r+c.t+c.d)
		}
	}
	if !(ior.R > 0 && ior.G > 0 && ior.B > 0) {
		return nil, fmt.Errorf("IOR must be > 0 per channel; got %+v", ior)
	}

	canon, dirs, err := starCellDirections(kind)
	if err != nil {
		return nil, err
	}

	R := rotFromAngles(angles)
	s := &StarCell{
		Kind:        canon,
		Center:      center,
		Scale:       scale,
		CoreRadius:  coreRadius,
		SpikeLength: spikeLength,
		Sharpness:   sharpness,
		OuterRadius: coreRadius + spikeLength,
		R:           R,
		RT:          R.Transpose(),
		Dirs:        dirs,
		Color:       color,
		Diffuse:     diffuse,
		Reflect:     reflectivity,
		Refract:     refractivity,
		IOR:         ior,
	}

	minCos := Real(1)
	for _, h := range starCellHoles(canon, dirs) {
		best := Real(-1)
		for _, d := range dirs {
			if c := d.Dot(h); c > best {
				best = c
			}
		}
		if best < minCos {
			minCos = best
		}
	}
	s.slope = starCellSlope(spikeLength, sharpness, minCos)

	// Conservative AABB from the rotated local box [-outer*scale, +outer*scale]^4.
	abs := func(x Real) Real {
		if x < 0 {
			return -x
		}
		return x
	}
	hx := s.OuterRadius * scale.X
	hy := s.OuterRadius * scale.Y
	hz := s.OuterRadius * scale.Z
	hw := s.OuterRadius * scale.W
	extent := func(row int) (minV, maxV Real) {
		off := abs(R.M[row][0])*hx + abs(R.M[row][1])*hy + abs(R.M[row][2])*hz + abs(R.M[row][3])*hw
		var c Real
		switch row {
		case 0:
			c = center.X
		case 1:
			c = center.Y
		case 2:
			c = center.Z
		default:
			c = center.W
		}
		return c - off, c + off
	}
	minX, maxX := extent(0)
	minY, maxY := extent(1)
	minZ, maxZ := extent(2)
	minW, maxW := extent(3)
	s.AABBMin = Point4{minX, minY, minZ, minW}
	s.AABBMax = Point4{maxX, maxY, maxZ, maxW}

	s.refl = [3]Real{reflectivity.R, reflectivity.G, reflectivity.B}
	s.refr = [3]Real{refractivity.R, refractivity.G, refractivity.B}
	s.colorArr = [3]Real{color.R, color.G, color.B}
	s.diff = [3]Real{diffuse.R, diffuse.G, diffuse.B}
	s.iorArr = [3]Real{ior.R, ior.G, ior.B}
	for i := 0; i < 3; i++ {
		p := 1 - s.refl[i] - s.refr[i] - s.diff[i]
		if p < 0 {
			p = 0
		}
		s.pAbs[i] = p
		s.iorInv[i] = 1 / s.iorArr[i]
		n := s.iorArr[i]
		r0 := (n - 1) / (n + 1)
		s.f0[i] = r0 * r0
	}

	DebugLog("Created %s: %+v", s.Kind, s)
	return s, nil
}

func (s *StarCell) evalNormalizedLocal(q Vector4) (f Real, bestDir Vector4, bestCos Real) {
	const eps = 1e-30
	r2 := q.Dot(q)
	if r2 <= eps {
		return -s.CoreRadius, Vector4{1, 0, 0, 0}, 0
	}
	r := math.Sqrt(r2)
	invR := 1 / r
	best := math.Inf(-1)
	for _, d := range s.Dirs {
		sc := d.Dot(q)
		if sc > best {
			best = sc
			bestDir = d
		}
	}
	if best < 0 {
		best = 0
		bestCos = 0
	} else {
		bestCos = best * invR
		if bestCos > 1 {
			bestCos = 1
		}
	}
	rho := s.CoreRadius
	if bestCos > 0 {
		rho += s.SpikeLength * math.Pow(bestCos, s.Sharpness)
	}
	return r - rho, bestDir, bestCos
}

func (s *StarCell) localGradient(q Vector4, bestDir Vector4, bestCos Real) Vector4 {
	const eps = 1e-30
	r2 := q.Dot(q)
	if r2 <= eps {
		return Vector4{1, 0, 0, 0}
	}
	r := math.Sqrt(r2)
	u := q.Mul(1 / r)
	if bestCos <= 0 {
		return u
	}
	coeff := s.SpikeLength * s.Sharpness * math.Pow(bestCos, s.Sharpness-1) / r
	return u.Sub(bestDir.Sub(u.Mul(bestCos)).Mul(coeff))
}

// boundDistance returns a lower bound of the distance from q to the surface
// in the normalized local space, negative inside. f = |q| - rho changes by at
// most sqrt(1 + (slope/x)²) per unit of length at the radius x, and the ball
// of the radius |q|/2 around q stays beyond the radius |q|/2.
func (s *StarCell) boundDistance(q Vector4) Real {
	f, _, _ := s.evalNormalizedLocal(q)
	r := q.Len()
	h := 0.5 * r
	b := Real(0)
	if h > 0 {
		b = math.Abs(f) * h / math.Hypot(h, s.slope)
		if b > h {
			b = h
		}
	}
	if g := s.CoreRadius - r; g > b {
		b = g
	}
	if g := r - s.OuterRadius; g > b {
		b = g
	}
	if f < 0 {
		return -b
	}
	return b
}

func intersectRayStarCell(O Point4, D Vector4, s *StarCell) (hit objectHit, ok bool) {
	Op := Vector4{O.X - s.Center.X, O.Y - s.Center.Y, O.Z - s.Center.Z, O.W - s.Center.W}
	Ol := s.RT.MulVec(Op)
	Dl := s.RT.MulVec(D)
	Oq := Vector4{Ol.X / s.Scale.X, Ol.Y / s.Scale.Y, Ol.Z / s.Scale.Z, Ol.W / s.Scale.W}
	Dq := Vector4{Dl.X / s.Scale.X, Dl.Y / s.Scale.Y, Dl.Z / s.Scale.Z, Dl.W / s.Scale.W}

	// Part of the ray inside the ball of the outer radius.
	a := Dq.Dot(Dq)
	b := Oq.Dot(Dq)
	c := Oq.Dot(Oq) - s.OuterRadius*s.OuterRadius
	disc := b*b - a*c
	if !(a > 0) || disc < 0 {
		return objectHit{}, false
	}
	sq := math.Sqrt(disc)

	dist := func(t Real) rayDist {
		return rayDist{d: s.boundDistance(Oq.Add(Dq.Mul(t)))}
	}
	tHit, inv, ok := traceSolid(dist, math.Sqrt(a), (-b-sq)/a, (-b+sq)/a)
	if !ok {
		return objectHit{}, false
	}

	qHit := Oq.Add(Dq.Mul(tHit))
	_, bestDir, bestCos := s.evalNormalizedLocal(qHit)
	gradQ := s.localGradient(qHit, bestDir, bestCos)
	gradScaled := Vector4{gradQ.X / s.Scale.X, gradQ.Y / s.Scale.Y, gradQ.Z / s.Scale.Z, gradQ.W / s.Scale.W}
	Nw := s.R.MulVec(gradScaled).Norm()
	if Nw.Len() == 0 {
		return objectHit{}, false
	}
	return objectHit{t: tHit, Nw: Nw, mat: s, inv: inv}, true
}

func (s *StarCell) PAbsCh(i int) Real   { return s.pAbs[i] }
func (s *StarCell) DiffCh(i int) Real   { return s.diff[i] }
func (s *StarCell) ColorCh(i int) Real  { return s.colorArr[i] }
func (s *StarCell) F0Ch(i int) Real     { return s.f0[i] }
func (s *StarCell) ReflCh(i int) Real   { return s.refl[i] }
func (s *StarCell) RefrCh(i int) Real   { return s.refr[i] }
func (s *StarCell) IORCh(i int) Real    { return s.iorArr[i] }
func (s *StarCell) IORInvCh(i int) Real { return s.iorInv[i] }
