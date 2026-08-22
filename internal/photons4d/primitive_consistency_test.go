package photons4d

import (
	"math"
	"math/rand"
	"testing"
)

// Universal black-box consistency checks for every primitive:
//   - marching a ray through the object must produce alternating
//     enter (inv=false, D·N<0) / exit (inv=true, D·N>0) crossings,
//   - normals must be unit length,
//   - all crossing points must lie inside the (conservative) AABB,
//   - the total number of crossings must be even (ray starts and ends outside).
type primUT struct {
	name      string
	min, max  Point4
	intersect func(Point4, Vector4) (objectHit, bool)
}

func buildAllPrimitives(t *testing.T, rng *rand.Rand) []primUT {
	t.Helper()
	rot := Rot4{
		XY: rng.Float64() - 0.5, XZ: rng.Float64() - 0.5, XW: rng.Float64() - 0.5,
		YZ: rng.Float64() - 0.5, YW: rng.Float64() - 0.5, ZW: rng.Float64() - 0.5,
	}
	ctr := Point4{rng.Float64()*0.5 - 0.25, rng.Float64()*0.5 - 0.25, rng.Float64()*0.5 - 0.25, rng.Float64()*0.5 - 0.25}
	sc := Vector4{0.8 + rng.Float64()*0.4, 0.8 + rng.Float64()*0.4, 0.8 + rng.Float64()*0.4, 0.8 + rng.Float64()*0.4}
	white := RGB{1, 1, 1}
	half := RGB{0.3, 0.3, 0.3}
	ior := RGB{1.5, 1.5, 1.5}

	var out []primUT
	check := func(name string, err error) {
		if err != nil {
			t.Fatalf("build %s: %v", name, err)
		}
	}
	add := func(name string, min, max Point4, fn func(Point4, Vector4) (objectHit, bool)) {
		out = append(out, primUT{name: name, min: min, max: max, intersect: fn})
	}

	hs, err := NewHyperSphere(ctr, sc, rot, white, half, half, half, ior)
	check("hypersphere", err)
	add("hypersphere", hs.AABBMin, hs.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, hs) })

	c5, err := NewCell5(ctr, sc, rot, white, half, half, half, ior)
	check("cell5", err)
	add("cell5", c5.AABBMin, c5.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell5(O, D, c5) })

	c8, err := NewCell8(ctr, sc, rot, white, half, half, half, ior)
	check("cell8", err)
	add("cell8", c8.AABBMin, c8.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell8(O, D, c8) })

	c16, err := NewCell16(ctr, sc, rot, white, half, half, half, ior)
	check("cell16", err)
	add("cell16", c16.AABBMin, c16.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell16(O, D, c16) })

	c24, err := NewCell24(ctr, sc, rot, white, half, half, half, ior)
	check("cell24", err)
	add("cell24", c24.AABBMin, c24.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell24(O, D, c24) })

	c120, err := NewCell120(ctr, sc, rot, white, half, half, half, ior)
	check("cell120", err)
	add("cell120", c120.AABBMin, c120.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &c120.cellPoly) })

	c600, err := NewCell600(ctr, sc, rot, white, half, half, half, ior)
	check("cell600", err)
	add("cell600", c600.AABBMin, c600.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &c600.cellPoly) })

	st, err := NewStarCell("cell16", ctr, sc, defaultStarCellCoreRadius, defaultStarCellSpikeLength, defaultStarCellSharpness, rot, white, half, half, half, ior)
	check("starcell", err)
	add("starcell", st.AABBMin, st.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayStarCell(O, D, st) })

	sp, err := NewSpherinder(ctr, sc, rot, white, half, half, half, ior)
	check("spherinder", err)
	add("spherinder", sp.AABBMin, sp.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRaySpherinder(O, D, sp) })

	hc, err := NewHyperCone(ctr, sc, rot, white, half, half, half, ior)
	check("hypercone", err)
	add("hypercone", hc.AABBMin, hc.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperCone(O, D, hc) })

	cap4, err := NewHyperCapsule(ctr, sc, 0.7, rot, white, half, half, half, ior)
	check("hypercapsule", err)
	add("hypercapsule", cap4.AABBMin, cap4.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperCapsule(O, D, cap4) })

	sto, err := NewSpheritorus(ctr, sc, 0.8, 0.3, rot, white, half, half, half, ior)
	check("spheritorus", err)
	add("spheritorus", sto.AABBMin, sto.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRaySpheritorus(O, D, sto) })

	dt, err := NewDuotorus(ctr, sc, 0.8, 0.7, 0.25, rot, white, half, half, half, ior)
	check("duotorus", err)
	add("duotorus", dt.AABBMin, dt.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayDuotorus(O, D, dt) })

	dc, err := NewDuocylinder(ctr, sc, rot, white, half, half, half, ior)
	check("duocylinder", err)
	add("duocylinder", dc.AABBMin, dc.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayDuocylinder(O, D, dc) })

	ts, err := NewTorisphere(ctr, sc, 0.8, 0.3, rot, white, half, half, half, ior)
	check("torisphere", err)
	add("torisphere", ts.AABBMin, ts.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayTorisphere(O, D, ts) })

	sq, err := NewSuperquadric(ctr, sc, 3.5, rot, white, half, half, half, ior)
	check("superquadric", err)
	add("superquadric", sq.AABBMin, sq.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRaySuperquadric(O, D, sq) })

	hf, err := NewHyperFrustum(ctr, sc, 0.5, 1.0, rot, white, half, half, half, ior)
	check("hyperfrustum", err)
	add("hyperfrustum", hf.AABBMin, hf.AABBMax, func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperFrustum(O, D, hf) })

	return out
}

func TestPrimitiveCrossingConsistency(t *testing.T) {
	rng := rand.New(rand.NewSource(12345))
	const raysPerObj = 400
	const marchEps = 1e-7
	const maxCrossings = 64

	for trial := 0; trial < 3; trial++ {
		prims := buildAllPrimitives(t, rng)
		for _, p := range prims {
			ctr := Point4{
				(p.min.X + p.max.X) / 2, (p.min.Y + p.max.Y) / 2,
				(p.min.Z + p.max.Z) / 2, (p.min.W + p.max.W) / 2,
			}
			diag := math.Sqrt((p.max.X-p.min.X)*(p.max.X-p.min.X) + (p.max.Y-p.min.Y)*(p.max.Y-p.min.Y) +
				(p.max.Z-p.min.Z)*(p.max.Z-p.min.Z) + (p.max.W-p.min.W)*(p.max.W-p.min.W))
			// AABB tolerance: scaled to object size (root-finding on curved
			// surfaces is only accurate to ~1e-6 of the object scale).
			tol := 1e-5 * (1 + diag)

			enterBad, exitBad, oddBad, unitBad, aabbBad, hits := 0, 0, 0, 0, 0, 0
			for r := 0; r < raysPerObj; r++ {
				// Random start outside the AABB, aimed at a random point inside it.
				dir := unitS3(rng)
				start := ctr.Add(dir.Mul(-(diag + 1.0)))
				target := Point4{
					p.min.X + rng.Float64()*(p.max.X-p.min.X),
					p.min.Y + rng.Float64()*(p.max.Y-p.min.Y),
					p.min.Z + rng.Float64()*(p.max.Z-p.min.Z),
					p.min.W + rng.Float64()*(p.max.W-p.min.W),
				}
				D := Vector4{target.X - start.X, target.Y - start.Y, target.Z - start.Z, target.W - start.W}.Norm()

				O := start
				inside := false
				crossings := 0
				for c := 0; c < maxCrossings; c++ {
					h, ok := p.intersect(O, D)
					if !ok {
						break
					}
					crossings++
					hits++
					P := O.Add(D.Mul(h.t))

					if l := h.Nw.Len(); math.Abs(l-1) > 1e-6 {
						unitBad++
					}
					if P.X < p.min.X-tol || P.X > p.max.X+tol ||
						P.Y < p.min.Y-tol || P.Y > p.max.Y+tol ||
						P.Z < p.min.Z-tol || P.Z > p.max.Z+tol ||
						P.W < p.min.W-tol || P.W > p.max.W+tol {
						aabbBad++
					}
					dn := D.Dot(h.Nw)
					if !inside {
						// entering: must not be flagged as inverse, normal opposes D
						if h.inv || dn > 1e-9 {
							enterBad++
						}
					} else {
						// exiting: inv flag set, normal along D
						if !h.inv || dn < -1e-9 {
							exitBad++
						}
					}
					inside = !inside
					O = P.Add(D.Mul(marchEps))
				}
				if inside {
					oddBad++
				}
			}
			if hits == 0 {
				t.Errorf("%s: no hits at all out of %d rays (intersection or AABB broken)", p.name, raysPerObj)
			}
			// Allow a tiny fraction of tangent/grazing numerical misclassifications.
			frac := func(n int) float64 { return float64(n) / float64(hits+1) }
			if frac(enterBad) > 0.01 || frac(exitBad) > 0.01 || frac(unitBad) > 0.001 || frac(aabbBad) > 0.001 || float64(oddBad)/raysPerObj > 0.01 {
				t.Errorf("%s: inconsistencies: enterBad=%d exitBad=%d unitBad=%d aabbBad=%d oddParity=%d (of %d hits, %d rays)",
					p.name, enterBad, exitBad, unitBad, aabbBad, oddBad, hits, raysPerObj)
			}
		}
	}
}
