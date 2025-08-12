package photons4d

type rayRecips struct {
	invX, invY, invZ, invW Real
	parX, parY, parZ, parW bool // parallel flags (|D| < eps)
}

func rayAABB(O Point4, minP, maxP Point4, rr rayRecips) (bool, Real) {
	tmin, tmax := -1e300, 1e300

	// X
	if !rr.parX {
		t1 := (minP.X - O.X) * rr.invX
		t2 := (maxP.X - O.X) * rr.invX
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
	} else if O.X < minP.X || O.X > maxP.X {
		return false, 0
	}

	// Y
	if !rr.parY {
		t1 := (minP.Y - O.Y) * rr.invY
		t2 := (maxP.Y - O.Y) * rr.invY
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
	} else if O.Y < minP.Y || O.Y > maxP.Y {
		return false, 0
	}

	// Z
	if !rr.parZ {
		t1 := (minP.Z - O.Z) * rr.invZ
		t2 := (maxP.Z - O.Z) * rr.invZ
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
	} else if O.Z < minP.Z || O.Z > maxP.Z {
		return false, 0
	}

	// W
	if !rr.parW {
		t1 := (minP.W - O.W) * rr.invW
		t2 := (maxP.W - O.W) * rr.invW
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
	} else if O.W < minP.W || O.W > maxP.W {
		return false, 0
	}

	if tmax < 0 || tmin > tmax {
		return false, 0
	}
	return true, tmin
}
