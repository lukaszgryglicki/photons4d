package photons4d

import (
	"math"
	"math/rand"
	"os"
	"testing"
)

// copy of nearest (exact part) with the window factor, the step of the search and the shortcut as parameters
func zzNearest(k *KleinSurface, q Vector4, factor, du, shortcut Real) Real {
	const most = 4096
	rho := math.Hypot(q.X, q.Y)
	phi := math.Atan2(q.Y, q.X)
	g0, g1, g2, f0 := k.at(q, phi)
	width := Real(math.Pi)
	if d0 := math.Sqrt(g0); rho > d0 {
		width = math.Min(width, factor*math.Asin(d0/rho)+1e-9)
	}
	if width <= shortcut {
		g, _ := k.descendFrom(q, phi, width, g0, g1, g2, f0)
		return math.Sqrt(g)
	}
	n := 2 * int(math.Ceil(width/du))
	if n < 4 {
		n = 4
	}
	if n > most {
		n = most
	}
	h := 2 * width / Real(n)
	var gs [most + 1]Real
	su, cu := math.Sincos(phi - width)
	sh, ch := math.Sincos(k.half * (phi - width))
	sd, cd := math.Sincos(h)
	shd, chd := math.Sincos(k.half * h)
	for i := 0; i <= n; i++ {
		gs[i] = k.sq(q, su, cu, sh, ch)
		su, cu = su*cd+cu*sd, cu*cd-su*sd
		sh, ch = sh*chd+ch*shd, ch*chd-sh*shd
	}
	slack := k.sweep * h / 2
	var cand [16]int
	nc := 0
	for i := 0; i <= n; i++ {
		if (i > 0 && gs[i] > gs[i-1]) || (i < n && gs[i] >= gs[i+1]) {
			continue
		}
		j := nc
		if nc < len(cand) {
			nc++
		} else if j--; gs[i] >= gs[cand[j]] {
			continue
		}
		for ; j > 0 && gs[cand[j-1]] > gs[i]; j-- {
			cand[j] = cand[j-1]
		}
		cand[j] = i
	}
	dist := math.Inf(1)
	for _, i := range cand[:nc] {
		if math.Sqrt(gs[i])-2*slack >= dist {
			break
		}
		g, _ := k.descend(q, phi-width+Real(i)*h, h)
		if d := math.Sqrt(g); d < dist {
			dist = d
		}
	}
	return dist
}

func TestZZWindow(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	one := Vector4{1, 1, 1, 1}
	var ks []*KleinSurface
	add := func(k *KleinSurface, err error) {
		if err != nil {
			t.Fatal(err)
		}
		ks = append(ks, k)
	}
	for _, tw := range []int{0, 1, 2, 3, 4, 6, 9, 16} {
		add(NewKleinBottle(Point4{}, one, 1, 0.4, 0.1, tw, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
		add(NewKleinBottle(Point4{}, one, 0.7, 0.65, 0.1, tw, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
		add(NewMobiusBand(Point4{}, one, 1, 0.45, 0.1, tw, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
		add(NewMobiusBand(Point4{}, one, 0.5, 1.5, 0.1, tw, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	}
	add(NewProjectivePlane(Point4{}, one, 0.5, 0.1, Rot4{}, testWhite, testNone, testHalf, testHalf, testIOR))
	points := 200000
	type variant struct {
		name                 string
		factor, du, shortcut Real // du and shortcut in units of k.du
	}
	vars := []variant{
		{"current", 1.01, 1, 1},
		{"halfwindow", 0.5, 1, 1},
		{"du x2", 1.01, 2, 2},
		{"du x3", 1.01, 3, 3},
		{"du x4", 1.01, 4, 4},
		{"du x2 short x1", 1.01, 2, 1},
		{"du x1 short x2", 1.01, 1, 2},
		{"du x1 short x4", 1.01, 1, 4},
	}
	dump := os.Getenv("ZZ_DUMP") != ""
	for _, k := range ks {
		reach := k.MajorRadius + k.MinorRadius
		miss := make([]int, len(vars))
		worst := make([]Real, len(vars))
		shown := 0
		for n := 0; n < points; n++ {
			var q Vector4
			onSurface := func() Vector4 {
				s := (rng.Float64()*2 - 1) * math.Pi
				if k.Band {
					s = (rng.Float64()*2 - 1) * k.MinorRadius
				}
				return kleinPoint(k, rng.Float64()*2*math.Pi, s)
			}
			switch n % 4 {
			case 0:
				q = onSurface().Add(unitS3(rng).Mul(0.3 * rng.Float64()))
			case 1:
				q = Vector4{(rng.Float64()*2 - 1) * reach * 1.2, (rng.Float64()*2 - 1) * reach * 1.2, (rng.Float64()*2 - 1) * k.MinorRadius * 1.5, (rng.Float64()*2 - 1) * k.MinorRadius * 1.5}
			case 2:
				a, b := onSurface(), onSurface()
				q = a.Add(b.Sub(a).Mul(0.3 + 0.4*rng.Float64()))
			default:
				q = onSurface().Add(unitS3(rng).Mul(k.tube * 2 * rng.Float64()))
			}
			ref := zzNearest(k, q, 1.01, k.du/4, 0)
			for i, v := range vars {
				d := zzNearest(k, q, v.factor, v.du*k.du, v.shortcut*k.du)
				if d > ref+1e-9 {
					miss[i]++
					worst[i] = math.Max(worst[i], d-ref)
					if dump && i == 1 && shown < 4 && k.Twist >= 3 && d-ref > 0.02 {
						shown++
						t.Logf("  halfwindow miss: %s tw%d R=%g M=%g q=Vector4{%.6f, %.6f, %.6f, %.6f} ref=%.6f got=%.6f", k.Name, k.Twist, k.MajorRadius, k.MinorRadius, q.X, q.Y, q.Z, q.W, ref, d)
					}
				} else if d < ref-1e-9 {
					t.Logf("  REFERENCE WORSE: %s tw%d q=%v ref=%.9f %s=%.9f", k.Name, k.Twist, q, ref, v.name, d)
				}
			}
		}
		line := ""
		for i, v := range vars {
			line += " | " + v.name + ": "
			if miss[i] == 0 {
				line += "0"
			} else {
				line += itoa(miss[i]) + " (max " + ftoa(worst[i]) + ")"
			}
		}
		t.Logf("%-16s tw%-2d R=%.2f M=%.2f%s", k.Name, k.Twist, k.MajorRadius, k.MinorRadius, line)
	}
}

func itoa(i int) string {
	s := ""
	if i == 0 {
		return "0"
	}
	for i > 0 {
		s = string(rune('0'+i%10)) + s
		i /= 10
	}
	return s
}
func ftoa(f Real) string {
	return itoa(int(f*1000)) + "e-3"
}
