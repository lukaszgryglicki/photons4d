package photons4d

import (
	"math/rand"
	"testing"
)

// A light inside a closed, purely diffuse box: every scatter happens on the
// inner side of a wall, so no photon may get out.
func TestDiffuseScatterStaysOnTheIncidentSide(t *testing.T) {
	const bounces = 12
	scene := NewScene(Point4{0, 0, 0, 0}, 2, 2, 2, 4, 4, 4, bounces, true)
	box, err := NewCell8(Point4{0, 0, 0, 0}, Vector4{1, 1.2, 0.8, 1.1}, Rot4{XY: 0.4, XW: -0.3, ZW: 0.9},
		RGB{1, 1, 1}, RGB{1, 1, 1}, RGB{0, 0, 0}, RGB{0, 0, 0}, RGB{1.5, 1.5, 1.5})
	if err != nil {
		t.Fatal(err)
	}
	scene.AddCell8(box)
	light, err := NewLight(Point4{0.1, -0.05, 0.02, 0.07}, Vector4{0, 0, 0, 1}, RGB{1, 1, 1}, 3.14159, 1)
	if err != nil {
		t.Fatal(err)
	}
	rng := rand.New(rand.NewSource(7))
	for i := 0; i < 2000; i++ {
		escaped, depth := castSingleRay(light, scene, rng, nil, false)
		if escaped || depth != bounces {
			t.Fatalf("ray %d left the diffuse box: escaped=%v after %d interactions", i, escaped, depth)
		}
	}

	// From outside the same box sends everything back: one scatter, then the photon escapes.
	outside, err := NewLight(Point4{0, 0, 0, 3}, Vector4{0, 0, 0, -1}, RGB{1, 1, 1}, 0.05, 1)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 2000; i++ {
		escaped, depth := castSingleRay(outside, scene, rng, nil, false)
		if !escaped || depth != 1 {
			t.Fatalf("ray %d from outside: escaped=%v after %d interactions, want one scatter", i, escaped, depth)
		}
	}
}
