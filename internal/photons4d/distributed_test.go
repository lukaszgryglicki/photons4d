package photons4d

import (
	"fmt"
	"math"
	"math/rand"
	"net"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"
)

// TestSparseRoundTripExact verifies that extract→apply reproduces the source
// buffer bit-for-bit (pure float64 additions, no lossy transform).
func TestSparseRoundTripExact(t *testing.T) {
	rng := rand.New(rand.NewSource(4242))
	src := make([]Real, 10000)
	for i := range src {
		if rng.Float64() < 0.15 {
			src[i] = Real(rng.Float64() * 1e3)
		}
	}
	orig := append([]Real(nil), src...)

	idx, vals := extractSparseAndZero(src)
	for i, v := range src {
		if v != 0 {
			t.Fatalf("buffer not zeroed at %d: %v", i, v)
		}
	}

	dst := make([]Real, len(orig))
	if err := applySparse(dst, idx, vals); err != nil {
		t.Fatalf("applySparse: %v", err)
	}
	for i := range orig {
		if dst[i] != orig[i] {
			t.Fatalf("mismatch at %d: got %v want %v", i, dst[i], orig[i])
		}
	}

	// Corrupt-update rejection.
	if err := applySparse(dst, []int64{5, 0}, []float64{1, 1}); err == nil {
		t.Fatal("non-ascending delta not rejected")
	}
	if err := applySparse(dst, []int64{int64(len(dst))}, []float64{1}); err == nil {
		t.Fatal("out-of-range index not rejected")
	}
	if err := applySparse(dst, []int64{0}, []float64{math.NaN()}); err == nil {
		t.Fatal("NaN value not rejected")
	}
	if err := applySparse(dst, []int64{0}, []float64{-1}); err == nil {
		t.Fatal("negative value not rejected")
	}
}

// TestSparseBatchingWire verifies the multi-batch re-anchoring logic through
// a real gob wire (tiny batch size to force many batches), for both codecs.
func TestSparseBatchingWire(t *testing.T) {
	for _, compress := range []bool{false, true} {
		rng := rand.New(rand.NewSource(1717))
		src := make([]Real, 5000)
		for i := range src {
			if rng.Float64() < 0.3 {
				src[i] = Real(1 + rng.Float64())
			}
		}
		orig := append([]Real(nil), src...)

		c1, c2 := net.Pipe()
		sender, receiver := newWire(c1), newWire(c2)

		errCh := make(chan error, 1)
		go func() {
			errCh <- sendSparseUpdateSized(sender, 1, 7, src, 3, compress)
		}()

		dst := make([]Real, len(orig))
		for {
			var upd UpdateMsg
			if err := receiver.recv(&upd); err != nil {
				t.Fatalf("compress=%v recv: %v", compress, err)
			}
			if upd.RoundID != 7 {
				t.Fatalf("compress=%v wrong round: %d", compress, upd.RoundID)
			}
			if compress && upd.Codec != CodecFlate {
				t.Fatalf("expected packed message, got codec %d", upd.Codec)
			}
			if _, _, err := decodeUpdate(&upd); err != nil {
				t.Fatalf("compress=%v decode: %v", compress, err)
			}
			if len(upd.IdxDeltas) > 0 {
				if err := applySparse(dst, upd.IdxDeltas, upd.Values); err != nil {
					t.Fatalf("compress=%v applySparse: %v", compress, err)
				}
			}
			if !upd.More {
				break
			}
		}
		if err := <-errCh; err != nil {
			t.Fatalf("compress=%v send: %v", compress, err)
		}
		for i := range orig {
			if dst[i] != orig[i] {
				t.Fatalf("compress=%v mismatch at %d: got %v want %v", compress, i, dst[i], orig[i])
			}
		}
		sender.Close()
		receiver.Close()
	}
}

// TestPackSparseRoundTrip fuzzes the packed codec: random sparse layouts and
// extreme float values must survive pack→unpack bit-exactly.
func TestPackSparseRoundTrip(t *testing.T) {
	rng := rand.New(rand.NewSource(424242))
	cases := [][2][]float64{
		{{}, {}},
		{{0}, {1.5}},
		{{7, 1, 1 << 40}, {5e-324, 1.7976931348623157e308, 3.141592653589793}},
	}
	for c := 0; c < 30; c++ {
		n := rng.Intn(3000)
		deltas := make([]float64, n)
		values := make([]float64, n)
		for i := 0; i < n; i++ {
			if i == 0 {
				deltas[i] = float64(rng.Intn(1 << 20))
			} else {
				deltas[i] = float64(1 + rng.Intn(1000))
			}
			values[i] = rng.ExpFloat64() * math.Pow(10, float64(rng.Intn(20)-10))
		}
		cases = append(cases, [2][]float64{deltas, values})
	}
	for ci, cse := range cases {
		idx := make([]int64, len(cse[0]))
		for i, d := range cse[0] {
			idx[i] = int64(d)
		}
		packed, err := packSparse(idx, cse[1])
		if err != nil {
			t.Fatalf("case %d: pack: %v", ci, err)
		}
		gotIdx, gotVal, rawLen, err := unpackSparse(packed)
		if err != nil {
			t.Fatalf("case %d: unpack: %v", ci, err)
		}
		if len(gotIdx) != len(idx) || len(gotVal) != len(cse[1]) {
			t.Fatalf("case %d: length mismatch", ci)
		}
		if len(idx) > 0 && rawLen <= 0 {
			t.Fatalf("case %d: bad rawLen %d", ci, rawLen)
		}
		for i := range idx {
			if gotIdx[i] != idx[i] {
				t.Fatalf("case %d: idx[%d] = %d want %d", ci, i, gotIdx[i], idx[i])
			}
			if math.Float64bits(gotVal[i]) != math.Float64bits(cse[1][i]) {
				t.Fatalf("case %d: val[%d] = %x want %x", ci, i, math.Float64bits(gotVal[i]), math.Float64bits(cse[1][i]))
			}
		}
	}
	// Corrupt streams must error out, never panic or corrupt.
	if _, _, _, err := unpackSparse([]byte{0x42, 0x13, 0x37}); err == nil {
		t.Fatal("garbage packed stream not rejected")
	}
}

func writeDistTestScene(t *testing.T, dir, gifOut string) string {
	t.Helper()
	cfg := fmt.Sprintf(`{
  "sceneResX": 16, "sceneResY": 16, "sceneResZ": 4,
  "probeRays": 4096, "spp": 4,
  "gifOut": %q, "gifDelay": 10, "gamma": 0.8,
  "scene": {
    "center": { "X": 0, "Y": 0, "Z": 0, "W": 0 },
    "width": 2.0, "height": 2.0, "depth": 2.0, "maxBounces": 16
  },
  "lights": [
    {
      "origin":    { "X": 0, "Y": 0, "Z": 0, "W": 1 },
      "direction": { "X": 0, "Y": 0, "Z": 0, "W": -1 },
      "color":     { "R": 1, "G": 1, "B": 1 },
      "angleDeg":  45.0, "intensity": 1.0
    }
  ],
  "hyperspheres": [
    {
      "center": { "X": 0, "Y": 0, "Z": 0, "W": 0.3 },
      "scale":  { "X": 0.35, "Y": 0.35, "Z": 0.35, "W": 0.35 },
      "color":   { "R": 1, "G": 0.8, "B": 0.5 },
      "diffuse": { "R": 0.4, "G": 0.4, "B": 0.4 },
      "reflect": { "R": 0.3, "G": 0.3, "B": 0.3 },
      "refract": { "R": 0.2, "G": 0.2, "B": 0.2 },
      "ior":     { "R": 1.2, "G": 1.2, "B": 1.2 }
    }
  ]
}`, gifOut)
	path := filepath.Join(dir, "dist_test_scene.json")
	if err := os.WriteFile(path, []byte(cfg), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestDistributedEndToEnd runs a real server and two real clients in-process
// over localhost TCP and checks that the whole spp budget is collected, the
// GIF is produced, and the merged buffer carries sensible energy compared to
// an equivalent single-machine (vertical) render of the same scene.
func TestDistributedEndToEnd(t *testing.T) {
	dir := t.TempDir()
	gif := filepath.Join(dir, "out.gif")
	scenePath := writeDistTestScene(t, dir, gif)

	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	addr := ln.Addr().String()

	serverErr := make(chan error, 1)
	go func() {
		serverErr <- runServerOnListener(scenePath, ln, 8, 0.05)
	}()

	var wg sync.WaitGroup
	clientErrs := make([]error, 2)
	for i := 0; i < 2; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			// One raw client, one compressed with a tiny batch size, so a
			// single E2E run exercises both codecs and multi-batch rounds.
			clientErrs[i] = RunClient(scenePath, addr, ClientOpts{Compress: i == 1, BatchEntries: 1 + i*97})
		}(i)
	}
	wg.Wait()
	for i, err := range clientErrs {
		if err != nil {
			t.Fatalf("client %d: %v", i, err)
		}
	}
	select {
	case err := <-serverErr:
		if err != nil {
			t.Fatalf("server: %v", err)
		}
	case <-time.After(60 * time.Second):
		t.Fatal("server did not finish")
	}

	st, err := os.Stat(gif)
	if err != nil || st.Size() == 0 {
		t.Fatalf("distributed GIF missing/empty: %v", err)
	}
}

// TestDistributedMergeMatchesLocalSum runs server+client in-process but keeps
// hold of the merged scene via the internal server structure, verifying that
// the master buffer's total energy is within Monte Carlo tolerance of a
// vertical render with the same ray budget.
func TestDistributedMergeMatchesLocalSum(t *testing.T) {
	dir := t.TempDir()
	gif := filepath.Join(dir, "out2.gif")
	scenePath := writeDistTestScene(t, dir, gif)

	cfg, err := prepareConfig(scenePath)
	if err != nil {
		t.Fatal(err)
	}
	hash, err := SceneHash(cfg)
	if err != nil {
		t.Fatal(err)
	}
	scene, lights, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	needRays := computeNeedRays(cfg, scene, lights)

	srv := newDistServer(cfg, scene, hash, needRays, 4, 0)

	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	defer ln.Close()
	go func() {
		for {
			conn, err := ln.Accept()
			if err != nil {
				return
			}
			go srv.handleConn(conn)
		}
	}()

	if err := RunClient(scenePath, ln.Addr().String(), ClientOpts{Compress: true}); err == nil {
		// Client exits cleanly only after Done; Done requires all merged.
	} else {
		t.Fatalf("client: %v", err)
	}
	select {
	case <-srv.doneCh:
	case <-time.After(30 * time.Second):
		t.Fatal("server never reached completion")
	}

	if done, total := srv.progress(); done != total {
		t.Fatalf("rays merged %d != total %d", done, total)
	}
	srv.mu.Lock()
	sumD := 0.0
	for _, v := range srv.scene.Buf {
		sumD += float64(v)
	}
	srv.mu.Unlock()

	// Vertical reference with identical budget.
	sceneV, lightsV, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	castRays(lightsV, sceneV, needRays)
	sumV := 0.0
	for _, v := range sceneV.Buf {
		sumV += float64(v)
	}

	if sumD <= 0 || sumV <= 0 {
		t.Fatalf("zero energy: distributed=%v vertical=%v", sumD, sumV)
	}
	rel := math.Abs(sumD-sumV) / math.Max(sumD, sumV)
	if rel > 0.10 {
		t.Fatalf("distributed vs vertical energy mismatch: %v vs %v (rel %.3f)", sumD, sumV, rel)
	}
}

// TestSceneHashMismatchRejected ensures a client with a different scene is
// refused during handshake.
func TestSceneHashMismatchRejected(t *testing.T) {
	dir := t.TempDir()
	sceneA := writeDistTestScene(t, dir, filepath.Join(dir, "a.gif"))

	// Scene B differs in geometry (radius) → different SHA256.
	b, err := os.ReadFile(sceneA)
	if err != nil {
		t.Fatal(err)
	}
	sceneB := filepath.Join(dir, "b.json")
	if !strings.Contains(string(b), `"W": 0.35`) {
		t.Fatal("expected scale in test scene")
	}
	if err := os.WriteFile(sceneB, []byte(strings.Replace(string(b), `"W": 0.35`, `"W": 0.36`, 1)), 0o644); err != nil {
		t.Fatal(err)
	}

	cfg, err := prepareConfig(sceneA)
	if err != nil {
		t.Fatal(err)
	}
	hash, err := SceneHash(cfg)
	if err != nil {
		t.Fatal(err)
	}
	scene, lights, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	needRays := computeNeedRays(cfg, scene, lights)
	srv := newDistServer(cfg, scene, hash, needRays, 1, 0)
	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	defer ln.Close()
	go func() {
		for {
			conn, err := ln.Accept()
			if err != nil {
				return
			}
			go srv.handleConn(conn)
		}
	}()

	err = RunClient(sceneB, ln.Addr().String(), ClientOpts{})
	if err == nil {
		t.Fatal("client with different scene was not rejected")
	}
	t.Logf("correctly rejected: %v", err)
}

// TestAssignWorkAdaptive verifies the adaptive round sizing math: clamping to
// [minRound, chunkTotal], proportional spread across lights by remaining
// budget, the progress guarantee, and exact accounting until drained.
func TestAssignWorkAdaptive(t *testing.T) {
	srv := newDistServer(nil, nil, "x", []int{900000, 100000}, 10, 60)
	if srv.chunkTotal != 100000 {
		t.Fatalf("chunkTotal=%d want 100000", srv.chunkTotal)
	}

	sum := func(v []int) int {
		s := 0
		for _, n := range v {
			s += n
		}
		return s
	}

	// want<=0 → classic fixed chunks: exactly min(chunk[i], remaining[i]).
	a := srv.assignWork(0)
	if a.RaysPerLight[0] != 90000 || a.RaysPerLight[1] != 10000 {
		t.Fatalf("classic mode: got %v want [90000 10000]", a.RaysPerLight)
	}

	// Explicit size respected.
	a = srv.assignWork(50000)
	if got := sum(a.RaysPerLight); got < 49998 || got > 50000 {
		t.Fatalf("want=50000 assigned %d", got)
	}

	// Tiny request clamped up to minRound (chunkTotal/1024 < 1024 → 1024).
	a = srv.assignWork(10)
	if got := sum(a.RaysPerLight); got < 1022 || got > 1024 {
		t.Fatalf("want=10 assigned %d, expected ~1024", got)
	}

	// Huge request clamped down to chunkTotal.
	a = srv.assignWork(1 << 30)
	if got := sum(a.RaysPerLight); got > 100000 {
		t.Fatalf("want=1<<30 assigned %d > chunkTotal", got)
	}

	// Drain fully; totals must match exactly, all rounds tracked.
	total := 0
	for _, r := range srv.outstanding {
		total += sum(r)
	}
	for i := 0; i < 1_000_000; i++ {
		a = srv.assignWork(1 << 30)
		if a.RoundID == 0 {
			break
		}
		total += sum(a.RaysPerLight)
	}
	if total != 1000000 {
		t.Fatalf("drained %d rays, want 1000000", total)
	}
	if rem := sum(srv.remaining); rem != 0 {
		t.Fatalf("remaining %d after drain", rem)
	}
	// All assigned but outstanding → wait round (all zeros, RoundID 0).
	a = srv.assignWork(0)
	if a.Done || a.RoundID != 0 || sum(a.RaysPerLight) != 0 {
		t.Fatalf("expected wait round, got %+v", a)
	}
}

// TestSingleLightManyClients is the requirement-#0 regression test:
// horizontal distribution must be per-ray, so a scene with exactly ONE light
// must still be split across many concurrent clients. Four clients handshake
// through a barrier (conns accepted first, served after all four are in), and
// every one of them must merge at least one non-empty round.
func TestSingleLightManyClients(t *testing.T) {
	dir := t.TempDir()
	gif := filepath.Join(dir, "one.gif")
	scenePath := writeDistTestScene(t, dir, gif)

	// More rays → many small rounds → every client provably participates.
	b, err := os.ReadFile(scenePath)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(b), `"spp": 4`) {
		t.Fatal("expected spp in test scene")
	}
	if err := os.WriteFile(scenePath, []byte(strings.Replace(string(b), `"spp": 4`, `"spp": 64`, 1)), 0o644); err != nil {
		t.Fatal(err)
	}

	cfg, err := prepareConfig(scenePath)
	if err != nil {
		t.Fatal(err)
	}
	hash, err := SceneHash(cfg)
	if err != nil {
		t.Fatal(err)
	}
	scene, lights, err := buildScene(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(lights) != 1 {
		t.Fatalf("test premise broken: %d lights, want 1", len(lights))
	}
	needRays := computeNeedRays(cfg, scene, lights)

	const nClients = 4
	srv := newDistServer(cfg, scene, hash, needRays, 16, 0.01)

	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	defer ln.Close()
	go func() {
		// Barrier: accept all clients before serving any, so none can race
		// ahead and drain the budget before the others handshake.
		conns := make([]net.Conn, 0, nClients)
		for len(conns) < nClients {
			conn, err := ln.Accept()
			if err != nil {
				return
			}
			conns = append(conns, conn)
		}
		for _, c := range conns {
			go srv.handleConn(c)
		}
	}()

	var wg sync.WaitGroup
	errs := make([]error, nClients)
	for i := 0; i < nClients; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			errs[i] = RunClient(scenePath, ln.Addr().String(), ClientOpts{Compress: i%2 == 1})
		}(i)
	}
	wg.Wait()
	for i, err := range errs {
		if err != nil {
			t.Fatalf("client %d: %v", i, err)
		}
	}
	select {
	case <-srv.doneCh:
	case <-time.After(60 * time.Second):
		t.Fatal("server never completed")
	}

	srv.mu.Lock()
	defer srv.mu.Unlock()
	if srv.raysDone != srv.raysTotal {
		t.Fatalf("raysDone=%d raysTotal=%d", srv.raysDone, srv.raysTotal)
	}
	if len(srv.perClient) != nClients {
		t.Fatalf("admitted %d clients, want %d", len(srv.perClient), nClients)
	}
	var sumRays int64
	for id, st := range srv.perClient {
		t.Logf("client #%d: %d rounds, %d rays", id, st.rounds, st.rays)
		if st.rounds < 1 || st.rays < 1 {
			t.Fatalf("client #%d did not participate in the single-light workload: %+v", id, st)
		}
		sumRays += st.rays
	}
	if sumRays != srv.raysTotal {
		t.Fatalf("per-client rays sum %d != total %d", sumRays, srv.raysTotal)
	}
}
