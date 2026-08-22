package photons4d

import (
	"fmt"
	"net"
	"sort"
	"sync"
	"time"
)

// statsLoop periodically prints cluster progress and every client's
// contribution until casting completes. Runs on the server only.
func (s *distServer) statsLoop(every time.Duration) {
	tick := time.NewTicker(every)
	defer tick.Stop()
	last := time.Now()
	var lastDone int64
	for {
		select {
		case <-s.doneCh:
			return
		case now := <-tick.C:
			type row struct {
				id           int
				host         string
				rounds, rays int64
			}
			s.mu.Lock()
			done, total := s.raysDone, s.raysTotal
			connected := s.clients
			outst := len(s.outstanding)
			rows := make([]row, 0, len(s.perClient))
			for id, st := range s.perClient {
				rows = append(rows, row{id, st.host, st.rounds, st.rays})
			}
			s.mu.Unlock()
			sort.Slice(rows, func(i, j int) bool { return rows[i].id < rows[j].id })
			rate := float64(done-lastDone) / now.Sub(last).Seconds()
			eta := "n/a"
			if rate > 0 && done < total {
				eta = time.Duration(float64(total-done) / rate * float64(time.Second)).Round(time.Second).String()
			}
			logf("[STATS] merged %.2f%% (%d/%d rays) | window rate %.2fM rays/s | ETA %s | clients connected %d | outstanding rounds %d\n",
				100*float64(done)/float64(total), done, total, rate/1e6, eta, connected, outst)
			for _, r := range rows {
				logf("[STATS]   client #%d (%s): %d rounds, %d rays (%.2f%% of total)\n",
					r.id, r.host, r.rounds, r.rays, 100*float64(r.rays)/float64(total))
			}
			last, lastDone = now, done
		}
	}
}

// clientStat tracks one admitted client's contribution (for the final
// summary and the adaptive round sizing).
type clientStat struct {
	host   string
	rounds int64
	rays   int64
}

// distServer owns the master scene buffer and doles out ray-casting rounds.
type distServer struct {
	cfg       *Config
	scene     *Scene
	sceneHash string
	needRays  []int // authoritative per-light totals (spp saturation target)

	mu          sync.Mutex // guards everything below + scene.Buf merges
	remaining   []int      // per-light rays not yet assigned
	outstanding map[int][]int
	nextRound   int
	nextClient  int
	raysDone    int64
	raysTotal   int64
	entriesIn   int64
	packedIn    int64 // compressed bytes received (CodecFlate batches)
	unpackedIn  int64 // same batches after decompression
	clients     int
	perClient   map[int]*clientStat
	doneCh      chan struct{}
	doneOnce    sync.Once
	chunk       []int // per-light max rays per assignment
	chunkTotal  int   // sum(chunk)
	roundTarget float64
}

// RunServer starts distributed-mode collector: it builds the scene, computes
// the spp-saturation ray budget exactly like vertical mode, serves rounds to
// clients until the budget is exhausted and every round is merged, then runs
// the unchanged single-node output pipeline (GIF/PNG/RAW).
// roundSec > 0 sizes rounds adaptively so each client's round takes about
// that long (measured per client); roundSec <= 0 uses fixed max-size chunks.
// statsSec sets the progress-report interval (<= 0 → 300s default).
func RunServer(cfgPath, listenAddr string, chunks int, roundSec, statsSec float64) error {
	ln, err := net.Listen("tcp", listenAddr)
	if err != nil {
		return fmt.Errorf("server listen %s: %w", listenAddr, err)
	}
	return runServerOnListener(cfgPath, ln, chunks, roundSec, statsSec)
}

// newDistServer wires up server state for a prepared scene. chunks controls
// the per-light max round size (budget/chunks); roundSec > 0 enables
// adaptive per-client round sizing.
func newDistServer(cfg *Config, scene *Scene, hash string, needRays []int, chunks int, roundSec float64) *distServer {
	if chunks < 1 {
		chunks = 1
	}
	srv := &distServer{
		cfg:         cfg,
		scene:       scene,
		sceneHash:   hash,
		needRays:    needRays,
		remaining:   append([]int(nil), needRays...),
		outstanding: make(map[int][]int),
		nextRound:   1,
		doneCh:      make(chan struct{}),
		chunk:       make([]int, len(needRays)),
		perClient:   make(map[int]*clientStat),
		roundTarget: roundSec,
	}
	for i, n := range needRays {
		c := n / chunks
		if c < 1 {
			c = 1
		}
		srv.chunk[i] = c
		srv.chunkTotal += c
		srv.raysTotal += int64(n)
	}
	return srv
}

func runServerOnListener(cfgPath string, ln net.Listener, chunks int, roundSec, statsSec float64) error {
	defer ln.Close()

	cfg, err := prepareConfig(cfgPath)
	if err != nil {
		return err
	}
	hash, err := SceneHash(cfg)
	if err != nil {
		return err
	}
	scene, lights, err := buildScene(cfg)
	if err != nil {
		return err
	}
	needRays := computeNeedRays(cfg, scene, lights)

	srv := newDistServer(cfg, scene, hash, needRays, chunks, roundSec)

	logf("[SERVER] listening on %s | scene %s | sha256 %s\n", ln.Addr(), cfgPath, hash)
	logf("[SERVER] total rays needed: %d (per light: %v), max chunk per light: %v, round target: %.0fs\n",
		srv.raysTotal, needRays, srv.chunk, roundSec)

	if statsSec <= 0 {
		statsSec = 300
	}
	go srv.statsLoop(time.Duration(statsSec * float64(time.Second)))

	var wg sync.WaitGroup
	acceptDone := make(chan struct{})
	go func() {
		defer close(acceptDone)
		for {
			conn, err := ln.Accept()
			if err != nil {
				return // listener closed on completion
			}
			wg.Add(1)
			go func() {
				defer wg.Done()
				srv.handleConn(conn)
			}()
		}
	}()

	<-srv.doneCh
	srv.mu.Lock()
	raysDone, entriesIn, packedIn, unpackedIn := srv.raysDone, srv.entriesIn, srv.packedIn, srv.unpackedIn
	srv.mu.Unlock()
	logf("[SERVER] all %d rays merged (%d sparse entries received), finalizing\n", raysDone, entriesIn)
	if packedIn > 0 {
		logf("[SERVER] compressed updates: %d bytes on the wire for %d raw bytes (%.2fx saved)\n",
			packedIn, unpackedIn, float64(unpackedIn)/float64(packedIn))
	}
	srv.mu.Lock()
	ids := make([]int, 0, len(srv.perClient))
	for id := range srv.perClient {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	for _, id := range ids {
		st := srv.perClient[id]
		logf("[SERVER] client #%d (%s): %d rounds, %d rays (%.2f%% of total)\n",
			id, st.host, st.rounds, st.rays, 100.0*float64(st.rays)/float64(srv.raysTotal))
	}
	srv.mu.Unlock()

	// Stop accepting, give connected clients a grace period to pick up their
	// Done notice, then force-close whatever is left.
	ln.Close()
	graceDone := make(chan struct{})
	go func() {
		wg.Wait()
		close(graceDone)
	}()
	select {
	case <-graceDone:
	case <-time.After(30 * time.Second):
		logf("[SERVER] grace period expired, proceeding with output\n")
	}
	<-acceptDone

	return saveOutputs(cfg, scene)
}

// assignWork returns the next round for a client, an all-zero wait assignment
// (stragglers still hold outstanding rounds), or Done.
//
// want <= 0 (default, -round 0): classic fixed-chunk mode — every light
// contributes min(chunk[i], remaining[i]), byte-for-byte the original
// behavior. want > 0 (optional adaptive mode, -round > 0): the
// client-specific desired total ray count (EWMA rays/s x round target),
// clamped to [minRound, chunkTotal] and spread across lights proportionally
// to their remaining budgets. Either way any light mix — including a single
// light — is split across however many clients are pulling work.
func (s *distServer) assignWork(want int) WorkAssign {
	s.mu.Lock()
	defer s.mu.Unlock()

	remTotal := 0
	for _, r := range s.remaining {
		remTotal += r
	}
	if remTotal == 0 {
		if len(s.outstanding) == 0 {
			return WorkAssign{Done: true}
		}
		// Everything assigned but not all merged yet: ask client to stand by
		// in case a straggler dies and its rays come back to the pool.
		return WorkAssign{RoundID: 0, RaysPerLight: make([]int, len(s.remaining))}
	}

	rays := make([]int, len(s.remaining))
	if want <= 0 {
		// Classic fixed-chunk distribution (default).
		for i, rem := range s.remaining {
			n := s.chunk[i]
			if n > rem {
				n = rem
			}
			rays[i] = n
		}
	} else {
		if want > s.chunkTotal {
			want = s.chunkTotal
		}
		minRound := s.chunkTotal / 1024
		if minRound < 1024 {
			minRound = 1024
		}
		if want < minRound {
			want = minRound
		}
		assigned := 0
		for i, rem := range s.remaining {
			n := int(float64(want) * float64(rem) / float64(remTotal))
			if n > rem {
				n = rem
			}
			rays[i] = n
			assigned += n
		}
		if assigned == 0 {
			// Rounding starved the round; guarantee progress on the largest
			// remaining light.
			bi := 0
			for i, r := range s.remaining {
				if r > s.remaining[bi] {
					bi = i
				}
			}
			n := s.remaining[bi]
			if n > minRound {
				n = minRound
			}
			rays[bi] = n
		}
	}
	for i, n := range rays {
		s.remaining[i] -= n
	}
	id := s.nextRound
	s.nextRound++
	s.outstanding[id] = rays
	return WorkAssign{RoundID: id, RaysPerLight: rays}
}

// returnWork puts a failed client's outstanding rounds back into the pool.
func (s *distServer) returnWork(rounds map[int]bool) {
	if len(rounds) == 0 {
		return
	}
	s.mu.Lock()
	defer s.mu.Unlock()
	for id := range rounds {
		rays, ok := s.outstanding[id]
		if !ok {
			continue
		}
		delete(s.outstanding, id)
		for i, n := range rays {
			s.remaining[i] += n
		}
		logf("[SERVER] round %d returned to pool (client lost): %v rays\n", id, rays)
	}
}

// completeRound transactionally merges all buffered batches of a round into
// the master buffer and updates progress accounting.
func (s *distServer) completeRound(clientID, roundID int, batches []UpdateMsg) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	rays, ok := s.outstanding[roundID]
	if !ok {
		return fmt.Errorf("unknown or already-completed round %d", roundID)
	}
	for _, b := range batches {
		if err := applySparse(s.scene.Buf, b.IdxDeltas, b.Values); err != nil {
			return err
		}
		s.entriesIn += int64(len(b.Values))
	}
	delete(s.outstanding, roundID)
	roundRays := int64(0)
	for _, n := range rays {
		roundRays += int64(n)
	}
	s.raysDone += roundRays
	if st, ok := s.perClient[clientID]; ok {
		st.rounds++
		st.rays += roundRays
	}

	pct := 100.0 * float64(s.raysDone) / float64(s.raysTotal)
	logf("[SERVER] round %d merged | rays %d/%d (%.2f%%) | outstanding rounds %d\n",
		roundID, s.raysDone, s.raysTotal, pct, len(s.outstanding))

	if s.raysDone >= s.raysTotal && len(s.outstanding) == 0 {
		remTotal := 0
		for _, r := range s.remaining {
			remTotal += r
		}
		if remTotal == 0 {
			s.doneOnce.Do(func() { close(s.doneCh) })
		}
	}
	return nil
}

func (s *distServer) progress() (done, total int64) {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.raysDone, s.raysTotal
}

// handleConn drives one client connection through handshake and rounds.
func (s *distServer) handleConn(conn net.Conn) {
	w := newWire(conn)
	defer w.Close()

	owned := make(map[int]bool) // rounds assigned to this conn, not yet merged
	defer func() { s.returnWork(owned) }()

	var hello HelloMsg
	if err := w.recv(&hello); err != nil {
		logf("[SERVER] %s handshake read failed: %v\n", conn.RemoteAddr(), err)
		return
	}
	if hello.Version != DistProtocolVersion {
		_ = w.send(HelloAck{OK: false, Reason: fmt.Sprintf("protocol version mismatch: server=%d client=%d", DistProtocolVersion, hello.Version)})
		return
	}
	if hello.SceneHash != s.sceneHash {
		logf("[SERVER] rejecting %s (%s): scene hash mismatch (server=%s client=%s)\n", conn.RemoteAddr(), hello.Host, s.sceneHash, hello.SceneHash)
		_ = w.send(HelloAck{OK: false, Reason: fmt.Sprintf("scene SHA256 mismatch: server=%s client=%s — make sure both run the identical scene file, binary version and env (DEBUG/FORCE_ESCAPE/SPP_ADJUST)", s.sceneHash, hello.SceneHash)})
		return
	}
	if hello.BufLen != len(s.scene.Buf) {
		_ = w.send(HelloAck{OK: false, Reason: fmt.Sprintf("buffer size mismatch: server=%d client=%d", len(s.scene.Buf), hello.BufLen)})
		return
	}
	if hello.NLights != len(s.needRays) {
		_ = w.send(HelloAck{OK: false, Reason: fmt.Sprintf("light count mismatch: server=%d client=%d", len(s.needRays), hello.NLights)})
		return
	}

	s.mu.Lock()
	s.nextClient++
	clientID := s.nextClient
	s.clients++
	s.perClient[clientID] = &clientStat{host: hello.Host}
	s.mu.Unlock()
	if err := w.send(HelloAck{OK: true, ClientID: clientID}); err != nil {
		return
	}
	logf("[SERVER] client #%d connected from %s (%s, %d workers)\n", clientID, conn.RemoteAddr(), hello.Host, hello.Workers)

	rate := 0.0 // EWMA rays/s for this client, 0 = not measured yet
	for {
		var req WorkRequest
		if err := w.recv(&req); err != nil {
			logf("[SERVER] client #%d disconnected: %v\n", clientID, err)
			return
		}
		want := 0
		if s.roundTarget > 0 {
			if rate > 0 {
				want = int(rate * s.roundTarget)
			} else {
				// Unknown speed: start with a small probe round and ramp up.
				want = s.chunkTotal / 8
			}
		}
		assign := s.assignWork(want)
		if err := w.send(assign); err != nil {
			return
		}
		if assign.Done {
			logf("[SERVER] client #%d released (all work done)\n", clientID)
			return
		}
		if assign.RoundID != 0 {
			owned[assign.RoundID] = true
		}
		roundRays := 0
		for _, n := range assign.RaysPerLight {
			roundRays += n
		}
		tAssign := time.Now()

		// Receive the round's batches (possibly zero entries for a wait
		// round). Packed batches are decompressed here, outside the merge
		// lock, so inflation from many clients runs in parallel.
		var batches []UpdateMsg
		var wireB, rawB int64
		for {
			var upd UpdateMsg
			if err := w.recv(&upd); err != nil {
				logf("[SERVER] client #%d died mid-round %d: %v\n", clientID, assign.RoundID, err)
				return
			}
			if upd.RoundID != assign.RoundID {
				_ = w.send(UpdateAck{OK: false, Reason: fmt.Sprintf("unexpected round %d (assigned %d)", upd.RoundID, assign.RoundID)})
				return
			}
			wl, rl, err := decodeUpdate(&upd)
			if err != nil {
				logf("[SERVER] client #%d round %d undecodable: %v\n", clientID, assign.RoundID, err)
				_ = w.send(UpdateAck{OK: false, Reason: err.Error()})
				return
			}
			wireB += int64(wl)
			rawB += int64(rl)
			if len(upd.IdxDeltas) > 0 {
				batches = append(batches, upd)
			}
			if !upd.More {
				break
			}
		}

		if assign.RoundID != 0 {
			// Update this client's measured throughput (cast + extract +
			// upload) for adaptive sizing of its next round.
			if roundRays > 0 {
				if el := time.Since(tAssign).Seconds(); el > 0.005 {
					r := float64(roundRays) / el
					if rate <= 0 {
						rate = r
					} else {
						rate = 0.5*rate + 0.5*r
					}
				}
			}
			if wireB > 0 {
				s.mu.Lock()
				s.packedIn += wireB
				s.unpackedIn += rawB
				s.mu.Unlock()
			}
			if err := s.completeRound(clientID, assign.RoundID, batches); err != nil {
				logf("[SERVER] client #%d round %d rejected: %v\n", clientID, assign.RoundID, err)
				_ = w.send(UpdateAck{OK: false, Reason: err.Error()})
				return
			}
			delete(owned, assign.RoundID)
		}
		done, total := s.progress()
		if err := w.send(UpdateAck{OK: true, RaysDone: done, RaysTotal: total}); err != nil {
			return
		}
	}
}
