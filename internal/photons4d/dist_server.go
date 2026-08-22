package photons4d

import (
	"fmt"
	"net"
	"sync"
	"time"
)

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
	doneCh      chan struct{}
	doneOnce    sync.Once
	chunk       []int // per-light rays per assignment
}

// RunServer starts distributed-mode collector: it builds the scene, computes
// the spp-saturation ray budget exactly like vertical mode, serves rounds to
// clients until the budget is exhausted and every round is merged, then runs
// the unchanged single-node output pipeline (GIF/PNG/RAW).
func RunServer(cfgPath, listenAddr string, chunks int) error {
	ln, err := net.Listen("tcp", listenAddr)
	if err != nil {
		return fmt.Errorf("server listen %s: %w", listenAddr, err)
	}
	return runServerOnListener(cfgPath, ln, chunks)
}

func runServerOnListener(cfgPath string, ln net.Listener, chunks int) error {
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
	}
	for i, n := range needRays {
		c := n / chunks
		if c < 1 {
			c = 1
		}
		srv.chunk[i] = c
		srv.raysTotal += int64(n)
	}

	fmt.Printf("[SERVER] listening on %s | scene %s | sha256 %s\n", ln.Addr(), cfgPath, hash)
	fmt.Printf("[SERVER] total rays needed: %d (per light: %v), chunks per light: %v\n", srv.raysTotal, needRays, srv.chunk)

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
	fmt.Printf("[SERVER] all %d rays merged (%d sparse entries received), finalizing\n", raysDone, entriesIn)
	if packedIn > 0 {
		fmt.Printf("[SERVER] compressed updates: %d bytes on the wire for %d raw bytes (%.2fx saved)\n",
			packedIn, unpackedIn, float64(unpackedIn)/float64(packedIn))
	}

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
		fmt.Printf("[SERVER] grace period expired, proceeding with output\n")
	}
	<-acceptDone

	return saveOutputs(cfg, scene)
}

// assignWork returns the next round for a client, an all-zero wait assignment
// (stragglers still hold outstanding rounds), or Done.
func (s *distServer) assignWork() WorkAssign {
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
	for i := range s.remaining {
		n := s.chunk[i]
		if n > s.remaining[i] {
			n = s.remaining[i]
		}
		rays[i] = n
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
		fmt.Printf("[SERVER] round %d returned to pool (client lost): %v rays\n", id, rays)
	}
}

// completeRound transactionally merges all buffered batches of a round into
// the master buffer and updates progress accounting.
func (s *distServer) completeRound(roundID int, batches []UpdateMsg) error {
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
	for _, n := range rays {
		s.raysDone += int64(n)
	}

	pct := 100.0 * float64(s.raysDone) / float64(s.raysTotal)
	fmt.Printf("[SERVER] round %d merged | rays %d/%d (%.2f%%) | outstanding rounds %d\n",
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
		fmt.Printf("[SERVER] %s handshake read failed: %v\n", conn.RemoteAddr(), err)
		return
	}
	if hello.Version != DistProtocolVersion {
		_ = w.send(HelloAck{OK: false, Reason: fmt.Sprintf("protocol version mismatch: server=%d client=%d", DistProtocolVersion, hello.Version)})
		return
	}
	if hello.SceneHash != s.sceneHash {
		fmt.Printf("[SERVER] rejecting %s (%s): scene hash mismatch (server=%s client=%s)\n", conn.RemoteAddr(), hello.Host, s.sceneHash, hello.SceneHash)
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
	s.mu.Unlock()
	if err := w.send(HelloAck{OK: true, ClientID: clientID}); err != nil {
		return
	}
	fmt.Printf("[SERVER] client #%d connected from %s (%s, %d workers)\n", clientID, conn.RemoteAddr(), hello.Host, hello.Workers)

	for {
		var req WorkRequest
		if err := w.recv(&req); err != nil {
			fmt.Printf("[SERVER] client #%d disconnected: %v\n", clientID, err)
			return
		}
		assign := s.assignWork()
		if err := w.send(assign); err != nil {
			return
		}
		if assign.Done {
			fmt.Printf("[SERVER] client #%d released (all work done)\n", clientID)
			return
		}
		if assign.RoundID != 0 {
			owned[assign.RoundID] = true
		}

		// Receive the round's batches (possibly zero entries for a wait
		// round). Packed batches are decompressed here, outside the merge
		// lock, so inflation from many clients runs in parallel.
		var batches []UpdateMsg
		var wireB, rawB int64
		for {
			var upd UpdateMsg
			if err := w.recv(&upd); err != nil {
				fmt.Printf("[SERVER] client #%d died mid-round %d: %v\n", clientID, assign.RoundID, err)
				return
			}
			if upd.RoundID != assign.RoundID {
				_ = w.send(UpdateAck{OK: false, Reason: fmt.Sprintf("unexpected round %d (assigned %d)", upd.RoundID, assign.RoundID)})
				return
			}
			wl, rl, err := decodeUpdate(&upd)
			if err != nil {
				fmt.Printf("[SERVER] client #%d round %d undecodable: %v\n", clientID, assign.RoundID, err)
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
			if wireB > 0 {
				s.mu.Lock()
				s.packedIn += wireB
				s.unpackedIn += rawB
				s.mu.Unlock()
			}
			if err := s.completeRound(assign.RoundID, batches); err != nil {
				fmt.Printf("[SERVER] client #%d round %d rejected: %v\n", clientID, assign.RoundID, err)
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
