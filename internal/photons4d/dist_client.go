package photons4d

import (
	"fmt"
	"net"
	"os"
	"runtime"
	"time"
)

// ClientOpts tunes the distributed worker's upload behavior.
type ClientOpts struct {
	// Compress enables DEFLATE compression of voxel updates (CodecFlate).
	// Bit-exact; trades a little CPU for a large network-traffic cut.
	Compress bool
	// BatchEntries caps sparse entries per UpdateMsg (<=0: default 4M ≈ 64 MB
	// raw). Bigger batches mean fewer messages and better compression at the
	// cost of per-message memory on both sides.
	BatchEntries int
}

// RunClient runs distributed-mode worker: it loads and builds the very same
// scene as the server (verified via SHA256 handshake), then loops requesting
// per-light ray budgets, casting them with the full local worker pool
// (identical castRays code path as vertical mode) and shipping sparse voxel
// deltas back. It performs no output of its own.
func RunClient(cfgPath, serverAddr string, opts ClientOpts) error {
	if opts.BatchEntries <= 0 {
		opts.BatchEntries = updateBatchEntries
	}
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

	host, _ := os.Hostname()

	// Retry the initial connect so clients can be started before the server
	// (convenient when spraying jobs across a fleet).
	var conn net.Conn
	const attempts = 60
	for i := 1; ; i++ {
		conn, err = net.DialTimeout("tcp", serverAddr, 10*time.Second)
		if err == nil {
			break
		}
		if i >= attempts {
			return fmt.Errorf("client: cannot reach server %s after %d attempts: %w", serverAddr, attempts, err)
		}
		fmt.Printf("[CLIENT] server %s not reachable (attempt %d/%d): %v — retrying in 2s\n", serverAddr, i, attempts, err)
		time.Sleep(2 * time.Second)
	}
	w := newWire(conn)
	defer w.Close()

	if err := w.send(HelloMsg{
		Version:   DistProtocolVersion,
		SceneHash: hash,
		BufLen:    len(scene.Buf),
		NLights:   len(lights),
		Workers:   runtime.GOMAXPROCS(0),
		Host:      host,
	}); err != nil {
		return fmt.Errorf("client handshake send: %w", err)
	}
	var ack HelloAck
	if err := w.recv(&ack); err != nil {
		return fmt.Errorf("client handshake recv: %w", err)
	}
	if !ack.OK {
		return fmt.Errorf("client rejected by server: %s", ack.Reason)
	}
	fmt.Printf("[CLIENT] connected to %s as client #%d | scene sha256 %s | compress=%v batch=%d\n", serverAddr, ack.ClientID, hash, opts.Compress, opts.BatchEntries)

	rounds := 0
	raysCast := int64(0)
	start := time.Now()
	for {
		if err := w.send(WorkRequest{ClientID: ack.ClientID}); err != nil {
			return fmt.Errorf("client: work request: %w", err)
		}
		var assign WorkAssign
		if err := w.recv(&assign); err != nil {
			return fmt.Errorf("client: work assignment: %w", err)
		}
		if assign.Done {
			fmt.Printf("[CLIENT] all work done: %d rounds, %d rays in %s\n", rounds, raysCast, time.Since(start).Truncate(time.Millisecond))
			return nil
		}
		if len(assign.RaysPerLight) != len(lights) {
			return fmt.Errorf("client: server assigned %d light budgets for %d lights", len(assign.RaysPerLight), len(lights))
		}

		roundRays := 0
		for _, n := range assign.RaysPerLight {
			roundRays += n
		}
		if roundRays > 0 {
			fmt.Printf("[CLIENT] round %d: casting %d rays %v\n", assign.RoundID, roundRays, assign.RaysPerLight)
			castRays(lights, scene, assign.RaysPerLight)
			rounds++
			raysCast += int64(roundRays)
		}

		if err := sendSparseUpdate(w, ack.ClientID, assign.RoundID, scene.Buf, opts); err != nil {
			return fmt.Errorf("client: sending round %d update: %w", assign.RoundID, err)
		}
		var uack UpdateAck
		if err := w.recv(&uack); err != nil {
			return fmt.Errorf("client: update ack: %w", err)
		}
		if !uack.OK {
			return fmt.Errorf("client: server rejected round %d: %s", assign.RoundID, uack.Reason)
		}
		if roundRays > 0 {
			pct := 0.0
			if uack.RaysTotal > 0 {
				pct = 100.0 * float64(uack.RaysDone) / float64(uack.RaysTotal)
			}
			fmt.Printf("[CLIENT] round %d merged | cluster progress %d/%d rays (%.2f%%)\n", assign.RoundID, uack.RaysDone, uack.RaysTotal, pct)
		} else {
			// Wait round: no work available right now, stand by.
			time.Sleep(2 * time.Second)
		}
	}
}

// sendSparseUpdate extracts all non-zero voxel deltas (zeroing the local
// buffer for the next round) and streams them in bounded batches, optionally
// DEFLATE-compressed.
func sendSparseUpdate(w *wire, clientID, roundID int, buf []Real, opts ClientOpts) error {
	return sendSparseUpdateSized(w, clientID, roundID, buf, opts.BatchEntries, opts.Compress)
}

func sendSparseUpdateSized(w *wire, clientID, roundID int, buf []Real, batchEntries int, compress bool) error {
	idxDeltas, values := extractSparseAndZero(buf)
	n := len(values)
	if n == 0 {
		return w.send(UpdateMsg{ClientID: clientID, RoundID: roundID, More: false})
	}
	prevAbs := int64(0)
	for off := 0; off < n; off += batchEntries {
		end := off + batchEntries
		if end > n {
			end = n
		}
		if off > 0 {
			// Re-anchor: first delta of a batch must be the absolute index.
			// prevAbs is the absolute index of the last entry of the previous
			// batch; idxDeltas is scratch, so rewriting in place is safe.
			idxDeltas[off] += prevAbs
		}
		prevAbs = idxDeltas[off]
		for _, d := range idxDeltas[off+1 : end] {
			prevAbs += d
		}
		msg := UpdateMsg{
			ClientID: clientID,
			RoundID:  roundID,
			More:     end < n,
		}
		if compress {
			packed, err := packSparse(idxDeltas[off:end], values[off:end])
			if err != nil {
				return err
			}
			msg.Codec = CodecFlate
			msg.Packed = packed
		} else {
			msg.IdxDeltas = idxDeltas[off:end]
			msg.Values = values[off:end]
		}
		if err := w.send(msg); err != nil {
			return err
		}
	}
	return nil
}
