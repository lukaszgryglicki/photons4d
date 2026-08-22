package photons4d

// Distributed (horizontal) rendering protocol.
//
// Topology: one server owns the master voxel buffer; N clients connect over
// TCP, prove (via SHA256 of the *effective* config) that they loaded exactly
// the same scene, and then repeatedly ask for work. A unit of work is a
// per-light ray-count vector. The client casts those rays with its full local
// worker pool (identical code path to vertical mode), then ships back the
// sparse non-zero voxel deltas and zeroes its local buffer. Because photon
// deposits are purely additive, merging sparse deltas on the server is
// mathematically exact — the master buffer is bit-for-bit a valid sum of all
// casts, indistinguishable from a single-machine run with the same rays.
//
// The wire protocol is gob over TCP with strict request/response alternation:
//
//	client                         server
//	  HelloMsg          ------>
//	                    <------    HelloAck        (hash verified here)
//	  WorkRequest       ------>
//	                    <------    WorkAssign      (Done=true terminates)
//	  UpdateMsg(More=t) ------>                    (0..n chunk batches)
//	  UpdateMsg(More=f) ------>
//	                    <------    UpdateAck
//	  WorkRequest       ------>    ...
//
// A round's batches are buffered server-side and applied to the master
// buffer transactionally only once the final batch arrives, so a client
// dying mid-upload never double-counts (its round is simply re-assigned).

import (
	"bufio"
	"crypto/sha256"
	"encoding/gob"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"net"
)

// DistProtocolVersion guards against mixing incompatible binaries.
// v2: optional per-message compression (UpdateMsg.Codec/Packed) and
// client-configurable update batch size.
const DistProtocolVersion = 2

// updateBatchEntries is the default cap on entries per UpdateMsg so neither
// side ever has to hold a whole huge round in one gob message
// (4M entries ≈ 64 MB raw). Clients may override it via -batch/BATCH.
const updateBatchEntries = 4 << 20

// HelloMsg is sent by a client immediately after connecting.
type HelloMsg struct {
	Version   int    // DistProtocolVersion
	SceneHash string // hex SHA256 of the canonical effective config
	BufLen    int    // len(scene.Buf) — belt-and-braces dimension check
	NLights   int
	Workers   int    // client worker count (informational)
	Host      string // client hostname (informational)
}

// HelloAck is the server's verdict on a HelloMsg.
type HelloAck struct {
	OK       bool
	Reason   string // set when !OK
	ClientID int
}

// WorkRequest asks the server for the next ray-casting round.
type WorkRequest struct {
	ClientID int
}

// WorkAssign hands a client its next round.
// Done=true: all work finished, client should exit.
// All-zero RaysPerLight with Done=false: no work available right now
// (stragglers still out) — client should wait briefly and re-ask.
type WorkAssign struct {
	Done         bool
	RoundID      int
	RaysPerLight []int
}

// UpdateMsg carries one batch of sparse voxel deltas for a round.
// Indices are delta-encoded in strictly ascending order: IdxDeltas[0] is the
// absolute flat buffer index of the first entry *of this batch*, every
// following element is the positive distance to the previous index.
//
// Codec selects the transport encoding: CodecRaw sends IdxDeltas/Values as
// plain gob fields; CodecFlate sends the same data DEFLATE-compressed in
// Packed (see dist_pack.go) and leaves IdxDeltas/Values empty. Both encodings
// are bit-exact — the server materializes packed batches back into
// IdxDeltas/Values before validation and merging.
type UpdateMsg struct {
	ClientID  int
	RoundID   int
	More      bool // true: further batches of this round follow
	Codec     uint8
	IdxDeltas []int64
	Values    []float64
	Packed    []byte
}

// UpdateAck confirms a fully-merged round.
type UpdateAck struct {
	OK     bool
	Reason string
	// Progress info so clients can log cluster-wide status.
	RaysDone  int64
	RaysTotal int64
}

// wire wraps a TCP connection with buffered gob encoding/decoding.
type wire struct {
	conn net.Conn
	bw   *bufio.Writer
	enc  *gob.Encoder
	dec  *gob.Decoder
}

func newWire(conn net.Conn) *wire {
	bw := bufio.NewWriterSize(conn, 1<<20)
	return &wire{
		conn: conn,
		bw:   bw,
		enc:  gob.NewEncoder(bw),
		dec:  gob.NewDecoder(bufio.NewReaderSize(conn, 1<<20)),
	}
}

func (w *wire) send(v interface{}) error {
	if err := w.enc.Encode(v); err != nil {
		return err
	}
	return w.bw.Flush()
}

func (w *wire) recv(v interface{}) error {
	return w.dec.Decode(v)
}

func (w *wire) Close() error {
	return w.conn.Close()
}

// SceneHash returns the hex SHA256 of the canonical JSON serialization of the
// *effective* config, i.e. after defaults, env overrides (FORCE_ESCAPE,
// SPP_ADJUST) and debug reductions have been applied by prepareConfig. Both
// sides compute it the same way, so any divergence in scene geometry,
// resolution, lights, spp or materials yields a different hash and the client
// is rejected before any voxel data flows.
func SceneHash(cfg *Config) (string, error) {
	b, err := json.Marshal(cfg)
	if err != nil {
		return "", fmt.Errorf("scene hash: %w", err)
	}
	sum := sha256.Sum256(b)
	return hex.EncodeToString(sum[:]), nil
}

// extractSparseAndZero scans buf, appends every non-zero entry as a
// delta-encoded (index, value) pair and zeroes it, leaving buf ready to
// accumulate the next round. Returns nil, nil when the buffer is all zero.
func extractSparseAndZero(buf []Real) ([]int64, []float64) {
	var idxDeltas []int64
	var values []float64
	prev := int64(-1)
	for i := range buf {
		v := buf[i]
		if v == 0 {
			continue
		}
		if prev < 0 {
			idxDeltas = append(idxDeltas, int64(i))
		} else {
			idxDeltas = append(idxDeltas, int64(i)-prev)
		}
		prev = int64(i)
		values = append(values, float64(v))
		buf[i] = 0
	}
	return idxDeltas, values
}

// applySparse adds one delta-encoded batch into buf, validating bounds,
// monotonicity and value sanity (finite, non-negative — photon deposits can
// never be negative). Any violation aborts the merge with an error before
// buf is touched.
func applySparse(buf []Real, idxDeltas []int64, values []float64) error {
	if len(idxDeltas) != len(values) {
		return fmt.Errorf("sparse update: %d indices vs %d values", len(idxDeltas), len(values))
	}
	// Validate first (don't partially apply a corrupt batch).
	idx := int64(-1)
	for n, d := range idxDeltas {
		if n == 0 {
			idx = d
		} else {
			if d <= 0 {
				return fmt.Errorf("sparse update: non-ascending index delta %d at %d", d, n)
			}
			idx += d
		}
		if idx < 0 || idx >= int64(len(buf)) {
			return fmt.Errorf("sparse update: index %d out of range [0,%d)", idx, len(buf))
		}
		v := values[n]
		if math.IsNaN(v) || math.IsInf(v, 0) || v < 0 {
			return fmt.Errorf("sparse update: invalid value %v at %d", v, n)
		}
	}
	idx = -1
	for n, d := range idxDeltas {
		if n == 0 {
			idx = d
		} else {
			idx += d
		}
		buf[idx] += Real(values[n])
	}
	return nil
}
