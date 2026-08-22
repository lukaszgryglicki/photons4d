package photons4d

// Optional wire compression for sparse voxel updates.
//
// A packed update batch is a DEFLATE stream of:
//
//	uvarint n                     — entry count
//	n × uvarint                   — index deltas (first is the absolute index)
//	n × 8 bytes                   — float64 values, little-endian, byte-transposed
//
// Index deltas are small repetitive integers, which DEFLATE shrinks heavily.
// Values are Monte-Carlo sums whose mantissa bytes are near-random, but the
// sign/exponent bytes of neighbouring voxels are strongly correlated;
// byte-transposing (all 0th bytes, then all 1st bytes, …) groups those
// correlated bytes into runs DEFLATE can exploit. Decoding reverses both
// steps bit-exactly — compression never changes any merged value.

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"fmt"
	"io"
	"math"
)

// Codec identifiers for UpdateMsg.Codec.
const (
	CodecRaw   = 0 // IdxDeltas/Values travel as plain gob fields
	CodecFlate = 1 // Packed carries the DEFLATE stream described above
)

// packSparse serializes one (idxDeltas, values) batch into the packed format.
func packSparse(idxDeltas []int64, values []float64) ([]byte, error) {
	if len(idxDeltas) != len(values) {
		return nil, fmt.Errorf("packSparse: %d indices vs %d values", len(idxDeltas), len(values))
	}
	n := len(values)
	raw := make([]byte, 0, binary.MaxVarintLen64*(n+1)+8*n)
	var tmp [binary.MaxVarintLen64]byte
	raw = append(raw, tmp[:binary.PutUvarint(tmp[:], uint64(n))]...)
	for _, d := range idxDeltas {
		if d < 0 {
			return nil, fmt.Errorf("packSparse: negative index delta %d", d)
		}
		raw = append(raw, tmp[:binary.PutUvarint(tmp[:], uint64(d))]...)
	}
	// Byte-transpose the values.
	vOff := len(raw)
	raw = append(raw, make([]byte, 8*n)...)
	tr := raw[vOff:]
	for i, v := range values {
		bits := math.Float64bits(v)
		for b := 0; b < 8; b++ {
			tr[b*n+i] = byte(bits >> (8 * b))
		}
	}
	var out bytes.Buffer
	out.Grow(len(raw) / 2)
	fw, err := flate.NewWriter(&out, flate.BestSpeed)
	if err != nil {
		return nil, err
	}
	if _, err := fw.Write(raw); err != nil {
		return nil, err
	}
	if err := fw.Close(); err != nil {
		return nil, err
	}
	return out.Bytes(), nil
}

// unpackSparse reverses packSparse. rawLen is the uncompressed payload size
// (for compression-ratio accounting). The result is bit-exact.
func unpackSparse(packed []byte) (idxDeltas []int64, values []float64, rawLen int, err error) {
	fr := flate.NewReader(bytes.NewReader(packed))
	raw, err := io.ReadAll(fr)
	_ = fr.Close()
	if err != nil {
		return nil, nil, 0, fmt.Errorf("unpackSparse: inflate: %w", err)
	}
	off := 0
	un, w := binary.Uvarint(raw[off:])
	if w <= 0 {
		return nil, nil, 0, fmt.Errorf("unpackSparse: bad entry count")
	}
	off += w
	if un > uint64(len(raw)) { // each entry needs ≥1 delta byte + 8 value bytes
		return nil, nil, 0, fmt.Errorf("unpackSparse: entry count %d exceeds payload", un)
	}
	n := int(un)
	idxDeltas = make([]int64, n)
	values = make([]float64, n)
	for i := 0; i < n; i++ {
		d, w := binary.Uvarint(raw[off:])
		if w <= 0 {
			return nil, nil, 0, fmt.Errorf("unpackSparse: truncated delta %d/%d", i, n)
		}
		if d > math.MaxInt64 {
			return nil, nil, 0, fmt.Errorf("unpackSparse: delta %d overflows int64", d)
		}
		idxDeltas[i] = int64(d)
		off += w
	}
	if len(raw)-off != 8*n {
		return nil, nil, 0, fmt.Errorf("unpackSparse: value block is %d bytes, want %d", len(raw)-off, 8*n)
	}
	tr := raw[off:]
	for i := 0; i < n; i++ {
		var bits uint64
		for b := 0; b < 8; b++ {
			bits |= uint64(tr[b*n+i]) << (8 * b)
		}
		values[i] = math.Float64frombits(bits)
	}
	return idxDeltas, values, len(raw), nil
}

// decodeUpdate materializes a packed UpdateMsg in place so downstream code
// (buffering, applySparse validation, merging) only ever sees plain
// IdxDeltas/Values. Returns wire/raw payload sizes for ratio accounting
// (0,0 for CodecRaw messages).
func decodeUpdate(upd *UpdateMsg) (wireLen, rawLen int, err error) {
	switch upd.Codec {
	case CodecRaw:
		return 0, 0, nil
	case CodecFlate:
		if len(upd.Packed) == 0 {
			return 0, 0, nil
		}
		wireLen = len(upd.Packed)
		upd.IdxDeltas, upd.Values, rawLen, err = unpackSparse(upd.Packed)
		upd.Packed = nil
		if err != nil {
			return 0, 0, err
		}
		return wireLen, rawLen, nil
	default:
		return 0, 0, fmt.Errorf("unknown update codec %d", upd.Codec)
	}
}
