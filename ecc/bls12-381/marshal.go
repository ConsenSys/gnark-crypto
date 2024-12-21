// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package bls12381

import (
	"encoding/binary"
	"errors"
	"io"
	"reflect"
	"sync/atomic"

	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/internal/fptower"
	"github.com/consensys/gnark-crypto/internal/parallel"
)

// To encode G1Affine and G2Affine points, we mask the most significant bits with these bits to specify without ambiguity
// metadata needed for point (de)compression
// we follow the BLS12-381 style encoding as specified in ZCash and now IETF
// see https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/11/
// Appendix C.  ZCash serialization format for BLS12_381
const (
	mMask                 byte = 0b111 << 5
	mUncompressed         byte = 0b000 << 5
	_                     byte = 0b001 << 5 // invalid
	mUncompressedInfinity byte = 0b010 << 5
	_                     byte = 0b011 << 5 // invalid
	mCompressedSmallest   byte = 0b100 << 5
	mCompressedLargest    byte = 0b101 << 5
	mCompressedInfinity   byte = 0b110 << 5
	_                     byte = 0b111 << 5 // invalid
)

// SizeOfGT represents the size in bytes that a GT element need in binary form
const SizeOfGT = fptower.SizeOfGT

var (
	ErrInvalidInfinityEncoding = errors.New("invalid infinity point encoding")
	ErrInvalidEncoding         = errors.New("invalid point encoding")
)

// Encoder writes bls12-381 object values to an output stream
type Encoder struct {
	w   io.Writer
	n   int64 // written bytes
	raw bool  // raw vs compressed encoding
}

// Decoder reads bls12-381 object values from an inbound stream
type Decoder struct {
	r             io.Reader
	n             int64 // read bytes
	subGroupCheck bool  // default to true
}

// NewDecoder returns a binary decoder supporting curve bls12-381 objects in both
// compressed and uncompressed (raw) forms
func NewDecoder(r io.Reader, options ...func(*Decoder)) *Decoder {
	d := &Decoder{r: r, subGroupCheck: true}

	for _, o := range options {
		o(d)
	}

	return d
}

// Decode reads the binary encoding of v from the stream
// type must be *uint64, *fr.Element, *fp.Element, *G1Affine, *G2Affine, *[]G1Affine or *[]G2Affine
func (dec *Decoder) Decode(v interface{}) (err error) {
	rv := reflect.ValueOf(v)
	if v == nil || rv.Kind() != reflect.Ptr || rv.IsNil() || !rv.Elem().CanSet() {
		return errors.New("bls12-381 decoder: unsupported type, need pointer")
	}

	// implementation note: code is a bit verbose (abusing code generation), but minimize allocations on the heap
	// in particular, careful attention must be given to usage of Bytes() method on Elements and Points
	// that return an array (not a slice) of bytes. Using this is beneficial to minimize memory allocations
	// in very large (de)serialization upstream in gnark.
	// (but detrimental to code readability here)

	var read64 int64
	if vf, ok := v.(io.ReaderFrom); ok {
		read64, err = vf.ReadFrom(dec.r)
		dec.n += read64
		return
	}

	var buf [SizeOfG2AffineUncompressed]byte
	var read int
	var sliceLen uint32

	switch t := v.(type) {
	case *[][]uint64:
		if sliceLen, err = dec.readUint32(); err != nil {
			return
		}
		*t = make([][]uint64, sliceLen)

		for i := range *t {
			if sliceLen, err = dec.readUint32(); err != nil {
				return
			}
			(*t)[i] = make([]uint64, sliceLen)
			for j := range (*t)[i] {
				if (*t)[i][j], err = dec.readUint64(); err != nil {
					return
				}
			}
		}
		return
	case *[]uint64:
		if sliceLen, err = dec.readUint32(); err != nil {
			return
		}
		*t = make([]uint64, sliceLen)
		for i := range *t {
			if (*t)[i], err = dec.readUint64(); err != nil {
				return
			}
		}
		return
	case *fr.Element:
		read, err = io.ReadFull(dec.r, buf[:fr.Bytes])
		dec.n += int64(read)
		if err != nil {
			return
		}
		err = t.SetBytesCanonical(buf[:fr.Bytes])
		return
	case *fp.Element:
		read, err = io.ReadFull(dec.r, buf[:fp.Bytes])
		dec.n += int64(read)
		if err != nil {
			return
		}
		err = t.SetBytesCanonical(buf[:fp.Bytes])
		return
	case *[]fr.Element:
		read64, err = (*fr.Vector)(t).ReadFrom(dec.r)
		dec.n += read64
		return
	case *[]fp.Element:
		read64, err = (*fp.Vector)(t).ReadFrom(dec.r)
		dec.n += read64
		return
	case *[][]fr.Element:
		if sliceLen, err = dec.readUint32(); err != nil {
			return
		}
		if len(*t) != int(sliceLen) {
			*t = make([][]fr.Element, sliceLen)
		}
		for i := range *t {
			read64, err = (*fr.Vector)(&(*t)[i]).ReadFrom(dec.r)
			dec.n += read64
		}
		return
	case *[][][]fr.Element:
		if sliceLen, err = dec.readUint32(); err != nil {
			return
		}
		if len(*t) != int(sliceLen) {
			*t = make([][][]fr.Element, sliceLen)
		}
		for i := range *t {
			if sliceLen, err = dec.readUint32(); err != nil {
				return
			}
			if len((*t)[i]) != int(sliceLen) {
				(*t)[i] = make([][]fr.Element, sliceLen)
			}
			for j := range (*t)[i] {
				read64, err = (*fr.Vector)(&(*t)[i][j]).ReadFrom(dec.r)
				dec.n += read64
			}
		}
		return
	case *G1Affine:
		// we start by reading compressed point size, if metadata tells us it is uncompressed, we read more.
		read, err = io.ReadFull(dec.r, buf[:SizeOfG1AffineCompressed])
		dec.n += int64(read)
		if err != nil {
			return
		}
		nbBytes := SizeOfG1AffineCompressed

		// 111, 011, 001  --> invalid mask
		if isMaskInvalid(buf[0]) {
			err = ErrInvalidEncoding
			return
		}

		// most significant byte contains metadata
		if !isCompressed(buf[0]) {
			nbBytes = SizeOfG1AffineUncompressed
			// we read more.
			read, err = io.ReadFull(dec.r, buf[SizeOfG1AffineCompressed:SizeOfG1AffineUncompressed])
			dec.n += int64(read)
			if err != nil {
				return
			}
		}
		_, err = t.setBytes(buf[:nbBytes], dec.subGroupCheck)
		return
	case *G2Affine:
		// we start by reading compressed point size, if metadata tells us it is uncompressed, we read more.
		read, err = io.ReadFull(dec.r, buf[:SizeOfG2AffineCompressed])
		dec.n += int64(read)
		if err != nil {
			return
		}
		nbBytes := SizeOfG2AffineCompressed

		// 111, 011, 001  --> invalid mask
		if isMaskInvalid(buf[0]) {
			err = ErrInvalidEncoding
			return
		}

		// most significant byte contains metadata
		if !isCompressed(buf[0]) {
			nbBytes = SizeOfG2AffineUncompressed
			// we read more.
			read, err = io.ReadFull(dec.r, buf[SizeOfG2AffineCompressed:SizeOfG2AffineUncompressed])
			dec.n += int64(read)
			if err != nil {
				return
			}
		}
		_, err = t.setBytes(buf[:nbBytes], dec.subGroupCheck)
		return
	case *[]G1Affine:
		sliceLen, err = dec.readUint32()
		if err != nil {
			return
		}
		if len(*t) != int(sliceLen) || *t == nil {
			*t = make([]G1Affine, sliceLen)
		}
		compressed := make([]bool, sliceLen)
		for i := 0; i < len(*t); i++ {

			// we start by reading compressed point size, if metadata tells us it is uncompressed, we read more.
			read, err = io.ReadFull(dec.r, buf[:SizeOfG1AffineCompressed])
			dec.n += int64(read)
			if err != nil {
				return
			}
			nbBytes := SizeOfG1AffineCompressed

			// 111, 011, 001  --> invalid mask
			if isMaskInvalid(buf[0]) {
				err = ErrInvalidEncoding
				return
			}

			// most significant byte contains metadata
			if !isCompressed(buf[0]) {
				nbBytes = SizeOfG1AffineUncompressed
				// we read more.
				read, err = io.ReadFull(dec.r, buf[SizeOfG1AffineCompressed:SizeOfG1AffineUncompressed])
				dec.n += int64(read)
				if err != nil {
					return
				}
				_, err = (*t)[i].setBytes(buf[:nbBytes], false)
				if err != nil {
					return
				}
			} else {
				var r bool
				if r, err = (*t)[i].unsafeSetCompressedBytes(buf[:nbBytes]); err != nil {
					return
				}
				compressed[i] = !r
			}
		}
		var nbErrs uint64
		parallel.Execute(len(compressed), func(start, end int) {
			for i := start; i < end; i++ {
				if compressed[i] {
					if err := (*t)[i].unsafeComputeY(dec.subGroupCheck); err != nil {
						atomic.AddUint64(&nbErrs, 1)
					}
				} else if dec.subGroupCheck {
					if !(*t)[i].IsInSubGroup() {
						atomic.AddUint64(&nbErrs, 1)
					}
				}
			}
		})
		if nbErrs != 0 {
			return errors.New("point decompression failed")
		}

		return nil
	case *[]G2Affine:
		sliceLen, err = dec.readUint32()
		if err != nil {
			return
		}
		if len(*t) != int(sliceLen) {
			*t = make([]G2Affine, sliceLen)
		}
		compressed := make([]bool, sliceLen)
		for i := 0; i < len(*t); i++ {

			// we start by reading compressed point size, if metadata tells us it is uncompressed, we read more.
			read, err = io.ReadFull(dec.r, buf[:SizeOfG2AffineCompressed])
			dec.n += int64(read)
			if err != nil {
				return
			}
			nbBytes := SizeOfG2AffineCompressed

			// 111, 011, 001  --> invalid mask
			if isMaskInvalid(buf[0]) {
				err = ErrInvalidEncoding
				return
			}

			// most significant byte contains metadata
			if !isCompressed(buf[0]) {
				nbBytes = SizeOfG2AffineUncompressed
				// we read more.
				read, err = io.ReadFull(dec.r, buf[SizeOfG2AffineCompressed:SizeOfG2AffineUncompressed])
				dec.n += int64(read)
				if err != nil {
					return
				}
				_, err = (*t)[i].setBytes(buf[:nbBytes], false)
				if err != nil {
					return
				}
			} else {
				var r bool
				if r, err = (*t)[i].unsafeSetCompressedBytes(buf[:nbBytes]); err != nil {
					return
				}
				compressed[i] = !r
			}
		}
		var nbErrs uint64
		parallel.Execute(len(compressed), func(start, end int) {
			for i := start; i < end; i++ {
				if compressed[i] {
					if err := (*t)[i].unsafeComputeY(dec.subGroupCheck); err != nil {
						atomic.AddUint64(&nbErrs, 1)
					}
				} else if dec.subGroupCheck {
					if !(*t)[i].IsInSubGroup() {
						atomic.AddUint64(&nbErrs, 1)
					}
				}
			}
		})
		if nbErrs != 0 {
			return errors.New("point decompression failed")
		}

		return nil
	default:
		n := binary.Size(t)
		if n == -1 {
			return errors.New("bls12-381 encoder: unsupported type")
		}
		err = binary.Read(dec.r, binary.BigEndian, t)
		if err == nil {
			dec.n += int64(n)
		}
		return
	}
}

// BytesRead return total bytes read from reader
func (dec *Decoder) BytesRead() int64 {
	return dec.n
}

func (dec *Decoder) readUint32() (r uint32, err error) {
	var read int
	var buf [4]byte
	read, err = io.ReadFull(dec.r, buf[:4])
	dec.n += int64(read)
	if err != nil {
		return
	}
	r = binary.BigEndian.Uint32(buf[:4])
	return
}

func (dec *Decoder) readUint64() (r uint64, err error) {
	var read int
	var buf [8]byte
	read, err = io.ReadFull(dec.r, buf[:])
	dec.n += int64(read)
	if err != nil {
		return
	}
	r = binary.BigEndian.Uint64(buf[:])
	return
}

// isMaskInvalid returns true if the mask is invalid
func isMaskInvalid(msb byte) bool {
	mData := msb & mMask
	return ((mData == (0b111 << 5)) || (mData == (0b011 << 5)) || (mData == (0b001 << 5)))
}

func isCompressed(msb byte) bool {
	mData := msb & mMask
	return !((mData == mUncompressed) || (mData == mUncompressedInfinity))
}

// NewEncoder returns a binary encoder supporting curve bls12-381 objects
func NewEncoder(w io.Writer, options ...func(*Encoder)) *Encoder {
	// default settings
	enc := &Encoder{
		w:   w,
		n:   0,
		raw: false,
	}

	// handle options
	for _, option := range options {
		option(enc)
	}

	return enc
}

// Encode writes the binary encoding of v to the stream
// type must be uint64, *fr.Element, *fp.Element, *G1Affine, *G2Affine, []G1Affine, []G2Affine, *[]G1Affine or *[]G2Affine
func (enc *Encoder) Encode(v interface{}) (err error) {
	if enc.raw {
		return enc.encodeRaw(v)
	}
	return enc.encode(v)
}

// BytesWritten return total bytes written on writer
func (enc *Encoder) BytesWritten() int64 {
	return enc.n
}

// RawEncoding returns an option to use in NewEncoder(...) which sets raw encoding mode to true
// points will not be compressed using this option
func RawEncoding() func(*Encoder) {
	return func(enc *Encoder) {
		enc.raw = true
	}
}

// NoSubgroupChecks returns an option to use in NewDecoder(...) which disable subgroup checks on the points
// the decoder will read. Use with caution, as crafted points from an untrusted source can lead to crypto-attacks.
func NoSubgroupChecks() func(*Decoder) {
	return func(dec *Decoder) {
		dec.subGroupCheck = false
	}
}

// isZeroed checks that the provided bytes are at 0
func isZeroed(firstByte byte, buf []byte) bool {
	if firstByte != 0 {
		return false
	}
	for _, b := range buf {
		if b != 0 {
			return false
		}
	}
	return true
}

func (enc *Encoder) encode(v interface{}) (err error) {
	rv := reflect.ValueOf(v)
	if v == nil || (rv.Kind() == reflect.Ptr && rv.IsNil()) {
		return errors.New("<no value> encoder: can't encode <nil>")
	}

	// implementation note: code is a bit verbose (abusing code generation), but minimize allocations on the heap

	var written64 int64
	if vw, ok := v.(io.WriterTo); ok {
		written64, err = vw.WriteTo(enc.w)
		enc.n += written64
		return
	}

	var written int

	switch t := v.(type) {
	case []uint64:
		return enc.writeUint64Slice(t)
	case [][]uint64:
		return enc.writeUint64SliceSlice(t)
	case *fr.Element:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *fp.Element:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *G1Affine:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *G2Affine:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case fr.Vector:
		written64, err = t.WriteTo(enc.w)
		enc.n += written64
		return
	case fp.Vector:
		written64, err = t.WriteTo(enc.w)
		enc.n += written64
		return
	case []fr.Element:
		written64, err = (*fr.Vector)(&t).WriteTo(enc.w)
		enc.n += written64
		return
	case []fp.Element:
		written64, err = (*fp.Vector)(&t).WriteTo(enc.w)
		enc.n += written64
		return
	case [][]fr.Element:
		// write slice length
		if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t))); err != nil {
			return
		}
		enc.n += 4
		for i := range t {
			written64, err = (*fr.Vector)(&t[i]).WriteTo(enc.w)
			enc.n += written64
		}
		return
	case [][][]fr.Element:
		// number of collections
		if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t))); err != nil {
			return
		}
		enc.n += 4
		for i := range t {
			// size of current collection
			if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t[i]))); err != nil {
				return
			}
			enc.n += 4
			// write each vector of the current collection
			for j := range t[i] {
				written64, err = (*fr.Vector)(&t[i][j]).WriteTo(enc.w)
				enc.n += written64
			}
		}
		return
	case *[]G1Affine:
		return enc.encode(*t)
	case []G1Affine:
		// write slice length
		err = binary.Write(enc.w, binary.BigEndian, uint32(len(t)))
		if err != nil {
			return
		}
		enc.n += 4

		var buf [SizeOfG1AffineCompressed]byte

		for i := 0; i < len(t); i++ {
			buf = t[i].Bytes()
			written, err = enc.w.Write(buf[:])
			enc.n += int64(written)
			if err != nil {
				return
			}
		}
		return nil
	case *[]G2Affine:
		return enc.encode(*t)
	case []G2Affine:
		// write slice length
		err = binary.Write(enc.w, binary.BigEndian, uint32(len(t)))
		if err != nil {
			return
		}
		enc.n += 4

		var buf [SizeOfG2AffineCompressed]byte

		for i := 0; i < len(t); i++ {
			buf = t[i].Bytes()
			written, err = enc.w.Write(buf[:])
			enc.n += int64(written)
			if err != nil {
				return
			}
		}
		return nil
	default:
		n := binary.Size(t)
		if n == -1 {
			return errors.New("<no value> encoder: unsupported type")
		}
		err = binary.Write(enc.w, binary.BigEndian, t)
		enc.n += int64(n)
		return
	}
}

func (enc *Encoder) encodeRaw(v interface{}) (err error) {
	rv := reflect.ValueOf(v)
	if v == nil || (rv.Kind() == reflect.Ptr && rv.IsNil()) {
		return errors.New("<no value> encoder: can't encode <nil>")
	}

	// implementation note: code is a bit verbose (abusing code generation), but minimize allocations on the heap

	var written64 int64
	if vw, ok := v.(io.WriterTo); ok {
		written64, err = vw.WriteTo(enc.w)
		enc.n += written64
		return
	}

	var written int

	switch t := v.(type) {
	case []uint64:
		return enc.writeUint64Slice(t)
	case [][]uint64:
		return enc.writeUint64SliceSlice(t)
	case *fr.Element:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *fp.Element:
		buf := t.Bytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *G1Affine:
		buf := t.RawBytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case *G2Affine:
		buf := t.RawBytes()
		written, err = enc.w.Write(buf[:])
		enc.n += int64(written)
		return
	case fr.Vector:
		written64, err = t.WriteTo(enc.w)
		enc.n += written64
		return
	case fp.Vector:
		written64, err = t.WriteTo(enc.w)
		enc.n += written64
		return
	case []fr.Element:
		written64, err = (*fr.Vector)(&t).WriteTo(enc.w)
		enc.n += written64
		return
	case []fp.Element:
		written64, err = (*fp.Vector)(&t).WriteTo(enc.w)
		enc.n += written64
		return
	case [][]fr.Element:
		// write slice length
		if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t))); err != nil {
			return
		}
		enc.n += 4
		for i := range t {
			written64, err = (*fr.Vector)(&t[i]).WriteTo(enc.w)
			enc.n += written64
		}
		return
	case [][][]fr.Element:
		// number of collections
		if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t))); err != nil {
			return
		}
		enc.n += 4
		for i := range t {
			// size of current collection
			if err = binary.Write(enc.w, binary.BigEndian, uint32(len(t[i]))); err != nil {
				return
			}
			enc.n += 4
			// write each vector of the current collection
			for j := range t[i] {
				written64, err = (*fr.Vector)(&t[i][j]).WriteTo(enc.w)
				enc.n += written64
			}
		}
		return
	case *[]G1Affine:
		return enc.encodeRaw(*t)
	case []G1Affine:
		// write slice length
		err = binary.Write(enc.w, binary.BigEndian, uint32(len(t)))
		if err != nil {
			return
		}
		enc.n += 4

		var buf [SizeOfG1AffineUncompressed]byte

		for i := 0; i < len(t); i++ {
			buf = t[i].RawBytes()
			written, err = enc.w.Write(buf[:])
			enc.n += int64(written)
			if err != nil {
				return
			}
		}
		return nil
	case *[]G2Affine:
		return enc.encodeRaw(*t)
	case []G2Affine:
		// write slice length
		err = binary.Write(enc.w, binary.BigEndian, uint32(len(t)))
		if err != nil {
			return
		}
		enc.n += 4

		var buf [SizeOfG2AffineUncompressed]byte

		for i := 0; i < len(t); i++ {
			buf = t[i].RawBytes()
			written, err = enc.w.Write(buf[:])
			enc.n += int64(written)
			if err != nil {
				return
			}
		}
		return nil
	default:
		n := binary.Size(t)
		if n == -1 {
			return errors.New("<no value> encoder: unsupported type")
		}
		err = binary.Write(enc.w, binary.BigEndian, t)
		enc.n += int64(n)
		return
	}
}

func (enc *Encoder) writeUint64Slice(t []uint64) (err error) {
	if err = enc.writeUint32(uint32(len(t))); err != nil {
		return
	}
	for i := range t {
		if err = enc.writeUint64(t[i]); err != nil {
			return
		}
	}
	return nil
}

func (enc *Encoder) writeUint64SliceSlice(t [][]uint64) (err error) {
	if err = enc.writeUint32(uint32(len(t))); err != nil {
		return
	}
	for i := range t {
		if err = enc.writeUint32(uint32(len(t[i]))); err != nil {
			return
		}
		for j := range t[i] {
			if err = enc.writeUint64(t[i][j]); err != nil {
				return
			}
		}
	}
	return nil
}

func (enc *Encoder) writeUint64(a uint64) error {
	var buff [64 / 8]byte
	binary.BigEndian.PutUint64(buff[:], a)
	written, err := enc.w.Write(buff[:])
	enc.n += int64(written)
	return err
}

func (enc *Encoder) writeUint32(a uint32) error {
	var buff [32 / 8]byte
	binary.BigEndian.PutUint32(buff[:], a)
	written, err := enc.w.Write(buff[:])
	enc.n += int64(written)
	return err
}

// SizeOfG1AffineCompressed represents the size in bytes that a G1Affine need in binary form, compressed
const SizeOfG1AffineCompressed = 48

// SizeOfG1AffineUncompressed represents the size in bytes that a G1Affine need in binary form, uncompressed
const SizeOfG1AffineUncompressed = SizeOfG1AffineCompressed * 2

// Marshal converts p to a byte slice (without point compression)
func (p *G1Affine) Marshal() []byte {
	b := p.RawBytes()
	return b[:]
}

// Unmarshal is an alias to SetBytes()
func (p *G1Affine) Unmarshal(buf []byte) error {
	_, err := p.SetBytes(buf)
	return err
}

// Bytes returns binary representation of p
// will store X coordinate in regular form and a parity bit
// we follow the BLS12-381 style encoding as specified in ZCash and now IETF
//
// The most significant bit, when set, indicates that the point is in compressed form. Otherwise, the point is in uncompressed form.
//
// The second-most significant bit indicates that the point is at infinity. If this bit is set, the remaining bits of the group element's encoding should be set to zero.
//
// The third-most significant bit is set if (and only if) this point is in compressed form and it is not the point at infinity and its y-coordinate is the lexicographically largest of the two associated with the encoded x-coordinate.
func (p *G1Affine) Bytes() (res [SizeOfG1AffineCompressed]byte) {

	// check if p is infinity point
	if p.X.IsZero() && p.Y.IsZero() {
		res[0] = mCompressedInfinity
		return
	}

	msbMask := mCompressedSmallest
	// compressed, we need to know if Y is lexicographically bigger than -Y
	// if p.Y ">" -p.Y
	if p.Y.LexicographicallyLargest() {
		msbMask = mCompressedLargest
	}

	// we store X  and mask the most significant word with our metadata mask
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[0:0+fp.Bytes]), p.X)

	res[0] |= msbMask

	return
}

// RawBytes returns binary representation of p (stores X and Y coordinate)
// see Bytes() for a compressed representation
func (p *G1Affine) RawBytes() (res [SizeOfG1AffineUncompressed]byte) {

	// check if p is infinity point
	if p.X.IsZero() && p.Y.IsZero() {

		res[0] = mUncompressedInfinity

		return
	}

	// not compressed
	// we store the Y coordinate
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[48:48+fp.Bytes]), p.Y)

	// we store X  and mask the most significant word with our metadata mask
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[0:0+fp.Bytes]), p.X)

	res[0] |= mUncompressed

	return
}

// SetBytes sets p from binary representation in buf and returns number of consumed bytes
//
// bytes in buf must match either RawBytes() or Bytes() output
//
// if buf is too short io.ErrShortBuffer is returned
//
// if buf contains compressed representation (output from Bytes()) and we're unable to compute
// the Y coordinate (i.e the square root doesn't exist) this function returns an error
//
// this check if the resulting point is on the curve and in the correct subgroup
func (p *G1Affine) SetBytes(buf []byte) (int, error) {
	return p.setBytes(buf, true)
}

func (p *G1Affine) setBytes(buf []byte, subGroupCheck bool) (int, error) {
	if len(buf) < SizeOfG1AffineCompressed {
		return 0, io.ErrShortBuffer
	}

	// most significant byte
	mData := buf[0] & mMask

	// 111, 011, 001  --> invalid mask
	if isMaskInvalid(mData) {
		return 0, ErrInvalidEncoding
	}

	// check buffer size
	if (mData == mUncompressed) || (mData == mUncompressedInfinity) {
		if len(buf) < SizeOfG1AffineUncompressed {
			return 0, io.ErrShortBuffer
		}
	}

	// infinity encoded, we still check that the buffer is full of zeroes.
	if mData == mCompressedInfinity {
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG1AffineCompressed]) {
			return 0, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return SizeOfG1AffineCompressed, nil
	}
	if mData == mUncompressedInfinity {
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG1AffineUncompressed]) {
			return 0, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return SizeOfG1AffineUncompressed, nil
	}

	// uncompressed point
	if mData == mUncompressed {
		// read X and Y coordinates
		if err := p.X.SetBytesCanonical(buf[:fp.Bytes]); err != nil {
			return 0, err
		}
		if err := p.Y.SetBytesCanonical(buf[fp.Bytes : fp.Bytes*2]); err != nil {
			return 0, err
		}

		// subgroup check
		if subGroupCheck && !p.IsInSubGroup() {
			return 0, errors.New("invalid point: subgroup check failed")
		}

		return SizeOfG1AffineUncompressed, nil
	}

	// we have a compressed coordinate
	// we need to
	// 	1. copy the buffer (to keep this method thread safe)
	// 	2. we need to solve the curve equation to compute Y

	var bufX [fp.Bytes]byte
	copy(bufX[:fp.Bytes], buf[:fp.Bytes])
	bufX[0] &= ^mMask

	// read X coordinate
	if err := p.X.SetBytesCanonical(bufX[:fp.Bytes]); err != nil {
		return 0, err
	}

	var YSquared, Y fp.Element

	YSquared.Square(&p.X).Mul(&YSquared, &p.X)
	YSquared.Add(&YSquared, &bCurveCoeff)
	if Y.Sqrt(&YSquared) == nil {
		return 0, errors.New("invalid compressed coordinate: square root doesn't exist")
	}

	if Y.LexicographicallyLargest() {
		// Y ">" -Y
		if mData == mCompressedSmallest {
			Y.Neg(&Y)
		}
	} else {
		// Y "<=" -Y
		if mData == mCompressedLargest {
			Y.Neg(&Y)
		}
	}

	p.Y = Y

	// subgroup check
	if subGroupCheck && !p.IsInSubGroup() {
		return 0, errors.New("invalid point: subgroup check failed")
	}

	return SizeOfG1AffineCompressed, nil
}

// unsafeComputeY called by Decoder when processing slices of compressed point in parallel (step 2)
// it computes the Y coordinate from the already set X coordinate and is compute intensive
func (p *G1Affine) unsafeComputeY(subGroupCheck bool) error {
	// stored in unsafeSetCompressedBytes

	mData := byte(p.Y[0])

	// we have a compressed coordinate, we need to solve the curve equation to compute Y
	var YSquared, Y fp.Element

	YSquared.Square(&p.X).Mul(&YSquared, &p.X)
	YSquared.Add(&YSquared, &bCurveCoeff)
	if Y.Sqrt(&YSquared) == nil {
		return errors.New("invalid compressed coordinate: square root doesn't exist")
	}

	if Y.LexicographicallyLargest() {
		// Y ">" -Y
		if mData == mCompressedSmallest {
			Y.Neg(&Y)
		}
	} else {
		// Y "<=" -Y
		if mData == mCompressedLargest {
			Y.Neg(&Y)
		}
	}

	p.Y = Y

	// subgroup check
	if subGroupCheck && !p.IsInSubGroup() {
		return errors.New("invalid point: subgroup check failed")
	}

	return nil
}

// unsafeSetCompressedBytes is called by Decoder when processing slices of compressed point in parallel (step 1)
// assumes buf[:8] mask is set to compressed
// returns true if point is infinity and need no further processing
// it sets X coordinate and uses Y for scratch space to store decompression metadata
func (p *G1Affine) unsafeSetCompressedBytes(buf []byte) (isInfinity bool, err error) {

	// read the most significant byte
	mData := buf[0] & mMask

	if mData == mCompressedInfinity {
		isInfinity = true
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG1AffineCompressed]) {
			return isInfinity, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return isInfinity, nil
	}

	// we need to copy the input buffer (to keep this method thread safe)
	var bufX [fp.Bytes]byte
	copy(bufX[:fp.Bytes], buf[:fp.Bytes])
	bufX[0] &= ^mMask

	// read X coordinate
	if err := p.X.SetBytesCanonical(bufX[:fp.Bytes]); err != nil {
		return false, err
	}
	// store mData in p.Y[0]
	p.Y[0] = uint64(mData)

	// recomputing Y will be done asynchronously
	return isInfinity, nil
}

// SizeOfG2AffineCompressed represents the size in bytes that a G2Affine need in binary form, compressed
const SizeOfG2AffineCompressed = 48 * 2

// SizeOfG2AffineUncompressed represents the size in bytes that a G2Affine need in binary form, uncompressed
const SizeOfG2AffineUncompressed = SizeOfG2AffineCompressed * 2

// Marshal converts p to a byte slice (without point compression)
func (p *G2Affine) Marshal() []byte {
	b := p.RawBytes()
	return b[:]
}

// Unmarshal is an alias to SetBytes()
func (p *G2Affine) Unmarshal(buf []byte) error {
	_, err := p.SetBytes(buf)
	return err
}

// Bytes returns binary representation of p
// will store X coordinate in regular form and a parity bit
// we follow the BLS12-381 style encoding as specified in ZCash and now IETF
//
// The most significant bit, when set, indicates that the point is in compressed form. Otherwise, the point is in uncompressed form.
//
// The second-most significant bit indicates that the point is at infinity. If this bit is set, the remaining bits of the group element's encoding should be set to zero.
//
// The third-most significant bit is set if (and only if) this point is in compressed form and it is not the point at infinity and its y-coordinate is the lexicographically largest of the two associated with the encoded x-coordinate.
func (p *G2Affine) Bytes() (res [SizeOfG2AffineCompressed]byte) {

	// check if p is infinity point
	if p.X.IsZero() && p.Y.IsZero() {
		res[0] = mCompressedInfinity
		return
	}

	msbMask := mCompressedSmallest
	// compressed, we need to know if Y is lexicographically bigger than -Y
	// if p.Y ">" -p.Y
	if p.Y.LexicographicallyLargest() {
		msbMask = mCompressedLargest
	}

	// we store X  and mask the most significant word with our metadata mask
	// p.X.A1 | p.X.A0
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[48:48+fp.Bytes]), p.X.A0)
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[0:0+fp.Bytes]), p.X.A1)

	res[0] |= msbMask

	return
}

// RawBytes returns binary representation of p (stores X and Y coordinate)
// see Bytes() for a compressed representation
func (p *G2Affine) RawBytes() (res [SizeOfG2AffineUncompressed]byte) {

	// check if p is infinity point
	if p.X.IsZero() && p.Y.IsZero() {

		res[0] = mUncompressedInfinity

		return
	}

	// not compressed
	// we store the Y coordinate
	// p.Y.A1 | p.Y.A0
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[144:144+fp.Bytes]), p.Y.A0)
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[96:96+fp.Bytes]), p.Y.A1)

	// we store X  and mask the most significant word with our metadata mask
	// p.X.A1 | p.X.A0
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[0:0+fp.Bytes]), p.X.A1)
	fp.BigEndian.PutElement((*[fp.Bytes]byte)(res[48:48+fp.Bytes]), p.X.A0)

	res[0] |= mUncompressed

	return
}

// SetBytes sets p from binary representation in buf and returns number of consumed bytes
//
// bytes in buf must match either RawBytes() or Bytes() output
//
// if buf is too short io.ErrShortBuffer is returned
//
// if buf contains compressed representation (output from Bytes()) and we're unable to compute
// the Y coordinate (i.e the square root doesn't exist) this function returns an error
//
// this check if the resulting point is on the curve and in the correct subgroup
func (p *G2Affine) SetBytes(buf []byte) (int, error) {
	return p.setBytes(buf, true)
}

func (p *G2Affine) setBytes(buf []byte, subGroupCheck bool) (int, error) {
	if len(buf) < SizeOfG2AffineCompressed {
		return 0, io.ErrShortBuffer
	}

	// most significant byte
	mData := buf[0] & mMask

	// 111, 011, 001  --> invalid mask
	if isMaskInvalid(mData) {
		return 0, ErrInvalidEncoding
	}

	// check buffer size
	if (mData == mUncompressed) || (mData == mUncompressedInfinity) {
		if len(buf) < SizeOfG2AffineUncompressed {
			return 0, io.ErrShortBuffer
		}
	}

	// infinity encoded, we still check that the buffer is full of zeroes.
	if mData == mCompressedInfinity {
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG2AffineCompressed]) {
			return 0, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return SizeOfG2AffineCompressed, nil
	}
	if mData == mUncompressedInfinity {
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG2AffineUncompressed]) {
			return 0, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return SizeOfG2AffineUncompressed, nil
	}

	// uncompressed point
	if mData == mUncompressed {
		// read X and Y coordinates
		// p.X.A1 | p.X.A0
		if err := p.X.A1.SetBytesCanonical(buf[:fp.Bytes]); err != nil {
			return 0, err
		}
		if err := p.X.A0.SetBytesCanonical(buf[fp.Bytes : fp.Bytes*2]); err != nil {
			return 0, err
		}
		// p.Y.A1 | p.Y.A0
		if err := p.Y.A1.SetBytesCanonical(buf[fp.Bytes*2 : fp.Bytes*3]); err != nil {
			return 0, err
		}
		if err := p.Y.A0.SetBytesCanonical(buf[fp.Bytes*3 : fp.Bytes*4]); err != nil {
			return 0, err
		}

		// subgroup check
		if subGroupCheck && !p.IsInSubGroup() {
			return 0, errors.New("invalid point: subgroup check failed")
		}

		return SizeOfG2AffineUncompressed, nil
	}

	// we have a compressed coordinate
	// we need to
	// 	1. copy the buffer (to keep this method thread safe)
	// 	2. we need to solve the curve equation to compute Y

	var bufX [fp.Bytes]byte
	copy(bufX[:fp.Bytes], buf[:fp.Bytes])
	bufX[0] &= ^mMask

	// read X coordinate
	// p.X.A1 | p.X.A0
	if err := p.X.A1.SetBytesCanonical(bufX[:fp.Bytes]); err != nil {
		return 0, err
	}
	if err := p.X.A0.SetBytesCanonical(buf[fp.Bytes : fp.Bytes*2]); err != nil {
		return 0, err
	}

	var YSquared, Y fptower.E2

	YSquared.Square(&p.X).Mul(&YSquared, &p.X)
	YSquared.Add(&YSquared, &bTwistCurveCoeff)
	if YSquared.Legendre() == -1 {
		return 0, errors.New("invalid compressed coordinate: square root doesn't exist")
	}
	Y.Sqrt(&YSquared)

	if Y.LexicographicallyLargest() {
		// Y ">" -Y
		if mData == mCompressedSmallest {
			Y.Neg(&Y)
		}
	} else {
		// Y "<=" -Y
		if mData == mCompressedLargest {
			Y.Neg(&Y)
		}
	}

	p.Y = Y

	// subgroup check
	if subGroupCheck && !p.IsInSubGroup() {
		return 0, errors.New("invalid point: subgroup check failed")
	}

	return SizeOfG2AffineCompressed, nil
}

// unsafeComputeY called by Decoder when processing slices of compressed point in parallel (step 2)
// it computes the Y coordinate from the already set X coordinate and is compute intensive
func (p *G2Affine) unsafeComputeY(subGroupCheck bool) error {
	// stored in unsafeSetCompressedBytes

	mData := byte(p.Y.A0[0])

	// we have a compressed coordinate, we need to solve the curve equation to compute Y
	var YSquared, Y fptower.E2

	YSquared.Square(&p.X).Mul(&YSquared, &p.X)
	YSquared.Add(&YSquared, &bTwistCurveCoeff)
	if YSquared.Legendre() == -1 {
		return errors.New("invalid compressed coordinate: square root doesn't exist")
	}
	Y.Sqrt(&YSquared)

	if Y.LexicographicallyLargest() {
		// Y ">" -Y
		if mData == mCompressedSmallest {
			Y.Neg(&Y)
		}
	} else {
		// Y "<=" -Y
		if mData == mCompressedLargest {
			Y.Neg(&Y)
		}
	}

	p.Y = Y

	// subgroup check
	if subGroupCheck && !p.IsInSubGroup() {
		return errors.New("invalid point: subgroup check failed")
	}

	return nil
}

// unsafeSetCompressedBytes is called by Decoder when processing slices of compressed point in parallel (step 1)
// assumes buf[:8] mask is set to compressed
// returns true if point is infinity and need no further processing
// it sets X coordinate and uses Y for scratch space to store decompression metadata
func (p *G2Affine) unsafeSetCompressedBytes(buf []byte) (isInfinity bool, err error) {

	// read the most significant byte
	mData := buf[0] & mMask

	if mData == mCompressedInfinity {
		isInfinity = true
		if !isZeroed(buf[0] & ^mMask, buf[1:SizeOfG2AffineCompressed]) {
			return isInfinity, ErrInvalidInfinityEncoding
		}
		p.X.SetZero()
		p.Y.SetZero()
		return isInfinity, nil
	}

	// we need to copy the input buffer (to keep this method thread safe)
	var bufX [fp.Bytes]byte
	copy(bufX[:fp.Bytes], buf[:fp.Bytes])
	bufX[0] &= ^mMask

	// read X coordinate
	// p.X.A1 | p.X.A0
	if err := p.X.A1.SetBytesCanonical(bufX[:fp.Bytes]); err != nil {
		return false, err
	}
	if err := p.X.A0.SetBytesCanonical(buf[fp.Bytes : fp.Bytes*2]); err != nil {
		return false, err
	}

	// store mData in p.Y.A0[0]
	p.Y.A0[0] = uint64(mData)

	// recomputing Y will be done asynchronously
	return isInfinity, nil
}
