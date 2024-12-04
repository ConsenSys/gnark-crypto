// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package mimc

import (
	"errors"
	stdhash "hash"
	"math/big"
	"sync"

	"github.com/consensys/gnark-crypto/ecc/bls24-315/fr"
	"github.com/consensys/gnark-crypto/hash"

	"golang.org/x/crypto/sha3"
)

func init() {
	hash.RegisterHash(hash.MIMC_BLS24_315, func() stdhash.Hash {
		return NewMiMC()
	})
}

const (
	mimcNbRounds = 109
	seed         = "seed"   // seed to derive the constants
	BlockSize    = fr.Bytes // BlockSize size that mimc consumes
)

// Params constants for the mimc hash function
var (
	mimcConstants [mimcNbRounds]fr.Element
	once          sync.Once
)

// digest represents the partial evaluation of the checksum
// along with the params of the mimc function
type digest struct {
	h         fr.Element
	data      []fr.Element // data to hash
	byteOrder fr.ByteOrder
}

// GetConstants exposed to be used in gnark
func GetConstants() []big.Int {
	once.Do(initConstants) // init constants
	res := make([]big.Int, mimcNbRounds)
	for i := 0; i < mimcNbRounds; i++ {
		mimcConstants[i].BigInt(&res[i])
	}
	return res
}

// NewMiMC returns a MiMC implementation, pure Go reference implementation. The
// returned instance also implements [hash.StateStorer] for recovering the internal
// state of the hasher.
func NewMiMC(opts ...Option) stdhash.Hash {
	d := new(digest)
	d.Reset()
	cfg := mimcOptions(opts...)
	d.byteOrder = cfg.byteOrder
	return d
}

// Reset resets the Hash to its initial state.
func (d *digest) Reset() {
	d.data = d.data[:0]
	d.h = fr.Element{0, 0, 0, 0}
}

// Sum appends the current hash to b and returns the resulting slice.
// It does not change the underlying hash state.
func (d *digest) Sum(b []byte) []byte {
	buffer := d.checksum()
	d.data = nil // flush the data already hashed
	hash := buffer.Bytes()
	b = append(b, hash[:]...)
	return b
}

// BlockSize returns the hash's underlying block size.
// The Write method must be able to accept any amount
// of data, but it may operate more efficiently if all writes
// are a multiple of the block size.
func (d *digest) Size() int {
	return BlockSize
}

// BlockSize returns the number of bytes Sum will return.
func (d *digest) BlockSize() int {
	return BlockSize
}

// Write (via the embedded io.Writer interface) adds more data to the running hash.
//
// Each []byte block of size BlockSize represents a big endian fr.Element.
//
// If len(p) is not a multiple of BlockSize and any of the []byte in p represent an integer
// larger than fr.Modulus, this function returns an error.
//
// To hash arbitrary data ([]byte not representing canonical field elements) use fr.Hash first
func (d *digest) Write(p []byte) (int, error) {
	// we usually expect multiple of block size. But sometimes we hash short
	// values (FS transcript). Instead of forcing to hash to field, we left-pad the
	// input here.
	if len(p) > 0 && len(p) < BlockSize {
		pp := make([]byte, BlockSize)
		copy(pp[len(pp)-len(p):], p)
		p = pp
	}

	var start int
	for start = 0; start < len(p); start += BlockSize {
		if elem, err := d.byteOrder.Element((*[BlockSize]byte)(p[start : start+BlockSize])); err == nil {
			d.data = append(d.data, elem)
		} else {
			return 0, err
		}
	}

	if start != len(p) {
		return 0, errors.New("invalid input length: must represent a list of field elements, expects a []byte of len m*BlockSize")
	}
	return len(p), nil
}

// Hash hash using Miyaguchi-Preneel:
// https://en.wikipedia.org/wiki/One-way_compression_function
// The XOR operation is replaced by field addition, data is in Montgomery form
func (d *digest) checksum() fr.Element {
	// Write guarantees len(data) % BlockSize == 0

	// TODO @ThomasPiellard shouldn't Sum() returns an error if there is no data?
	// TODO: @Tabaie, @Thomas Piellard Now sure what to make of this
	/*if len(d.data) == 0 {
		d.data = make([]byte, BlockSize)
	}*/

	for i := range d.data {
		r := d.encrypt(d.data[i])
		d.h.Add(&r, &d.h).Add(&d.h, &d.data[i])
	}

	return d.h
}

// plain execution of a mimc run
// m: message
// k: encryption key
func (d *digest) encrypt(m fr.Element) fr.Element {
	once.Do(initConstants) // init constants

	var tmp fr.Element
	for i := 0; i < mimcNbRounds; i++ {
		// m = (m+k+c)^5
		tmp.Add(&m, &d.h).Add(&tmp, &mimcConstants[i])
		m.Square(&tmp).
			Square(&m).
			Mul(&m, &tmp)
	}
	m.Add(&m, &d.h)
	return m
}

// Sum computes the mimc hash of msg from seed
func Sum(msg []byte) ([]byte, error) {
	var d digest
	if _, err := d.Write(msg); err != nil {
		return nil, err
	}
	h := d.checksum()
	bytes := h.Bytes()
	return bytes[:], nil
}

func initConstants() {
	bseed := ([]byte)(seed)

	hash := sha3.NewLegacyKeccak256()
	_, _ = hash.Write(bseed)
	rnd := hash.Sum(nil) // pre hash before use
	hash.Reset()
	_, _ = hash.Write(rnd)

	for i := 0; i < mimcNbRounds; i++ {
		rnd = hash.Sum(nil)
		mimcConstants[i].SetBytes(rnd)
		hash.Reset()
		_, _ = hash.Write(rnd)
	}
}

// WriteString writes a string that doesn't necessarily consist of field elements
func (d *digest) WriteString(rawBytes []byte) error {
	if elems, err := fr.Hash(rawBytes, []byte("string:"), 1); err != nil {
		return err
	} else {
		d.data = append(d.data, elems[0])
	}
	return nil
}

// SetState manually sets the state of the hasher to an user-provided value. In
// the context of MiMC, the method expects a byte slice of 32 elements.
func (d *digest) SetState(newState []byte) error {

	if len(newState) != 32 {
		return errors.New("the mimc state expects a state of 32 bytes")
	}

	if err := d.h.SetBytesCanonical(newState); err != nil {
		return errors.New("the provided newState does not represent a valid state")
	}

	d.data = nil

	return nil
}

// State returns the internal state of the hasher
func (d *digest) State() []byte {
	_ = d.Sum(nil) // this flushes the hasher
	b := d.h.Bytes()
	return b[:]
}
