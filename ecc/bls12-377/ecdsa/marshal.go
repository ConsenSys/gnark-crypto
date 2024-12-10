// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package ecdsa

import (
	"crypto/subtle"
	"errors"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
	"io"
	"math/big"
)

var errWrongSize = errors.New("wrong size buffer")
var errRBiggerThanRMod = errors.New("r >= r_mod")
var errSBiggerThanRMod = errors.New("s >= r_mod")
var errZero = errors.New("zero value")

// Bytes returns the binary representation of the public key
// follows https://tools.ietf.org/html/rfc8032#section-3.1
// and returns a compressed representation of the point (x,y)
//
// x, y are the coordinates of the point
// on the curve as big endian integers.
// compressed representation store x with a parity bit to recompute y
func (pk *PublicKey) Bytes() []byte {
	var res [sizePublicKey]byte
	pkBin := pk.A.Bytes()
	subtle.ConstantTimeCopy(1, res[:sizePublicKey], pkBin[:])
	return res[:]
}

// SetBytes sets p from binary representation in buf.
// buf represents a public key as x||y where x, y are
// interpreted as big endian binary numbers corresponding
// to the coordinates of a point on the curve.
// It returns the number of bytes read from the buffer.
func (pk *PublicKey) SetBytes(buf []byte) (int, error) {
	n := 0
	if len(buf) < sizePublicKey {
		return n, io.ErrShortBuffer
	}
	if _, err := pk.A.SetBytes(buf[:sizePublicKey]); err != nil {
		return 0, err
	}
	n += sizeFp
	return n, nil
}

// Bytes returns the binary representation of pk,
// as byte array publicKey||scalar
// where publicKey is as publicKey.Bytes(), and
// scalar is in big endian, of size sizeFr.
func (privKey *PrivateKey) Bytes() []byte {
	var res [sizePrivateKey]byte
	pubkBin := privKey.PublicKey.A.Bytes()
	subtle.ConstantTimeCopy(1, res[:sizePublicKey], pubkBin[:])
	subtle.ConstantTimeCopy(1, res[sizePublicKey:sizePrivateKey], privKey.scalar[:])
	return res[:]
}

// SetBytes sets pk from buf, where buf is interpreted
// as  publicKey||scalar
// where publicKey is as publicKey.Bytes(), and
// scalar is in big endian, of size sizeFr.
// It returns the number byte read.
func (privKey *PrivateKey) SetBytes(buf []byte) (int, error) {
	n := 0
	if len(buf) < sizePrivateKey {
		return n, io.ErrShortBuffer
	}
	if _, err := privKey.PublicKey.A.SetBytes(buf[:sizePublicKey]); err != nil {
		return 0, err
	}
	n += sizePublicKey
	subtle.ConstantTimeCopy(1, privKey.scalar[:], buf[sizePublicKey:sizePrivateKey])
	n += sizeFr
	return n, nil
}

// Bytes returns the binary representation of sig
// as a byte array of size 2*sizeFr r||s
func (sig *Signature) Bytes() []byte {
	var res [sizeSignature]byte
	subtle.ConstantTimeCopy(1, res[:sizeFr], sig.R[:])
	subtle.ConstantTimeCopy(1, res[sizeFr:], sig.S[:])
	return res[:]
}

// SetBytes sets sig from a buffer in binary.
// buf is read interpreted as r||s
// It returns the number of bytes read from buf.
func (sig *Signature) SetBytes(buf []byte) (int, error) {
	n := 0
	if len(buf) != sizeSignature {
		return n, errWrongSize
	}

	// S, R < R_mod (to avoid malleability)
	frMod := fr.Modulus()
	zero := big.NewInt(0)
	bufBigInt := new(big.Int)
	bufBigInt.SetBytes(buf[:sizeFr])
	if bufBigInt.Cmp(zero) == 0 {
		return 0, errZero
	}
	if bufBigInt.Cmp(frMod) != -1 {
		return 0, errRBiggerThanRMod
	}
	bufBigInt.SetBytes(buf[sizeFr : 2*sizeFr])
	if bufBigInt.Cmp(zero) == 0 {
		return 0, errZero
	}
	if bufBigInt.Cmp(frMod) != -1 {
		return 0, errSBiggerThanRMod
	}

	subtle.ConstantTimeCopy(1, sig.R[:], buf[:sizeFr])
	n += sizeFr
	subtle.ConstantTimeCopy(1, sig.S[:], buf[sizeFr:2*sizeFr])
	n += sizeFr
	return n, nil
}
