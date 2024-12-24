// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package mpcsetup

import (
	"bytes"
	curve "github.com/consensys/gnark-crypto/ecc/bls24-317"
	curve "github.com/consensys/gnark-crypto/ecc/bls24-317/fr"
	"github.com/stretchr/testify/require"
	"math/big"
	"testing"
)

// small tests for sub-functionalities of the mpc setup
// this file is not autogenerated, and not generified for other curves

func TestContributionPok(t *testing.T) {
	const (
		pokChallenge = "challenge"
		pokDst       = 1
	)
	x0, err := curve.HashToG1([]byte("contribution test"), nil)
	require.NoError(t, err)
	proof, d := newValueUpdate([]byte(pokChallenge), pokDst)
	var (
		x1 curve.G1Affine
		d  fr.Element
		dI big.Int
	)
	d.BigInt(&dI)
	x1.ScalarMultiplication(&x0, &dI)

	// verify proof - no G2
	require.NoError(t, proof.verify(pair{x0, nil}, pair{x1, nil}, []byte(pokChallenge), pokDst))

	// verify proof - with G2
	y0, err := curve.RandomOnG2()
	require.NoError(t, err)
	var y1 curve.G2Affine
	y1.ScalarMultiplication(&y0, &dI)

	require.NoError(t, proof.verify(pair{x0, &y0}, pair{x1, &y1}, []byte(pokChallenge), pokDst))

	// read/write round-trip
	var bb bytes.Buffer
	n0, err := proof.WriteTo(&bb)
	require.NoError(t, err)
	var proofBack valueUpdate
	n1, err := proofBack.ReadFrom(&bb)
	require.NoError(t, err)
	require.Equal(t, n0, n1)

	require.NoError(t, proofBack.verify(pair{x0, nil}, pair{x1, nil}, []byte(pokChallenge), pokDst))
	require.NoError(t, proofBack.verify(pair{x0, &y0}, pair{x1, &y1}, []byte(pokChallenge), pokDst))
}
