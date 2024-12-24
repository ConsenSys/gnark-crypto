// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package mpcsetup

import (
	"bytes"
	curve "github.com/consensys/gnark-crypto/ecc/bls24-317"
	"github.com/stretchr/testify/require"
	"testing"
)

func TestContributionPok(t *testing.T) {
	const (
		pokChallenge = "challenge"
		pokDst       = 1
	)
	x0, err := curve.HashToG1([]byte("contribution test"), nil)
	require.NoError(t, err)
	y0, err := curve.RandomOnG2()
	require.NoError(t, err)
	x1, y1 := x0, y0
	proof := UpdateValues(nil, []byte(pokChallenge), pokDst, &x1, &y1)

	representations := []ValueUpdate{
		{Previous: x0, Next: x1},
		{Previous: y0, Next: y1},
	}

	// verify proof - G1 only
	require.NoError(t, proof.Verify([]byte(pokChallenge), pokDst, representations[0]))

	// verify proof - G2 only
	require.NoError(t, proof.Verify([]byte(pokChallenge), pokDst, representations[1]))

	// verify proof - G1 and G2
	require.NoError(t, proof.Verify([]byte(pokChallenge), pokDst, representations...))

	// read/write round-trip
	var bb bytes.Buffer
	n0, err := proof.WriteTo(&bb)
	require.NoError(t, err)
	var proofBack UpdateProof
	n1, err := proofBack.ReadFrom(&bb)
	require.NoError(t, err)
	require.Equal(t, n0, n1)

	require.NoError(t, proofBack.Verify([]byte(pokChallenge), pokDst, representations[0]))
	require.NoError(t, proofBack.Verify([]byte(pokChallenge), pokDst, representations[1]))
	require.NoError(t, proofBack.Verify([]byte(pokChallenge), pokDst, representations...))
}

func TestSameRatioMany(t *testing.T) {
	_, _, g1, g2 := curve.Generators()
	g1Slice := []curve.G1Affine{g1, g1, g1}
	g2Slice := []curve.G2Affine{g2, g2}
	require.NoError(t, SameRatioMany(g1Slice, g2Slice, g1Slice, g1Slice))
}
