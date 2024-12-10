// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fflonk

import (
	"crypto/sha256"
	"math/big"
	"testing"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bls24-317/fr"
	"github.com/consensys/gnark-crypto/ecc/bls24-317/kzg"
	"github.com/stretchr/testify/require"
)

// Test SRS re-used across tests of the KZG scheme
var testSrs *kzg.SRS
var bAlpha *big.Int

func init() {
	const srsSize = 600
	bAlpha = new(big.Int).SetInt64(42) // randomise ?
	testSrs, _ = kzg.NewSRS(ecc.NextPowerOfTwo(srsSize), bAlpha)
}

func TestFflonk(t *testing.T) {

	assert := require.New(t)

	// sample random polynomials of various sizes
	nbSets := 5
	p := make([][][]fr.Element, nbSets)
	for i := 0; i < nbSets; i++ {
		nbPolysInSet := 9
		p[i] = make([][]fr.Element, nbPolysInSet)
		for j := 0; j < nbPolysInSet; j++ {
			curSizePoly := j + 10
			p[i][j] = make([]fr.Element, curSizePoly)
			for k := 0; k < curSizePoly; k++ {
				p[i][j][k].SetRandom()
			}
		}
	}

	// sample random sets Sᵢ
	x := make([][]fr.Element, nbSets)
	for i := 0; i < nbSets; i++ {
		curSetSize := i + 4
		x[i] = make([]fr.Element, curSetSize)
		for j := 0; j < curSetSize; j++ {
			x[i][j].SetRandom()
		}
	}

	// commit to the folded polynomials
	digests := make([]kzg.Digest, nbSets)
	var err error
	for i := 0; i < nbSets; i++ {
		digests[i], err = FoldAndCommit(p[i], testSrs.Pk)
		assert.NoError(err)
	}

	// compute flonk opening proof
	hf := sha256.New()
	proof, err := BatchOpen(p, digests, x, hf, testSrs.Pk)
	assert.NoError(err)

	// check opening proof
	err = BatchVerify(proof, digests, x, hf, testSrs.Vk)
	assert.NoError(err)

	// tamper the proof
	proof.ClaimedValues[0][0][0].SetRandom()
	err = BatchVerify(proof, digests, x, hf, testSrs.Vk)
	assert.Error(err)

}

func TestCommit(t *testing.T) {

	assert := require.New(t)

	// sample polynomials
	nbPolys := 2
	p := make([][]fr.Element, nbPolys)
	for i := 0; i < nbPolys; i++ {
		p[i] = make([]fr.Element, i+10)
		for j := 0; j < i+10; j++ {
			p[i][j].SetRandom()
		}
	}

	// fflonk commit to them
	var x fr.Element
	x.SetRandom()
	proof, err := kzg.Open(Fold(p), x, testSrs.Pk)
	assert.NoError(err)

	// check that Open(C, x) = ∑_{i<t}Pᵢ(xᵗ)xⁱ
	var xt fr.Element
	var expo big.Int
	expo.SetUint64(uint64(nbPolys))
	xt.Exp(x, &expo)
	px := make([]fr.Element, nbPolys)
	for i := 0; i < nbPolys; i++ {
		px[i] = eval(p[i], xt)
	}
	y := eval(px, x)
	assert.True(y.Equal(&proof.ClaimedValue))
}

func TestGetIthRootOne(t *testing.T) {
	assert := require.New(t)

	order := getNextDivisorRMinusOne(9)
	omega, err := getIthRootOne(order)
	assert.NoError(err)
	var orderBigInt big.Int
	orderBigInt.SetUint64(uint64(order))
	omega.Exp(omega, &orderBigInt)
	assert.True(omega.IsOne())
}
