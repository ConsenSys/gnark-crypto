// Copyright 2020 ConsenSys Software Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package pedersen

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls24-315"
	"github.com/consensys/gnark-crypto/ecc/bls24-315/fr"
	"github.com/consensys/gnark-crypto/utils"
	"github.com/stretchr/testify/assert"
	"testing"
)

func interfaceSliceToFrSlice(t *testing.T, values ...interface{}) []fr.Element {
	res := make([]fr.Element, len(values))
	for i, v := range values {
		_, err := res[i].SetInterface(v)
		assert.NoError(t, err)
	}
	return res
}

func randomFrSlice(t *testing.T, size int) []interface{} {
	res := make([]interface{}, size)
	var err error
	for i := range res {
		var v fr.Element
		res[i], err = v.SetRandom()
		assert.NoError(t, err)
	}
	return res
}

func randomOnG1() (curve.G1Affine, error) { // TODO: Add to G1.go?
	if gBytes, err := randomFrSizedBytes(); err == nil {
		return curve.HashToG1(gBytes, []byte("random on g1"))
	} else {
		return curve.G1Affine{}, err
	}
}

func randomG1Slice(t *testing.T, size int) []curve.G1Affine {
	res := make([]curve.G1Affine, size)
	for i := range res {
		var err error
		res[i], err = randomOnG1()
		assert.NoError(t, err)
	}
	return res
}

func testCommit(t *testing.T, values ...interface{}) {

	basis := randomG1Slice(t, len(values))

	var (
		pk              ProvingKey
		vk              VerifyingKey
		err             error
		commitment, pok curve.G1Affine
	)

	pk, vk, err = Setup(basis)
	assert.NoError(t, err)
	commitment, pok, err = pk.Commit(interfaceSliceToFrSlice(t, values...))
	assert.NoError(t, err)
	assert.NoError(t, vk.Verify(commitment, pok))

	pok.Neg(&pok)
	assert.NotNil(t, vk.Verify(commitment, pok))
}

func TestCommitToOne(t *testing.T) {
	testCommit(t, 1)
}

func TestCommitSingle(t *testing.T) {
	testCommit(t, randomFrSlice(t, 1)...)
}

func TestCommitFiveElements(t *testing.T) {
	testCommit(t, randomFrSlice(t, 5)...)
}

func TestMarshal(t *testing.T) {
	var pk ProvingKey
	pk.basisExpSigma = randomG1Slice(t, 5)
	pk.basis = randomG1Slice(t, 5)

	var (
		vk  VerifyingKey
		err error
	)
	vk.g, err = randomOnG2()
	assert.NoError(t, err)
	vk.gRootSigmaNeg, err = randomOnG2()
	assert.NoError(t, err)

	t.Run("ProvingKey -> Bytes -> ProvingKey must remain identical.", utils.SerializationRoundTrip(&pk))
	t.Run("VerifyingKey -> Bytes -> VerifyingKey must remain identical.", utils.SerializationRoundTrip(&vk))
}
