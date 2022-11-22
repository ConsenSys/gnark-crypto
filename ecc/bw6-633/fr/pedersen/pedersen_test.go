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
	"github.com/consensys/gnark-crypto/ecc/bw6-633"
	"github.com/consensys/gnark-crypto/ecc/bw6-633/fr"
	"github.com/stretchr/testify/assert"
	"math/rand"
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

func randomOnG1() (bw6633.G1Affine, error) { // TODO: Add to G1.go?
	gBytes := make([]byte, fr.Bytes)
	if _, err := rand.Read(gBytes); err != nil {
		return bw6633.G1Affine{}, err
	}
	return bw6633.HashToG1(gBytes, []byte("random on g2"))
}

func testCommit(t *testing.T, values ...interface{}) {

	basis := make([]bw6633.G1Affine, len(values))
	for i := range basis {
		var err error
		basis[i], err = randomOnG1()
		assert.NoError(t, err)
	}

	var (
		key             Key
		err             error
		commitment, pok bw6633.G1Affine
	)

	key, err = Setup(basis)
	assert.NoError(t, err)
	commitment, pok, err = key.Commit(interfaceSliceToFrSlice(t, values...))
	assert.NoError(t, err)
	assert.NoError(t, key.VerifyKnowledgeProof(commitment, pok))

	pok.Neg(&pok)
	assert.NotNil(t, key.VerifyKnowledgeProof(commitment, pok))
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
