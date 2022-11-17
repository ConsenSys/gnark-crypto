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
	"crypto/rand"
	"fmt"
	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bls24-315"
	"github.com/consensys/gnark-crypto/ecc/bls24-315/fr"
	"math/big"
)

// Key for proof and verification
type Key struct {
	g             bls24315.G2Affine // TODO @tabaie: does this really have to be randomized?
	gRootSigmaNeg bls24315.G2Affine //gRootSigmaNeg = g^{-1/σ}
	basis         []*bls24315.G1Affine
	basisExpSigma []bls24315.G1Affine
}

func randomOnG2() (bls24315.G2Affine, error) { // TODO: Add to G2.go?
	gBytes := make([]byte, fr.Bytes)
	if _, err := rand.Read(gBytes); err != nil {
		return bls24315.G2Affine{}, err
	}
	return bls24315.HashToG2(gBytes, []byte("random on g2"))
}

func Setup(basis []*bls24315.G1Affine) (Key, error) {
	var (
		k   Key
		err error
	)

	if k.g, err = randomOnG2(); err != nil {
		return k, err
	}

	var modMinusOne big.Int
	modMinusOne.Sub(fr.Modulus(), big.NewInt(1))
	var sigma *big.Int
	if sigma, err = rand.Int(rand.Reader, &modMinusOne); err != nil {
		return k, err
	}
	sigma.Add(sigma, big.NewInt(1))

	var sigmaInvNeg big.Int
	sigmaInvNeg.ModInverse(sigma, fr.Modulus())
	sigmaInvNeg.Sub(fr.Modulus(), &sigmaInvNeg)
	k.gRootSigmaNeg.ScalarMultiplication(&k.g, &sigmaInvNeg)

	k.basisExpSigma = make([]bls24315.G1Affine, len(basis))
	for i, gᵢ := range basis {
		k.basisExpSigma[i].ScalarMultiplication(gᵢ, sigma)
	}

	k.basis = basis
	return k, err
}

// TODO: If this takes too long in practice, edit MultiExp to accept pointers too
func ptrSliceToSlice[T any](ptrSlice []*T) []T {
	slice := make([]T, len(ptrSlice))
	for i, p := range ptrSlice {
		slice[i] = *p
	}
	return slice
}

func (k *Key) Commit(values []*fr.Element) (commitment bls24315.G1Affine, knowledgeProof bls24315.G1Affine, err error) {

	if len(values) != len(k.basis) {
		err = fmt.Errorf("unexpected number of values")
		return
	}

	valuesNoPtr := ptrSliceToSlice(values)
	config := ecc.MultiExpConfig{
		NbTasks:     1, // TODO Experiment
		ScalarsMont: true,
	}

	if _, err = commitment.MultiExp(ptrSliceToSlice(k.basis), valuesNoPtr, config); err != nil {
		return
	}

	_, err = knowledgeProof.MultiExp(k.basisExpSigma, valuesNoPtr, config)

	return
}

// VerifyKnowledgeProof checks if the proof of knowledge is valid
func (k *Key) VerifyKnowledgeProof(commitment bls24315.G1Affine, knowledgeProof bls24315.G1Affine) error {

	if !commitment.IsInSubGroup() || !knowledgeProof.IsInSubGroup() {
		return fmt.Errorf("subgroup check failed")
	}

	product, err := bls24315.Pair([]bls24315.G1Affine{commitment, knowledgeProof}, []bls24315.G2Affine{k.g, k.gRootSigmaNeg})
	if err != nil {
		return err
	}
	if product.IsOne() {
		return nil
	}
	return fmt.Errorf("proof rejected")
}
