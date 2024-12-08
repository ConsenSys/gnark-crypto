//go:build !purego

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

package koalabear

// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Add: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call addVecGeneric
		addVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 16
	addVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n%blockSize != 0 {
		// call addVecGeneric on the rest
		start := n - n%blockSize
		addVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

//go:noescape
func addVec(res, a, b *Element, n uint64)

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	subVecGeneric(*vector, a, b)
}

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *Element) {
	scalarMulVecGeneric(*vector, a, b)
}

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res Element) {
	sumVecGeneric(&res, *vector)
	return
}

// InnerProduct computes the inner product of two vectors.
// It panics if the vectors don't have the same length.
func (vector *Vector) InnerProduct(other Vector) (res Element) {
	innerProductVecGeneric(&res, *vector, other)
	return
}

// Mul multiplies two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Mul(a, b Vector) {
	mulVecGeneric(*vector, a, b)
}
