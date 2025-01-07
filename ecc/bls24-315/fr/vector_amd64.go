//go:build !purego

// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fr

// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Add: vectors don't have the same length")
	}
	n := uint64(len(a))
	addVec(&(*vector)[0], &a[0], &b[0], n)
}

//go:noescape
func addVec(res, a, b *Element, n uint64)

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Sub: vectors don't have the same length")
	}
	subVec(&(*vector)[0], &a[0], &b[0], uint64(len(a)))
}

//go:noescape
func subVec(res, a, b *Element, n uint64)

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *Element) {
	if len(a) != len(*vector) {
		panic("vector.ScalarMul: vectors don't have the same length")
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || uint64(len(a)) >= maxN {
		// call scalarMulVecGeneric
		scalarMulVecGeneric(*vector, a, b)
		return
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	// the code for scalarMul is identical to mulVec; and it expects at least
	// 2 elements in the vector to fill the Z registers
	var bb [2]Element
	bb[0] = *b
	bb[1] = *b
	const blockSize = 16
	scalarMulVec(&(*vector)[0], &a[0], &bb[0], n/blockSize, qInvNeg)
	if n%blockSize != 0 {
		// call scalarMulVecGeneric on the rest
		start := n - n%blockSize
		scalarMulVecGeneric((*vector)[start:], a[start:], b)
	}
}

//go:noescape
func scalarMulVec(res, a, b *Element, n uint64, qInvNeg uint64)

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res Element) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	const minN = 16 * 7 // AVX512 slower than generic for small n
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n <= minN || n >= maxN {
		// call sumVecGeneric
		sumVecGeneric(&res, *vector)
		return
	}
	sumVec(&res, &(*vector)[0], uint64(len(*vector)))
	return
}

//go:noescape
func sumVec(res *Element, a *Element, n uint64)

// InnerProduct computes the inner product of two vectors.
// It panics if the vectors don't have the same length.
func (vector *Vector) InnerProduct(other Vector) (res Element) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	if n != uint64(len(other)) {
		panic("vector.InnerProduct: vectors don't have the same length")
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n >= maxN {
		// call innerProductVecGeneric
		// note; we could split the vector into smaller chunks and call innerProductVec
		innerProductVecGeneric(&res, *vector, other)
		return
	}
	innerProdVec(&res[0], &(*vector)[0], &other[0], uint64(len(*vector)))

	return
}

//go:noescape
func innerProdVec(res *uint64, a, b *Element, n uint64)

// Mul multiplies two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Mul(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Mul: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n >= maxN {
		// call mulVecGeneric
		mulVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 16
	mulVec(&(*vector)[0], &a[0], &b[0], n/blockSize, qInvNeg)
	if n%blockSize != 0 {
		// call mulVecGeneric on the rest
		start := n - n%blockSize
		mulVecGeneric((*vector)[start:], a[start:], b[start:])
	}

}

// Patterns use for transposing the vectors in mulVec
var (
	pattern1 = [8]uint64{0, 8, 1, 9, 2, 10, 3, 11}
	pattern2 = [8]uint64{12, 4, 13, 5, 14, 6, 15, 7}
	pattern3 = [8]uint64{0, 1, 8, 9, 2, 3, 10, 11}
	pattern4 = [8]uint64{12, 13, 4, 5, 14, 15, 6, 7}
)

//go:noescape
func mulVec(res, a, b *Element, n uint64, qInvNeg uint64)
