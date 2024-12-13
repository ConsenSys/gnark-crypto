//go:build !purego

// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package babybear

//go:noescape
func addVec(res, a, b *Element, n uint64)

//go:noescape
func subVec(res, a, b *Element, n uint64)

//go:noescape
func sumVec(t *uint64, a *Element, n uint64)

//go:noescape
func mulVec(res, a, b *Element, n uint64)

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

	const blockSize = 4
	addVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n%blockSize != 0 {
		// call addVecGeneric on the rest
		start := n - n%blockSize
		addVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Sub: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}

	const blockSize = 4
	subVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n%blockSize != 0 {
		// call subVecGeneric on the rest
		start := n - n%blockSize
		subVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *Element) {
	scalarMulVecGeneric(*vector, a, b)
}

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res Element) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}

	const blockSize = 16
	var t [4]uint64 // stores the accumulators (not reduced mod q)
	sumVec(&t[0], &(*vector)[0], n/blockSize)
	// we reduce the accumulators mod q and add to res
	var v Element
	for i := 0; i < 4; i++ {
		v[0] = uint32(t[i] % q)
		res.Add(&res, &v)
	}
	if n%blockSize != 0 {
		// call sumVecGeneric on the rest
		start := n - n%blockSize
		sumVecGeneric(&res, (*vector)[start:])
	}

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
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Mul: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}

	const blockSize = 4
	mulVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n%blockSize != 0 {
		// call mulVecGeneric on the rest
		start := n - n%blockSize
		mulVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}
