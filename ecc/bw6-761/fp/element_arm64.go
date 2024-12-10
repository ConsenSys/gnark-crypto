//go:build !purego

// Copyright 2020 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fp

// Butterfly sets
//
//	a = a + b (mod q)
//	b = a - b (mod q)
func Butterfly(a, b *Element) {
	_butterflyGeneric(a, b)
}

//go:noescape
func mul(res, x, y *Element)

// Mul z = x * y (mod q)
//
// x and y must be less than q
func (z *Element) Mul(x, y *Element) *Element {
	mul(z, x, y)
	return z
}

// Square z = x * x (mod q)
//
// x must be less than q
func (z *Element) Square(x *Element) *Element {
	// see Mul for doc.
	mul(z, x, x)
	return z
}

// MulBy3 x *= 3 (mod q)
func MulBy3(x *Element) {
	_x := *x
	x.Double(x).Add(x, &_x)
}

// MulBy5 x *= 5 (mod q)
func MulBy5(x *Element) {
	_x := *x
	x.Double(x).Double(x).Add(x, &_x)
}

// MulBy13 x *= 13 (mod q)
func MulBy13(x *Element) {
	var y = Element{
		4345973640412121648,
		16340807117537158706,
		14673764841507373218,
		5587754667198343811,
		12846753860245084942,
		4041391838244625385,
		8324122986343791677,
		8773809490091176420,
		5465994123296109449,
		6649773564661156048,
		9147430723089113754,
		54281803719730243,
	}
	x.Mul(x, &y)
}

func fromMont(z *Element) {
	_fromMontGeneric(z)
}

//go:noescape
func reduce(res *Element)
