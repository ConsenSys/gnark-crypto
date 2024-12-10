// Copyright 2020 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fp

// expBySqrtExp is equivalent to z.Exp(x, 3fffffffffffffffffffffffffffffffffffffffffffffffffffffffbfffff0c)
//
// uses github.com/mmcloughlin/addchain v0.4.0 to generate a shorter addition chain
func (z *Element) expBySqrtExp(x Element) *Element {
	// addition chain:
	//
	//	_10      = 2*1
	//	_11      = 1 + _10
	//	_1100    = _11 << 2
	//	_1111    = _11 + _1100
	//	_11110   = 2*_1111
	//	_11111   = 1 + _11110
	//	_1111100 = _11111 << 2
	//	_1111111 = _11 + _1111100
	//	x11      = _1111111 << 4 + _1111
	//	x22      = x11 << 11 + x11
	//	x27      = x22 << 5 + _11111
	//	x54      = x27 << 27 + x27
	//	x108     = x54 << 54 + x54
	//	x216     = x108 << 108 + x108
	//	x223     = x216 << 7 + _1111111
	//	return     ((x223 << 23 + x22) << 6 + _11) << 2
	//
	// Operations: 253 squares 13 multiplies

	// Allocate Temporaries.
	var (
		t0 = new(Element)
		t1 = new(Element)
		t2 = new(Element)
		t3 = new(Element)
	)

	// var t0,t1,t2,t3 Element
	// Step 1: z = x^0x2
	z.Square(&x)

	// Step 2: z = x^0x3
	z.Mul(&x, z)

	// Step 4: t0 = x^0xc
	t0.Square(z)
	for s := 1; s < 2; s++ {
		t0.Square(t0)
	}

	// Step 5: t0 = x^0xf
	t0.Mul(z, t0)

	// Step 6: t1 = x^0x1e
	t1.Square(t0)

	// Step 7: t2 = x^0x1f
	t2.Mul(&x, t1)

	// Step 9: t1 = x^0x7c
	t1.Square(t2)
	for s := 1; s < 2; s++ {
		t1.Square(t1)
	}

	// Step 10: t1 = x^0x7f
	t1.Mul(z, t1)

	// Step 14: t3 = x^0x7f0
	t3.Square(t1)
	for s := 1; s < 4; s++ {
		t3.Square(t3)
	}

	// Step 15: t0 = x^0x7ff
	t0.Mul(t0, t3)

	// Step 26: t3 = x^0x3ff800
	t3.Square(t0)
	for s := 1; s < 11; s++ {
		t3.Square(t3)
	}

	// Step 27: t0 = x^0x3fffff
	t0.Mul(t0, t3)

	// Step 32: t3 = x^0x7ffffe0
	t3.Square(t0)
	for s := 1; s < 5; s++ {
		t3.Square(t3)
	}

	// Step 33: t2 = x^0x7ffffff
	t2.Mul(t2, t3)

	// Step 60: t3 = x^0x3ffffff8000000
	t3.Square(t2)
	for s := 1; s < 27; s++ {
		t3.Square(t3)
	}

	// Step 61: t2 = x^0x3fffffffffffff
	t2.Mul(t2, t3)

	// Step 115: t3 = x^0xfffffffffffffc0000000000000
	t3.Square(t2)
	for s := 1; s < 54; s++ {
		t3.Square(t3)
	}

	// Step 116: t2 = x^0xfffffffffffffffffffffffffff
	t2.Mul(t2, t3)

	// Step 224: t3 = x^0xfffffffffffffffffffffffffff000000000000000000000000000
	t3.Square(t2)
	for s := 1; s < 108; s++ {
		t3.Square(t3)
	}

	// Step 225: t2 = x^0xffffffffffffffffffffffffffffffffffffffffffffffffffffff
	t2.Mul(t2, t3)

	// Step 232: t2 = x^0x7fffffffffffffffffffffffffffffffffffffffffffffffffffff80
	for s := 0; s < 7; s++ {
		t2.Square(t2)
	}

	// Step 233: t1 = x^0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffff
	t1.Mul(t1, t2)

	// Step 256: t1 = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff800000
	for s := 0; s < 23; s++ {
		t1.Square(t1)
	}

	// Step 257: t0 = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffffbfffff
	t0.Mul(t0, t1)

	// Step 263: t0 = x^0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc0
	for s := 0; s < 6; s++ {
		t0.Square(t0)
	}

	// Step 264: z = x^0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc3
	z.Mul(z, t0)

	// Step 266: z = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffffbfffff0c
	for s := 0; s < 2; s++ {
		z.Square(z)
	}

	return z
}

// expByLegendreExp is equivalent to z.Exp(x, 7fffffffffffffffffffffffffffffffffffffffffffffffffffffff7ffffe17)
//
// uses github.com/mmcloughlin/addchain v0.4.0 to generate a shorter addition chain
func (z *Element) expByLegendreExp(x Element) *Element {
	// addition chain:
	//
	//	_10       = 2*1
	//	_100      = 2*_10
	//	_110      = _10 + _100
	//	_111      = 1 + _110
	//	_1110     = 2*_111
	//	_10101    = _111 + _1110
	//	_10111    = _10 + _10101
	//	_101110   = 2*_10111
	//	_10111000 = _101110 << 2
	//	_11100110 = _101110 + _10111000
	//	_11111101 = _10111 + _11100110
	//	x11       = _11111101 << 3 + _10111
	//	x22       = x11 << 11 + x11
	//	i29       = 2*x22
	//	i31       = i29 << 2
	//	i54       = i31 << 22 + i31
	//	i122      = (i54 << 20 + i29) << 46 + i54
	//	x223      = i122 << 110 + i122 + _111
	//	return      (x223 << 23 + x22) << 9 + _10111
	//
	// Operations: 253 squares 15 multiplies

	// Allocate Temporaries.
	var (
		t0 = new(Element)
		t1 = new(Element)
		t2 = new(Element)
		t3 = new(Element)
		t4 = new(Element)
	)

	// var t0,t1,t2,t3,t4 Element
	// Step 1: z = x^0x2
	z.Square(&x)

	// Step 2: t0 = x^0x4
	t0.Square(z)

	// Step 3: t0 = x^0x6
	t0.Mul(z, t0)

	// Step 4: t1 = x^0x7
	t1.Mul(&x, t0)

	// Step 5: t0 = x^0xe
	t0.Square(t1)

	// Step 6: t0 = x^0x15
	t0.Mul(t1, t0)

	// Step 7: z = x^0x17
	z.Mul(z, t0)

	// Step 8: t0 = x^0x2e
	t0.Square(z)

	// Step 10: t2 = x^0xb8
	t2.Square(t0)
	for s := 1; s < 2; s++ {
		t2.Square(t2)
	}

	// Step 11: t0 = x^0xe6
	t0.Mul(t0, t2)

	// Step 12: t0 = x^0xfd
	t0.Mul(z, t0)

	// Step 15: t0 = x^0x7e8
	for s := 0; s < 3; s++ {
		t0.Square(t0)
	}

	// Step 16: t0 = x^0x7ff
	t0.Mul(z, t0)

	// Step 27: t2 = x^0x3ff800
	t2.Square(t0)
	for s := 1; s < 11; s++ {
		t2.Square(t2)
	}

	// Step 28: t0 = x^0x3fffff
	t0.Mul(t0, t2)

	// Step 29: t3 = x^0x7ffffe
	t3.Square(t0)

	// Step 31: t2 = x^0x1fffff8
	t2.Square(t3)
	for s := 1; s < 2; s++ {
		t2.Square(t2)
	}

	// Step 53: t4 = x^0x7ffffe000000
	t4.Square(t2)
	for s := 1; s < 22; s++ {
		t4.Square(t4)
	}

	// Step 54: t2 = x^0x7ffffffffff8
	t2.Mul(t2, t4)

	// Step 74: t4 = x^0x7ffffffffff800000
	t4.Square(t2)
	for s := 1; s < 20; s++ {
		t4.Square(t4)
	}

	// Step 75: t3 = x^0x7fffffffffffffffe
	t3.Mul(t3, t4)

	// Step 121: t3 = x^0x1ffffffffffffffff800000000000
	for s := 0; s < 46; s++ {
		t3.Square(t3)
	}

	// Step 122: t2 = x^0x1fffffffffffffffffffffffffff8
	t2.Mul(t2, t3)

	// Step 232: t3 = x^0x7ffffffffffffffffffffffffffe0000000000000000000000000000
	t3.Square(t2)
	for s := 1; s < 110; s++ {
		t3.Square(t3)
	}

	// Step 233: t2 = x^0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffff8
	t2.Mul(t2, t3)

	// Step 234: t1 = x^0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffff
	t1.Mul(t1, t2)

	// Step 257: t1 = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff800000
	for s := 0; s < 23; s++ {
		t1.Square(t1)
	}

	// Step 258: t0 = x^0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffffbfffff
	t0.Mul(t0, t1)

	// Step 267: t0 = x^0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffff7ffffe00
	for s := 0; s < 9; s++ {
		t0.Square(t0)
	}

	// Step 268: z = x^0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffff7ffffe17
	z.Mul(z, t0)

	return z
}
