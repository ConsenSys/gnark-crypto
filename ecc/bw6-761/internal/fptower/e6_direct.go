// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

package fptower

import (
	"github.com/consensys/gnark-crypto/ecc/bw6-761/fp"
)

// E6D is a degree 6 finite field extension of fp
type E6D struct {
	A0, A1, A2, A3, A4, A5 fp.Element
}

// Equal returns true if z equals x, false otherwise
func (z *E6D) Equal(x *E6D) bool {
	return z.A0.Equal(&x.A0) &&
		z.A1.Equal(&x.A1) &&
		z.A2.Equal(&x.A2) &&
		z.A3.Equal(&x.A3) &&
		z.A4.Equal(&x.A4) &&
		z.A5.Equal(&x.A5)
}

// String puts E6D elmt in string form
func (z *E6D) String() string {
	return (z.A0.String() + "+(" + z.A1.String() + ")*u+(" + z.A2.String() + ")*u**2+(" + z.A3.String() + ")*u**3+(" + z.A4.String() + ")*u**4+(" + z.A5.String() + ")*u**5")
}

// SetString sets a E6D elmt from string
func (z *E6D) SetString(s1, s2, s3, s4, s5, s6 string) *E6D {
	z.A0.SetString(s1)
	z.A1.SetString(s2)
	z.A2.SetString(s3)
	z.A3.SetString(s4)
	z.A4.SetString(s5)
	z.A5.SetString(s6)
	return z
}

// Set copies x into z and returns z
func (z *E6D) Set(x *E6D) *E6D {
	*z = *x
	return z
}

// SetOne sets z to 1 in Montgomery form and returns z
func (z *E6D) SetOne() *E6D {
	z.A0.SetOne()
	z.A1.SetZero()
	z.A2.SetZero()
	z.A3.SetZero()
	z.A4.SetZero()
	z.A5.SetZero()
	return z
}

// Add sets z=x+y in E6D and returns z
func (z *E6D) Add(x, y *E6D) *E6D {
	z.A0.Add(&x.A0, &y.A0)
	z.A1.Add(&x.A1, &y.A1)
	z.A2.Add(&x.A2, &y.A2)
	z.A3.Add(&x.A3, &y.A3)
	z.A4.Add(&x.A4, &y.A4)
	z.A5.Add(&x.A5, &y.A5)
	return z
}

// Sub sets z to x-y and returns z
func (z *E6D) Sub(x, y *E6D) *E6D {
	z.A0.Sub(&x.A0, &y.A0)
	z.A1.Sub(&x.A1, &y.A1)
	z.A2.Sub(&x.A2, &y.A2)
	z.A3.Sub(&x.A3, &y.A3)
	z.A4.Sub(&x.A4, &y.A4)
	z.A5.Sub(&x.A5, &y.A5)
	return z
}

// Double sets z=2*x and returns z
func (z *E6D) Double(x *E6D) *E6D {
	z.A0.Double(&x.A0)
	z.A1.Double(&x.A1)
	z.A2.Double(&x.A2)
	z.A3.Double(&x.A3)
	z.A4.Double(&x.A4)
	z.A5.Double(&x.A5)
	return z
}

// SetRandom used only in tests
func (z *E6D) SetRandom() (*E6D, error) {
	if _, err := z.A0.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A1.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A2.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A3.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A4.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A5.SetRandom(); err != nil {
		return nil, err
	}
	return z, nil
}

// IsZero returns true if z is zero, false otherwise
func (z *E6D) IsZero() bool {
	return z.A0.IsZero() && z.A1.IsZero() && z.A2.IsZero() && z.A3.IsZero() && z.A4.IsZero() && z.A5.IsZero()
}

// IsOne returns true if z is one, false otherwise
func (z *E6D) IsOne() bool {
	return z.A0.IsOne() && z.A1.IsZero() && z.A2.IsZero() && z.A3.IsZero() && z.A4.IsZero() && z.A5.IsZero()
}

// Mul sets z=x*y in E6D and returns z
func (z *E6D) Mul(x, y *E6D) *E6D {
	return z.mulTower(x, y)
}

func (z *E6D) mulTower(x, y *E6D) *E6D {
	_x := ToTower(x)
	_y := ToTower(y)
	_x.Mul(_x, _y)
	_z := FromTower(_x)
	return z.Set(_z)
}

func (z *E6D) mulMontgomery6(x, y *E6D) *E6D {
	// Ref.: Peter L. Montgomery. Five, six, and seven-term Karatsuba-like formulae. IEEE
	// Transactions on Computers, 54(3):362–369, 2005.
	//
	// Fixing the polynomial C to X^6 we first compute the interpolation points
	// vi = x(Xi)*y(Xi) at Xi={0, ±1, ±2, ±3, ±4, 5,∞}:
	//
	//		v0 = (a0 + a1 + a2 + a3 + a4 + a5)(b0 + b1 + b2 + b3 + b4 + b5)
	//		v2 = (a0 + a1 + a3 + a4)(b0 + b1 + b3 + b4)
	//		v3 = (a0 − a2 − a3 + a5)(b0 − b2 − b3 + b5)
	//		v4 = (a0 − a2 − a5)(b0 − b2 − b5)
	//		v5 = (a0 + a3 − a5)(b0 + b3 − b5)
	//		v6 = (a0 + a1 + a2)(b0 + b1 + b2)
	//		v7 = (a3 + a4 + a5)(b3 + b4 + b5)
	//		v8 = (a2 + a3)(b2 + b3)
	//		v9 = (a1 − a4)(b1 − b4)
	//		v10 = (a1 + a2)(b1 + b2)
	//		v11 = (a3 + a4)(b3 + b4)
	//		v12 = (a0 + a1)(b0 + b1)
	//		v13 = (a4 + a5)(b4 + b5)
	//		v14 = a0b0
	//		v15 = a1b1
	//		v16 = a4b4
	//		v17 = a5b5

	var _t0, _s0, t0, t1, t2, t3, t4, s0, s1, s2, s3, s4 fp.Element
	_t0.Add(&x.A0, &x.A1)
	t0.Add(&_t0, &x.A2)
	t1.Add(&x.A3, &x.A4)
	t2.Add(&_t0, &t1)
	t3.Add(&t2, &x.A5)
	t3.Add(&t3, &x.A2)

	_s0.Add(&y.A0, &y.A1)
	s0.Add(&_s0, &y.A2)
	s1.Add(&y.A3, &y.A4)
	s2.Add(&_s0, &s1)
	s3.Add(&s2, &y.A5)
	s3.Add(&s3, &y.A2)

	var v0, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17 fp.Element
	v0.Mul(&t3, &s3)
	v2.Mul(&t2, &s2)
	v6.Mul(&t0, &s0)
	t4.Add(&t1, &x.A5)
	s4.Add(&s1, &y.A5)
	v7.Mul(&t4, &s4)
	v12.Mul(&_t0, &_s0)
	v11.Mul(&t1, &s1)
	t0.Add(&x.A2, &x.A3)
	s0.Add(&y.A2, &y.A3)
	v8.Mul(&t0, &s0)
	_t0.Sub(&x.A1, &x.A4)
	_s0.Sub(&y.A1, &y.A4)
	v9.Mul(&_t0, &_s0)
	t1.Add(&x.A1, &x.A2)
	s1.Add(&y.A1, &y.A2)
	v10.Mul(&t1, &s1)
	t1.Add(&x.A4, &x.A5)
	s1.Add(&y.A4, &y.A5)
	v13.Mul(&t1, &s1)
	v3.Add(&x.A0, &x.A5)
	v3.Sub(&v3, &t0)
	s1.Add(&y.A0, &y.A5)
	s1.Sub(&s1, &s0)
	v3.Mul(&v3, &s1)
	t1.Add(&x.A2, &x.A5)
	t2.Sub(&x.A0, &t1)
	s1.Add(&y.A2, &y.A5)
	s2.Sub(&y.A0, &s1)
	v4.Mul(&t2, &s2)
	t1.Add(&x.A0, &x.A3)
	t1.Sub(&t1, &x.A5)
	s1.Add(&y.A0, &y.A3)
	s1.Sub(&s1, &y.A5)
	v5.Mul(&t1, &s1)
	v14.Mul(&x.A0, &y.A0)
	v15.Mul(&x.A1, &y.A1)
	v16.Mul(&x.A4, &y.A4)
	v17.Mul(&x.A5, &y.A5)

	// Then we compute the coefficients c0,c1,c3,c4 and c5 in the direct sextic
	// extension of the product x*y as follows:
	//
	// 	c0 = v14 + β(v0 − v2 + v4 + 2(v3+v5+v6-v12) + 3(v7+v15-v8-v10-v11) +
	// 	4(v16-v13) − 5(v14+v17))
	//
	//  c1 = v12 − (v14 + v15) + β(v8 + v10 + v12 − (v3 + v5 + v6 + v15) +
	//  2(v14 + v17 + v13 - v7) + 3(v11 - v16))
	//
	// 	c2 = 2v15 + v6 − (v10 + v12) + β(2v16 + v7 − (v11 + v13))
	//
	// 	c3 = v8 + v11 + v13 − (v3 + v4 + v7 + v16) + 3(v10 - v15) + 2(v12 + v14
	// 	+ v17 - v6) + β(v13 − (v16 + v17))
	//
	// 	c4 = v2 + v3 + v4 + v7 + v15 + v9 − (v8 + v13) − 3v12 + 2(v6 − (v17 +
	// 	v10 + v11 + v14)) + βv17
	//
	//  c5 = −(v3 + v4 + v5 + v9 + v15 + v16) + 2(v8 + v10 + v11 + v12 + v13 −
	//  (v6 + v7)) + 3(v14 + v17)

	var c0, c1, c2, c3, c4, c5, s811, s81110, s35, s1012, s34 fp.Element
	var twelve, twenty, twentyone fp.Element
	twelve.SetUint64(12)
	twenty.SetUint64(20)
	twentyone.SetUint64(21)
	c0.Double(&v2).Double(&c0)
	s811.Add(&v8, &v11)
	s81110.Add(&s811, &v10)
	s1.Mul(&s81110, &twelve)
	c0.Add(&c0, &s1)
	s1.Double(&v12).Double(&s1).Double(&s1)
	c0.Add(&c0, &s1)
	s1.Double(&v13).Double(&s1).Double(&s1).Double(&s1)
	c0.Add(&c0, &s1)
	s1.Mul(&v14, &twentyone)
	c0.Add(&c0, &s1)
	s1.Mul(&v17, &twenty)
	c0.Add(&c0, &s1)
	s1.Mul(&v15, &twelve)
	s2.Double(&v16).Double(&s2).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v0).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v3).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v4).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v5).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v6).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Mul(&v7, &twelve)
	s1.Add(&s1, &s2)
	c0.Sub(&c0, &s1)

	s35.Add(&v3, &v5)
	c1.Add(&s35, &v6)
	c1.Double(&c1).Double(&c1)
	s1.Double(&v7).Double(&s1).Double(&s1)
	c1.Add(&c1, &s1)
	s1.Mul(&v16, &twelve)
	c1.Add(&c1, &s1)
	s1 = v15
	fp.MulBy3(&s1)
	c1.Add(&c1, &s1)
	s1 = v12
	fp.MulBy3(&s1)
	s2 = v14
	fp.MulBy3(&s2)
	fp.MulBy3(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v8).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v10).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Mul(&v11, &twelve)
	s1.Add(&s1, &s2)
	s2.Double(&v13).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v17).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	c1.Sub(&c1, &s1)

	c2.Double(&v15)
	c2.Add(&c2, &v6)
	s1.Double(&v11).Double(&s1)
	c2.Add(&c2, &s1)
	s1.Double(&v13).Double(&s1)
	c2.Add(&c2, &s1)
	s1012.Add(&v10, &v12)
	s2.Double(&v7).Double(&s2)
	s1.Add(&s1012, &s2)
	s2.Double(&v16).Double(&s2).Double(&s2)
	s1.Add(&s1, &s2)
	c2.Sub(&c2, &s1)

	s1 = v10
	fp.MulBy3(&s1)
	c3.Add(&s811, &s1)
	s1.Double(&v12)
	c3.Add(&c3, &s1)
	s1.Double(&v14)
	c3.Add(&c3, &s1)
	s1 = v16
	fp.MulBy3(&s1)
	c3.Add(&c3, &s1)
	s1.Double(&v17)
	fp.MulBy3(&s1)
	c3.Add(&c3, &s1)
	s34.Add(&v3, &v4)
	s1.Add(&s34, &v7)
	s2.Double(&v6)
	s1.Add(&s1, &s2)
	s2 = v13
	fp.MulBy3(&s2)
	s1.Add(&s1, &s2)
	s2 = v15
	fp.MulBy3(&s2)
	s1.Add(&s1, &s2)
	c3.Sub(&c3, &s1)

	c4.Add(&v2, &v15)
	c4.Add(&c4, &v9)
	c4.Add(&c4, &v7)
	c4.Add(&c4, &s34)
	s1.Double(&v6)
	c4.Add(&c4, &s1)
	s1.Add(&v13, &v8)
	s2.Double(&v10)
	s1.Add(&s1, &s2)
	s2.Double(&v11)
	s1.Add(&s1, &s2)
	s2 = v12
	fp.MulBy3(&s2)
	s1.Add(&s1, &s2)
	s2.Double(&v14)
	s1.Add(&s1, &s2)
	s2.Double(&v17)
	fp.MulBy3(&s2)
	s1.Add(&s1, &s2)
	c4.Sub(&c4, &s1)

	c5.Add(&s81110, &v12)
	c5.Add(&c5, &v13)
	c5.Double(&c5)
	s1 = v14
	fp.MulBy3(&s1)
	c5.Add(&c5, &s1)
	s1 = v17
	fp.MulBy3(&s1)
	c5.Add(&c5, &s1)
	s1.Add(&v15, &v16)
	s1.Add(&s1, &s34)
	s1.Add(&s1, &v5)
	s1.Add(&s1, &v9)
	s2.Double(&v6)
	s1.Add(&s1, &s2)
	s2.Double(&v7)
	s1.Add(&s1, &s2)
	c5.Sub(&c5, &s1)

	z.A0.Set(&c0)
	z.A1.Set(&c1)
	z.A2.Set(&c2)
	z.A3.Set(&c3)
	z.A4.Set(&c4)
	z.A5.Set(&c5)

	return z
}

// MulBy023 multiplies z by an E6 sparse element 023
func (z *E6D) MulBy023(c0, c1 *fp.Element) *E6D {
	var a, b, tmp, a0, a1, a2, b0, b1, b2, d, zC10, zC11, zC12, one, t0, t1, t2 fp.Element
	a.Mul(&z.A0, c0)
	b.Mul(&z.A2, c1)
	tmp.Add(&z.A2, &z.A4)
	a0.Mul(c1, &tmp)
	a0.Sub(&b, &a0)
	a0.Double(&a0).Double(&a0)
	a0.Add(&a0, &a)
	a2.Mul(&z.A4, c0)
	a2.Add(&a2, &b)
	a1.Add(c0, c1)
	tmp.Add(&z.A0, &z.A2)
	a1.Mul(&a1, &tmp)
	a1.Sub(&a1, &a)
	a1.Sub(&a1, &b)

	b0.Double(&z.A5).Double(&b0)
	b2.Neg(&z.A3)
	b1.Neg(&z.A1)

	one.SetOne()
	d.Add(c1, &one)

	zC10.Add(&z.A1, &z.A0)
	zC11.Add(&z.A3, &z.A2)
	zC12.Add(&z.A5, &z.A4)

	a.Mul(&zC10, c0)
	b.Mul(&zC11, &d)
	tmp.Add(&zC11, &zC12)
	t0.Mul(&d, &tmp)
	t0.Sub(&b, &t0)
	t0.Double(&t0).Double(&t0)
	t0.Add(&t0, &a)
	t2.Mul(&zC12, c0)
	t2.Add(&t2, &b)
	t1.Add(c0, &d)
	tmp.Add(&zC10, &zC11)
	t1.Mul(&t1, &tmp)
	t1.Sub(&t1, &a)
	t1.Sub(&t1, &b)

	zC10.Sub(&t0, &a0)
	zC11.Sub(&t1, &a1)
	zC12.Sub(&t2, &a2)

	zC10.Add(&zC10, &b0)
	zC11.Add(&zC11, &b1)
	zC12.Add(&zC12, &b2)

	var zC00, zC01, zC02 fp.Element
	zC00.Double(&b2).Double(&zC00)
	zC00.Add(&zC00, &a0)
	zC01.Sub(&a1, &b0)
	zC02.Sub(&a2, &b1)

	z.A0.Set(&zC00)
	z.A1.Set(&zC10)
	z.A2.Set(&zC01)
	z.A3.Set(&zC11)
	z.A4.Set(&zC02)
	z.A5.Set(&zC12)

	return z
}

// Square sets z=x*x in E6D and returns z
func (z *E6D) Square(x *E6D) *E6D {
	_x := ToTower(x)
	_x.Square(_x)
	_z := FromTower(_x)
	return z.Set(_z)
}

// Inverse sets z to the inverse of x in E6D and returns z
//
// if x == 0, sets and returns z = x
func (z *E6D) Inverse(x *E6D) *E6D {
	_x := ToTower(x)
	_x.Inverse(_x)
	_z := FromTower(_x)
	return z.Set(_z)
}

// InverseUnitary inverses a unitary element
func (z *E6D) InverseUnitary(x *E6D) *E6D {
	return z.Conjugate(x)
}

// Conjugate sets z to x conjugated and returns z
func (z *E6D) Conjugate(x *E6D) *E6D {
	z.A0.Set(&x.A0)
	z.A1.Neg(&x.A1)
	z.A2.Set(&x.A2)
	z.A3.Neg(&x.A3)
	z.A4.Set(&x.A4)
	z.A5.Neg(&x.A5)
	return z
}

// FromTower
func FromTower(x *E6) *E6D {
	// gnark-crypto uses a quadratic over cubic sextic extension of Fp.
	// The two towers are isomorphic and the coefficients are permuted as follows:
	// 		a00 a01 a02 a10 a11 a12
	// 		A0  A2  A4  A1  A3  A5
	var z E6D
	z.A0.Set(&x.B0.A0)
	z.A1.Set(&x.B1.A0)
	z.A2.Set(&x.B0.A1)
	z.A3.Set(&x.B1.A1)
	z.A4.Set(&x.B0.A2)
	z.A5.Set(&x.B1.A2)
	return &z
}

// ToTower
func ToTower(x *E6D) *E6 {
	// gnark-crypto uses a quadratic over cubic sextic extension of Fp.
	// The two towers are isomorphic and the coefficients are permuted as follows:
	// 		a00 a01 a02 a10 a11 a12
	// 		A0  A2  A4  A1  A3  A5
	var z E6
	z.B0.A0.Set(&x.A0)
	z.B1.A0.Set(&x.A1)
	z.B0.A1.Set(&x.A2)
	z.B1.A1.Set(&x.A3)
	z.B0.A2.Set(&x.A4)
	z.B1.A2.Set(&x.A5)
	return &z
}
