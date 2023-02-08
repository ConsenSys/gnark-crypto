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

package iop

import (
	"math/big"
	"math/bits"

	"github.com/consensys/gnark-crypto/ecc/bls24-317/fr"
	"github.com/consensys/gnark-crypto/ecc/bls24-317/fr/fft"
)

// Basis indicates the basis in which a polynomial is represented.
type Basis int64

const (
	Canonical Basis = iota
	Lagrange
	LagrangeCoset
)

// Layout indicates if a polynomial has a BitReverse or a Regular layout
type Layout int64

const (
	Regular Layout = iota
	BitReverse
)

// Form describes the form of a polynomial.
type Form struct {
	Basis  Basis
	Layout Layout
}

// enum of the possible Form values for type-safe switches
// in this package
var (
	canonicalRegular        = Form{Canonical, Regular}
	canonicalBitReverse     = Form{Canonical, BitReverse}
	lagrangeRegular         = Form{Lagrange, Regular}
	lagrangeBitReverse      = Form{Lagrange, BitReverse}
	lagrangeCosetRegular    = Form{LagrangeCoset, Regular}
	lagrangeCosetBitReverse = Form{LagrangeCoset, BitReverse}
)

// Polynomial represents a polynomial, the vector of coefficients
// along with the basis and the layout.
type Polynomial struct {
	Coefficients []fr.Element
	Form
}

// NewPolynomial creates a new polynomial. The slice coeff NOT copied
// but directly assigned to the new polynomial.
func NewPolynomial(coeffs []fr.Element, form Form) *Polynomial {
	return &Polynomial{Coefficients: coeffs, Form: form}
}

// Clone returns a deep copy of the underlying data structure.
func (p *Polynomial) Clone() *Polynomial {
	r := &Polynomial{
		Coefficients: make([]fr.Element, len(p.Coefficients)),
		Form:         p.Form,
	}
	copy(r.Coefficients, p.Coefficients)
	return r
}

// WrappedPolynomial wraps a polynomial so that it is
// interpreted as P'(X)=P(\omega^{s}X).
// Size is the real size of the polynomial (seen as a vector).
// For instance if len(P)=32 but P.Size=8, it means that P has been
// extended (e.g. it is evaluated on a larger set) but P is a polynomial
// of degree 7.
// BlindedSize is the size of the polynomial when it is blinded. By
// default BlindedSize=Size, until the polynomial is blinded.
type WrappedPolynomial struct {
	*Polynomial
	Shift       int
	Size        int
	BlindedSize int
}

// NewWrappedPolynomial returned a WrappedPolynomial from p.
// ! Warning this does not do a deep copy of p, and modifications on the wrapped
// ! polynomial will modify the underlying coefficients of p.
func NewWrappedPolynomial(p *Polynomial) *WrappedPolynomial {
	return &WrappedPolynomial{
		Polynomial:  p,
		Size:        len(p.Coefficients),
		BlindedSize: len(p.Coefficients),
	}
}

// TODO @gbotrel rename to Shift
// ShiftMe the wrapped polynomial; it doesn't modify the underlying data structure,
// but flag the WrappedPolynomial such that it will be interpreted as p(\omega^shift X)
func (wp *WrappedPolynomial) ShiftMe(shift int) *WrappedPolynomial {
	wp.Shift = shift
	return wp
}

//----------------------------------------------------
// Blind a polynomial

// blindPoly blinds a polynomial q by adding Q(X)*(X^{n}-1),
// where deg Q = blindingOrder and Q is random, and n is the
// size of q. Sets the result to p and returns it.
//
// * bo blinding order,  it's the degree of Q, where the blinding is Q(X)*(X^{n}-1)
// where n is the size of wp. The size of wp is modified since the underlying
// polynomial is of bigger degree now. The new size is wp.Size+1+blindingOrder.
//
// /!\ The code panics if wq is not in canonical, regular layout
func (wp *WrappedPolynomial) Blind(wq *WrappedPolynomial, blindingOrder int) *WrappedPolynomial {

	// check that q is in canonical basis
	if wq.Polynomial.Basis != Canonical || wq.Polynomial.Layout != Regular {
		panic("the input must be in canonical basis, regular layout")
	}

	// ensure we don't mutate input
	// TODO @gbotrel this is inconsistent with other APIs
	if wp != wq {
		wp.Polynomial = wq.Polynomial.Clone()
		wp.Shift = wq.Shift
		wp.Size = wq.Size
	}

	// we add Q*(x^{n}-1) so the new size is deg(Q)+n+1
	// where n is the size of wq.
	newSize := wp.Size + blindingOrder + 1

	// Resize wp. The size of wq might has already been increased
	// (e.g. when the polynomial is evaluated on a larger domain),
	// if that's the case we don't resize the polynomial.
	offset := newSize - len(wp.Polynomial.Coefficients)
	if offset > 0 {
		z := make([]fr.Element, offset)
		wp.Polynomial.Coefficients = append(wp.Polynomial.Coefficients, z...)
	}

	// blinding: we add Q(X)(X^{n}-1) to P, where deg(Q)=blindingOrder
	var r fr.Element

	for i := 0; i <= blindingOrder; i++ {
		r.SetRandom()
		wp.Polynomial.Coefficients[i].Sub(&wp.Polynomial.Coefficients[i], &r)
		wp.Polynomial.Coefficients[i+wp.Size].Add(&wp.Polynomial.Coefficients[i+wp.Size], &r)
	}
	wp.BlindedSize = newSize

	return wp
}

//----------------------------------------------------
// Evaluation

// Evaluate evaluates p at x.
// The code panics if the function is not in canonical form.
func (p *Polynomial) Evaluate(x fr.Element) fr.Element {

	var r fr.Element
	if p.Basis != Canonical {
		panic("p must be in canonical basis")
	}

	if p.Layout == Regular {
		for i := len(p.Coefficients) - 1; i >= 0; i-- {
			r.Mul(&r, &x).Add(&r, &p.Coefficients[i])
		}
	} else {
		nn := uint64(64 - bits.TrailingZeros(uint(len(p.Coefficients))))
		for i := len(p.Coefficients) - 1; i >= 0; i-- {
			iRev := bits.Reverse64(uint64(i)) >> nn
			r.Mul(&r, &x).Add(&r, &p.Coefficients[iRev])
		}
	}

	return r

}

// Evaluate evaluates p at x.
// The code panics if the function is not in canonical form.
func (wp *WrappedPolynomial) Evaluate(x fr.Element) fr.Element {

	if wp.Shift == 0 {
		return wp.Polynomial.Evaluate(x)
	}

	// TODO find a way to retrieve the root properly instead of re generating the fft domain
	d := fft.NewDomain(uint64(wp.Size))
	var g fr.Element
	if wp.Shift <= 5 {
		g = smallExp(d.Generator, wp.Shift)
		x.Mul(&x, &g)
		return wp.Polynomial.Evaluate(x)
	}

	bs := big.NewInt(int64(wp.Shift))
	g = *g.Exp(g, bs)
	x.Mul(&x, &g)
	return wp.Polynomial.Evaluate(x)
}

// Clone returns a deep copy of wp. The underlying polynomial is cloned;
// see also ShallowClone to perform a ShallowClone on the underlying polynomial.
func (wp *WrappedPolynomial) Clone() *WrappedPolynomial {
	var res WrappedPolynomial
	res.Polynomial = wp.Polynomial.Clone()
	res.Shift = wp.Shift
	res.Size = wp.Size
	return &res
}

// ShallowClone returns a shallow copy of wp. The underlying polynomial coefficient
// is NOT cloned and both objects will point to the same data structure.
func (wp *WrappedPolynomial) ShallowClone() *WrappedPolynomial {
	res := *wp
	return &res
}

// GetCoeff returns the i-th entry of wp, taking the layout in account.
func (wp *WrappedPolynomial) GetCoeff(i int) fr.Element {

	n := len(wp.Polynomial.Coefficients)
	rho := n / wp.Size
	if wp.Polynomial.Form.Layout == Regular {
		return wp.Polynomial.Coefficients[(i+rho*wp.Shift)%n]
	} else {
		nn := uint64(64 - bits.TrailingZeros(uint(n)))
		iRev := bits.Reverse64(uint64((i+rho*wp.Shift)%n)) >> nn
		return wp.Polynomial.Coefficients[iRev]
	}

}

//----------------------------------------------------
// ToRegular

func (p *Polynomial) ToRegular() *Polynomial {
	if p.Layout == Regular {
		return p
	}
	fft.BitReverse(p.Coefficients)
	p.Layout = Regular
	return p
}

//----------------------------------------------------
// ToBitReverse

func (p *Polynomial) ToBitReverse() *Polynomial {
	if p.Layout == BitReverse {
		return p
	}
	fft.BitReverse(p.Coefficients)
	p.Layout = BitReverse
	return p
}

//----------------------------------------------------
// toLagrange

// the numeration corresponds to the following formatting:
// num = int(p.Basis)*2 + int(p.Layout)

// CANONICAL REGULAR
func (p *Polynomial) toLagrange0(d *fft.Domain) *Polynomial {
	p.Basis = Lagrange
	p.Layout = BitReverse
	d.FFT(p.Coefficients, fft.DIF)
	return p
}

// CANONICAL BITREVERSE
func (p *Polynomial) toLagrange1(d *fft.Domain) *Polynomial {
	p.Basis = Lagrange
	p.Layout = Regular
	d.FFT(p.Coefficients, fft.DIT)
	return p
}

// LAGRANGE_COSET REGULAR
func (p *Polynomial) toLagrange4(d *fft.Domain) *Polynomial {
	p.Basis = Lagrange
	p.Layout = Regular
	d.FFTInverse(p.Coefficients, fft.DIF, true)
	d.FFT(p.Coefficients, fft.DIT)
	return p
}

// LAGRANGE_COSET BITREVERSE
func (p *Polynomial) toLagrange5(d *fft.Domain) *Polynomial {
	p.Basis = Lagrange
	p.Layout = BitReverse
	d.FFTInverse(p.Coefficients, fft.DIT, true)
	d.FFT(p.Coefficients, fft.DIF)
	return p
}

// Set p to q in Lagrange form and returns it.
func (p *Polynomial) ToLagrange(d *fft.Domain) *Polynomial {
	id := p.Form
	resize(p, d.Cardinality)
	switch id {
	case canonicalRegular:
		return p.toLagrange0(d)
	case canonicalBitReverse:
		return p.toLagrange1(d)
	case lagrangeRegular, lagrangeBitReverse:
		return p
	case lagrangeCosetRegular:
		return p.toLagrange4(d)
	case lagrangeCosetBitReverse:
		return p.toLagrange5(d)
	default:
		panic("unknown ID")
	}
}

//----------------------------------------------------
// toCanonical

// LAGRANGE REGULAR
func (p *Polynomial) toCanonical2(d *fft.Domain) *Polynomial {
	p.Basis = Canonical
	p.Layout = BitReverse
	d.FFTInverse(p.Coefficients, fft.DIF)
	return p
}

// LAGRANGE BITREVERSE
func (p *Polynomial) toCanonical3(d *fft.Domain) *Polynomial {
	p.Basis = Canonical
	p.Layout = Regular
	d.FFTInverse(p.Coefficients, fft.DIT)
	return p
}

// LAGRANGE_COSET REGULAR
func (p *Polynomial) toCanonical4(d *fft.Domain) *Polynomial {
	p.Basis = Canonical
	p.Layout = BitReverse
	d.FFTInverse(p.Coefficients, fft.DIF, true)
	return p
}

// LAGRANGE_COSET BITREVERSE
func (p *Polynomial) toCanonical5(d *fft.Domain) *Polynomial {
	p.Basis = Canonical
	p.Layout = Regular
	d.FFTInverse(p.Coefficients, fft.DIT, true)
	return p
}

// ToCanonical Sets p to q, in canonical form and returns it.
func (p *Polynomial) ToCanonical(d *fft.Domain) *Polynomial {
	id := p.Form
	resize(p, d.Cardinality)
	switch id {
	case canonicalRegular, canonicalBitReverse:
		return p
	case lagrangeRegular:
		return p.toCanonical2(d)
	case lagrangeBitReverse:
		return p.toCanonical3(d)
	case lagrangeCosetRegular:
		return p.toCanonical4(d)
	case lagrangeCosetBitReverse:
		return p.toCanonical5(d)
	default:
		panic("unknown ID")
	}
}

//-----------------------------------------------------
// ToLagrangeCoset

func resize(p *Polynomial, newSize uint64) {
	z := make([]fr.Element, int(newSize)-len(p.Coefficients))
	p.Coefficients = append(p.Coefficients, z...)
}

// CANONICAL REGULAR
func (p *Polynomial) toLagrangeCoset0(d *fft.Domain) *Polynomial {
	p.Basis = LagrangeCoset
	p.Layout = BitReverse
	d.FFT(p.Coefficients, fft.DIF, true)
	return p
}

// CANONICAL BITREVERSE
func (p *Polynomial) toLagrangeCoset1(d *fft.Domain) *Polynomial {
	p.Basis = LagrangeCoset
	p.Layout = Regular
	d.FFT(p.Coefficients, fft.DIT, true)
	return p
}

// LAGRANGE REGULAR
func (p *Polynomial) toLagrangeCoset2(d *fft.Domain) *Polynomial {
	p.Basis = LagrangeCoset
	p.Layout = Regular
	d.FFTInverse(p.Coefficients, fft.DIF)
	d.FFT(p.Coefficients, fft.DIT, true)
	return p
}

// LAGRANGE BITREVERSE
func (p *Polynomial) toLagrangeCoset3(d *fft.Domain) *Polynomial {
	p.Basis = LagrangeCoset
	p.Layout = BitReverse
	d.FFTInverse(p.Coefficients, fft.DIT)
	d.FFT(p.Coefficients, fft.DIF, true)
	return p
}

// ToLagrangeCoset Sets p to q, in LagrangeCoset form and returns it.
func (p *Polynomial) ToLagrangeCoset(d *fft.Domain) *Polynomial {
	id := p.Form
	resize(p, d.Cardinality)
	switch id {
	case canonicalRegular:
		return p.toLagrangeCoset0(d)
	case canonicalBitReverse:
		return p.toLagrangeCoset1(d)
	case lagrangeRegular:
		return p.toLagrangeCoset2(d)
	case lagrangeBitReverse:
		return p.toLagrangeCoset3(d)
	case lagrangeCosetRegular, lagrangeCosetBitReverse:
		return p
	default:
		panic("unknown ID")
	}
}
