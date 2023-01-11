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

package iop

import (
	"errors"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr/fft"
)

//-----------------------------------------------------
// univariate polynomials

// Enum to tell in which basis a polynomial is represented.
type Basis int64

const (
	Canonical Basis = iota
	Lagrange
	LagrangeCoset
)

// Enum to tell if a polynomial is in bit reverse form or
// in the regular form.
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

// Polynomial represents a polynomial, the vector of coefficients
// along with the basis and the layout.
type Polynomial struct {
	Coefficients []fr.Element
	Form
}

// return a copy of p
// return a copy of p
func (p *Polynomial) Copy() *Polynomial {
	size := len(p.Coefficients)
	var r Polynomial
	r.Coefficients = make([]fr.Element, size)
	copy(r.Coefficients, p.Coefficients)
	r.Form = p.Form
	return &r
}

// return an ID corresponding to the polynomial extra data
func getShapeID(p Polynomial) int {
	return int(p.Basis)*2 + int(p.Layout)
}

// WrappedPolynomial wrapps a polynomial so that it is
// interpreted as P'(X)=P(\omega^{s}X)
type WrappedPolynomial struct {
	*Polynomial
	Shift int
}

//----------------------------------------------------
// ToRegular

func (p *Polynomial) ToRegular(q *Polynomial) *Polynomial {

	if p != q {
		*p = *q.Copy()
	}
	if p.Layout == Regular {
		return p
	}
	fft.BitReverse(p.Coefficients)
	p.Layout = Regular
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

// LAGRANGE REGULAR
func (p *Polynomial) toLagrange2(d *fft.Domain) *Polynomial {
	return p
}

// LAGRANGE BITREVERSE
func (p *Polynomial) toLagrange3(d *fft.Domain) *Polynomial {
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
func (p *Polynomial) ToLagrange(q *Polynomial, d *fft.Domain) *Polynomial {
	id := getShapeID(*q)
	if q != p {
		*p = *q.Copy()
	}
	switch id {
	case 0:
		return p.toLagrange0(d)
	case 1:
		return p.toLagrange1(d)
	case 2:
		return p.toLagrange2(d)
	case 3:
		return p.toLagrange3(d)
	case 4:
		return p.toLagrange4(d)
	case 5:
		return p.toLagrange5(d)
	default:
		panic("unknown ID")
	}
}

//----------------------------------------------------
// toCanonical

// CANONICAL REGULAR
func (p *Polynomial) toCanonical0(d *fft.Domain) *Polynomial {
	return p
}

// CANONICAL BITREVERSE
func (p *Polynomial) toCanonical1(d *fft.Domain) *Polynomial {
	return p
}

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
func (p *Polynomial) ToCanonical(q *Polynomial, d *fft.Domain) *Polynomial {
	id := getShapeID(*q)
	if q != p {
		*p = *q.Copy()
	}
	switch id {
	case 0:
		return p.toCanonical0(d)
	case 1:
		return p.toCanonical1(d)
	case 2:
		return p.toCanonical2(d)
	case 3:
		return p.toCanonical3(d)
	case 4:
		return p.toCanonical4(d)
	case 5:
		return p.toCanonical5(d)
	default:
		panic("unknown ID")
	}
}

//-----------------------------------------------------
// toLagrangeCoset

// /!\ those functions are for internal use only. It is assumed
// that the polynomial is in canonical basis and regular layout.
// Usually those functions are used to evaluate a polynomial of
// degree n on a coset of size k*n to compute quotients by X^n-1
// (hence the evaluation on the coset, to avoid the zeroes of X^n-1
// when dividing).

// toLagrangeCoset sets p to q in Lagrange coset, possibly resized to
// fit the size of the domain d.
func (p *Polynomial) toLagrangeCoset(q *Polynomial, d *fft.Domain) *Polynomial {

	if p != q {
		*p = *q.Copy()
	}

	// When using a WrappedPolynomial, it's possible that q has already
	// been expressed in LagrangeCoset, BitReverse layout. In that case
	// we do nothing and return p directly.
	if p.Basis == LagrangeCoset && p.Layout == BitReverse {
		return p
	}

	// resizing if needed.
	if len(p.Coefficients) < int(d.Cardinality) {
		z := make([]fr.Element, int(d.Cardinality)-len(p.Coefficients))
		p.Coefficients = append(p.Coefficients, z...)
	}

	// computing the transform, without any checks: since this function
	// is used internally only, it is assumed that the expected form
	// provided by the caller is correct, that is it is in Canonical basis,
	// Regular layout.
	d.FFT(p.Coefficients, fft.DIF, true)
	p.Basis = LagrangeCoset
	p.Layout = BitReverse

	return p
}

//-----------------------------------------------------
// multivariate polynomials

// errors related to the polynomials.
var ErrInconsistantNumberOfVariable = errors.New("the number of variables is not consistant")

// monomial represents a monomial encoded as
// coeff*X₁^{i₁}*..*X_n^{i_n} if exponents = [i₁,..iₙ]
type monomial struct {
	coeff     fr.Element
	exponents []int
}

// it is supposed that the number of variables matches
func (m monomial) evaluate(x []fr.Element) fr.Element {

	var res, tmp fr.Element

	nbVars := len(x)
	res.SetOne()
	for i := 0; i < nbVars; i++ {
		if m.exponents[i] <= 5 {
			tmp = smallExp(x[i], m.exponents[i])
			res.Mul(&res, &tmp)
			continue
		}
		bi := big.NewInt(int64(i))
		tmp.Exp(x[i], bi)
		res.Mul(&res, &tmp)
	}
	res.Mul(&res, &m.coeff)

	return res

}

// reprensents a multivariate polynomial as a list of monomial,
// the multivariate polynomial being the sum of the monomials.
type MultivariatePolynomial []monomial

// degree returns the total degree
func (m *MultivariatePolynomial) Degree() uint64 {
	r := 0
	for i := 0; i < len(*m); i++ {
		t := 0
		for j := 0; j < len((*m)[i].exponents); j++ {
			t += (*m)[i].exponents[j]
		}
		if t > r {
			r = t
		}
	}
	return uint64(r)
}

// AddMonomial adds a monomial to m. If m is empty, the monomial is
// added no matter what. But if m is already populated, an error is
// returned if len(e)\neq size of the previous list of exponents. This
// ensure that the number of variables is given by the size of any of
// the slices of exponent in any monomial.
func (m *MultivariatePolynomial) AddMonomial(c fr.Element, e []int) error {

	// if m is empty, we add the first monomial.
	if len(*m) == 0 {
		r := monomial{c, e}
		*m = append(*m, r)
		return nil
	}

	// at this stage all of exponennt in m are supposed to be of
	// the same size.
	if len((*m)[0].exponents) != len(e) {
		return ErrInconsistantNumberOfVariable
	}
	r := monomial{c, e}
	*m = append(*m, r)
	return nil

}

// evaluate a multivariate polynomial in x
// /!\ It is assumed that the multivariate polynomial has been
// built correctly, that is the sizes of the slices in exponents
// are the same /!\
func (m *MultivariatePolynomial) Evaluate(x []fr.Element) fr.Element {

	var res fr.Element

	for i := 0; i < len(*m); i++ {
		tmp := (*m)[i].evaluate(x)
		res.Add(&res, &tmp)
	}
	return res
}
