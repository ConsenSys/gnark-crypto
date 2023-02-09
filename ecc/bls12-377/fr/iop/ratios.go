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
	"errors"
	"math/bits"

	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/fft"
)

// errors related to the computation of the quotient and the ratios.
var (
	ErrMustBeRegular              = errors.New("the layout must be Regular")
	ErrMustBeCanonical            = errors.New("the basis must be Canonical")
	ErrMustBeLagrangeCoset        = errors.New("the basis must be LagrangeCoset")
	ErrInconsistentFormat         = errors.New("the format of the polynomials must be the same")
	ErrInconsistentSize           = errors.New("the sizes of the polynomial must be the same as the size of the domain")
	ErrNumberPolynomials          = errors.New("the number of polynomials in the denominator and the numerator must be the same")
	ErrSizeNotPowerOfTwo          = errors.New("the size of the polynomials must be a power of two")
	ErrInconsistentSizeDomain     = errors.New("the size of the domain must be consistent with the size of the polynomials")
	ErrIncorrectNumberOfVariables = errors.New("the number of variables is incorrect")
)

// Build an 'accumulating ratio' polynomial.
// * numerator list of polynomials that will form the numerator of the ratio
// * denominator list of polynomials that will form the denominator of the ratio
// The polynomials in the denominator and the numerator are expected to be of
// the same size and the size must be a power of 2. The polynomials are given as
// pointers in case the caller wants to FFTInv the polynomials during the process.
// * beta variable at which the numerator and denominators are evaluated
// * expectedForm expected form of the resulting polynomial
// * Return: say beta=β, numerator = [P₁,...,P_m], denominator = [Q₁,..,Q_m]. The function
// returns a polynomial whose evaluation on the j-th root of unity is
// (Π_{k<j}Π_{i<m}(β-Pᵢ(ωᵏ)))/(β-Qᵢ(ωᵏ))
func BuildRatioShuffledVectors(numerator, denominator []*Polynomial, beta fr.Element, expectedForm Form, domain *fft.Domain) (Polynomial, error) {

	var res Polynomial

	// check that len(numerator)=len(denominator)
	if len(numerator) != len(denominator) {
		return res, ErrNumberPolynomials
	}
	nbPolynomials := len(numerator)

	// check that the sizes are consistent
	err := checkSize(numerator, denominator)
	if err != nil {
		return res, err
	}

	// create the domain + some checks on the sizes of the polynomials
	n := len(numerator[0].Coefficients)
	domain, err = buildDomain(n, domain)
	if err != nil {
		return res, err
	}

	// put every polynomials in Lagrange form. Also make sure
	// that we don't modify the slices numerator and denominator, but
	// only their entries. If the polynomials are unlocked, the
	// entries of the slices numerator and denominator will be
	// modified.
	for i := 0; i < nbPolynomials; i++ {
		numerator[i].ToLagrange(domain)
		denominator[i].ToLagrange(domain)
	}

	// build the ratio (careful with the indices of
	// the polynomials which are bit reversed)
	res.Coefficients = make([]fr.Element, n)
	t := make([]fr.Element, n)
	res.Coefficients[0].SetOne()
	t[0].SetOne()
	var a, b, c, d fr.Element

	nn := uint64(64 - bits.TrailingZeros(uint(n)))
	for i := 0; i < n-1; i++ {

		b.SetOne()
		d.SetOne()

		iRev := bits.Reverse64(uint64(i)) >> nn

		for j := 0; j < nbPolynomials; j++ {

			if numerator[j].Layout == BitReverse {
				a.Sub(&beta, &numerator[j].Coefficients[iRev])
			} else {
				a.Sub(&beta, &numerator[j].Coefficients[i])
			}
			b.Mul(&b, &a)

			if denominator[j].Layout == BitReverse {
				c.Sub(&beta, &denominator[j].Coefficients[iRev])
			} else {
				c.Sub(&beta, &denominator[j].Coefficients[i])
			}
			d.Mul(&d, &c)
		}
		// b = Πₖ (β-Pₖ(ωⁱ⁻¹))
		// d = Πₖ (β-Qₖ(ωⁱ⁻¹))

		res.Coefficients[i+1].Mul(&res.Coefficients[i], &b)
		t[i+1].Mul(&t[i], &d)

	}

	t = fr.BatchInvert(t)
	for i := 1; i < n; i++ {
		res.Coefficients[i].Mul(&res.Coefficients[i], &t[i])
	}

	res.Basis = expectedForm.Basis
	res.Layout = expectedForm.Layout

	// at this stage the result is in Lagrange form, Regular layout
	putInExpectedFormFromLagrangeRegular(&res, domain, expectedForm)

	return res, nil
}

// BuildRatioCopyConstraint builds the accumulating ratio polynomial to prove that
// [P₁ ∥ .. ∥ P_{n—1}] is invariant by the permutation \sigma.
// Namely it returns the polynomial Z whose evaluation on the j-th root of unity is
// Z(ω^j) = Π_{i<j}(Π_{k<n}(P_k(ω^i)+β*u^k+γ))/(P_k(ω^i)+σ(kn+i)+γ)))
// * entries list of polynomials whose evaluation are invariant under \sigma
// * beta, gamma challenges
// * expectedForm expected form of the resulting polynomial
func BuildRatioCopyConstraint(
	entries []*Polynomial,
	permutation []int64,
	beta, gamma fr.Element,
	expectedForm Form,
	domain *fft.Domain) (Polynomial, error) {

	var res Polynomial

	nbPolynomials := len(entries)

	// check that the sizes are consistent
	err := checkSize(entries)
	if err != nil {
		return res, err
	}

	// create the domain + some checks on the sizes of the polynomials
	n := len(entries[0].Coefficients)
	domain, err = buildDomain(n, domain)
	if err != nil {
		return res, err
	}

	// put every polynomials in Lagrange form. Also make sure
	// that we don't modify the slice entries
	for i := 0; i < nbPolynomials; i++ {
		entries[i].ToLagrange(domain)
	}

	// get the support for the permutation
	evaluationIDSmallDomain := getSupportIdentityPermutation(nbPolynomials, domain)

	// build the ratio (careful with the indices of
	// the polynomials which are bit reversed)
	res.Coefficients = make([]fr.Element, n)
	t := make([]fr.Element, n)
	res.Coefficients[0].SetOne()
	t[0].SetOne()
	var a, b, c, d fr.Element

	nn := uint64(64 - bits.TrailingZeros(uint(n)))
	for i := 0; i < n-1; i++ {

		b.SetOne()
		d.SetOne()

		iRev := int(bits.Reverse64(uint64(i)) >> nn)

		for j, p := range entries {
			idx := i
			if p.Layout == BitReverse {
				idx = iRev
			}

			a.Mul(&beta, &evaluationIDSmallDomain[i+j*n]).
				Add(&a, &gamma).
				Add(&a, &p.Coefficients[idx])

			b.Mul(&b, &a)

			c.Mul(&beta, &evaluationIDSmallDomain[permutation[i+j*n]]).
				Add(&c, &gamma).
				Add(&c, &p.Coefficients[idx])
			d.Mul(&d, &c)
		}

		// b = Πⱼ(Pⱼ(ωⁱ)+β*ωⁱνʲ+γ)
		// d = Πⱼ(Qⱼ(ωⁱ)+β*σ(j*n+i)+γ)
		res.Coefficients[i+1].Mul(&res.Coefficients[i], &b)
		t[i+1].Mul(&t[i], &d)
	}

	t = fr.BatchInvert(t)
	for i := 1; i < n; i++ {
		res.Coefficients[i].Mul(&res.Coefficients[i], &t[i])
	}

	// at this stage the result is in Lagrange form, Regular layout
	putInExpectedFormFromLagrangeRegular(&res, domain, expectedForm)

	return res, nil

}

func putInExpectedFormFromLagrangeRegular(p *Polynomial, domain *fft.Domain, expectedForm Form) {

	p.Basis = expectedForm.Basis
	p.Layout = expectedForm.Layout

	if expectedForm.Basis == Canonical {
		domain.FFTInverse(p.Coefficients, fft.DIF)
		if expectedForm.Layout == Regular {
			fft.BitReverse(p.Coefficients)
		}
		return
	}

	if expectedForm.Basis == LagrangeCoset {
		domain.FFTInverse(p.Coefficients, fft.DIF)
		domain.FFT(p.Coefficients, fft.DIT, true)
		if expectedForm.Layout == BitReverse {
			fft.BitReverse(p.Coefficients)
		}
		return
	}

	if expectedForm.Layout == BitReverse {
		fft.BitReverse(p.Coefficients)
	}

}

// check that the polynomials are of the same size.
// It assumes that pols contains slices of the same size.
func checkSize(pols ...[]*Polynomial) error {

	// check sizes between one another
	m := len(pols)
	n := len(pols[0][0].Coefficients)
	for i := 0; i < m; i++ {
		for j := 0; j < len(pols); j++ {
			if len(pols[i][j].Coefficients) != n {
				return ErrInconsistentSize
			}
		}
	}

	return nil
}

// buildDomain builds the fft domain necessary to do FFTs.
// n is the cardinality of the domain, it must be a power of 2.
func buildDomain(n int, domain *fft.Domain) (*fft.Domain, error) {

	// check if the sizes are a power of 2
	if n&(n-1) != 0 {
		return nil, ErrSizeNotPowerOfTwo
	}

	// if the domain doesn't exist we create it.
	if domain == nil {
		domain = fft.NewDomain(uint64(n))
	}

	// in case domain was not nil, it must match the size of the polynomials.
	if domain.Cardinality != uint64(n) {
		return nil, ErrInconsistentSizeDomain
	}

	return domain, nil
}

// getSupportIdentityPermutation returns the support on which the permutation acts.
// Concretely it's X evaluated on
// [1,ω,..,ωˢ⁻¹,g,g*ω,..,g*ωˢ⁻¹,..,gⁿ⁻¹,gⁿ⁻¹*ω,..,gⁿ⁻¹*ωˢ⁻¹]
// nbCopies is the number of cosets of the roots of unity that are needed, including the set of
// roots of unity itself.
func getSupportIdentityPermutation(nbCopies int, domain *fft.Domain) []fr.Element {

	res := make([]fr.Element, uint64(nbCopies)*domain.Cardinality)
	sizePoly := int(domain.Cardinality)

	res[0].SetOne()
	for i := 0; i < sizePoly-1; i++ {
		res[i+1].Mul(&res[i], &domain.Generator)
	}
	for i := 1; i < nbCopies; i++ {
		copy(res[i*sizePoly:], res[(i-1)*sizePoly:i*int(domain.Cardinality)])
		for j := 0; j < sizePoly; j++ {
			res[i*sizePoly+j].Mul(&res[i*sizePoly+j], &domain.FrMultiplicativeGen)
		}
	}

	return res
}
