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
	"math/big"
	"math/bits"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr/fft"
)

// ComputeQuotient returns h(f₁,..,fₙ)/Xⁿ-1 where n=len(f_i).
func ComputeQuotient(entries []WrappedPolynomial, h MultivariatePolynomial, expectedForm Form, domains [2]*fft.Domain) (Polynomial, error) {

	var quotientLagrangeCosetBitReverse Polynomial

	// check that the sizes are consistant
	nbPolynomials := len(entries)
	n := len(entries[0].Coefficients)
	for i := 0; i < nbPolynomials; i++ {
		if len(entries[i].Coefficients) != n {
			return quotientLagrangeCosetBitReverse, ErrInconsistantSize
		}

	}

	// create the domains for the individual polynomials + for the quotient
	var err error
	nbElmts := len(entries[0].Coefficients)
	domains[0], err = buildDomain(nbElmts, domains[0])
	if err != nil {
		return quotientLagrangeCosetBitReverse, err
	}
	nbElmtsExtended := ecc.NextPowerOfTwo(h.Degree() * domains[0].Cardinality)
	domains[1], err = buildDomain(int(nbElmtsExtended), domains[1])
	if err != nil {
		return quotientLagrangeCosetBitReverse, err
	}

	// Note: we will need to interpret the obtained polynomials in
	// canonical form but of degree the size of the big domain. So
	// we will padd the obtained polynomials with zeroes, but this
	// works only if the obtained polynomials are in regular form.
	// So we call bitReverse here if the polynomials are in bit reverse
	// layout.
	for i := 0; i < nbPolynomials; i++ {
		entries[i].ToCanonical(entries[i].Polynomial, domains[0])
		entries[i].ToRegular(entries[i].Polynomial)
	}

	// compute h(f₁,..,fₙ) on a coset
	// entriesLagrangeBigDomain := make([]Polynomial, nbPolynomials)
	for i := 0; i < nbPolynomials; i++ {

		entries[i].toLagrangeCoset(entries[i].Polynomial, domains[1])

	}

	// prepare the evaluations of x^n-1 on the big domain's coset
	xnMinusOneInverseLagrangeCoset := evaluateXnMinusOneDomainBigCoset(domains)
	ratio := int(domains[1].Cardinality / domains[0].Cardinality)

	// compute the division. We take care of the indices of the
	// polnyomials which are bit reversed.
	// The result is temporarily stored in bit reversed Lagrange form,
	// before it is actually transformed into the expected format.
	nbEntries := len(entries)
	x := make([]fr.Element, nbEntries)

	quotientLagrangeCosetBitReverse.Coefficients = make([]fr.Element, nbElmtsExtended)

	rho := int(domains[1].Cardinality / domains[0].Cardinality)

	nn := uint64(64 - bits.TrailingZeros(uint(nbElmtsExtended)))
	// TODO use one goroutines per monomials here
	for i := 0; i < int(nbElmtsExtended); i++ {

		iRev := bits.Reverse64(uint64(i)) >> nn

		for j := 0; j < nbEntries; j++ {

			// take in account the fact that the polynomial mght be shifted...
			iRev := bits.Reverse64(uint64((i+rho*entries[j].Shift))%domains[1].Cardinality) >> nn

			// set the variable. The polynomials in entriesLagrangeBigDomain
			// are in bit reverse.
			x[j].Set(&entries[j].Coefficients[iRev])

		}

		// evaluate h on x
		quotientLagrangeCosetBitReverse.Coefficients[iRev] = h.Evaluate(x)

		// divide by x^n-1 evaluated on the correct point.
		quotientLagrangeCosetBitReverse.Coefficients[iRev].
			Mul(&quotientLagrangeCosetBitReverse.Coefficients[iRev], &xnMinusOneInverseLagrangeCoset[i%ratio])
	}

	// at this stage the result is in Lagrange, bitreversed format.
	// We put it in the expected format.
	putInExpectedFormFromLagrangeCosetBitReversed(&quotientLagrangeCosetBitReverse, domains[1], expectedForm)

	return quotientLagrangeCosetBitReverse, nil
}

// evaluateXnMinusOneDomainBigCoset evalutes Xᵐ-1 on DomainBig coset
func evaluateXnMinusOneDomainBigCoset(domains [2]*fft.Domain) []fr.Element {

	ratio := domains[1].Cardinality / domains[0].Cardinality

	res := make([]fr.Element, ratio)

	expo := big.NewInt(int64(domains[0].Cardinality))
	res[0].Exp(domains[1].FrMultiplicativeGen, expo)

	var t fr.Element
	t.Exp(domains[1].Generator, big.NewInt(int64(domains[0].Cardinality)))

	for i := 1; i < int(ratio); i++ {
		res[i].Mul(&res[i-1], &t)
	}

	var one fr.Element
	one.SetOne()
	for i := 0; i < int(ratio); i++ {
		res[i].Sub(&res[i], &one)
	}

	res = fr.BatchInvert(res)

	return res
}

func putInExpectedFormFromLagrangeCosetBitReversed(p *Polynomial, domain *fft.Domain, expectedForm Form) {

	p.Basis = expectedForm.Basis
	p.Layout = expectedForm.Layout

	if expectedForm.Basis == Canonical {
		domain.FFTInverse(p.Coefficients, fft.DIT, true)
		if expectedForm.Layout == BitReverse {
			fft.BitReverse(p.Coefficients)
		}
		return
	}

	if expectedForm.Basis == Lagrange {
		domain.FFTInverse(p.Coefficients, fft.DIT, true)
		domain.FFT(p.Coefficients, fft.DIF)
		if expectedForm.Layout == Regular {
			fft.BitReverse(p.Coefficients)
		}
		return
	}

	if expectedForm.Layout == Regular {
		fft.BitReverse(p.Coefficients)
	}

}
