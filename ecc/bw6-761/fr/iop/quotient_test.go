// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package iop

import (
	"testing"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-761/fr"
	"github.com/consensys/gnark-crypto/ecc/bw6-761/fr/fft"
)

// computes x₃ in h(x₁,x₂,x₃) = x₁^{2}*x₂ + x₃ - x₁^{3}
// from x₁ and x₂.
func computex3(x []fr.Element) fr.Element {

	var a, b fr.Element
	a.Square(&x[0]).Mul(&a, &x[1])
	b.Square(&x[0]).Mul(&b, &x[0])
	a.Sub(&b, &a)
	return a

}

func buildPoly(size int, form Form) *Polynomial {
	v := make([]fr.Element, size)
	return NewPolynomial(&v, form)
}

func TestDivideByXMinusOne(t *testing.T) {

	f := func(_ int, x ...fr.Element) fr.Element {
		var a, b fr.Element
		a.Square(&x[0]).Mul(&a, &x[1]).Add(&a, &x[2])
		b.Square(&x[0]).Mul(&b, &x[0])
		a.Sub(&a, &b)
		return a
	}

	// create the multivariate polynomial h
	// h(x₁,x₂,x₃) = x₁^{2}*x₂ + x₃ - x₁^{3}
	nbEntries := 3

	// create an instance (f_i) where h holds
	sizeSystem := 8

	form := Form{Basis: Lagrange, Layout: Regular}

	entries := make([]*Polynomial, nbEntries)
	entries[0] = buildPoly(sizeSystem, form)
	entries[1] = buildPoly(sizeSystem, form)
	entries[2] = buildPoly(sizeSystem, form)

	for i := 0; i < sizeSystem; i++ {

		entries[0].Coefficients()[i].SetRandom()
		entries[1].Coefficients()[i].SetRandom()
		tmp := computex3(
			[]fr.Element{entries[0].Coefficients()[i],
				entries[1].Coefficients()[i]})
		entries[2].Coefficients()[i].Set(&tmp)

		x := []fr.Element{
			entries[0].GetCoeff(i),
			entries[1].GetCoeff(i),
			entries[2].GetCoeff(i),
		}
		a := f(0, x...)
		if !a.IsZero() {
			t.Fatal("system does not vanish on x^n-1")
		}
	}

	// compute the quotient where the entries are in Regular layout
	var domains [2]*fft.Domain
	domains[0] = fft.NewDomain(uint64(sizeSystem))
	domains[1] = fft.NewDomain(ecc.NextPowerOfTwo(uint64(3 * sizeSystem)))

	entries[0].ToCanonical(domains[0]).
		ToRegular().
		ToLagrangeCoset(domains[1]).
		ToRegular()

	entries[1].ToCanonical(domains[0]).
		ToRegular().
		ToLagrangeCoset(domains[1]).
		ToRegular()

	entries[2].ToCanonical(domains[0]).
		ToRegular().
		ToLagrangeCoset(domains[1]).
		ToRegular()

	expectedForm := Form{Layout: BitReverse, Basis: LagrangeCoset}
	h, err := Evaluate(f, nil, expectedForm, entries...)
	if err != nil {
		t.Fatal(err)
	}

	q, err := DivideByXMinusOne(h, domains)
	if err != nil {
		t.Fatal(err)
	}

	// evaluate the quotient at a random point and check that
	// the relation holds.
	var x fr.Element
	x.SetRandom()
	qx := q.Evaluate(x)
	entries[0].ToCanonical(domains[1])
	entries[1].ToCanonical(domains[1])
	entries[2].ToCanonical(domains[1])
	ax := entries[0].Evaluate(x)
	bx := entries[1].Evaluate(x)
	cx := entries[2].Evaluate(x)
	hx := f(0, ax, bx, cx)

	var xnminusone, one fr.Element
	one.SetOne()
	xnminusone.Set(&x).
		Square(&xnminusone).
		Square(&xnminusone).
		Square(&xnminusone).
		Sub(&xnminusone, &one)
	qx.Mul(&qx, &xnminusone)
	if !qx.Equal(&hx) {
		t.Fatal("error computing quotient")
	}
}
