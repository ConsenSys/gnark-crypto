// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package iop

import (
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bls24-317/fr"
)

func TestEvaluate(t *testing.T) {

	f := func(_ int, x ...fr.Element) fr.Element {
		var a fr.Element
		a.Add(&x[0], &x[1]).Add(&a, &x[2])
		return a
	}

	size := 64
	u := make([]fr.Element, size)
	v := make([]fr.Element, size)
	w := make([]fr.Element, size)
	for i := 0; i < size; i++ {
		u[i].SetUint64(uint64(i))
		v[i].SetUint64(uint64(i + 1))
		w[i].SetUint64(uint64(i + 2))
	}
	r := make([]fr.Element, size)
	for i := 0; i < size; i++ {
		r[i].SetUint64(uint64(3 * (i + 1)))
	}
	form := Form{Layout: Regular, Basis: Canonical}
	wu := NewPolynomial(&u, form)
	wv := NewPolynomial(&v, form)
	ww := NewPolynomial(&w, form)

	rr, err := Evaluate(f, nil, form, wu, wv, ww)
	if err != nil {
		t.Fatal(err)
	}

	wu.ToBitReverse()
	rrb, err := Evaluate(f, nil, form, wu, wv, ww)
	if err != nil {
		t.Fatal(err)
	}

	wv.ToBitReverse()
	ww.ToBitReverse()
	rrc, err := Evaluate(f, nil, form, wu, wv, ww)
	if err != nil {
		t.Fatal(err)
	}

	// compare with the expected result
	for i := 0; i < size; i++ {
		if !rr.Coefficients()[i].Equal(&r[i]) {
			t.Fatal("error evaluation")
		}
		if !rrb.Coefficients()[i].Equal(&r[i]) {
			t.Fatal("error evaluation")
		}
		if !rrc.Coefficients()[i].Equal(&r[i]) {
			t.Fatal("error evaluation")
		}

	}
}
