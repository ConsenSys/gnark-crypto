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
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/fft"
)

func randomVector(size int) []fr.Element {

	r := make([]fr.Element, size)
	for i := 0; i < size; i++ {
		r[i].SetRandom()
	}
	return r
}

// list of functions to turn a polynomial in Lagrange-regular form
// to all different forms in ordered using this encoding:
// int(p.Basis)*4 + int(p.Layout)*2 + int(p.Status)
// p is in Lagrange/Regular here. This function is for testing purpose
// only.
type TransfoTest func(p Polynomial, d *fft.Domain) Polynomial

// CANONICAL REGULAR LOCKED
func fromLagrange0(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Canonical
	r.Layout = Regular
	r.Status = Locked
	d.FFTInverse(r.Coefficients, fft.DIF)
	fft.BitReverse(r.Coefficients)
	return &r
}

// CANONICAL REGULAR UNLOCKED
func fromLagrange1(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Canonical
	r.Layout = Regular
	r.Status = Unlocked
	d.FFTInverse(r.Coefficients, fft.DIF)
	fft.BitReverse(r.Coefficients)
	return &r
}

// CANONICAL BITREVERSE LOCKED
func fromLagrange2(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Canonical
	r.Layout = BitReverse
	r.Status = Locked
	d.FFTInverse(r.Coefficients, fft.DIF)
	return &r
}

// CANONICAL BITREVERSE UNLOCKED
func fromLagrange3(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Canonical
	r.Layout = BitReverse
	r.Status = Unlocked
	d.FFTInverse(r.Coefficients, fft.DIF)
	return &r
}

// LAGRANGE REGULAR LOCKED
func fromLagrange4(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Lagrange
	r.Layout = Regular
	r.Status = Locked
	return &r
}

// LAGRANGE REGULAR UNLOCKED
func fromLagrange5(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Lagrange
	r.Layout = Regular
	r.Status = Unlocked

	return &r
}

// LAGRANGE BITREVERSE LOCKED
func fromLagrange6(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Lagrange
	r.Layout = BitReverse
	r.Status = Locked
	fft.BitReverse(r.Coefficients)
	return &r
}

// LAGRANGE BITREVERSE UNLOCKED
func fromLagrange7(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = Lagrange
	r.Layout = BitReverse
	r.Status = Unlocked
	fft.BitReverse(r.Coefficients)
	return &r
}

// LAGRANGE_COSET REGULAR LOCKED
func fromLagrange8(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = LagrangeCoset
	r.Layout = Regular
	r.Status = Locked
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	return &r
}

// LAGRANGE_COSET REGULAR UNLOCKED
func fromLagrange9(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = LagrangeCoset
	r.Layout = Regular
	r.Status = Unlocked
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	return &r
}

// LAGRANGE_COSET BITREVERSE LOCKED
func fromLagrange10(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = LagrangeCoset
	r.Layout = BitReverse
	r.Status = Locked
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	fft.BitReverse(r.Coefficients)
	return &r
}

// LAGRANGE_COSET BITREVERSE UNLOCKED
func fromLagrange11(p *Polynomial, d *fft.Domain) *Polynomial {
	r := copyPoly(*p)
	r.Basis = LagrangeCoset
	r.Layout = BitReverse
	r.Status = Unlocked
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	fft.BitReverse(r.Coefficients)
	return &r
}

func fromLagrange(p *Polynomial, d *fft.Domain) *Polynomial {
	id := getShapeID(*p)
	switch id {
	case 0:
		return fromLagrange0(p, d)
	case 1:
		return fromLagrange1(p, d)
	case 2:
		return fromLagrange2(p, d)
	case 3:
		return fromLagrange3(p, d)
	case 4:
		return fromLagrange4(p, d)
	case 5:
		return fromLagrange5(p, d)
	case 6:
		return fromLagrange6(p, d)
	case 7:
		return fromLagrange7(p, d)
	case 8:
		return fromLagrange8(p, d)
	case 9:
		return fromLagrange9(p, d)
	case 10:
		return fromLagrange10(p, d)
	case 11:
		return fromLagrange11(p, d)
	default:
		panic("unknown id")
	}
}

func getCopy(l []Polynomial) []Polynomial {
	r := make([]Polynomial, len(l))
	for i := 0; i < len(l); i++ {
		r[i].Coefficients = make([]fr.Element, len(l[i].Coefficients))
		copy(r[i].Coefficients, l[i].Coefficients)
		r[i] = l[i]
	}
	return r
}

func cmpCoefficents(p, q []fr.Element) bool {
	if len(p) != len(q) {
		return false
	}
	res := true
	for i := 0; i < len(p); i++ {
		res = res && (p[i].Equal(&q[i]))
	}
	return res
}

func TestPutInLagrangeForm(t *testing.T) {

	size := 64
	domain := fft.NewDomain(uint64(size))

	// reference vector in Lagrange-regular form
	c := randomVector(size)
	var p Polynomial
	p.Coefficients = c
	p.Basis = Canonical
	p.Layout = Regular
	p.Status = Locked

	// CANONICAL REGULAR LOCKED
	{
		_p := fromLagrange(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL REGULAR UNLOCKED
	{
		_p := fromLagrange1(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE LOCKED
	{
		_p := fromLagrange2(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE UNLOCKED
	{
		_p := fromLagrange3(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR LOCKED
	{
		_p := fromLagrange4(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Locked {
			t.Fatal("expected status is Locked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR UNLOCKED
	{
		_p := fromLagrange5(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is UnLocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE LOCKED
	{
		_p := fromLagrange6(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		if q.Status != Locked {
			t.Fatal("expected status is Locked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE UNLOCKED
	{
		_p := fromLagrange7(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET REGULAR LOCKED
	{
		_p := fromLagrange8(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET REGULAR UNLOCKED
	{
		_p := fromLagrange9(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET BITREVERSE LOCKED
	{
		_p := fromLagrange10(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitRervese")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET BITREVERSE UNLOCKED
	{
		_p := fromLagrange11(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToLagrange(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitRervese")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is Unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

}

// CANONICAL REGULAR LOCKED
func fromCanonical0(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Canonical
	_p.Layout = Regular
	_p.Status = Locked
	return &_p
}

// CANONICAL REGULAR UNLOCKED
func fromCanonical1(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Canonical
	_p.Layout = Regular
	_p.Status = Unlocked
	return &_p
}

// CANONICAL BITREVERSE LOCKED
func fromCanonical2(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Canonical
	_p.Layout = BitReverse
	_p.Status = Locked
	return &_p
}

// CANONICAL BITREVERSE UNLOCKED
func fromCanonical3(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Canonical
	_p.Layout = BitReverse
	_p.Status = Unlocked
	return &_p
}

// LAGRANGE REGULAR LOCKED
func fromCanonical4(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Lagrange
	_p.Layout = Regular
	_p.Status = Locked
	d.FFT(_p.Coefficients, fft.DIF)
	fft.BitReverse(_p.Coefficients)
	return &_p
}

// LAGRANGE REGULAR UNLOCKED
func fromCanonical5(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Lagrange
	_p.Layout = Regular
	_p.Status = Unlocked
	d.FFT(_p.Coefficients, fft.DIF)
	fft.BitReverse(_p.Coefficients)
	return &_p
}

// LAGRANGE BITREVERSE LOCKED
func fromCanonical6(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Lagrange
	_p.Layout = BitReverse
	_p.Status = Locked
	d.FFT(_p.Coefficients, fft.DIF)
	return &_p
}

// LAGRANGE BITREVERSE UNLOCKED
func fromCanonical7(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = Lagrange
	_p.Layout = BitReverse
	_p.Status = Unlocked
	d.FFT(_p.Coefficients, fft.DIF)
	return &_p
}

// LAGRANGE_COSET REGULAR LOCKED
func fromCanonical8(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = LagrangeCoset
	_p.Layout = Regular
	_p.Status = Locked
	d.FFT(_p.Coefficients, fft.DIF, true)
	fft.BitReverse(_p.Coefficients)
	return &_p
}

// LAGRANGE_COSET REGULAR UNLOCKED
func fromCanonical9(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = LagrangeCoset
	_p.Layout = Regular
	_p.Status = Unlocked
	d.FFT(_p.Coefficients, fft.DIF, true)
	fft.BitReverse(_p.Coefficients)
	return &_p
}

// LAGRANGE_COSET BITREVERSE LOCKED
func fromCanonical10(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = LagrangeCoset
	_p.Layout = BitReverse
	_p.Status = Unlocked
	d.FFT(_p.Coefficients, fft.DIF, true)
	return &_p
}

// LAGRANGE_COSET BITREVERSE UNLOCKED
func fromCanonical11(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := copyPoly(*p)
	_p.Basis = LagrangeCoset
	_p.Layout = BitReverse
	_p.Status = Unlocked
	d.FFT(_p.Coefficients, fft.DIF, true)
	return &_p
}

func TestPutInCanonicalForm(t *testing.T) {

	size := 64
	domain := fft.NewDomain(uint64(size))

	// reference vector in canonical-regular form
	c := randomVector(size)
	var p Polynomial
	p.Coefficients = c
	p.Basis = Canonical
	p.Layout = Regular
	p.Status = Locked

	// CANONICAL REGULAR LOCKED
	{
		_p := fromCanonical0(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Locked {
			t.Fatal("expected status is locked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL REGULAR UNLOCKED
	{
		_p := fromCanonical1(&p, domain)
		q := _p.ToCanonical(domain)
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is locked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE LOCKED
	{
		_p := fromCanonical2(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		if q.Status != Locked {
			t.Fatal("expected status is locked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE UNLOCKED
	{
		_p := fromCanonical3(&p, domain)
		q := _p.ToCanonical(domain)
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unLocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR LOCKED
	{
		_p := fromCanonical4(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(p.Coefficients, q.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR UNLOCKED
	{
		_p := fromCanonical5(&p, domain)
		q := _p.ToCanonical(domain)
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE LOCKED
	{
		_p := fromCanonical6(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE UNLOCKED
	{
		_p := fromCanonical7(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET REGULAR LOCKED
	{
		_p := fromCanonical8(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitreverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		// fft.BitReverse(q.Coefficients)
		// if !cmpCoefficents(q.Coefficients, p.Coefficients) {
		// 	t.Fatal("wrong coefficients")
		// }
	}

	// LAGRANGE_COSET REGULAR UNLOCKED
	{
		_p := fromCanonical9(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitreverse")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET BITREVERSE LOCKED
	{
		_p := fromCanonical10(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET BITREVERSE UNLOCKED
	{
		_p := fromCanonical11(&p, domain)
		// backup := copyPoly(*_p)
		q := _p.ToCanonical(domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if q.Status != Unlocked {
			t.Fatal("expected status is unlocked")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

}
