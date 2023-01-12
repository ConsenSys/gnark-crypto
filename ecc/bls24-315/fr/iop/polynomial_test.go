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

	"github.com/consensys/gnark-crypto/ecc/bls24-315/fr"
	"github.com/consensys/gnark-crypto/ecc/bls24-315/fr/fft"
)

func TestEvaluateSinglePoint(t *testing.T) {

	// monomial
	{
		var m monomial
		m.coeff.SetInt64(3)
		m.exponents = make([]int, 3)
		m.exponents[0] = 1
		m.exponents[1] = 2
		m.exponents[2] = 3

		x := make([]fr.Element, 3)
		x[0].SetUint64(2)
		x[1].SetUint64(2)
		x[2].SetUint64(2)

		var res fr.Element
		res.SetUint64(192)

		r := m.evaluate(x)

		if !r.Equal(&res) {
			t.Fatal("evaluation monomial failed")
		}
	}

	// multivariate polynomial
	{
		p := make(MultivariatePolynomial, 3)

		p[0].coeff.SetUint64(1)
		p[0].exponents = []int{1, 1, 1}

		p[1].coeff.SetUint64(2)
		p[1].exponents = []int{2, 2, 2}

		p[2].coeff.SetUint64(0)
		p[2].exponents = []int{3, 2, 1}

		x := make([]fr.Element, 3)
		x[0].SetUint64(2)
		x[1].SetUint64(2)
		x[2].SetUint64(2)

		var res fr.Element
		res.SetUint64(136)

		r := p.EvaluateSinglePoint(x)
		if !r.Equal(&res) {
			t.Fatal("evaluation multivariate polynomial failed")
		}

	}

}

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

// CANONICAL REGULAR
func fromLagrange0(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = Canonical
	r.Layout = Regular
	d.FFTInverse(r.Coefficients, fft.DIF)
	fft.BitReverse(r.Coefficients)
	return r
}

// CANONICAL BITREVERSE
func fromLagrange1(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = Canonical
	r.Layout = BitReverse
	d.FFTInverse(r.Coefficients, fft.DIF)
	return r
}

// LAGRANGE REGULAR
func fromLagrange2(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = Lagrange
	r.Layout = Regular
	return r
}

// LAGRANGE BITREVERSE
func fromLagrange3(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = Lagrange
	r.Layout = BitReverse
	fft.BitReverse(r.Coefficients)
	return r
}

// LAGRANGE_COSET REGULAR
func fromLagrange4(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = LagrangeCoset
	r.Layout = Regular
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	return r
}

// LAGRANGE_COSET BITREVERSE
func fromLagrange5(p *Polynomial, d *fft.Domain) *Polynomial {
	r := p.Copy()
	r.Basis = LagrangeCoset
	r.Layout = BitReverse
	d.FFTInverse(r.Coefficients, fft.DIF)
	d.FFT(r.Coefficients, fft.DIT, true)
	fft.BitReverse(r.Coefficients)
	return r
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
	default:
		panic("unknown id")
	}
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

	// CANONICAL REGULAR
	{
		_p := fromLagrange(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE
	{
		_p := fromLagrange1(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR
	{
		_p := fromLagrange2(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE
	{
		_p := fromLagrange3(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitReverse")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET REGULAR
	{
		_p := fromLagrange4(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is Regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET BITREVERSE
	{
		_p := fromLagrange5(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToLagrange(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Lagrange {
			t.Fatal("expected basis is Lagrange")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is BitRervese")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

}

// CANONICAL REGULAR
func fromCanonical0(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = Canonical
	_p.Layout = Regular
	return _p
}

// CANONICAL BITREVERSE
func fromCanonical1(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = Canonical
	_p.Layout = BitReverse
	return _p
}

// LAGRANGE REGULAR
func fromCanonical2(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = Lagrange
	_p.Layout = Regular
	d.FFT(_p.Coefficients, fft.DIF)
	fft.BitReverse(_p.Coefficients)
	return _p
}

// LAGRANGE BITREVERSE
func fromCanonical3(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = Lagrange
	_p.Layout = BitReverse
	d.FFT(_p.Coefficients, fft.DIF)
	return _p
}

// LAGRANGE_COSET REGULAR
func fromCanonical4(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = LagrangeCoset
	_p.Layout = Regular
	d.FFT(_p.Coefficients, fft.DIF, true)
	fft.BitReverse(_p.Coefficients)
	return _p
}

// LAGRANGE_COSET BITREVERSE
func fromCanonical5(p *Polynomial, d *fft.Domain) *Polynomial {
	_p := p.Copy()
	_p.Basis = LagrangeCoset
	_p.Layout = BitReverse
	d.FFT(_p.Coefficients, fft.DIF, true)
	return _p
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

	// CANONICAL REGULAR
	{
		_p := fromCanonical0(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToCanonical(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// CANONICAL BITREVERSE
	{
		_p := fromCanonical1(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToCanonical(_p, domain)
		// if !reflect.DeepEqual(_p, backup) {
		// 	t.Fatal("locked polynomial should not be modified")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE REGULAR
	{
		_p := fromCanonical2(&p, domain)
		var q Polynomial
		q.ToCanonical(_p, domain)
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitReverse")
		}
		fft.BitReverse(q.Coefficients)
		if !cmpCoefficents(p.Coefficients, q.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE BITREVERSE
	{
		_p := fromCanonical3(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToCanonical(_p, domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

	// LAGRANGE_COSET REGULAR
	{
		_p := fromCanonical4(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToCanonical(_p, domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != BitReverse {
			t.Fatal("epxected layout is bitreverse")
		}
		// fft.BitReverse(q.Coefficients)
		// if !cmpCoefficents(q.Coefficients, p.Coefficients) {
		// 	t.Fatal("wrong coefficients")
		// }
	}

	// LAGRANGE_COSET BITREVERSE
	{
		_p := fromCanonical5(&p, domain)
		// backup := copyPoly(*_p)
		var q Polynomial
		q.ToCanonical(_p, domain)
		// if !reflect.DeepEqual(backup, *_p){
		// 	t.Fatal("")
		// }
		if q.Basis != Canonical {
			t.Fatal("expected basis is canonical")
		}
		if q.Layout != Regular {
			t.Fatal("epxected layout is regular")
		}
		if !cmpCoefficents(q.Coefficients, p.Coefficients) {
			t.Fatal("wrong coefficients")
		}
	}

}
