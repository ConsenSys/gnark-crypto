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
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark-crypto/internal/parallel"
	"math/bits"
)

// Expression represents a multivariate polynomial.
type Expression func(x ...fr.Element) fr.Element

// Evaluate evaluates f on each entry of x. The returned value is
// the vector of evaluations of e on x.
// The form of the result is form.
// The Size field of the result is the same as the one of x[0].
// The blindedSize field of the result is the same as Size.
// The Shift field of the result is 0.
func Evaluate(f Expression, form Form, x ...*Polynomial) (Polynomial, error) {
	var res Polynomial

	if len(x) == 0 {
		return res, errors.New("need at lest one input")
	}

	// check that the sizes are consistent
	n := x[0].coefficients.Len()
	m := len(x)
	for i := 1; i < m; i++ {
		if n != x[i].coefficients.Len() {
			return res, ErrInconsistentSize
		}
	}

	// result coefficients
	r := make([]fr.Element, n)
	idx := func(i int) int {
		return i
	}
	if form.Layout != Regular {
		nn := uint64(64 - bits.TrailingZeros(uint(n)))
		idx = func(i int) int {
			return int(bits.Reverse64(uint64(i)) >> nn)
		}
	}

	parallel.Execute(n, func(start, end int) {
		vx := make([]fr.Element, m)
		for i := start; i < end; i++ {
			for j := 0; j < m; j++ {
				vx[j] = x[j].GetCoeff(i)
			}
			r[idx(i)] = f(vx...)
		}
	})

	res.polynomial = newPolynomial(&r, form)
	res.size = x[0].size
	res.blindedSize = x[0].size
	res.shift = 0

	return res, nil
}
