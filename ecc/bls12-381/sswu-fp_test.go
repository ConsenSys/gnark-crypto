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

package bls12381

import (
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"testing"
)

func TestSqrtRatio(t *testing.T) {
	testSqrtRatio(&fp.Element{0}, &fp.Element{1}, t)
	testSqrtRatio(&fp.Element{1}, &fp.Element{1}, t)

	for i := 0; i < 1000; i++ {
		var u fp.Element
		var v fp.Element
		u.SetRandom()
		v.SetRandom()
		testSqrtRatio(&u, &v, t)
	}
}

func testSqrtRatio(u *fp.Element, v *fp.Element, t *testing.T) {
	var ref fp.Element
	ref.Div(u, v)
	var qrRef bool
	if ref.Legendre() == -1 {
		fp.MulBy11(&ref)
		qrRef = false
	} else {
		qrRef = true
	}

	var seen fp.Element
	qr := sqrtRatio(&seen, u, v)
	seen.Square(&seen)

	if qr != qrRef || seen != ref {
		t.Error(*u, *v)
	}
}
