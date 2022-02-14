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
	"github.com/consensys/gnark-crypto/ecc/bls12-381/internal/fptower"
	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	"testing"
)

func TestG2SqrtRatio(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	properties := gopter.NewProperties(parameters)
	gen := genCoordPairG2(t)

	properties.Property("G2SqrtRatio must square back to the right value", prop.ForAll(
		func(uv []fptower.E2) bool {
			//uv = []fptower.E2{{fp.Element{1}, fp.Element{1}}, {fp.Element{1}, fp.Element{1}}}
			u := &uv[0]
			v := &uv[1]

			var seen fptower.E2
			qr := g2SqrtRatio(&seen, u, v) == 0

			if u.IsZero() {
				return seen.IsZero()
			}

			seen.
				Square(&seen).
				Mul(&seen, v)

			var ref fptower.E2
			if qr {
				ref = *u
			} else {
				g2MulByZ(&ref, u)
			}

			return seen.Equal(&ref)

		}, gen))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func genCoordPairG2(t *testing.T) gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {

		genRandomPair := func() (fptower.E2, fptower.E2) {
			var a, b fptower.E2

			if _, err := a.SetRandom(); err != nil {
				t.Error(err)
			}

			if _, err := b.SetRandom(); err != nil {
				t.Error(err)
			}

			return a, b
		}
		a, b := genRandomPair()

		genResult := gopter.NewGenResult([]fptower.E2{a, b}, gopter.NoShrinker)
		return genResult
	}
}
