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

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	"math/rand"
	"testing"
)

func TestG1SqrtRatio(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	gen := GenFp()

	properties.Property("G1SqrtRatio must square back to the right value", prop.ForAll(
		func(u fp.Element, v fp.Element) bool {

			var seen fp.Element
			qr := g1SqrtRatio(&seen, &u, &v) == 0

			seen.
				Square(&seen).
				Mul(&seen, &v)

			var ref fp.Element
			if qr {
				ref = u
			} else {
				g1MulByZ(&ref, &u)
			}

			return seen.Equal(&ref)
		}, gen, gen))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

//TODO: Crude. Do something clever in Jacobian
func isOnEPrimeG1(p G1Affine) bool {

	var A, B fp.Element

	// TODO: Value already in Mont form, set string without mont conversion

	A.SetString(
		"3282900303753830146348712611671767197635358960998995512158339553859262788125660792597406590695608535632558761028177",
	)

	B.SetString(
		"1001795216334152487852409076051911202261191374951378131193782595415765734425112784164810052880664366634483522380256",
	)

	A.FromMont()
	B.FromMont()

	var LHS fp.Element
	LHS.
		Square(&p.Y).
		Sub(&LHS, &B)

	var RHS fp.Element
	RHS.
		Square(&p.X).
		Add(&RHS, &A).
		Mul(&RHS, &p.X)

	return LHS.Equal(&RHS)
}

func TestG1SSWU(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[G1] hash outputs must be in appropriate groups", prop.ForAll(
		func(a fp.Element) bool {

			g := sswuMapG1(&a)

			if !isOnEPrimeG1(g) {
				t.Log("SSWU output not on E' curve")
				return false
			}

			g1Isogeny(&g)

			if !g.IsOnCurve() {
				t.Log("Isogeny/SSWU output not on curve")
				return false
			}

			return true
		},
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func g1TestMatchCoord(t *testing.T, coordName string, msg string, expectedStr string, seen *fp.Element) {
	var expected fp.Element

	expected.SetString(expectedStr)

	if !expected.Equal(seen) {
		t.Errorf("mismatch on \"%s\", %s:\n\texpected %s\n\tsaw      %s", msg, coordName, expected.String(), seen)
	}
}

func g1TestMatch(t *testing.T, c hashTestCase, seen *G1Affine) {
	g1TestMatchCoord(t, "x", c.msg, c.x, &seen.X)
	g1TestMatchCoord(t, "y", c.msg, c.y, &seen.Y)
}

func TestEncodeToG1(t *testing.T) {
	t.Parallel()
	for _, c := range g1EncodeToCurveSSWUVector.cases {
		seen, err := EncodeToG1([]byte(c.msg), g1EncodeToCurveSSWUVector.dst)
		if err != nil {
			t.Fatal(err)
		}
		g1TestMatch(t, c, &seen)
	}
}

func TestHashToG1(t *testing.T) {
	t.Parallel()
	for _, c := range g1HashToCurveSSWUVector.cases {
		seen, err := HashToG1([]byte(c.msg), g1HashToCurveSSWUVector.dst)
		if err != nil {
			t.Fatal(err)
		}
		g1TestMatch(t, c, &seen)
	}
	t.Log(len(g1HashToCurveSSWUVector.cases), "cases verified")
}

func BenchmarkEncodeToG1(b *testing.B) {
	const size = 54
	bytes := make([]byte, size)
	dst := g1EncodeToCurveSSWUVector.dst
	b.ResetTimer()

	for i := 0; i < b.N; i++ {

		bytes[rand.Int()%size] = byte(rand.Int())

		if _, err := EncodeToG1(bytes, dst); err != nil {
			b.Fail()
		}
	}
}

func BenchmarkHashToG1(b *testing.B) {
	const size = 54
	bytes := make([]byte, size)
	dst := g1HashToCurveSSWUVector.dst
	b.ResetTimer()

	for i := 0; i < b.N; i++ {

		bytes[rand.Int()%size] = byte(rand.Int())

		if _, err := HashToG1(bytes, dst); err != nil {
			b.Fail()
		}
	}
}

type hashTestVector struct {
	dst   []byte
	cases []hashTestCase
}

type hashTestCase struct {
	msg string
	x   string
	y   string
}

var g1HashToCurveSSWUVector hashTestVector
var g1EncodeToCurveSSWUVector hashTestVector
