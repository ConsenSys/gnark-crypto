// Copyright 2020 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fri

import (
	"crypto/sha256"
	"fmt"
	"math/big"
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/gen"
	"github.com/leanovate/gopter/prop"
)

// logFiber returns u, v such that {g^u, g^v} = f⁻¹((g²)^{_p})
func logFiber(_p, _n int) (_u, _v big.Int) {
	if _p%2 == 0 {
		_u.SetInt64(int64(_p / 2))
		_v.SetInt64(int64(_p/2 + _n/2))
	} else {
		l := (_n - 1 - _p) / 2
		_u.SetInt64(int64(_n - 1 - l))
		_v.SetInt64(int64(_n - 1 - l - _n/2))
	}
	return
}

func randomPolynomial(size uint64, seed int32) []fr.Element {
	p := make([]fr.Element, size)
	p[0].SetUint64(uint64(seed))
	for i := 1; i < len(p); i++ {
		p[i].Square(&p[i-1])
	}
	return p
}

// convertOrderCanonical convert the index i, an entry in a
// sorted polynomial, to the corresponding entry in canonical
// representation. n is the size of the polynomial.
func convertSortedCanonical(i, n int) int {
	if i%2 == 0 {
		return i / 2
	} else {
		l := (n - 1 - i) / 2
		return n - 1 - l
	}
}

func TestFRI(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10

	properties := gopter.NewProperties(parameters)

	size := 4096

	properties.Property("verifying wrong opening should fail", prop.ForAll(

		func(m int32) bool {

			_s := RADIX_2_FRI.New(uint64(size), sha256.New())
			s := _s.(radixTwoFri)

			p := randomPolynomial(uint64(size), m)

			pos := int64(m % 4096)
			pp, _ := s.BuildProofOfProximity(p)

			openingProof, err := s.Open(p, uint64(pos))
			if err != nil {
				t.Fatal(err)
			}

			// check the Merkle path
			tamperedPosition := pos + 1
			err = s.VerifyOpening(uint64(tamperedPosition), openingProof, pp)

			return err != nil

		},
		gen.Int32Range(1, int32(rho*size)),
	))

	properties.Property("verifying correct opening should succeed", prop.ForAll(

		func(m int32) bool {

			_s := RADIX_2_FRI.New(uint64(size), sha256.New())
			s := _s.(radixTwoFri)

			p := randomPolynomial(uint64(size), m)

			pos := uint64(m % int32(size))
			pp, _ := s.BuildProofOfProximity(p)

			openingProof, err := s.Open(p, uint64(pos))
			if err != nil {
				t.Fatal(err)
			}

			// check the Merkle path
			err = s.VerifyOpening(uint64(pos), openingProof, pp)

			return err == nil

		},
		gen.Int32Range(0, int32(rho*size)),
	))

	properties.Property("The claimed value of a polynomial should match P(x)", prop.ForAll(
		func(m int32) bool {

			_s := RADIX_2_FRI.New(uint64(size), sha256.New())
			s := _s.(radixTwoFri)

			p := randomPolynomial(uint64(size), m)

			// check the opening value
			var g fr.Element
			pos := int64(m % 4096)
			g.Set(&s.domain.Generator)
			g.Exp(g, big.NewInt(pos))

			var val fr.Element
			for i := len(p) - 1; i >= 0; i-- {
				val.Mul(&val, &g)
				val.Add(&p[i], &val)
			}

			openingProof, err := s.Open(p, uint64(pos))
			if err != nil {
				t.Fatal(err)
			}

			return openingProof.ClaimedValue.Equal(&val)

		},
		gen.Int32Range(0, int32(rho*size)),
	))

	properties.Property("Derive queries position: points should belong the correct fiber", prop.ForAll(

		func(m int32) bool {

			_s := RADIX_2_FRI.New(uint64(size), sha256.New())
			s := _s.(radixTwoFri)

			var g fr.Element

			_m := int(m) % size
			pos := s.deriveQueriesPositions(_m, int(s.domain.Cardinality))
			g.Set(&s.domain.Generator)
			n := int(s.domain.Cardinality)

			for i := 0; i < len(pos)-1; i++ {

				u, v := logFiber(pos[i], n)

				var g1, g2, g3 fr.Element
				g1.Exp(g, &u).Square(&g1)
				g2.Exp(g, &v).Square(&g2)
				nextPos := convertSortedCanonical(pos[i+1], n/2)
				g3.Square(&g).Exp(g3, big.NewInt(int64(nextPos)))

				if !g1.Equal(&g2) || !g1.Equal(&g3) {
					return false
				}
				g.Square(&g)
				n = n >> 1
			}
			return true
		},
		gen.Int32Range(0, int32(rho*size)),
	))

	properties.Property("verifying a correctly formed proof should succeed", prop.ForAll(

		func(s int32) bool {

			p := randomPolynomial(uint64(size), s)

			iop := RADIX_2_FRI.New(uint64(size), sha256.New())
			proof, err := iop.BuildProofOfProximity(p)
			if err != nil {
				t.Fatal(err)
			}

			err = iop.VerifyProofOfProximity(proof)
			return err == nil
		},
		gen.Int32Range(0, int32(rho*size)),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))

}

// Benchmarks

func BenchmarkProximityVerification(b *testing.B) {

	baseSize := 16

	for i := 0; i < 10; i++ {

		size := baseSize << i
		p := make([]fr.Element, size)
		for k := 0; k < size; k++ {
			p[k].SetRandom()
		}

		iop := RADIX_2_FRI.New(uint64(size), sha256.New())
		proof, _ := iop.BuildProofOfProximity(p)

		b.Run(fmt.Sprintf("Polynomial size %d", size), func(b *testing.B) {
			b.ResetTimer()
			for l := 0; l < b.N; l++ {
				iop.VerifyProofOfProximity(proof)
			}
		})

	}
}
