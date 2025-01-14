// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package plookup

import (
	"math/big"
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/kzg"
)

func TestLookupVector(t *testing.T) {

	lookupVector := make(fr.Vector, 8)
	fvector := make(fr.Vector, 7)
	for i := 0; i < 8; i++ {
		lookupVector[i].SetUint64(uint64(2 * i))
	}
	for i := 0; i < 7; i++ {
		fvector[i].Set(&lookupVector[(4*i+1)%8])
	}

	kzgSrs, err := kzg.NewSRS(64, big.NewInt(13))
	if err != nil {
		t.Fatal(err)
	}

	// correct proof vector
	{
		proof, err := ProveLookupVector(kzgSrs.Pk, fvector, lookupVector)
		if err != nil {
			t.Fatal(err)
		}

		err = VerifyLookupVector(kzgSrs.Vk, proof)
		if err != nil {
			t.Fatal(err)
		}
	}

	// wrong proofs vector
	{
		fvector[0].SetRandom()

		proof, err := ProveLookupVector(kzgSrs.Pk, fvector, lookupVector)
		if err != nil {
			t.Fatal(err)
		}

		err = VerifyLookupVector(kzgSrs.Vk, proof)
		if err == nil {
			t.Fatal(err)
		}
	}

}

func TestLookupTable(t *testing.T) {

	kzgSrs, err := kzg.NewSRS(64, big.NewInt(13))
	if err != nil {
		t.Fatal(err)
	}

	lookupTable := make([]fr.Vector, 3)
	fTable := make([]fr.Vector, 3)
	for i := 0; i < 3; i++ {
		lookupTable[i] = make(fr.Vector, 8)
		fTable[i] = make(fr.Vector, 7)
		for j := 0; j < 8; j++ {
			lookupTable[i][j].SetUint64(uint64(2*i + j))
		}
		for j := 0; j < 7; j++ {
			fTable[i][j].Set(&lookupTable[i][(4*j+1)%8])
		}
	}

	// correct proof
	{
		proof, err := ProveLookupTables(kzgSrs.Pk, fTable, lookupTable)
		if err != nil {
			t.Fatal(err)
		}

		err = VerifyLookupTables(kzgSrs.Vk, proof)
		if err != nil {
			t.Fatal(err)
		}
	}

	// wrong proof
	{
		fTable[0][0].SetRandom()
		proof, err := ProveLookupTables(kzgSrs.Pk, fTable, lookupTable)
		if err != nil {
			t.Fatal(err)
		}

		err = VerifyLookupTables(kzgSrs.Vk, proof)
		if err == nil {
			t.Fatal(err)
		}
	}

}

func BenchmarkPlookup(b *testing.B) {

	srsSize := 1 << 15
	polySize := 1 << 14

	kzgSrs, _ := kzg.NewSRS(uint64(srsSize), big.NewInt(13))
	a := make(fr.Vector, polySize)
	c := make(fr.Vector, polySize)

	for i := 0; i < 1<<14; i++ {
		a[i].SetUint64(uint64(i))
		c[i].SetUint64(uint64((8 * i) % polySize))
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ProveLookupVector(kzgSrs.Pk, a, c)
	}
}
