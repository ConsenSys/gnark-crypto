// Copyright 2020-2024 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package ecdsa

import (
	"crypto/rand"
	"crypto/subtle"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

const (
	nbFuzzShort = 10
	nbFuzz      = 100
)

func TestSerialization(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS24-317] ECDSA serialization: SetBytes(Bytes()) should stay the same", prop.ForAll(
		func() bool {
			privKey, _ := GenerateKey(rand.Reader)

			var end PrivateKey
			buf := privKey.Bytes()
			n, err := end.SetBytes(buf[:])
			if err != nil {
				return false
			}
			if n != sizePrivateKey {
				return false
			}

			return end.PublicKey.Equal(&privKey.PublicKey) && subtle.ConstantTimeCompare(end.scalar[:], privKey.scalar[:]) == 1

		},
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}
