
	// Copyright 2020 Consensys Software Inc.
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

package fr

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/consensys/gnark-crypto/ecc"
)

// Generator returns a generator for Z/2^(log(m))Z
// or an error if m is too big (required root of unity doesn't exist)
func Generator(m uint64) (Element, error) {
	x := ecc.NextPowerOfTwo(m)

	var rootOfUnity Element
    rootOfUnity.SetString("19103219067921713944291392827692070036145651957329286315305642004821462161904")
    const maxOrderRoot uint64 = 28

	// find generator for Z/2^(log(m))Z
	logx := uint64(bits.TrailingZeros64(x))
	if logx > maxOrderRoot {
		return Element{}, fmt.Errorf("m (%d) is too big: the required root of unity does not exist", m)
	}

	expo := uint64(1 << (maxOrderRoot - logx))
	var generator Element
	generator.Exp(rootOfUnity, big.NewInt(int64(expo))) // order x
	return generator, nil
}