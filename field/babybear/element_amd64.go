//go:build !purego

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

package babybear

//go:noescape
func MulBy3(x *Element)

//go:noescape
func MulBy5(x *Element)

//go:noescape
func MulBy13(x *Element)

//go:noescape
func mul(res, x, y *Element)

//go:noescape
func fromMont(res *Element)

//go:noescape
func reduce(res *Element)

// Butterfly sets
//
//	a = a + b (mod q)
//	b = a - b (mod q)
//
//go:noescape
func Butterfly(a, b *Element)

// Mul z = x * y (mod q)
//
// x and y must be less than q
func (z *Element) Mul(x, y *Element) *Element {

	// Algorithm 2 of "Faster Montgomery Multiplication and Multi-Scalar-Multiplication for SNARKS"
	// by Y. El Housni and G. Botrel https://doi.org/10.46586/tches.v2023.i3.504-521

	mul(z, x, y)
	return z
}

// Square z = x * x (mod q)
//
// x must be less than q
func (z *Element) Square(x *Element) *Element {
	// see Mul for doc.
	mul(z, x, x)
	return z
}
