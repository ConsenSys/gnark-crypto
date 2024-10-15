// Copyright 2022 ConsenSys Software Inc.
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

package arm64

import (
	"github.com/consensys/bavard/arm64"
)

func (f *FFArm64) generateAdd() {
	f.Comment("add(res, x, y *Element)")
	//stackSize := f.StackSize(f.NbWords*2, 0, 0)
	registers := f.FnHeader("add", 0, 24)
	defer f.AssertCleanStack(0, 0)

	// registers
	t := registers.PopN(f.NbWords)
	z := registers.PopN(f.NbWords)
	x := registers.PopN(f.NbWords)
	xPtr := registers.Pop()
	yPtr := registers.Pop()
	zPtr := registers.Pop()

	f.LDP("x+8(FP)", xPtr, yPtr)
	f.Comment("load operands and add mod 2^r")

	for i := 0; i < f.NbWords-1; i += 2 {
		f.LDP(xPtr.At(i), x[i], x[i+1])
		f.LDP(yPtr.At(i), z[i], z[i+1])
	}

	f.ADDS(x[0], z[0], z[0])
	for i := 1; i < f.NbWords; i++ {
		f.ADCS(x[i], z[i], z[i])
	}

	f.reduce(z, t)

	f.Comment("store")

	f.MOVD("res+0(FP)", zPtr)

	for i := 0; i < f.NbWords-1; i += 2 {
		f.STP(z[i], z[i+1], zPtr.At(i))
	}

	f.RET()

}

func (f *FFArm64) generateDouble() {
	f.Comment("double(res, x *Element)")
	registers := f.FnHeader("double", 0, 16)
	defer f.AssertCleanStack(0, 0)

	// registers
	z := registers.PopN(f.NbWords)
	xPtr := registers.Pop()
	zPtr := registers.Pop()
	//ops := registers.PopN(2)

	f.LDP("res+0(FP)", zPtr, xPtr)
	f.Comment("load operands and add mod 2^r")

	op0 := f.ADDS
	for i := 0; i < f.NbWords-1; i += 2 {
		f.LDP(xPtr.At(i), z[i], z[i+1])

		op0(z[i], z[i], z[i])
		op0 = f.ADCS

		f.ADCS(z[i+1], z[i+1], z[i+1])
	}

	if f.NbWords%2 == 1 {
		i := f.NbWords - 1
		f.MOVD(xPtr.At(i), z[i])
		op0(z[i], z[i], z[i])
	}
	registers.Push(xPtr)

	t := registers.PopN(f.NbWords)
	f.reduce(z, t)
	registers.Push(t...)

	f.Comment("store")
	f.storeVector(z, zPtr)

	f.RET()

}

// generateSub NO LONGER uses one more register than generateAdd, but that's okay since we have 29 registers available.
func (f *FFArm64) generateSub() {
	f.Comment("sub(res, x, y *Element)")

	registers := f.FnHeader("sub", 0, 24)
	defer f.AssertCleanStack(0, 0)

	// registers
	z := registers.PopN(f.NbWords)
	xPtr := registers.Pop()
	yPtr := registers.Pop()
	ops := registers.PopN(2)

	f.LDP("x+8(FP)", xPtr, yPtr)
	f.Comment("load operands and subtract mod 2^r")

	op0 := f.SUBS
	for i := 0; i < f.NbWords-1; i += 2 {
		f.LDP(xPtr.At(i), z[i], ops[0])
		f.LDP(yPtr.At(i), z[i+1], ops[1])

		op0(z[i+1], z[i], z[i])
		op0 = f.SBCS

		f.SBCS(ops[1], ops[0], z[i+1])
	}

	registers.Push(xPtr, yPtr)
	registers.Push(ops...)

	f.Comment("load modulus and select")

	t := registers.PopN(f.NbWords)
	zero := registers.Pop()
	f.MOVD(0, zero)

	for i := 0; i < f.NbWords-1; i += 2 {
		f.LDP(f.qAt(i), t[i], t[i+1])

		f.CSEL("CS", zero, t[i], t[i])
		f.CSEL("CS", zero, t[i+1], t[i+1])
	}

	if f.NbWords%2 == 1 {
		i := f.NbWords - 1
		f.MOVD(f.qAt(i), t[i])

		f.CSEL("CS", zero, t[i], t[i])
	}

	registers.Push(zero)

	f.Comment("augment (or not)")

	op0 = f.ADDS
	for i := 0; i < f.NbWords; i++ {
		op0(z[i], t[i], z[i])
		op0 = f.ADCS
	}

	registers.Push(t...)

	f.Comment("store")
	zPtr := registers.Pop()
	f.MOVD("res+0(FP)", zPtr)
	f.storeVector(z, zPtr)

	f.RET()

}

func (f *FFArm64) generateNeg() {
	f.Comment("neg(res, x *Element)")
	registers := f.FnHeader("neg", 0, 16)
	defer f.AssertCleanStack(0, 0)

	// registers
	z := registers.PopN(f.NbWords)
	xPtr := registers.Pop()
	zPtr := registers.Pop()
	ops := registers.PopN(2)
	xNotZero := registers.Pop()

	f.LDP("res+0(FP)", zPtr, xPtr)
	f.Comment("load operands and subtract")

	f.MOVD(0, xNotZero)
	op0 := f.SUBS
	for i := 0; i < f.NbWords-1; i += 2 {
		f.LDP(xPtr.At(i), z[i], z[i+1])
		f.LDP(f.qAt(i), ops[0], ops[1])

		f.ORR(z[i], xNotZero, xNotZero, "has x been 0 so far?")
		f.ORR(z[i+1], xNotZero, xNotZero)

		op0(z[i], ops[0], z[i])
		op0 = f.SBCS

		f.SBCS(z[i+1], ops[1], z[i+1])
	}

	registers.Push(xPtr)
	registers.Push(ops...)

	f.TST(-1, xNotZero)
	for i := 0; i < f.NbWords; i++ {
		f.CSEL("EQ", xNotZero, z[i], z[i])
	}

	f.Comment("store")
	f.storeVector(z, zPtr)

	f.RET()

}

func (f *FFArm64) reduce(z, t []arm64.Register) {

	if len(z) != f.NbWords || len(t) != f.NbWords {
		panic("need 2*nbWords registers")
	}

	f.Comment("load modulus and subtract")

	for i := 0; i < f.NbWords; i++ {
		f.MOVD(f.qi(i), t[i])
	}
	f.SUBS(t[0], z[0], t[0])
	for i := 1; i < f.NbWords; i++ {
		f.SBCS(t[i], z[i], t[i])
	}

	f.Comment("reduce if necessary")

	for i := 0; i < f.NbWords; i++ {
		f.CSEL("CS", t[i], z[i], z[i])
	}
}

func (f *FFArm64) storeVector(vector interface{}, baseAddress arm64.Register) {
	// f.callTemplate("storeVector", toInterfaceSlice(baseAddress, vector)...)
}
