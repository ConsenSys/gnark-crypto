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

package amd64

import "github.com/consensys/bavard/amd64"

// addVec res = a + b
// func addVec(res, a, b *{{.ElementName}}, n uint64)
func (f *FFAmd64) generateAddVec() {
	f.Comment("addVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] + b[0...n]")

	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader("addVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)

	a := f.PopN(&registers)
	t := f.PopN(&registers)

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	// load arguments
	f.MOVQ("res+0(FP)", addrRes)
	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("b+16(FP)", addrB)
	f.MOVQ("n+24(FP)", len)

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	// a = a + b
	f.LabelRegisters("a", a...)
	f.Mov(addrA, a)
	f.Add(addrB, a)

	// reduce a
	f.ReduceElement(a, t)

	// save a into res
	f.Mov(a, addrRes)

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrB)
	f.ADDQ("$32", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, a...)
	f.Push(&registers, t...)
	f.Push(&registers, addrA, addrB, addrRes, len)

}

// subVec res = a - b
// func subVec(res, a, b *{{.ElementName}}, n uint64)
func (f *FFAmd64) generateSubVec() {
	f.Comment("subVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] - b[0...n]")

	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+5, 0, 0)
	registers := f.FnHeader("subVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)
	zero := f.Pop(&registers)

	a := f.PopN(&registers)
	q := f.PopN(&registers)

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	// load arguments
	f.MOVQ("res+0(FP)", addrRes)
	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("b+16(FP)", addrB)
	f.MOVQ("n+24(FP)", len)

	f.XORQ(zero, zero)

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	// a = a - b
	f.LabelRegisters("a", a...)
	f.Mov(addrA, a)
	f.Sub(addrB, a)

	// reduce a
	f.Comment("reduce (a-b) mod q")
	f.LabelRegisters("q", q...)
	f.Mov(f.Q, q)
	for i := 0; i < f.NbWords; i++ {
		f.CMOVQCC(zero, q[i])
	}
	// add registers (q or 0) to a, and set to result
	f.Comment("add registers (q or 0) to a, and set to result")
	f.Add(q, a)

	// save a into res
	f.Mov(a, addrRes)

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrB)
	f.ADDQ("$32", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, a...)
	f.Push(&registers, q...)
	f.Push(&registers, addrA, addrB, addrRes, len, zero)

}

// scalarMulVec res = a * b
// func scalarMulVec(res, a, b *{{.ElementName}}, n uint64)
func (f *FFAmd64) generateScalarMulVec() {
	f.Comment("scalarMulVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] * b")

	const argSize = 4 * 8
	const minStackSize = 7 * 8 // 2 slices (3 words each) + pointer to the scalar
	stackSize := f.StackSize(f.NbWords*2+3, 2, minStackSize)
	reserved := []amd64.Register{amd64.DX, amd64.AX}
	registers := f.FnHeader("scalarMulVec", stackSize, argSize, reserved...)
	defer f.AssertCleanStack(stackSize, minStackSize)

	// labels & registers we need
	noAdx := f.NewLabel("noAdx")
	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	t := registers.PopN(f.NbWords)
	scalar := registers.PopN(f.NbWords)

	addrB := registers.Pop()
	addrA := registers.Pop()
	addrRes := addrB
	len := registers.Pop()

	// check ADX instruction support
	f.CMPB("·supportAdx(SB)", 1)
	f.JNE(noAdx)

	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("b+16(FP)", addrB)
	f.MOVQ("n+24(FP)", len)

	// we store b, the scalar, fully in registers
	f.LabelRegisters("scalar", scalar...)
	f.Mov(addrB, scalar)

	xat := func(i int) string {
		return string(scalar[i])
	}

	f.MOVQ("res+0(FP)", addrRes)

	f.LABEL(loop)
	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	yat := func(i int) string {
		return addrA.At(i)
	}

	f.Comment("TODO @gbotrel this is generated from the same macro as the unit mul, we should refactor this in a single asm function")

	f.MulADX(&registers, xat, yat, t)

	// registers.Push(addrA)

	// reduce; we need at least 4 extra registers
	registers.Push(amd64.AX, amd64.DX)
	f.Comment("reduce t mod q")
	f.Reduce(&registers, t)
	f.Mov(t, addrRes)

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)
	f.RET()

	// no ADX support
	f.LABEL(noAdx)

	f.MOVQ("n+24(FP)", amd64.DX)

	f.MOVQ("res+0(FP)", amd64.AX)
	f.MOVQ(amd64.AX, "(SP)")
	f.MOVQ(amd64.DX, "8(SP)")  // len
	f.MOVQ(amd64.DX, "16(SP)") // cap
	f.MOVQ("a+8(FP)", amd64.AX)
	f.MOVQ(amd64.AX, "24(SP)")
	f.MOVQ(amd64.DX, "32(SP)") // len
	f.MOVQ(amd64.DX, "40(SP)") // cap
	f.MOVQ("b+16(FP)", amd64.AX)
	f.MOVQ(amd64.AX, "48(SP)")
	f.WriteLn("CALL ·scalarMulVecGeneric(SB)")
	f.RET()

}
