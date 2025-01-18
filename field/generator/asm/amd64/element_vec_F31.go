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

import (
	"fmt"

	"github.com/consensys/bavard/amd64"
)

func (f *FFAmd64) generateAddVecF31() {
	f.Comment("addVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] + b[0...n]")
	f.Comment("n is the number of blocks of 16 elements to process")

	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader("addVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	b := amd64.Register("Z1")
	t := amd64.Register("Z2")
	q := amd64.Register("Z3")

	// load q in Z3
	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTD("AX", q)

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
	f.VMOVDQU32(addrA.At(0), a)
	f.VMOVDQU32(addrB.At(0), b)
	f.VPADDD(a, b, a, "a = a + b")
	// t = a - q
	f.VPSUBD(q, a, t, "t = a - q")
	// b = min(t, a)
	f.VPMINUD(a, t, b, "b = min(t, a)")

	// move b to res
	f.VMOVDQU32(b, addrRes.At(0), "res = b")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$64", addrA)
	f.ADDQ("$64", addrB)
	f.ADDQ("$64", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, addrA, addrB, addrRes, len)

}

func (f *FFAmd64) generateSubVecF31() {
	f.Comment("subVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] - b[0...n]")
	f.Comment("n is the number of blocks of 16 elements to process")

	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader("subVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	b := amd64.Register("Z1")
	t := amd64.Register("Z2")
	q := amd64.Register("Z3")

	// load q in Z3
	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTD("AX", q)

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

	// a = a - b
	f.VMOVDQU32(addrA.At(0), a)
	f.VMOVDQU32(addrB.At(0), b)

	f.VPSUBD(b, a, a, "a = a - b")

	// t = a + q
	f.VPADDD(q, a, t, "t = a + q")

	// b = min(t, a)
	f.VPMINUD(a, t, b, "b = min(t, a)")

	// move b to res
	f.VMOVDQU32(b, addrRes.At(0), "res = b")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$64", addrA)
	f.ADDQ("$64", addrB)
	f.ADDQ("$64", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, addrA, addrB, addrRes, len)

}

// sumVec res = sum(a[0...n])
func (f *FFAmd64) generateSumVecF31() {
	f.Comment("sumVec(res *uint64, a *[]uint32, n uint64) res = sum(a[0...n])")
	f.Comment("n is the number of blocks of 16 elements to process")
	const argSize = 3 * 8
	stackSize := f.StackSize(f.NbWords*3+2, 0, 0)
	registers := f.FnHeader("sumVec", stackSize, argSize, amd64.DX, amd64.AX)
	defer f.AssertCleanStack(stackSize, 0)

	f.WriteLn(`
	// We load 8 31bits values at a time and accumulate them into an accumulator of
	// 8 quadwords (64bits). The caller then needs to reduce the result mod q.
	// We can safely accumulate ~2**33 31bits values into a single accumulator.
	// That gives us a maximum of 2**33 * 8 = 2**36 31bits values to sum safely.
	`)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrT := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a1 := amd64.Register("Z0")
	a2 := amd64.Register("Z1")
	acc1 := amd64.Register("Z2")
	acc2 := amd64.Register("Z3")

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	// load arguments
	f.MOVQ("t+0(FP)", addrT)
	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("n+16(FP)", len)

	// zeroize the accumulators
	f.VXORPS(acc1, acc1, acc1, "acc1 = 0")
	f.VMOVDQA64(acc1, acc2, "acc2 = 0")

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	// 1 cache line is typically 64 bytes, so we maintain 2 accumulators
	f.VPMOVZXDQ(addrA.At(0), a1, "load 8 31bits values in a1")
	f.VPMOVZXDQ(addrA.At(4), a2, "load 8 31bits values in a2")

	f.VPADDQ(a1, acc1, acc1, "acc1 += a1")
	f.VPADDQ(a2, acc2, acc2, "acc2 += a2")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$64", addrA)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	// store t into res
	f.VPADDQ(acc1, acc2, acc1, "acc1 += acc2")
	f.VMOVDQU64(acc1, addrT.At(0), "res = acc1")

	f.RET()

	f.Push(&registers, addrA, addrT, len)
}

// mulVec res = a * b
func (f *FFAmd64) generateMulVecF31() {
	f.Comment("mulVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] * b[0...n]")
	f.Comment("n is the number of blocks of 8 elements to process")
	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader("mulVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	b := amd64.Register("Z1")
	P := amd64.Register("Z2")
	q := amd64.Register("Z3")
	qInvNeg := amd64.Register("Z4")
	PL := amd64.Register("Z5")
	LSW := amd64.Register("Z6")

	// load q in Z3
	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTQ("AX", q)
	f.WriteLn("MOVD $const_qInvNeg, AX")
	f.VPBROADCASTQ("AX", qInvNeg)

	f.Comment("Create mask for low dword in each qword")
	f.VPCMPEQB("Y0", "Y0", "Y0")
	f.VPMOVZXDQ("Y0", LSW)

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

	// a = a * b
	f.VPMOVZXDQ(addrA.At(0), a)
	f.VPMOVZXDQ(addrB.At(0), b)
	f.VPMULUDQ(a, b, P, "P = a * b")
	f.VPANDQ(LSW, P, PL, "m = uint32(P)")
	f.VPMULUDQ(PL, qInvNeg, PL, "m = m * qInvNeg")
	f.VPANDQ(LSW, PL, PL, "m = uint32(m)")
	f.VPMULUDQ(PL, q, PL, "m = m * q")
	f.VPADDQ(P, PL, P, "P = P + m")
	f.VPSRLQ("$32", P, P, "P = P >> 32")

	f.VPSUBD(q, P, PL, "PL = P - q")
	f.VPMINUD(P, PL, P, "P = min(P, PL)")

	// move P to res
	f.VPMOVQD(P, addrRes.At(0), "res = P")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrB)
	f.ADDQ("$32", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, addrA, addrB, addrRes, len)

}

// scalarMulVec res = a * b
func (f *FFAmd64) generateScalarMulVecF31() {
	f.Comment("scalarMulVec(res, a, b *Element, n uint64) res[0...n] = a[0...n] * b")
	f.Comment("n is the number of blocks of 8 elements to process")
	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader("scalarMulVec", stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrRes := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	b := amd64.Register("Z1")
	P := amd64.Register("Z2")
	q := amd64.Register("Z3")
	qInvNeg := amd64.Register("Z4")
	PL := amd64.Register("Z5")
	LSW := amd64.Register("Z6")

	// load q in Z3
	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTQ("AX", q)
	f.WriteLn("MOVD $const_qInvNeg, AX")
	f.VPBROADCASTQ("AX", qInvNeg)

	f.Comment("Create mask for low dword in each qword")
	f.VPCMPEQB("Y0", "Y0", "Y0")
	f.VPMOVZXDQ("Y0", LSW)

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	// load arguments
	f.MOVQ("res+0(FP)", addrRes)
	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("b+16(FP)", addrB)
	f.MOVQ("n+24(FP)", len)

	f.VPBROADCASTD(addrB.At(0), b)

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	// a = a * b
	f.VPMOVZXDQ(addrA.At(0), a)

	f.VPMULUDQ(a, b, P, "P = a * b")
	f.VPANDQ(LSW, P, PL, "m = uint32(P)")
	f.VPMULUDQ(PL, qInvNeg, PL, "m = m * qInvNeg")
	f.VPANDQ(LSW, PL, PL, "m = uint32(m)")
	f.VPMULUDQ(PL, q, PL, "m = m * q")
	f.VPADDQ(P, PL, P, "P = P + m")
	f.VPSRLQ("$32", P, P, "P = P >> 32")

	f.VPSUBD(q, P, PL, "PL = P - q")
	f.VPMINUD(P, PL, P, "P = min(P, PL)")

	// move P to res
	f.VPMOVQD(P, addrRes.At(0), "res = P")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrRes)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, addrA, addrB, addrRes, len)
}

// innerProdVec res = sum(a * b)
func (f *FFAmd64) generateInnerProdVecF31() {
	f.Comment("innerProdVec(t *uint64, a,b *[]uint32, n uint64) res = sum(a[0...n] * b[0...n])")
	f.Comment("n is the number of blocks of 8 elements to process")
	const argSize = 4 * 8
	stackSize := f.StackSize(f.NbWords*4+2, 0, 0)
	registers := f.FnHeader("innerProdVec", stackSize, argSize, amd64.DX, amd64.AX)
	defer f.AssertCleanStack(stackSize, 0)

	f.WriteLn(`
	// Similar to mulVec; we do most of the montgomery multiplication but don't do
	// the final reduction. We accumulate the result like in sumVec and let the caller
	// reduce mod q.
	`)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrB := f.Pop(&registers)
	addrT := f.Pop(&registers)
	len := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	b := amd64.Register("Z1")
	acc := amd64.Register("Z2")
	q := amd64.Register("Z3")
	qInvNeg := amd64.Register("Z4")
	PL := amd64.Register("Z5")
	LSW := amd64.Register("Z6")
	P := amd64.Register("Z7")

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTQ("AX", q)
	f.WriteLn("MOVD $const_qInvNeg, AX")
	f.VPBROADCASTQ("AX", qInvNeg)

	f.Comment("Create mask for low dword in each qword")
	f.VPCMPEQB("Y0", "Y0", "Y0")
	f.VPMOVZXDQ("Y0", LSW)

	// zeroize the accumulators
	f.VXORPS(acc, acc, acc, "acc = 0")

	// load arguments
	f.MOVQ("t+0(FP)", addrT)
	f.MOVQ("a+8(FP)", addrA)
	f.MOVQ("b+16(FP)", addrB)
	f.MOVQ("n+24(FP)", len)

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	f.VPMOVZXDQ(addrA.At(0), a)
	f.VPMOVZXDQ(addrB.At(0), b)

	f.VPMULUDQ(a, b, P, "P = a * b")
	f.VPANDQ(LSW, P, PL, "m = uint32(P)")
	f.VPMULUDQ(PL, qInvNeg, PL, "m = m * qInvNeg")
	f.VPANDQ(LSW, PL, PL, "m = uint32(m)")
	f.VPMULUDQ(PL, q, PL, "m = m * q")
	f.VPADDQ(P, PL, P, "P = P + m")
	f.VPSRLQ("$32", P, P, "P = P >> 32")

	// we can accumulate ~2**32 32bits values into a single accumulator without overflow;
	// that gives us a maximum of 2**32 * 8 = 2**35 32bits values to sum safely.
	f.Comment("accumulate P into acc, P is in [0, 2q] on 32bits max")
	f.VPADDQ(P, acc, acc, "acc += P")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrB)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	// store t into res
	f.VMOVDQU64(acc, addrT.At(0), "res = acc")

	f.RET()

	f.Push(&registers, addrA, addrT, len)
}

// func kerDIFNP_{{.sizeKernel}}generic(a []{{ .FF }}.Element, twiddles [][]{{ .FF }}.Element, stage int) {
// 	n := 1 << {{.sizeKernelLog2}}
// 	m := n >> 1
// 	split := 1

// 	for step := 0; step < {{.sizeKernelLog2}}; step++ {
// 		bound := split * n
// 		if m > 8 {
// 			for offset :=0; offset  < bound; offset += n {
// 				innerDIFWithTwiddles(a[offset:offset + n], twiddles[stage + step], 0, m, m)
// 			}
// 		} else {
// 			for offset :=0; offset  < bound; offset += n {
// 				// innerDIFWithTwiddles(a[offset:offset + n], twiddles[stage + step], 0, m, m)
// 				aa := a[offset:offset + n]
// 				t := twiddles[stage + step]
// 				for i := 0; i < m; i++ {
// 					{{ .FF }}.Butterfly(&aa[i], &aa[i+m])
// 					aa[i+m].Mul(&aa[i+m], &t[i])
// 				}
// 			}
// 		}
// 		n = n >> 1
// 		m = n >> 1
// 		split = split << 1
// 	}
// }

func (f *FFAmd64) generateFFTKernelF31(klog2 int) {
	ksize := 1 << klog2
	f.Comment(fmt.Sprintf("kerDIFNP_%d_avx512(a, twiddles [][]{{ .FF }}.Element, stage int)", ksize))
	// f.Comment("kerDIFNP_256_avx512(a, twiddles *Element, m uint64)")
	// f.Comment("n is the number of blocks of 8 elements to process")
	const argSize = 3 * 8
	stackSize := f.StackSize(f.NbWords*2+4, 0, 0)
	registers := f.FnHeader(fmt.Sprintf("kerDIFNP_%d_avx512", ksize), stackSize, argSize)
	defer f.AssertCleanStack(stackSize, 0)

	// registers & labels we need
	addrA := f.Pop(&registers)
	addrAPlusM := f.Pop(&registers)
	addrTwiddles := f.Pop(&registers)
	len := f.Pop(&registers)
	m := f.Pop(&registers)

	// AVX512 registers
	a := amd64.Register("Z0")
	am := amd64.Register("Z1")
	twiddles := amd64.Register("Z2")
	P := amd64.Register("Z3")
	q := amd64.Register("Z4")
	qInvNeg := amd64.Register("Z5")
	PL := amd64.Register("Z6")
	LSW := amd64.Register("Z7")
	b0 := amd64.Register("Z8")
	b1 := amd64.Register("Z9")

	// load q in Z3
	f.WriteLn("MOVD $const_q, AX")
	f.VPBROADCASTQ("AX", q)
	f.WriteLn("MOVD $const_qInvNeg, AX")
	f.VPBROADCASTQ("AX", qInvNeg)

	f.Comment("Create mask for low dword in each qword")
	f.VPCMPEQB("Y0", "Y0", "Y0")
	f.VPMOVZXDQ("Y0", LSW)

	loop := f.NewLabel("loop")
	done := f.NewLabel("done")

	// load arguments
	f.MOVQ("a+0(FP)", addrA)
	f.MOVQ("twiddles+8(FP)", addrTwiddles)
	f.MOVQ("m+16(FP)", m)

	f.MOVQ(m, len)

	// divide len by 8 (nb of iterations of 8 elements)
	f.SHRQ("$3", len)

	// multiply m by 4 (number of bytes in one element)
	f.SHLQ("$2", m)

	f.MOVQ(addrA, addrAPlusM)
	f.ADDQ(m, addrAPlusM)

	f.LABEL(loop)

	f.TESTQ(len, len)
	f.JEQ(done, "n == 0, we are done")

	// a[i] = a[i] + a[i+m]
	// a[m] = (a[i] - a[i+m])
	// Butterfly(&a[0], &a[m])
	f.VPMOVZXDQ(addrA.At(0), a)
	f.VPMOVZXDQ(addrAPlusM.At(0), am)

	// b0 = a + am
	f.VPADDD(a, am, b0, "b0 = a + am")
	// b1 = a - am
	f.VPSUBD(am, a, b1, "b1 = a - am")

	// reduce b0 mod q
	f.VPSUBD(q, b0, PL, "PL = b0 - q")
	f.VPMINUD(b0, PL, b0, "b0 = min(b0, PL)")

	// store it
	f.VPMOVQD(b0, addrA.At(0), "a = b0")

	// reduce b1 mod q
	f.VPADDD(q, b1, b1, "PL = b1 + q")
	// f.VPMINUD(b1, PL, b1, "b1 = min(b1, PL)")

	// now we just mul b1 by twiddles
	// a[i+m].Mul(&a[i+m], &twiddles[i])
	f.VPMOVZXDQ(addrTwiddles.At(0), twiddles)

	f.VPMULUDQ(b1, twiddles, P, "P = b1 * twiddles")
	f.VPANDQ(LSW, P, PL, "m = uint32(P)")
	f.VPMULUDQ(PL, qInvNeg, PL, "m = m * qInvNeg")
	f.VPANDQ(LSW, PL, PL, "m = uint32(m)")
	f.VPMULUDQ(PL, q, PL, "m = m * q")
	f.VPADDQ(P, PL, P, "P = P + m")
	f.VPSRLQ("$32", P, P, "P = P >> 32")

	f.VPSUBD(q, P, PL, "PL = P - q")
	f.VPMINUD(P, PL, P, "P = min(P, PL)")

	// move P to res
	f.VPMOVQD(P, addrAPlusM.At(0), "res = P")

	f.Comment("increment pointers to visit next element")
	f.ADDQ("$32", addrA)
	f.ADDQ("$32", addrTwiddles)
	f.ADDQ("$32", addrAPlusM)
	f.DECQ(len, "decrement n")
	f.JMP(loop)

	f.LABEL(done)

	f.RET()

	f.Push(&registers, addrA, addrTwiddles, addrAPlusM, len)

}
