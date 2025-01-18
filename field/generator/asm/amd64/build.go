// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Package amd64 contains syntactic sugar to generate amd64 assembly code
package amd64

import (
	"fmt"
	"io"
	"strings"

	"github.com/consensys/bavard/amd64"
)

const SmallModulus = 6

func NewFFAmd64(w io.Writer, nbWords int) *FFAmd64 {
	F := &FFAmd64{
		amd64.NewAmd64(w),
		0,
		0,
		nbWords,
		nbWords - 1,
		make([]int, nbWords),
		make([]int, nbWords-1),
		make(map[string]defineFn),
	}

	// indexes (template helpers)
	for i := 0; i < F.NbWords; i++ {
		F.NbWordsIndexesFull[i] = i
		if i > 0 {
			F.NbWordsIndexesNoZero[i-1] = i
		}
	}

	return F
}

type FFAmd64 struct {
	// *config.FieldConfig
	*amd64.Amd64
	nbElementsOnStack    int
	maxOnStack           int
	NbWords              int
	NbWordsLastIndex     int
	NbWordsIndexesFull   []int
	NbWordsIndexesNoZero []int
	mDefines             map[string]defineFn
}

type defineFn func(args ...amd64.Register)

func (f *FFAmd64) StackSize(maxNbRegistersNeeded, nbRegistersReserved, minStackSize int) int {
	got := amd64.NbRegisters - nbRegistersReserved
	r := got - maxNbRegistersNeeded
	if r >= 0 {
		return minStackSize
	}
	r *= -8
	return max(r, minStackSize)
}

func (f *FFAmd64) DefineFn(name string) (fn defineFn, err error) {
	fn, ok := f.mDefines[name]
	if !ok {
		return nil, fmt.Errorf("function %s not defined", name)
	}
	return fn, nil
}

func (f *FFAmd64) Define(name string, nbInputs int, fn defineFn) defineFn {

	inputs := make([]string, nbInputs)
	for i := 0; i < nbInputs; i++ {
		inputs[i] = fmt.Sprintf("in%d", i)
	}
	name = strings.ToUpper(name)

	for _, ok := f.mDefines[name]; ok; {
		// name already exist, for code generation purpose we add a suffix
		// should happen only with e2 deprecated functions
		fmt.Println("WARNING: function name already defined, adding suffix")
		i := 0
		for {
			newName := fmt.Sprintf("%s_%d", name, i)
			if _, ok := f.mDefines[newName]; !ok {
				name = newName
				goto startDefine
			}
			i++
		}
	}
startDefine:

	f.StartDefine()
	f.WriteLn("#define " + name + "(" + strings.Join(inputs, ", ") + ")")
	inputsRegisters := make([]amd64.Register, nbInputs)
	for i := 0; i < nbInputs; i++ {
		inputsRegisters[i] = amd64.Register(inputs[i])
	}
	fn(inputsRegisters...)
	f.EndDefine()
	f.WriteLn("")

	toReturn := func(args ...amd64.Register) {
		if len(args) != nbInputs {
			panic("invalid number of arguments")
		}
		inputsStr := make([]string, len(args))
		for i := 0; i < len(args); i++ {
			inputsStr[i] = string(args[i])
		}
		f.WriteLn(name + "(" + strings.Join(inputsStr, ", ") + ")")
	}

	f.mDefines[name] = toReturn

	return toReturn
}

func (f *FFAmd64) AssertCleanStack(reservedStackSize, minStackSize int) {
	if f.nbElementsOnStack != 0 {
		panic("missing f.Push stack elements")
	}
	if reservedStackSize < minStackSize {
		panic("invalid minStackSize or reservedStackSize")
	}
	usedStackSize := f.maxOnStack * 8
	if usedStackSize > reservedStackSize {
		panic("using more stack size than reserved")
	} else if max(usedStackSize, minStackSize) < reservedStackSize {
		// this panic is for dev purposes as this may be by design for aligment
		panic("reserved more stack size than needed")
	}

	f.maxOnStack = 0
}

func (f *FFAmd64) Push(registers *amd64.Registers, rIn ...amd64.Register) {
	for _, r := range rIn {
		if strings.HasPrefix(string(r), "s") {
			// it's on the stack, decrease the offset
			f.nbElementsOnStack--
			continue
		}
		registers.Push(r)
	}
}

func (f *FFAmd64) Pop(registers *amd64.Registers, forceStack ...bool) amd64.Register {
	if registers.Available() >= 1 && !(len(forceStack) > 0 && forceStack[0]) {
		return registers.Pop()
	}
	r := amd64.Register(fmt.Sprintf("s%d-%d(SP)", f.nbElementsOnStack, 8+f.nbElementsOnStack*8))
	f.nbElementsOnStack++
	if f.nbElementsOnStack > f.maxOnStack {
		f.maxOnStack = f.nbElementsOnStack
	}
	return r
}

func (f *FFAmd64) PopN(registers *amd64.Registers, forceStack ...bool) []amd64.Register {
	if len(forceStack) > 0 && forceStack[0] {
		nbStack := f.NbWords
		var u []amd64.Register

		for i := f.nbElementsOnStack; i < nbStack+f.nbElementsOnStack; i++ {
			u = append(u, amd64.Register(fmt.Sprintf("s%d-%d(SP)", i, 8+i*8)))
		}
		f.nbElementsOnStack += nbStack
		if f.nbElementsOnStack > f.maxOnStack {
			f.maxOnStack = f.nbElementsOnStack
		}
		return u
	}
	if registers.Available() >= f.NbWords {
		return registers.PopN(f.NbWords)
	}
	nbStack := f.NbWords - registers.Available()
	u := registers.PopN(registers.Available())

	for i := f.nbElementsOnStack; i < nbStack+f.nbElementsOnStack; i++ {
		u = append(u, amd64.Register(fmt.Sprintf("s%d-%d(SP)", i, 8+i*8)))
	}
	f.nbElementsOnStack += nbStack
	if f.nbElementsOnStack > f.maxOnStack {
		f.maxOnStack = f.nbElementsOnStack
	}
	return u
}

func (f *FFAmd64) qAt(index int) string {
	return fmt.Sprintf("·qElement+%d(SB)", index*8)
}

func (f *FFAmd64) qAt_bcst(index int) string {
	return fmt.Sprintf("·qElement+%d(SB)", index*4)
}

func (f *FFAmd64) qInv0() string {
	return "$const_qInvNeg"
}

func (f *FFAmd64) mu() string {
	return "$const_mu"
}

// GenerateCommonASM generates assembly code for the base field provided to goff
// see internal/templates/ops*
func GenerateCommonASM(w io.Writer, nbWords, nbBits int, hasVector bool) error {
	f := NewFFAmd64(w, nbWords)
	f.Comment("Code generated by gnark-crypto/generator. DO NOT EDIT.")

	f.WriteLn("#include \"textflag.h\"")
	f.WriteLn("#include \"funcdata.h\"")
	f.WriteLn("#include \"go_asm.h\"")
	f.WriteLn("")

	if nbWords == 1 {
		if nbBits == 31 {
			return GenerateF31ASM(f, hasVector)
		} else {
			panic("not implemented")
		}
	}

	f.GenerateReduceDefine()

	// reduce
	f.generateReduce()

	// mul by constants
	f.generateMulBy3()
	f.generateMulBy5()
	f.generateMulBy13()

	// fft butterflies
	f.generateButterfly()

	// mul
	f.generateMul(false)

	// from mont
	f.generateFromMont(false)

	if hasVector {
		f.WriteLn("")
		f.Comment("Vector operations are partially derived from Dag Arne Osvik's work in github.com/a16z/vectorized-fields")
		f.WriteLn("")

		f.generateAddVecW4()
		f.generateSubVecW4()
		f.generateSumVecW4()
		f.generateInnerProductW4()
		f.generateMulVecW4("scalarMulVec")
		f.generateMulVecW4("mulVec")
	}

	return nil
}

func GenerateF31FFTKernels(w io.Writer, nbBits int, kernels []int) error {
	if nbBits != 31 {
		return fmt.Errorf("only 31 bits supported for now")
	}
	f := NewFFAmd64(w, 1)
	f.Comment("Code generated by gnark-crypto/generator. DO NOT EDIT.")

	f.WriteLn("#include \"textflag.h\"")
	f.WriteLn("#include \"funcdata.h\"")
	f.WriteLn("#include \"go_asm.h\"")
	f.WriteLn("")

	f.generateFFTKernelF31(kernels[0])

	return nil
}

func GenerateF31ASM(f *FFAmd64, hasVector bool) error {
	if !hasVector {
		return nil // nothing for now.
	}

	f.generateAddVecF31()
	f.generateSubVecF31()
	f.generateSumVecF31()
	f.generateMulVecF31()
	f.generateScalarMulVecF31()
	f.generateInnerProdVecF31()

	return nil
}

func ElementASMFileName(nbWords, nbBits int) string {
	const nameW1 = "element_%db_amd64.s"
	const nameWN = "element_%dw_amd64.s"
	if nbWords == 1 {
		if nbBits >= 32 {
			panic("not implemented")
		}
		return fmt.Sprintf(nameW1, 31)
	}
	return fmt.Sprintf(nameWN, nbWords)
}
