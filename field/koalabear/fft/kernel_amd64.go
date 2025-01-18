//go:build !purego

// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fft

import (
	"github.com/consensys/gnark-crypto/field/koalabear"
	"golang.org/x/sys/cpu"
)

var (
	supportAVX512 = cpu.X86.HasAVX512 && cpu.X86.HasAVX512DQ
)

// q + r'.r = 1, i.e., qInvNeg = - q⁻¹ mod r
// used for Montgomery reduction
const qInvNeg = 2130706431
const q = 2130706433

//go:noescape
func kerDIFNP_32_avx512(a []koalabear.Element, twiddles [][]koalabear.Element, stage int)

func kerDIFNP_32(a []koalabear.Element, twiddles [][]koalabear.Element, stage int) {
	if !supportAVX512 {
		kerDIFNP_32generic(a, twiddles, stage)
		return
	}
	kerDIFNP_32_avx512(a, twiddles, stage)

	for offset := 0; offset < 32; offset += 8 {
		innerDIFWithTwiddles(a[offset:offset+8], twiddles[stage+2], 0, 4, 4)
	}
	for offset := 0; offset < 32; offset += 4 {
		innerDIFWithTwiddles(a[offset:offset+4], twiddles[stage+3], 0, 2, 2)
	}
	for offset := 0; offset < 32; offset += 2 {
		koalabear.Butterfly(&a[offset], &a[offset+1])
	}
}
func kerDITNP_32(a []koalabear.Element, twiddles [][]koalabear.Element, stage int) {
	kerDITNP_32generic(a, twiddles, stage)
}

//go:noescape
func kerDIFNP_256_avx512(a []koalabear.Element, twiddles [][]koalabear.Element, stage int)

func kerDIFNP_256(a []koalabear.Element, twiddles [][]koalabear.Element, stage int) {
	if !supportAVX512 {
		kerDIFNP_256generic(a, twiddles, stage)
		return
	}
	kerDIFNP_256_avx512(a, twiddles, stage)

	for offset := 0; offset < 256; offset += 8 {
		innerDIFWithTwiddles(a[offset:offset+8], twiddles[stage+5], 0, 4, 4)
	}
	for offset := 0; offset < 256; offset += 4 {
		innerDIFWithTwiddles(a[offset:offset+4], twiddles[stage+6], 0, 2, 2)
	}
	for offset := 0; offset < 256; offset += 2 {
		koalabear.Butterfly(&a[offset], &a[offset+1])
	}
}
func kerDITNP_256(a []koalabear.Element, twiddles [][]koalabear.Element, stage int) {
	kerDITNP_256generic(a, twiddles, stage)
}