// Copyright 2020 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package kzg

import (
	"fmt"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/consensys/gnark-crypto/ecc"
	curve "github.com/consensys/gnark-crypto/ecc/bw6-761"
	"github.com/consensys/gnark-crypto/ecc/bw6-761/fr"
	"github.com/consensys/gnark-crypto/internal/parallel"
)

// ToLagrangeG1 in place transform of coeffs canonical form into Lagrange form.
// From the formula Lᵢ(τ) = 1/n∑_{j<n}(τ/ωⁱ)ʲ we
// see that [L₁(τ),..,Lₙ(τ)] = FFT_inv(∑_{j<n}τʲXʲ), so it suffices to apply the inverse
// fft on the vector consisting of the original SRS.
// Size of coeffs must be a power of 2.
func ToLagrangeG1(coeffs []curve.G1Affine) ([]curve.G1Affine, error) {
	if bits.OnesCount64(uint64(len(coeffs))) != 1 {
		return nil, fmt.Errorf("len(coeffs) must be a power of 2")
	}
	size := len(coeffs)

	numCPU := uint64(runtime.NumCPU())
	maxSplits := bits.TrailingZeros64(ecc.NextPowerOfTwo(numCPU)) << 1

	twiddlesInv, err := computeTwiddlesInv(size)
	if err != nil {
		return nil, err
	}

	// batch convert to Jacobian
	jCoeffs := make([]curve.G1Jac, len(coeffs))
	for i := 0; i < len(coeffs); i++ {
		jCoeffs[i].FromAffine(&coeffs[i])
	}

	difFFTG1(jCoeffs, twiddlesInv, 0, maxSplits, nil)

	// TODO @gbotrel generify the cobra bitreverse function, benchmark it and use it everywhere
	bitReverse(jCoeffs)

	var invBigint big.Int
	var frCardinality fr.Element
	frCardinality.SetUint64(uint64(size))
	frCardinality.Inverse(&frCardinality)
	frCardinality.BigInt(&invBigint)

	parallel.Execute(size, func(start, end int) {
		for i := start; i < end; i++ {
			jCoeffs[i].ScalarMultiplication(&jCoeffs[i], &invBigint)
		}
	})

	// batch convert to affine
	return curve.BatchJacobianToAffineG1(jCoeffs), nil
}

func computeTwiddlesInv(cardinality int) ([]*big.Int, error) {
	generator, err := fr.Generator(uint64(cardinality))
	if err != nil {
		return nil, err
	}

	// inverse the generator
	generator.Inverse(&generator)

	// nb fft stages
	nbStages := uint64(bits.TrailingZeros64(uint64(cardinality)))

	r := make([]*big.Int, 1+(1<<(nbStages-1)))

	w := generator
	r[0] = new(big.Int).SetUint64(1)
	if len(r) == 1 {
		return r, nil
	}
	r[1] = new(big.Int)
	w.BigInt(r[1])
	for j := 2; j < len(r); j++ {
		w.Mul(&w, &generator)
		r[j] = new(big.Int)
		w.BigInt(r[j])
	}

	return r, nil
}

func bitReverse[T any](a []T) {
	n := uint64(len(a))
	nn := uint64(64 - bits.TrailingZeros64(n))

	for i := uint64(0); i < n; i++ {
		irev := bits.Reverse64(i) >> nn
		if irev > i {
			a[i], a[irev] = a[irev], a[i]
		}
	}
}

func butterflyG1(a *curve.G1Jac, b *curve.G1Jac) {
	t := *a
	a.AddAssign(b)
	t.SubAssign(b)
	b.Set(&t)
}

func difFFTG1(a []curve.G1Jac, twiddles []*big.Int, stage, maxSplits int, chDone chan struct{}) {
	if chDone != nil {
		defer close(chDone)
	}

	n := len(a)
	if n == 1 {
		return
	}
	m := n >> 1

	butterflyG1(&a[0], &a[m])
	// stage determines the stride
	// if stage == 0, then we use 1, w, w**2, w**3, w**4, w**5, w**6, ...
	// if stage == 1, then we use 1, w**2, w**4, w**6, ... that is, indexes 0, 2, 4, 6, ... of stage 0
	// if stage == 2, then we use 1, w**4, w**8, w**12, ... that is indexes 0, 4, 8, 12, ... of stage 0
	stride := 1 << stage

	const butterflyThreshold = 8
	if m >= butterflyThreshold {
		// 1 << stage == estimated used CPUs
		numCPU := runtime.NumCPU() / (1 << (stage))
		parallel.Execute(m, func(start, end int) {
			if start == 0 {
				start = 1
			}
			j := start * stride
			for i := start; i < end; i++ {
				butterflyG1(&a[i], &a[i+m])
				a[i+m].ScalarMultiplication(&a[i+m], twiddles[j])
				j += stride
			}
		}, numCPU)
	} else {
		j := stride
		for i := 1; i < m; i++ {
			butterflyG1(&a[i], &a[i+m])
			a[i+m].ScalarMultiplication(&a[i+m], twiddles[j])
			j += stride
		}
	}

	if m == 1 {
		return
	}

	nextStage := stage + 1
	if stage < maxSplits {
		chDone := make(chan struct{}, 1)
		go difFFTG1(a[m:n], twiddles, nextStage, maxSplits, chDone)
		difFFTG1(a[0:m], twiddles, nextStage, maxSplits, nil)
		<-chDone
	} else {
		difFFTG1(a[0:m], twiddles, nextStage, maxSplits, nil)
		difFFTG1(a[m:n], twiddles, nextStage, maxSplits, nil)
	}
}
