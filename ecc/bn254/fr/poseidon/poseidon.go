package poseidon

import (
	"errors"
	"hash"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
)

const (
	BlockSize = fr.Bytes // BlockSize size that poseidon consumes
)

func zeroElement() *fr.Element {
	return &fr.Element{0, 0, 0, 0}
}

func deepCopy(dst, src []*fr.Element) {
	if len(src) > len(dst) {
		panic("Cannot copy to a smaller destination")
	}
	for i := 0; i < len(src); i++ {
		v := *src[i]
		dst[i] = &v
	}
}

// Add round constants
func arc(state []*fr.Element, C []*fr.Element, t, offset int) {
	for i := 0; i < t; i++ {
		state[i].Add(state[i], C[offset+i])
	}
}

// power 5 as s-box for full state
func sbox(state []*fr.Element, t int) {
	for i := 0; i < t; i++ {
		state[i].Exp(*state[i], alpha)
	}
}

// Matrix vector multiplication
func mix(state []*fr.Element, M [][]*fr.Element, t int) []*fr.Element {
	newState := make([]*fr.Element, t)

	for i := 0; i < t; i++ {
		newState[i] = zeroElement()
		for j := 0; j < t; j++ {
			newState[i].Add(newState[i], zeroElement().Mul(M[j][i], state[j]))
		}
	}
	return newState
}

func permutation(state []*fr.Element) []*fr.Element {
	// Minimum length of state = nInput + nOutput = 2
	t := len(state)
	index := t - 2
	RP := rp[index]
	C := c[index]
	M := m[index]
	S := s[index]
	P := p[index]

	// 1. Pre-step to the first-half of full rounds: add round constant for round=0
	arc(state, C, t, 0)

	// 2. First-half of full rounds starting at roundNumber = 1 except last round
	for i := 0; i < rf/2-1; i++ {
		sbox(state, t)
		arc(state, C, t, (i+1)*t)
		state = mix(state, M, t)
	}

	// 3. Last round of first-half of full rounds
	sbox(state, t)
	arc(state, C, t, (rf/2)*t)
	state = mix(state, P, t)

	// 4. Partial rounds
	for i := 0; i < RP; i++ {
		state[0].Exp(*state[0], alpha)
		state[0].Add(state[0], C[(rf/2+1)*t+i])
		// S[i] is a vector of [t*2-1] elements where first t elements are used to compute state[0]
		// and the remaining elements starting at [t] are used to compute state[1,..,t-1]
		offset := (t*2 - 1) * i
		newState0 := zeroElement()
		for j := 0; j < len(state); j++ {
			newState0.Add(newState0, zeroElement().Mul(state[j], S[offset+j]))
		}
		offset += t - 1
		for k := 1; k < t; k++ {
			state[k].Add(state[k], zeroElement().Mul(state[0], S[offset+k]))
		}
		state[0] = newState0
	}

	// 5. Second-half of full rounds except last round
	for i := 0; i < rf/2-1; i++ {
		sbox(state, t)
		arc(state, C, t, (rf/2+1)*t+RP+i*t)
		state = mix(state, M, t)
	}

	// 6. Last round of the second-half of full rounds
	sbox(state, t)
	state = mix(state, M, t)
	return state
}

func Poseidon(input ...*fr.Element) *fr.Element {
	inputLength := len(input)
	if inputLength == 0 {
		panic("No support for dummy input")
	}

	const maxLength = 12
	state := make([]*fr.Element, maxLength+1)
	state[0] = zeroElement()
	startIndex := 0
	lastIndex := 0

	// Make a hash chain of the input if its length > maxLength
	if inputLength > maxLength {
		count := inputLength / maxLength
		for i := 0; i < count; i++ {
			lastIndex = (i + 1) * maxLength
			deepCopy(state[1:], input[startIndex:lastIndex])
			state = permutation(state)
			startIndex = lastIndex
		}
	}

	// For the remaining part of the input OR if 2 <= inputLength <= 12
	if lastIndex < inputLength {
		lastIndex = inputLength
		remainigLength := lastIndex - startIndex
		deepCopy(state[1:], input[startIndex:lastIndex])
		state = permutation(state[:remainigLength+1])
	}
	return state[1]
}

func PoseidonBytes(input ...[]byte) []byte {
	inputElements := make([]*fr.Element, len(input))
	for i, ele := range input {
		num := new(big.Int).SetBytes(ele)
		if num.Cmp(fr.Modulus()) >= 0 {
			panic("not support bytes bigger than modulus")
		}
		e := fr.Element{0, 0, 0, 0}
		e.SetBigInt(num)
		inputElements[i] = &e
	}
	res := Poseidon(inputElements...).Bytes()
	return res[:]
}

type digest struct {
	h    fr.Element
	data [][]byte // data to hash
}

func NewPoseidon() hash.Hash {
	d := new(digest)
	d.Reset()
	return d
}

// Reset resets the Hash to its initial state.
func (d *digest) Reset() {
	d.data = nil
	d.h = fr.Element{0, 0, 0, 0}
}

// Only receive byte slice less than fr.Modulus()
func (d *digest) Write(p []byte) (n int, err error) {
	n = len(p)
	num := new(big.Int).SetBytes(p)
	if num.Cmp(fr.Modulus()) >= 0 {
		return 0, errors.New("not support bytes bigger than modulus")
	}
	d.data = append(d.data, p)
	return n, nil
}

func (d *digest) Size() int {
	return BlockSize
}

// BlockSize returns the number of bytes Sum will return.
func (d *digest) BlockSize() int {
	return BlockSize
}

// Sum appends the current hash to b and returns the resulting slice.
// It does not change the underlying hash state.
func (d *digest) Sum(b []byte) []byte {
	e := fr.Element{0, 0, 0, 0}
	e.SetBigInt(new(big.Int).SetBytes(PoseidonBytes(d.data...)))
	d.h = e
	d.data = nil // flush the data already hashed
	hash := d.h.Bytes()
	b = append(b, hash[:]...)
	return b
}
