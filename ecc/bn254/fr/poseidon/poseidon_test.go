package poseidon

import (
	"math/big"
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/stretchr/testify/assert"
)

func elementFromHexString(v string) *fr.Element {
	n, success := new(big.Int).SetString(v, 16)
	if !success {
		panic("Error parsing hex number")
	}
	e := fr.Element{0, 0, 0, 0}
	e.SetBigInt(n)
	return &e
}

func TestPoseidon1(t *testing.T) {
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	inputsStr := []string{"1"}
	expectedHash := elementFromHexString("112A4F9241E384B0EDE4655E6D2BBF7EBD9595775DE9E7536DF87CD487852FC4")
	inputs := make([]*fr.Element, len(inputsStr))
	for i := 0; i < len(inputsStr); i++ {
		inputs[i] = elementFromHexString(inputsStr[i])
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon2(t *testing.T) {
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	inputsStr := []string{"1", "2"}
	expectedHash := elementFromHexString("FCA49B798923AB0239DE1C9E7A4A9A2210312B6A2F616D18B5A87F9B628AE29")
	inputs := make([]*fr.Element, len(inputsStr))
	for i := 0; i < len(inputsStr); i++ {
		inputs[i] = elementFromHexString(inputsStr[i])
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon4(t *testing.T) {
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	inputsStr := []string{"1", "2", "3", "4"}
	expectedHash := elementFromHexString("1148AAEF609AA338B27DAFD89BB98862D8BB2B429ACEAC47D86206154FFE053D")
	inputs := make([]*fr.Element, len(inputsStr))
	for i := 0; i < len(inputsStr); i++ {
		inputs[i] = elementFromHexString(inputsStr[i])
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon13(t *testing.T) {
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	expectedHash := elementFromHexString("4EA9108A1551C780A7408570BDEAA3A0294B01F21198B72FF01545A60DA677F")
	length := 13
	inputs := make([]*fr.Element, length)
	for i := 0; i < length; i++ {
		e := fr.NewElement((uint64)(i + 1))
		inputs[i] = &e
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon24(t *testing.T) {
	// WARNING: No test vector to compare with
	expectedHash := elementFromHexString("6C7676E83EF8CB9EF6C25746A5F6B2D39FBA4548B4C29B3D41490BBF3C1108D")
	length := 24
	inputs := make([]*fr.Element, length)
	for i := 0; i < length; i++ {
		e := fr.NewElement((uint64)(i + 1))
		inputs[i] = &e
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon30(t *testing.T) {
	// WARNING: No test vector to compare with
	expectedHash := elementFromHexString("2FF47AB8E9E9F6134600A8DE8B8E99596E573620A7D8D39ED7B2C7CEF9F105F1")
	length := 30
	inputs := make([]*fr.Element, length)
	for i := 0; i < length; i++ {
		e := fr.NewElement((uint64)(i + 1))
		inputs[i] = &e
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestPoseidon256(t *testing.T) {
	// WARNING: No test vector to compare with
	expectedHash := elementFromHexString("182AF1C3FFD14FA66CDF5FE5D5199473678F221CA3BAB09B44758EF80641C1E0")
	length := 256
	inputs := make([]*fr.Element, length)
	for i := 0; i < length; i++ {
		e := fr.NewElement((uint64)(i + 1))
		inputs[i] = &e
	}
	actualHash := Poseidon(inputs...)
	assert.True(t, actualHash.Equal(expectedHash), "%s != %s", actualHash, expectedHash)
}

func TestConsistency(t *testing.T) {
	// Check whether Poseidon returns the same value for the same input
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	inputsStr := []string{"1", "2", "3", "4"}
	inputs := make([]*fr.Element, len(inputsStr))
	for i := 0; i < len(inputsStr); i++ {
		inputs[i] = elementFromHexString(inputsStr[i])
	}
	actualHash1 := Poseidon(inputs...)
	actualHash2 := Poseidon(inputs...)
	assert.True(t, actualHash1.Equal(actualHash2), "%s != %s", actualHash1, actualHash2)
}

func TestPoseidonBytes(t *testing.T) {
	// Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
	expectedHash := elementFromHexString("FCA49B798923AB0239DE1C9E7A4A9A2210312B6A2F616D18B5A87F9B628AE29")
	inputs := make([][]byte, 2)
	inputs[0] = make([]byte, 1)
	inputs[0][0] = 1
	inputs[1] = make([]byte, 1)
	inputs[1][0] = 2
	actualHash := PoseidonBytes(inputs...)
	actualHashEle := fr.Element{0, 0, 0, 0}
	actualHashEle.SetBytes(actualHash)
	assert.True(t, actualHashEle.Equal(expectedHash), "%s != %s", actualHashEle, expectedHash)
}

func TestDigest(t *testing.T) {
	expectedHash := elementFromHexString("FCA49B798923AB0239DE1C9E7A4A9A2210312B6A2F616D18B5A87F9B628AE29")
	hFunc := NewPoseidon()
	inputs := make([][]byte, 2)
	inputs[0] = make([]byte, 1)
	inputs[0][0] = 1
	inputs[1] = make([]byte, 1)
	inputs[1][0] = 2
	hFunc.Write(inputs[0])
	hFunc.Write(inputs[1])
	actualHash := hFunc.Sum(nil)
	actualHashEle := fr.Element{0, 0, 0, 0}
	actualHashEle.SetBytes(actualHash)
	assert.True(t, actualHashEle.Equal(expectedHash), "%s != %s", actualHashEle, expectedHash)

	hFunc.Reset()
	bigNumber, _ := new(big.Int).SetString("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10)
	inputs[0] = bigNumber.Bytes()
	n, err := hFunc.Write(inputs[0])
	assert.EqualError(t, err, "not support bytes bigger than modulus")
	assert.Equal(t, n, 0)
}
