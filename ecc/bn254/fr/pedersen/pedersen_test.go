package pedersen

import (
	"github.com/consensys/gnark-crypto/ecc/bn254"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/stretchr/testify/assert"
	"math/rand"
	"testing"
)

func interfaceSliceToFrPtrSlice(t *testing.T, values ...interface{}) []*fr.Element {
	res := make([]*fr.Element, len(values))
	for i, v := range values {
		var V fr.Element
		_, err := V.SetInterface(v)
		assert.NoError(t, err)
		res[i] = &V
	}
	return res
}

func randomFrSlice(t *testing.T, size int) []interface{} {
	res := make([]interface{}, size)
	var err error
	for i := range res {
		var v fr.Element
		res[i], err = v.SetRandom()
		assert.NoError(t, err)
	}
	return res
}

func randomOnG1() (bn254.G1Affine, error) { // TODO: Add to G1.go?
	gBytes := make([]byte, fr.Bytes)
	if _, err := rand.Read(gBytes); err != nil {
		return bn254.G1Affine{}, err
	}
	return bn254.HashToG1(gBytes, []byte("random on g2"))
}

func testCommit(t *testing.T, values ...interface{}) {
	// g₁ g₂

	basis := make([]*bn254.G1Affine, len(values))
	for i := range basis {
		basisI, err := randomOnG1()
		assert.NoError(t, err)
		basis[i] = &basisI
	}

	key, err := Setup(basis)
	commitment, pok, err := key.Commit(interfaceSliceToFrPtrSlice(t, values...))
	assert.NoError(t, err)
	assert.NoError(t, key.VerifyKnowledgeProof(commitment, pok))

	pok.Neg(&pok)
	assert.NotNil(t, key.VerifyKnowledgeProof(commitment, pok))
}

func TestCommitToOne(t *testing.T) {
	testCommit(t, 1)
}

func TestCommitSingle(t *testing.T) {
	testCommit(t, randomFrSlice(t, 1)...)
}

func TestCommitFiveElements(t *testing.T) {
	testCommit(t, randomFrSlice(t, 5)...)
}
