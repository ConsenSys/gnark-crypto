package sumcheck

import (
	"fmt"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr/polynomial"
	"math/bits"
	"testing"
)

type singleMultilinClaim struct {
	g polynomial.MultiLin
}

func (c singleMultilinClaim) ProveFinalEval([]fr.Element) interface{} {
	return nil // verifier can compute the final eval itself
}

func (c singleMultilinClaim) VarsNum() int {
	return bits.TrailingZeros(uint(len(c.g)))
}

func (c singleMultilinClaim) ClaimsNum() int {
	return 1
}

func sumForX1One(g polynomial.MultiLin) polynomial.Polynomial {
	sum := g[len(g)/2]
	for i := len(g)/2 + 1; i < len(g); i++ {
		sum.Add(&sum, &g[i])
	}
	return []fr.Element{sum}
}

func (c singleMultilinClaim) Combine(fr.Element) polynomial.Polynomial {
	return sumForX1One(c.g)
}

func (c *singleMultilinClaim) Next(r fr.Element) polynomial.Polynomial {
	c.g.Fold(r)
	return sumForX1One(c.g)
}

type singleMultilinLazyClaim struct {
	g          polynomial.MultiLin
	claimedSum fr.Element
}

func (c singleMultilinLazyClaim) VerifyFinalEval(r []fr.Element, _ fr.Element, purportedValue fr.Element, _ interface{}) bool {
	val := c.g.Evaluate(r)
	return val.Equal(&purportedValue)
}

func (c singleMultilinLazyClaim) CombinedSum(fr.Element) fr.Element {
	return c.claimedSum
}

func (c singleMultilinLazyClaim) Degree(int) int {
	return 1
}

func (c singleMultilinLazyClaim) ClaimsNum() int {
	return 1
}

func (c singleMultilinLazyClaim) VarsNum() int {
	return bits.TrailingZeros(uint(len(c.g)))
}

func testSumcheckSingleClaimMultilin(polyInt []uint64, hashGenerator func() ArithmeticTranscript) bool {
	poly := make(polynomial.MultiLin, len(polyInt))
	for i, n := range polyInt {
		poly[i].SetUint64(n)
	}

	claim := singleMultilinClaim{g: poly.Clone()}

	proof := Prove(&claim, hashGenerator())

	fmt.Print("For hash ", hashGenerator(), " and poly ")
	printPoly(poly)
	fmt.Print("Proof = ")
	printProof(proof)

	lazyClaim := singleMultilinLazyClaim{g: poly, claimedSum: poly.Sum()}

	return Verify(lazyClaim, proof, hashGenerator())
}

// For debugging TODO Remove
func printMsws(limit int) {

	for i := -limit; i <= limit; i++ {
		if i == 0 {
			continue
		}
		var iElem fr.Element
		iElem.SetInt64(int64(i))
		fmt.Printf("%d: %d\n", i, iElem[fr.Limbs-1])
	}
}

func printPoly(poly []fr.Element) {
	text := make([]string, len(poly))
	for i, element := range poly {
		text[i] = element.Text(10)
	}
	fmt.Println(text)
}

func printProof(proof Proof) {
	fmt.Println("[")

	for _, line := range proof.PartialSumPolys {
		fmt.Print("\t")
		printPoly(line)
	}

	fmt.Println("],", proof.FinalEvalProof)
}

// Autogen proof for SNARK circuit to verify (the others are not usable since they require "interpolation" on points where the polynomial is already specified)
func TestSumcheckDeterministicHashSingleClaimMultilinSnark(t *testing.T) {
	if !testSumcheckSingleClaimMultilin([]uint64{1, 2, 3, 4}, NewMessageCounterGenerator(0, 1)) {
		t.Error()
	}
}

func TestSumcheckDeterministicHashSingleClaimMultilin(t *testing.T) {
	//printMsws(100)

	polys := [][]uint64{
		{1, 2, 3, 4},             // 1 + 2X₁ + X₂
		{1, 2, 3, 4, 5, 6, 7, 8}, // 1 + 4X₁ + 2X₂ + X₃
		{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, // 1 + 8X₁ + 4X₂ + 2X₃ + X₄
	}

	const MaxStep = 4
	const MaxStart = 4
	hashGens := make([]func() ArithmeticTranscript, 0, MaxStart*MaxStep)

	for step := 0; step < MaxStep; step++ {
		for startState := 0; startState < MaxStart; startState++ {
			hashGens = append(hashGens, NewMessageCounterGenerator(startState, step))
		}
	}

	for _, poly := range polys {
		for _, hashGen := range hashGens {
			if !testSumcheckSingleClaimMultilin(poly, hashGen) {
				t.Error(poly, hashGen())
			}
		}
	}
}