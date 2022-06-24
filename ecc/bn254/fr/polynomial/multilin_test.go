package polynomial

import (
	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/gen"
	"github.com/leanovate/gopter/prop"
	"testing"
)

//TODO: Property based tests?
func TestFoldBilinear(t *testing.T) {

	for i := 0; i < 100; i++ {

		// f = c₀ + c₁ X₁ + c₂ X₂ + c₃ X₁ X₂
		var coefficients [4]fr.Element
		for i := 0; i < 4; i++ {
			if _, err := coefficients[i].SetRandom(); err != nil {
				t.Error(err)
			}
		}

		var r fr.Element
		if _, err := r.SetRandom(); err != nil {
			t.Error(err)
		}

		// interpolate at {0,1}²:
		m := make(MultiLin, 4)
		m[0] = coefficients[0]
		m[1].Add(&coefficients[0], &coefficients[2])
		m[2].Add(&coefficients[0], &coefficients[1])
		m[3].
			Add(&m[1], &coefficients[1]).
			Add(&m[3], &coefficients[3])

		m.Fold(r)

		// interpolate at {r}×{0,1}:
		var expected0, expected1 fr.Element
		expected0.
			Mul(&r, &coefficients[1]).
			Add(&expected0, &coefficients[0])

		expected1.
			Mul(&r, &coefficients[3]).
			Add(&expected1, &coefficients[2]).
			Add(&expected0, &expected1)

		if !m[0].Equal(&expected0) || !m[1].Equal(&expected1) {
			t.Fail()
		}
	}
}

func TestPrecomputeLagrange(t *testing.T) {

	testForDomainSize := func(domainSize uint8) bool {
		polys := precomputeLagrangeCoefficients(domainSize)

		for l := uint8(0); l < domainSize; l++ {
			for i := uint8(0); i < domainSize; i++ {
				var I fr.Element
				I.SetUint64(uint64(i))
				y := polys[l].Eval(&I)

				if i == l && !y.IsOne() || i != l && !y.IsZero() {
					t.Errorf("domainSize = %d: p_%d(%d) = %s", domainSize, l, i, signedText(&y, 10))
					return false
				}
			}
		}
		return true
	}

	//testForDomainSize(4)

	t.Parallel()
	parameters := gopter.DefaultTestParameters()

	parameters.MinSuccessfulTests = int(maxLagrangeDomainSize)

	properties := gopter.NewProperties(parameters)

	properties.Property("l'th lagrange polynomials must evaluate to 1 on l and 0 on other values in the domain", prop.ForAll(
		testForDomainSize,
		gen.UInt8Range(2, maxLagrangeDomainSize),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// TODO: Benchmark folding? Algorithms is pretty straightforward; unless we want to measure how well memory management is working

func signedText(v *fr.Element, base int) string {
	i := signedBigInt(v)
	return i.Text(base)
}