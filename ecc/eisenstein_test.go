package ecc

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

const (
	nbFuzzShort = 10
	nbFuzz      = 50
)

func TestEisensteinReceiverIsOperand(t *testing.T) {

	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	genE := GenComplexNumber()

	properties.Property("Having the receiver as operand (addition) should output the same result", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c, d ComplexNumber
			d.Set(a)
			c.Add(a, b)
			a.Add(a, b)
			b.Add(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genE,
		genE,
	))

	properties.Property("Having the receiver as operand (sub) should output the same result", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c, d ComplexNumber
			d.Set(a)
			c.Sub(a, b)
			a.Sub(a, b)
			b.Sub(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genE,
		genE,
	))

	properties.Property("Having the receiver as operand (mul) should output the same result", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c, d ComplexNumber
			d.Set(a)
			c.Mul(a, b)
			a.Mul(a, b)
			b.Mul(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genE,
		genE,
	))

	properties.Property("Having the receiver as operand (neg) should output the same result", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			b.Neg(a)
			a.Neg(a)
			return a.Equal(&b)
		},
		genE,
	))

	properties.Property("Having the receiver as operand (conjugate) should output the same result", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			b.Conjugate(a)
			a.Conjugate(a)
			return a.Equal(&b)
		},
		genE,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestEisensteinArithmetic(t *testing.T) {

	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	genE := GenComplexNumber()

	properties.Property("sub & add should leave an element invariant", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c ComplexNumber
			c.Set(a)
			c.Add(&c, b).Sub(&c, b)
			return c.Equal(a)
		},
		genE,
		genE,
	))

	properties.Property("neg twice should leave an element invariant", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			b.Neg(a).Neg(&b)
			return a.Equal(&b)
		},
		genE,
	))

	properties.Property("conj twice should leave an element invariant", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			b.Conjugate(a).Conjugate(&b)
			return a.Equal(&b)
		},
		genE,
	))

	properties.Property("add zero should leave element invariant", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			zero := new(ComplexNumber)
			b.Add(a, zero)
			return a.Equal(&b)
		},
		genE,
	))

	properties.Property("mul by one should leave element invariant", prop.ForAll(
		func(a *ComplexNumber) bool {
			var b ComplexNumber
			one := &ComplexNumber{
				*big.NewInt(1),
				*big.NewInt(0),
			}
			b.Mul(a, one)
			return a.Equal(&b)
		},
		genE,
	))

	properties.Property("add should be commutative", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c, d ComplexNumber
			c.Add(a, b)
			d.Add(b, a)
			return c.Equal(&d)
		},
		genE,
		genE,
	))

	properties.Property("add should be assiocative", prop.ForAll(
		func(a, b, c *ComplexNumber) bool {
			var d, e ComplexNumber
			d.Add(a, b).Add(&d, c)
			e.Add(c, b).Add(&e, a)
			return e.Equal(&d)
		},
		genE,
		genE,
		genE,
	))

	properties.Property("mul should be commutative", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			var c, d ComplexNumber
			c.Mul(a, b)
			d.Mul(b, a)
			return c.Equal(&d)
		},
		genE,
		genE,
	))

	properties.Property("mul should be assiocative", prop.ForAll(
		func(a, b, c *ComplexNumber) bool {
			var d, e ComplexNumber
			d.Mul(a, b).Mul(&d, c)
			e.Mul(c, b).Mul(&e, a)
			return e.Equal(&d)
		},
		genE,
		genE,
		genE,
	))

	properties.Property("norm should always be positive", prop.ForAll(
		func(a *ComplexNumber) bool {
			return a.Norm().Sign() > 0
		},
		genE,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestEisensteinHalfGCD(t *testing.T) {

	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = nbFuzzShort
	} else {
		parameters.MinSuccessfulTests = nbFuzz
	}

	properties := gopter.NewProperties(parameters)

	genE := GenComplexNumber()

	properties.Property("half-GCD", prop.ForAll(
		func(a, b *ComplexNumber) bool {
			res := HalfGCD(a, b)
			var c, d ComplexNumber
			c.Mul(b, res[1])
			d.Mul(a, res[2])
			d.Add(&c, &d)
			return d.Equal(res[0])
		},
		genE,
		genE,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// GenNumber generates a random integer
func GenNumber() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var prime, _ = new(big.Int).SetString("7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16) // 2^255 - 19
		elmt, _ := rand.Int(rand.Reader, prime)
		genResult := gopter.NewGenResult(*elmt, gopter.NoShrinker)
		return genResult
	}
}

// GenComplexNumber generates a random integer
func GenComplexNumber() gopter.Gen {
	return gopter.CombineGens(
		GenNumber(),
		GenNumber(),
	).Map(func(values []interface{}) *ComplexNumber {
		return &ComplexNumber{A0: values[0].(big.Int), A1: values[1].(big.Int)}
	})
}