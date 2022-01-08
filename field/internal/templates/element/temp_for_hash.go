package element

const TempForHash = `

// Sgn0 is an algebraic substitute for the notion of sign in ordered fields
// Namely, every non-zero quadratic residue in a finite field of characteristic =/= 2 has exactly two square roots, one of each sign
// Taken from https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/ section 4.1
// The sign of an element is not obviously related to that of its Montgomery form
func (z *{{.ElementName}}) Sgn0() bool {
	nonMont := *z
	nonMont.FromMont()

	return nonMont[0]%2 == 1
}

func (z *{{.ElementName}}) SetHex(hex string) {
	var i big.Int
	i.SetString(hex, 16)
	if _, b := i.SetString(hex, 16); !b {
		panic("SetString failed")
	}
	z.SetBigInt(&i)
}

func (z *{{.ElementName}}) EvalPolynomialHex(x *{{.ElementName}}, cHex []string) {
	c := make([]{{.ElementName}}, len(cHex))

	for i, hex := range cHex {
		c[i].SetHex(hex)
	}

	z.EvalPolynomial(x, c)
}

func (z *{{.ElementName}}) EvalPolynomial(x *{{.ElementName}}, c []{{.ElementName}}) {
	f := c[len(c)-1]

	for i := len(c) - 2; i >= 0; i-- {
		f.Mul(&f, x)
		f.Add(&f, &c[i])
	}

	*z = f
}
`