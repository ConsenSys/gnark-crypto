package element

const VectorOpsPureGo = `
// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	addVecGeneric(*vector, a, b)
}

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	subVecGeneric(*vector, a, b)
}

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *{{.ElementName}}) {
	scalarMulVecGeneric(*vector, a, b)
}

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res {{.ElementName}}) {
	sumVecGeneric(&res, *vector)
	return
}

// InnerProduct computes the inner product of two vectors.
// It panics if the vectors don't have the same length.
func (vector *Vector) InnerProduct(other Vector) (res {{.ElementName}}) {
	innerProductVecGeneric(&res, *vector, other)
	return
}

// Mul multiplies two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Mul(a, b Vector) {
	mulVecGeneric(*vector, a, b)
}

// ButterflyMul used in FFT; TODO complete.
func (vector *Vector) ButterflyMul(twiddles Vector, start, end, m int) {
	butterflyMulVecGeneric(*vector, twiddles, start, end, m)
}
`
