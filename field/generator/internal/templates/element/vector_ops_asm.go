package element

const VectorOpsAmd64 = `
// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Add: vectors don't have the same length")
	}
		n := uint64(len(a))
	addVec(&(*vector)[0], &a[0], &b[0], n)
}

//go:noescape
func addVec(res, a, b *{{.ElementName}}, n uint64)

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Sub: vectors don't have the same length")
	}
	subVec(&(*vector)[0], &a[0], &b[0], uint64(len(a)))
}

//go:noescape
func subVec(res, a, b *{{.ElementName}}, n uint64)

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *{{.ElementName}}) {
	if len(a) != len(*vector) {
		panic("vector.ScalarMul: vectors don't have the same length")
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || uint64(len(a)) >= maxN {
		// call scalarMulVecGeneric
		scalarMulVecGeneric(*vector, a, b)
		return
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	// the code for scalarMul is identical to mulVec; and it expects at least 
	// 2 elements in the vector to fill the Z registers
	var bb  [2]{{.ElementName}}
	bb[0] = *b
	bb[1] = *b
	const blockSize = 16
	scalarMulVec(&(*vector)[0], &a[0], &bb[0], n/blockSize, qInvNeg)
	if n % blockSize != 0 {
		// call scalarMulVecGeneric on the rest
		start := n - n % blockSize
		scalarMulVecGeneric((*vector)[start:], a[start:], b)
	}
}

//go:noescape
func scalarMulVec(res, a, b *{{.ElementName}}, n uint64, qInvNeg uint64)

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res {{.ElementName}}) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	const minN = 16*7 // AVX512 slower than generic for small n
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n <= minN || n >= maxN {
		// call sumVecGeneric
		sumVecGeneric(&res, *vector)
		return
	}
	sumVec(&res, &(*vector)[0], uint64(len(*vector)))
	return
}

//go:noescape
func sumVec(res *{{.ElementName}}, a *{{.ElementName}}, n uint64)

// InnerProduct computes the inner product of two vectors.
// It panics if the vectors don't have the same length.
func (vector *Vector) InnerProduct(other Vector) (res {{.ElementName}}) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	if n != uint64(len(other)) {
		panic("vector.InnerProduct: vectors don't have the same length")
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n >= maxN {
		// call innerProductVecGeneric
		// note; we could split the vector into smaller chunks and call innerProductVec
		innerProductVecGeneric(&res, *vector, other)
		return
	}
	innerProdVec(&res[0], &(*vector)[0], &other[0], uint64(len(*vector)))

	return
}

//go:noescape
func innerProdVec(res *uint64, a,b *{{.ElementName}}, n uint64)

// Mul multiplies two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Mul(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Mul: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	const maxN = (1 << 32) - 1
	if !supportAvx512 || n >= maxN {
		// call mulVecGeneric
		mulVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 16
	mulVec(&(*vector)[0], &a[0], &b[0], n/blockSize, qInvNeg)
	if n % blockSize != 0 {
		// call mulVecGeneric on the rest
		start := n - n % blockSize
		mulVecGeneric((*vector)[start:], a[start:], b[start:])
	}
	
}

// Patterns use for transposing the vectors in mulVec
var (
	pattern1 = [8]uint64{0, 8, 1, 9, 2, 10, 3, 11}
	pattern2 = [8]uint64{12, 4, 13, 5, 14, 6, 15, 7}
	pattern3 = [8]uint64{0, 1, 8, 9, 2, 3, 10, 11}
	pattern4 = [8]uint64{12, 13, 4, 5, 14, 15, 6, 7}
)

//go:noescape
func mulVec(res, a, b *{{.ElementName}}, n uint64, qInvNeg uint64)


// ButterflyMul used in FFT; TODO complete.
func (vector *Vector) ButterflyMul(twiddles Vector, start, end, m int) {
	butterflyMulVecGeneric(*vector, twiddles, start, end, m)
}

`

const VectorOpsArm64 = VectorOpsPureGo

const VectorOpsArm64F31 = `

//go:noescape
func addVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func subVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func sumVec(t *uint64, a *{{.ElementName}}, n uint64)

// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Add: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}

	const blockSize = 4
	addVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n % blockSize != 0 {
		// call addVecGeneric on the rest
		start := n - n % blockSize
		addVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Sub: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}

	const blockSize = 4
	subVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n % blockSize != 0 {
		// call subVecGeneric on the rest
		start := n - n % blockSize
		subVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res {{.ElementName}}) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}

	const blockSize = 16
	var t [4]uint64 // stores the accumulators (not reduced mod q)
	sumVec(&t[0], &(*vector)[0], n/blockSize)
	// we reduce the accumulators mod q and add to res
	var v {{.ElementName}}
	for i := 0; i < 4; i++ {
		v[0] = uint32(t[i] % q)
		res.Add(&res, &v)
	}
	if n % blockSize != 0 {
		// call sumVecGeneric on the rest
		start := n - n % blockSize
		sumVecGeneric(&res, (*vector)[start:])
	}

	return
}

// note: unfortunately, as of Dec. 2024, Golang doesn't support enough NEON instructions
// for these to be worth it in assembly. Will hopefully revisit in future versions.

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *{{.ElementName}}) {
	scalarMulVecGeneric(*vector, a, b)
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

const VectorOpsAmd64F31 = `

//go:noescape
func addVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func subVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func sumVec(t *uint64, a *{{.ElementName}}, n uint64)

//go:noescape
func mulVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func scalarMulVec(res, a, b *{{.ElementName}}, n uint64)

//go:noescape
func innerProdVec(t *uint64, a, b *{{.ElementName}}, n uint64)

//go:noescape
func butterflyMulVec(a, twiddles *{{.ElementName}}, m int)

// Add adds two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Add(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Add: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call addVecGeneric
		addVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 16
	addVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n % blockSize != 0 {
		// call addVecGeneric on the rest
		start := n - n % blockSize
		addVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// Sub subtracts two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Sub(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Sub: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call subVecGeneric
		subVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 16
	subVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n % blockSize != 0 {
		// call subVecGeneric on the rest
		start := n - n % blockSize
		subVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}

// ScalarMul multiplies a vector by a scalar element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) ScalarMul(a Vector, b *{{.ElementName}}) {
	if len(a) != len(*vector) {
		panic("vector.ScalarMul: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call scalarMulVecGeneric
		scalarMulVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 8
	scalarMulVec(&(*vector)[0], &a[0], b, n/blockSize)
	if n % blockSize != 0 {
		// call scalarMulVecGeneric on the rest
		start := n - n % blockSize
		scalarMulVecGeneric((*vector)[start:], a[start:], b)
	}
}

// Sum computes the sum of all elements in the vector.
func (vector *Vector) Sum() (res {{.ElementName}}) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call sumVecGeneric
		sumVecGeneric(&res, *vector)
		return
	}

	const blockSize = 16
	var t [8]uint64 // stores the accumulators (not reduced mod q)
	sumVec(&t[0], &(*vector)[0], n/blockSize)
	// we reduce the accumulators mod q and add to res
	var v {{.ElementName}}
	for i := 0; i < 8; i++ {
		v[0] = uint32(t[i] % q)
		res.Add(&res, &v)
	}
	if n % blockSize != 0 {
		// call sumVecGeneric on the rest
		start := n - n % blockSize
		sumVecGeneric(&res, (*vector)[start:])
	}

	return
}

// InnerProduct computes the inner product of two vectors.
// It panics if the vectors don't have the same length.
func (vector *Vector) InnerProduct(other Vector) (res {{.ElementName}}) {
	n := uint64(len(*vector))
	if n == 0 {
		return
	}
	if n != uint64(len(other)) {
		panic("vector.InnerProduct: vectors don't have the same length")
	}
	if !supportAvx512 {
		// call innerProductVecGeneric
		innerProductVecGeneric(&res, *vector, other)
		return
	}

	const blockSize = 8
	var t [8]uint64 // stores the accumulators (not reduced mod q)
	innerProdVec(&t[0], &(*vector)[0], &other[0], n/blockSize)
	// we reduce the accumulators mod q and add to res
	var v {{.ElementName}}
	for i := 0; i < 8; i++ {
		v[0] = uint32(t[i] % q)
		res.Add(&res, &v)
	}
	if n % blockSize != 0 {
		// call innerProductVecGeneric on the rest
		start := n - n % blockSize
		innerProductVecGeneric(&res, (*vector)[start:], other[start:])
	}

	return
}

// Mul multiplies two vectors element-wise and stores the result in self.
// It panics if the vectors don't have the same length.
func (vector *Vector) Mul(a, b Vector) {
	if len(a) != len(b) || len(a) != len(*vector) {
		panic("vector.Mul: vectors don't have the same length")
	}
	n := uint64(len(a))
	if n == 0 {
		return
	}
	if !supportAvx512 {
		// call mulVecGeneric
		mulVecGeneric(*vector, a, b)
		return
	}

	const blockSize = 8
	mulVec(&(*vector)[0], &a[0], &b[0], n/blockSize)
	if n % blockSize != 0 {
		// call mulVecGeneric on the rest
		start := n - n % blockSize
		mulVecGeneric((*vector)[start:], a[start:], b[start:])
	}
}


// ButterflyMul used in FFT; TODO complete.
func (vector *Vector) ButterflyMul(twiddles Vector, start, end, m int) {
	n := end - start
	if !supportAvx512 || start != 0 ||m % 8 != 0 || n !=m{
		// call butterflyMulVec
		butterflyMulVecGeneric(*vector, twiddles, start, end, m)
		return
	}

	butterflyMulVec(&(*vector)[0], &twiddles[0], m)
}


`
