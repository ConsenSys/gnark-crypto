// Copyright 2020 ConsenSys Software Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fp

// /!\ WARNING /!\
// this code has not been audited and is provided as-is. In particular,
// there is no security guarantees such as constant time implementation
// or side-channel attack resistance
// /!\ WARNING /!\

import (
	"crypto/rand"
	"encoding/binary"
	"errors"
	"io"
	"math/big"
	"math/bits"
	"reflect"
	"strconv"
	"strings"
	"sync"
)

// Element represents a field element stored on 1 words (uint64)
// Element are assumed to be in Montgomery form in all methods
// field modulus q =
//
// 326667333367
type Element [1]uint64

// Limbs number of 64 bits words needed to represent Element
const Limbs = 1

// Bits number bits needed to represent Element
const Bits = 39

// Bytes number bytes needed to represent Element
const Bytes = Limbs * 8

// field modulus stored as big.Int
var _modulus big.Int

// Modulus returns q as a big.Int
// q =
//
// 326667333367
func Modulus() *big.Int {
	return new(big.Int).Set(&_modulus)
}

// q (modulus)
const qElementWord0 uint64 = 326667333367
const q uint64 = qElementWord0

var qElement = Element{
	qElementWord0,
}

// rSquare
var rSquare = Element{
	242079875569,
}

var bigIntPool = sync.Pool{
	New: func() interface{} {
		return new(big.Int)
	},
}

func init() {
	// base10: 326667333367
	_modulus.SetString("4c0ee3eef7", 16)
}

// NewElement returns a new Element from a uint64 value
//
// it is equivalent to
// 		var v NewElement
// 		v.SetUint64(...)
func NewElement(v uint64) Element {
	z := Element{v}
	z.Mul(&z, &rSquare)
	return z
}

// SetUint64 sets z to v and returns z
func (z *Element) SetUint64(v uint64) *Element {
	//  sets z LSB to v (non-Montgomery form) and convert z to Montgomery form
	*z = Element{v}
	return z.Mul(z, &rSquare) // z.ToMont()
}

// SetInt64 sets z to v and returns z
func (z *Element) SetInt64(v int64) *Element {

	// absolute value of v
	m := v >> 63
	z.SetUint64(uint64((v ^ m) - m))

	if m != 0 {
		// v is negative
		z.Neg(z)
	}

	return z
}

// Set z = x
func (z *Element) Set(x *Element) *Element {
	z[0] = x[0]
	return z
}

// SetInterface converts provided interface into Element
// returns an error if provided type is not supported
// supported types: Element, *Element, uint64, int, string (interpreted as base10 integer),
// *big.Int, big.Int, []byte
func (z *Element) SetInterface(i1 interface{}) (*Element, error) {
	switch c1 := i1.(type) {
	case Element:
		return z.Set(&c1), nil
	case *Element:
		if c1 == nil {
			return nil, errors.New("can't set fp.Element with <nil>")
		}
		return z.Set(c1), nil
	case uint8:
		return z.SetUint64(uint64(c1)), nil
	case uint16:
		return z.SetUint64(uint64(c1)), nil
	case uint32:
		return z.SetUint64(uint64(c1)), nil
	case uint:
		return z.SetUint64(uint64(c1)), nil
	case uint64:
		return z.SetUint64(c1), nil
	case int8:
		return z.SetInt64(int64(c1)), nil
	case int16:
		return z.SetInt64(int64(c1)), nil
	case int32:
		return z.SetInt64(int64(c1)), nil
	case int64:
		return z.SetInt64(c1), nil
	case int:
		return z.SetInt64(int64(c1)), nil
	case string:
		return z.SetString(c1), nil
	case *big.Int:
		if c1 == nil {
			return nil, errors.New("can't set fp.Element with <nil>")
		}
		return z.SetBigInt(c1), nil
	case big.Int:
		return z.SetBigInt(&c1), nil
	case []byte:
		return z.SetBytes(c1), nil
	default:
		return nil, errors.New("can't set fp.Element from type " + reflect.TypeOf(i1).String())
	}
}

// SetZero z = 0
func (z *Element) SetZero() *Element {
	z[0] = 0
	return z
}

// SetOne z = 1 (in Montgomery form)
func (z *Element) SetOne() *Element {
	z[0] = 152135744813
	return z
}

// Div z = x*y^-1 mod q
func (z *Element) Div(x, y *Element) *Element {
	var yInv Element
	yInv.Inverse(y)
	z.Mul(x, &yInv)
	return z
}

// Bit returns the i'th bit, with lsb == bit 0.
// It is the responsibility of the caller to convert from Montgomery to Regular form if needed
func (z *Element) Bit(i uint64) uint64 {
	j := i / 64
	if j >= 1 {
		return 0
	}
	return uint64(z[j] >> (i % 64) & 1)
}

// Equal returns z == x; constant-time
func (z *Element) Equal(x *Element) bool {
	return z.NotEqual(x) == 0
}

// NotEqual returns 0 if and only if z == x; constant-time
func (z *Element) NotEqual(x *Element) uint64 {
	return (z[0] ^ x[0])
}

// IsZero returns z == 0
func (z *Element) IsZero() bool {
	return (z[0]) == 0
}

// IsOne returns z == 1
func (z *Element) IsOne() bool {
	return z[0] == 152135744813
}

// IsUint64 reports whether z can be represented as an uint64.
func (z *Element) IsUint64() bool {
	return true
}

// Uint64 returns the uint64 representation of x. If x cannot be represented in a uint64, the result is undefined.
func (z *Element) Uint64() uint64 {
	zz := *z
	zz.FromMont()
	return zz[0]
}

// FitsOnOneWord reports whether z words (except the least significant word) are 0
func (z *Element) FitsOnOneWord() bool {
	return true
}

// Cmp compares (lexicographic order) z and x and returns:
//
//   -1 if z <  x
//    0 if z == x
//   +1 if z >  x
//
func (z *Element) Cmp(x *Element) int {
	_z := *z
	_x := *x
	_z.FromMont()
	_x.FromMont()
	if _z[0] > _x[0] {
		return 1
	} else if _z[0] < _x[0] {
		return -1
	}
	return 0
}

// LexicographicallyLargest returns true if this element is strictly lexicographically
// larger than its negation, false otherwise
func (z *Element) LexicographicallyLargest() bool {
	// adapted from github.com/zkcrypto/bls12_381
	// we check if the element is larger than (q-1) / 2
	// if z - (((q -1) / 2) + 1) have no underflow, then z > (q-1) / 2

	_z := *z
	_z.FromMont()

	var b uint64
	_, b = bits.Sub64(_z[0], 163333666684, 0)

	return b == 0
}

// SetRandom sets z to a random element < q
func (z *Element) SetRandom() (*Element, error) {
	var bytes [8]byte
	if _, err := io.ReadFull(rand.Reader, bytes[:]); err != nil {
		return nil, err
	}
	z[0] = binary.BigEndian.Uint64(bytes[0:8])
	z[0] %= 326667333367

	return z, nil
}

// One returns 1 (in montgommery form)
func One() Element {
	var one Element
	one.SetOne()
	return one
}

// Halve sets z to z / 2 (mod p)
func (z *Element) Halve() {

	if z[0]&1 == 1 {
		// z = z + q
		z[0], _ = bits.Add64(z[0], 326667333367, 0)

	}
	// z = z >> 1
	z[0] >>= 1

}

// Mul z = x * y mod q
// see https://hackmd.io/@gnark/modular_multiplication
func (z *Element) Mul(x, y *Element) *Element {

	// CIOS multiplication
	// Used for Montgomery reduction. (qInvNeg) q + r'.r = 1, i.e., qInvNeg = - q⁻¹ mod r
	const qInvNegLsw uint64 = 14763286405836492089

	var r uint64
	hi, lo := bits.Mul64(x[0], y[0])
	m := lo * qInvNegLsw
	hi2, lo2 := bits.Mul64(m, q)
	_, carry := bits.Add64(lo2, lo, 0)
	r, carry = bits.Add64(hi2, hi, carry)

	if carry != 0 || r >= q {
		// we need to reduce
		r -= q

	}
	z[0] = r

	return z
}

// Square z = x * x mod q
// see https://hackmd.io/@gnark/modular_multiplication
func (z *Element) Square(x *Element) *Element {

	// CIOS multiplication
	// Used for Montgomery reduction. (qInvNeg) q + r'.r = 1, i.e., qInvNeg = - q⁻¹ mod r
	const qInvNegLsw uint64 = 14763286405836492089

	var r uint64
	hi, lo := bits.Mul64(x[0], x[0])
	m := lo * qInvNegLsw
	hi2, lo2 := bits.Mul64(m, q)
	_, carry := bits.Add64(lo2, lo, 0)
	r, carry = bits.Add64(hi2, hi, carry)

	if carry != 0 || r >= q {
		// we need to reduce
		r -= q

	}
	z[0] = r

	return z
}

// FromMont converts z in place (i.e. mutates) from Montgomery to regular representation
// sets and returns z = z * 1
func (z *Element) FromMont() *Element {
	fromMont(z)
	return z
}

// Add z = x + y mod q
func (z *Element) Add(x, y *Element) *Element {
	add(z, x, y)
	return z
}

// Double z = x + x mod q, aka Lsh 1
func (z *Element) Double(x *Element) *Element {
	double(z, x)
	return z
}

// Sub  z = x - y mod q
func (z *Element) Sub(x, y *Element) *Element {
	sub(z, x, y)
	return z
}

// Neg z = q - x
func (z *Element) Neg(x *Element) *Element {
	neg(z, x)
	return z
}

// Select is a constant-time conditional move.
// If c=0, z = x0. Else z = x1
func (z *Element) Select(c int, x0 *Element, x1 *Element) *Element {
	cC := uint64((int64(c) | -int64(c)) >> 63) // "canonicized" into: 0 if c=0, -1 otherwise
	z[0] = x0[0] ^ cC&(x0[0]^x1[0])
	return z
}

// Generic (no ADX instructions, no AMD64) versions of multiplication and squaring algorithms

func _mulGeneric(z, x, y *Element) {

	// CIOS multiplication
	// Used for Montgomery reduction. (qInvNeg) q + r'.r = 1, i.e., qInvNeg = - q⁻¹ mod r
	const qInvNegLsw uint64 = 14763286405836492089

	var r uint64
	hi, lo := bits.Mul64(x[0], y[0])
	m := lo * qInvNegLsw
	hi2, lo2 := bits.Mul64(m, q)
	_, carry := bits.Add64(lo2, lo, 0)
	r, carry = bits.Add64(hi2, hi, carry)

	if carry != 0 || r >= q {
		// we need to reduce
		r -= q

	}
	z[0] = r

}

func _fromMontGeneric(z *Element) {
	// the following lines implement z = z * 1
	// with a modified CIOS montgomery multiplication
	{
		// m = z[0]n'[0] mod W
		m := z[0] * 14763286405836492089
		C := madd0(m, 326667333367, z[0])
		z[0] = C
	}

	// if z >= q → z -= q
	// note: this is NOT constant time
	if z[0] >= q {
		z[0] -= q
	}
}

func _addGeneric(z, x, y *Element) {

	z[0], _ = bits.Add64(x[0], y[0], 0)
	if z[0] >= q {
		z[0] -= q
	}
}

func _doubleGeneric(z, x *Element) {
	if x[0]&(1<<63) == (1 << 63) {
		// if highest bit is set, then we have a carry to x + x, we shift and subtract q
		z[0] = (x[0] << 1) - q
	} else {
		// highest bit is not set, but x + x can still be >= q
		z[0] = (x[0] << 1)
		if z[0] >= q {
			z[0] -= q
		}
	}
}

func _subGeneric(z, x, y *Element) {
	var b uint64
	z[0], b = bits.Sub64(x[0], y[0], 0)
	if b != 0 {
		z[0] += q
	}
}

func _negGeneric(z, x *Element) {
	if x.IsZero() {
		z.SetZero()
		return
	}
	z[0] = q - x[0]
}

func _reduceGeneric(z *Element) {

	// if z >= q → z -= q
	// note: this is NOT constant time
	if z[0] >= q {
		z[0] -= q
	}
}

func mulByConstant(z *Element, c uint8) {
	switch c {
	case 0:
		z.SetZero()
		return
	case 1:
		return
	case 2:
		z.Double(z)
		return
	case 3:
		_z := *z
		z.Double(z).Add(z, &_z)
	case 5:
		_z := *z
		z.Double(z).Double(z).Add(z, &_z)
	case 11:
		_z := *z
		z.Double(z).Double(z).Add(z, &_z).Double(z).Add(z, &_z)
	default:
		var y Element
		y.SetUint64(uint64(c))
		z.Mul(z, &y)
	}
}

// BatchInvert returns a new slice with every element inverted.
// Uses Montgomery batch inversion trick
func BatchInvert(a []Element) []Element {
	res := make([]Element, len(a))
	if len(a) == 0 {
		return res
	}

	zeroes := make([]bool, len(a))
	accumulator := One()

	for i := 0; i < len(a); i++ {
		if a[i].IsZero() {
			zeroes[i] = true
			continue
		}
		res[i] = accumulator
		accumulator.Mul(&accumulator, &a[i])
	}

	accumulator.Inverse(&accumulator)

	for i := len(a) - 1; i >= 0; i-- {
		if zeroes[i] {
			continue
		}
		res[i].Mul(&res[i], &accumulator)
		accumulator.Mul(&accumulator, &a[i])
	}

	return res
}

func _butterflyGeneric(a, b *Element) {
	t := *a
	a.Add(a, b)
	b.Sub(&t, b)
}

// BitLen returns the minimum number of bits needed to represent z
// returns 0 if z == 0
func (z *Element) BitLen() int {
	return bits.Len64(z[0])
}

// Exp z = x^exponent mod q
func (z *Element) Exp(x Element, exponent *big.Int) *Element {
	var bZero big.Int
	if exponent.Cmp(&bZero) == 0 {
		return z.SetOne()
	}

	z.Set(&x)

	for i := exponent.BitLen() - 2; i >= 0; i-- {
		z.Square(z)
		if exponent.Bit(i) == 1 {
			z.Mul(z, &x)
		}
	}

	return z
}

// ToMont converts z to Montgomery form
// sets and returns z = z * r²
func (z *Element) ToMont() *Element {
	return z.Mul(z, &rSquare)
}

// ToRegular returns z in regular form (doesn't mutate z)
func (z Element) ToRegular() Element {
	return *z.FromMont()
}

// String returns the decimal representation of z as generated by
// z.Text(10).
func (z *Element) String() string {
	return z.Text(10)
}

// Text returns the string representation of z in the given base.
// Base must be between 2 and 36, inclusive. The result uses the
// lower-case letters 'a' to 'z' for digit values 10 to 35.
// No prefix (such as "0x") is added to the string. If z is a nil
// pointer it returns "<nil>".
// If base == 10 and -z fits in a uint16 prefix "-" is added to the string.
func (z *Element) Text(base int) string {
	if base < 2 || base > 36 {
		panic("invalid base")
	}
	if z == nil {
		return "<nil>"
	}

	const maxUint16 = 65535
	if base == 10 {
		var zzNeg Element
		zzNeg.Neg(z)
		zzNeg.FromMont()
		if zzNeg[0] <= maxUint16 && zzNeg[0] != 0 {
			return "-" + strconv.FormatUint(zzNeg[0], base)
		}
	}
	zz := *z
	zz.FromMont()
	return strconv.FormatUint(zz[0], base)
}

// ToBigInt returns z as a big.Int in Montgomery form
func (z *Element) ToBigInt(res *big.Int) *big.Int {
	var b [Limbs * 8]byte
	binary.BigEndian.PutUint64(b[0:8], z[0])

	return res.SetBytes(b[:])
}

// ToBigIntRegular returns z as a big.Int in regular form
func (z Element) ToBigIntRegular(res *big.Int) *big.Int {
	z.FromMont()
	return z.ToBigInt(res)
}

// Bytes returns the regular (non montgomery) value
// of z as a big-endian byte array.
func (z *Element) Bytes() (res [Limbs * 8]byte) {
	_z := z.ToRegular()
	binary.BigEndian.PutUint64(res[0:8], _z[0])

	return
}

// Marshal returns the regular (non montgomery) value
// of z as a big-endian byte slice.
func (z *Element) Marshal() []byte {
	b := z.Bytes()
	return b[:]
}

// SetBytes interprets e as the bytes of a big-endian unsigned integer,
// sets z to that value (in Montgomery form), and returns z.
func (z *Element) SetBytes(e []byte) *Element {
	if len(e) == 8 {
		// fast path
		z[0] = binary.BigEndian.Uint64(e)
		return z.ToMont()
	}
	// get a big int from our pool
	vv := bigIntPool.Get().(*big.Int)
	vv.SetBytes(e)

	// set big int
	z.SetBigInt(vv)

	// put temporary object back in pool
	bigIntPool.Put(vv)

	return z
}

// SetBigInt sets z to v (regular form) and returns z in Montgomery form
func (z *Element) SetBigInt(v *big.Int) *Element {
	z.SetZero()

	var zero big.Int

	// fast path
	c := v.Cmp(&_modulus)
	if c == 0 {
		// v == 0
		return z
	} else if c != 1 && v.Cmp(&zero) != -1 {
		// 0 < v < q
		return z.setBigInt(v)
	}

	// get temporary big int from the pool
	vv := bigIntPool.Get().(*big.Int)

	// copy input + modular reduction
	vv.Set(v)
	vv.Mod(v, &_modulus)

	// set big int byte value
	z.setBigInt(vv)

	// release object into pool
	bigIntPool.Put(vv)
	return z
}

// setBigInt assumes 0 ⩽ v < q
func (z *Element) setBigInt(v *big.Int) *Element {
	vBits := v.Bits()

	if bits.UintSize == 64 {
		for i := 0; i < len(vBits); i++ {
			z[i] = uint64(vBits[i])
		}
	} else {
		for i := 0; i < len(vBits); i++ {
			if i%2 == 0 {
				z[i/2] = uint64(vBits[i])
			} else {
				z[i/2] |= uint64(vBits[i]) << 32
			}
		}
	}

	return z.ToMont()
}

// SetString creates a big.Int with number and calls SetBigInt on z
//
// The number prefix determines the actual base: A prefix of
// ''0b'' or ''0B'' selects base 2, ''0'', ''0o'' or ''0O'' selects base 8,
// and ''0x'' or ''0X'' selects base 16. Otherwise, the selected base is 10
// and no prefix is accepted.
//
// For base 16, lower and upper case letters are considered the same:
// The letters 'a' to 'f' and 'A' to 'F' represent digit values 10 to 15.
//
// An underscore character ''_'' may appear between a base
// prefix and an adjacent digit, and between successive digits; such
// underscores do not change the value of the number.
// Incorrect placement of underscores is reported as a panic if there
// are no other errors.
//
func (z *Element) SetString(number string) *Element {
	// get temporary big int from the pool
	vv := bigIntPool.Get().(*big.Int)

	if _, ok := vv.SetString(number, 0); !ok {
		panic("Element.SetString failed -> can't parse number into a big.Int " + number)
	}

	z.SetBigInt(vv)

	// release object into pool
	bigIntPool.Put(vv)

	return z
}

// MarshalJSON returns json encoding of z (z.Text(10))
// If z == nil, returns null
func (z *Element) MarshalJSON() ([]byte, error) {
	if z == nil {
		return []byte("null"), nil
	}
	const maxSafeBound = 15 // we encode it as number if it's small
	s := z.Text(10)
	if len(s) <= maxSafeBound {
		return []byte(s), nil
	}
	var sbb strings.Builder
	sbb.WriteByte('"')
	sbb.WriteString(s)
	sbb.WriteByte('"')
	return []byte(sbb.String()), nil
}

// UnmarshalJSON accepts numbers and strings as input
// See Element.SetString for valid prefixes (0x, 0b, ...)
func (z *Element) UnmarshalJSON(data []byte) error {
	s := string(data)
	if len(s) > Bits*3 {
		return errors.New("value too large (max = Element.Bits * 3)")
	}

	// we accept numbers and strings, remove leading and trailing quotes if any
	if len(s) > 0 && s[0] == '"' {
		s = s[1:]
	}
	if len(s) > 0 && s[len(s)-1] == '"' {
		s = s[:len(s)-1]
	}

	// get temporary big int from the pool
	vv := bigIntPool.Get().(*big.Int)

	if _, ok := vv.SetString(s, 0); !ok {
		return errors.New("can't parse into a big.Int: " + s)
	}

	z.SetBigInt(vv)

	// release object into pool
	bigIntPool.Put(vv)
	return nil
}

// Legendre returns the Legendre symbol of z (either +1, -1, or 0.)
func (z *Element) Legendre() int {
	var l Element
	// z^((q-1)/2)
	l.expByLegendreExp(*z)

	if l.IsZero() {
		return 0
	}

	// if l == 1
	if l[0] == 152135744813 {
		return 1
	}
	return -1
}

// Sqrt z = √x mod q
// if the square root doesn't exist (x is not a square mod q)
// Sqrt leaves z unchanged and returns nil
func (z *Element) Sqrt(x *Element) *Element {
	// q ≡ 3 (mod 4)
	// using  z ≡ ± x^((p+1)/4) (mod q)
	var y, square Element
	y.expBySqrtExp(*x)
	// as we didn't compute the legendre symbol, ensure we found y such that y * y = x
	square.Square(&y)
	if square.Equal(x) {
		return z.Set(&y)
	}
	return nil
}

// Inverse z = x⁻¹ mod q
// Algorithm 16 in "Efficient Software-Implementation of Finite Fields with Applications to Cryptography"
// if x == 0, sets and returns z = x
func (z *Element) Inverse(x *Element) *Element {
	const q uint64 = qElementWord0
	if x.IsZero() {
		z.SetZero()
		return z
	}

	var r, s, u, v uint64
	u = q            // u = q
	s = 242079875569 // s = r^2
	r = 0
	v = x[0]

	var carry, borrow uint64

	for (u != 1) && (v != 1) {
		for v&1 == 0 {
			v >>= 1
			if s&1 == 0 {
				s >>= 1
			} else {
				s, carry = bits.Add64(s, q, 0)
				s >>= 1
				if carry != 0 {
					s |= (1 << 63)
				}
			}
		}
		for u&1 == 0 {
			u >>= 1
			if r&1 == 0 {
				r >>= 1
			} else {
				r, carry = bits.Add64(r, q, 0)
				r >>= 1
				if carry != 0 {
					r |= (1 << 63)
				}
			}
		}
		if v >= u {
			v -= u
			s, borrow = bits.Sub64(s, r, 0)
			if borrow == 1 {
				s += q
			}
		} else {
			u -= v
			r, borrow = bits.Sub64(r, s, 0)
			if borrow == 1 {
				r += q
			}
		}
	}

	if u == 1 {
		z[0] = r
	} else {
		z[0] = s
	}

	return z
}
