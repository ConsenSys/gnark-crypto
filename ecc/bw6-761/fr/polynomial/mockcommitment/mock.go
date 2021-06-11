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

package mockcommitment

import (
	"io"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bw6-761/fr"
	bw6761 "github.com/consensys/gnark-crypto/ecc/bw6-761/fr/polynomial"
	"github.com/consensys/gnark-crypto/polynomial"
)

// MockProof empty struct
type MockProof struct {
	Point        fr.Element
	ClaimedValue fr.Element
}

func (mp *MockProof) Marshal() []byte {
	panic("not implemented")
}

// GetClaimedValue returns the serialized claimed value.
func (mp *MockProof) GetClaimedValue() []byte {
	return mp.ClaimedValue.Marshal()
}

// MockBatchProofsSinglePoint empty struct
type MockBatchProofsSinglePoint struct {
	Point         fr.Element
	ClaimedValues []fr.Element
}

func (mp *MockBatchProofsSinglePoint) Marshal() []byte {
	panic("not implemented")
}

// GetClaimedValues returns a slice of the claimed values,
// serialized.
func (mp *MockBatchProofsSinglePoint) GetClaimedValues() [][]byte {
	res := make([][]byte, len(mp.ClaimedValues))
	for i := 0; i < len(mp.ClaimedValues); i++ {
		res[i] = mp.ClaimedValues[i].Marshal()
	}
	return res
}

// Scheme mock commitment, useful for testing polynomial based IOP
// like PLONK, where the scheme should not depend on which polynomial commitment scheme
// is used.
type Scheme struct{}

// Digest commitment of a polynomials
type Digest struct {
	data fr.Element
}

// Marshal serializes the point as in bw6761.G1Affine.
func (d *Digest) Marshal() []byte {
	return d.data.Marshal()
}

// Marshal serializes the digest.
func (d *Digest) Unmarshal(buf []byte) error {
	d.data.SetBytes(buf)
	return nil
}

// Add adds two digest.
func (d *Digest) Add(d1, d2 polynomial.Digest) polynomial.Digest {
	_d1 := d1.(*Digest)
	_d2 := d2.(*Digest)
	d.data.Add(&_d1.data, &_d2.data)
	return d
}

// Sub adds two digest.
func (d *Digest) Sub(d1, d2 polynomial.Digest) polynomial.Digest {
	_d1 := d1.(*Digest)
	_d2 := d2.(*Digest)
	d.data.Sub(&_d1.data, &_d2.data)
	return d
}

// Add adds two digest.
func (d *Digest) ScalarMul(d1 polynomial.Digest, s big.Int) polynomial.Digest {
	_d1 := d1.(*Digest)
	var _s fr.Element
	_s.SetBigInt(&s)
	d.data.Mul(&_s, &_d1.data)
	return d
}

// WriteTo panics
func (s *Scheme) WriteTo(w io.Writer) (n int64, err error) {
	panic("not implemented")
}

// ReadFrom panics
func (s *Scheme) ReadFrom(r io.Reader) (n int64, err error) {
	panic("not implemented")
}

// Commit returns the first coefficient of p
func (s *Scheme) Commit(p polynomial.Polynomial) (polynomial.Digest, error) {
	_p := p.(*bw6761.Polynomial)
	var res Digest
	res.data.SetInterface((*_p)[0])
	return &res, nil
}

// Open computes an opening proof of _p at _val.
// Returns a MockProof, which is an empty interface.
func (s *Scheme) Open(point interface{}, p polynomial.Polynomial) (polynomial.OpeningProof, error) {

	res := MockProof{}
	res.Point.SetInterface(point)
	res.ClaimedValue.SetInterface(p.Eval(point))

	return &res, nil
}

// Verify mock implementation of verify
func (s *Scheme) Verify(commitment polynomial.Digest, proof polynomial.OpeningProof) error {
	return nil
}

// BatchOpenSinglePoint computes a batch opening proof for _p at _val.
func (s *Scheme) BatchOpenSinglePoint(point interface{}, digests []polynomial.Digest, polynomials []polynomial.Polynomial) (polynomial.BatchOpeningProofSinglePoint, error) {

	var res MockBatchProofsSinglePoint
	res.ClaimedValues = make([]fr.Element, len(polynomials))
	res.Point.SetInterface(point)

	for i := 0; i < len(polynomials); i++ {
		res.ClaimedValues[i].SetInterface(polynomials[i].Eval(point))
	}

	return &res, nil
}

// BatchVerifySinglePoint computes a batch opening proof for
func (s *Scheme) BatchVerifySinglePoint(digests []polynomial.Digest, batchOpeningProof polynomial.BatchOpeningProofSinglePoint) error {

	return nil

}
