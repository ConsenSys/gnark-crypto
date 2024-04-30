// Copyright 2020 Consensys Software Inc.
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

package kzg

import (
	"bytes"
	"encoding/binary"
	"github.com/consensys/gnark-crypto/ecc/bls12-378"
	"github.com/consensys/gnark-crypto/ecc/bls12-378/fp"
	"io"
)

// WriteTo writes binary encoding of the ProvingKey
func (pk *ProvingKey) WriteTo(w io.Writer) (int64, error) {
	return pk.writeTo(w)
}

// WriteRawTo writes binary encoding of ProvingKey to w without point compression
func (pk *ProvingKey) WriteRawTo(w io.Writer) (int64, error) {
	return pk.writeTo(w, bls12378.RawEncoding())
}

func (pk *ProvingKey) writeTo(w io.Writer, options ...func(*bls12378.Encoder)) (int64, error) {
	// encode the ProvingKey
	enc := bls12378.NewEncoder(w, options...)
	if err := enc.Encode(pk.G1); err != nil {
		return enc.BytesWritten(), err
	}
	return enc.BytesWritten(), nil
}

// WriteRawTo writes binary encoding of VerifyingKey to w without point compression
func (vk *VerifyingKey) WriteRawTo(w io.Writer) (int64, error) {
	return vk.writeTo(w, bls12378.RawEncoding())
}

// WriteTo writes binary encoding of the VerifyingKey
func (vk *VerifyingKey) WriteTo(w io.Writer) (int64, error) {
	return vk.writeTo(w)
}

func (vk *VerifyingKey) writeTo(w io.Writer, options ...func(*bls12378.Encoder)) (int64, error) {
	// encode the VerifyingKey
	enc := bls12378.NewEncoder(w, options...)
	nLines := 63
	toEncode := make([]interface{}, 0, 4*nLines+3)
	toEncode = append(toEncode, &vk.G2[0])
	toEncode = append(toEncode, &vk.G2[1])
	toEncode = append(toEncode, &vk.G1)
	for k := 0; k < 2; k++ {
		for j := 0; j < 2; j++ {
			for i := nLines - 1; i >= 0; i-- {
				toEncode = append(toEncode, &vk.Lines[k][j][i].R0)
				toEncode = append(toEncode, &vk.Lines[k][j][i].R1)
			}
		}
	}

	for _, v := range toEncode {
		if err := enc.Encode(v); err != nil {
			return enc.BytesWritten(), err
		}
	}

	return enc.BytesWritten(), nil
}

// UnsafeToBytes returns the binary encoding of the entire SRS memory representation
// It is meant to be use to achieve fast serialization/deserialization and
// is not compatible with WriteTo / ReadFrom. It does not do any validation
// and doesn't encode points in a canonical form.
// @unstable: the format may change in the future
// If maxPkPoints is provided, the number of points in the ProvingKey will be limited to maxPkPoints
func (srs *SRS) UnsafeToBytes(maxPkPoints ...int) ([]byte, error) {
	maxG1 := len(srs.Pk.G1)
	if len(maxPkPoints) > 0 && maxPkPoints[0] < maxG1 && maxPkPoints[0] > 0 {
		maxG1 = maxPkPoints[0]
	}
	// first we write the VerifyingKey; it is small so we re-use WriteTo
	var buf bytes.Buffer

	if _, err := srs.Vk.writeTo(&buf, bls12378.RawEncoding()); err != nil {
		return nil, err
	}

	buf.Grow(2*maxG1*fp.Bytes + 8) // pre-allocate space for the ProvingKey

	// write nb points we encode.
	if err := binary.Write(&buf, binary.LittleEndian, uint64(maxG1)); err != nil {
		return nil, err
	}

	// write the limbs directly
	var bbuf [fp.Bytes * 2]byte
	for i := 0; i < maxG1; i++ {
		for j := 0; j < fp.Limbs; j++ {
			binary.LittleEndian.PutUint64(bbuf[j*8:j*8+8], srs.Pk.G1[i].X[j])
		}
		for j := 0; j < fp.Limbs; j++ {
			binary.LittleEndian.PutUint64(bbuf[fp.Bytes+j*8:fp.Bytes+j*8+8], srs.Pk.G1[i].Y[j])
		}
		if _, err := buf.Write(bbuf[:]); err != nil {
			return nil, err
		}
	}
	return buf.Bytes(), nil
}

// UnsafeFromBytes deserializes the SRS from a byte slice
// It is meant to be use to achieve fast serialization/deserialization and
// is not compatible with WriteTo / ReadFrom. It does not do any validation
// and doesn't encode points in a canonical form.
// @unstable: the format may change in the future
func (srs *SRS) UnsafeFromBytes(data []byte, maxPkPoints ...int) error {
	buf := bytes.NewReader(data)

	// first we read the VerifyingKey; it is small so we re-use ReadFrom
	if _, err := srs.Vk.ReadFrom(buf); err != nil {
		return err
	}

	// read nb points we encode.
	var nbPoints uint64
	if err := binary.Read(buf, binary.LittleEndian, &nbPoints); err != nil {
		return err
	}

	if len(maxPkPoints) == 1 && maxPkPoints[0] > 0 && int(nbPoints) > maxPkPoints[0] {
		nbPoints = uint64(maxPkPoints[0])
	}

	srs.Pk.G1 = make([]bls12378.G1Affine, nbPoints)

	// read the limbs directly
	var bbuf [fp.Bytes * 2]byte
	for i := 0; i < int(nbPoints); i++ {
		if _, err := io.ReadFull(buf, bbuf[:]); err != nil {
			return err
		}
		for j := 0; j < fp.Limbs; j++ {
			srs.Pk.G1[i].X[j] = binary.LittleEndian.Uint64(bbuf[j*8 : j*8+8])
		}
		for j := 0; j < fp.Limbs; j++ {
			srs.Pk.G1[i].Y[j] = binary.LittleEndian.Uint64(bbuf[fp.Bytes+j*8 : fp.Bytes+j*8+8])
		}
	}
	return nil
}

// WriteTo writes binary encoding of the entire SRS
func (srs *SRS) WriteTo(w io.Writer) (int64, error) {
	// encode the SRS
	var pn, vn int64
	var err error
	if pn, err = srs.Pk.WriteTo(w); err != nil {
		return pn, err
	}
	vn, err = srs.Vk.WriteTo(w)
	return pn + vn, err
}

// WriteRawTo writes binary encoding of the entire SRS without point compression
func (srs *SRS) WriteRawTo(w io.Writer) (int64, error) {
	// encode the SRS
	var pn, vn int64
	var err error
	if pn, err = srs.Pk.WriteRawTo(w); err != nil {
		return pn, err
	}
	vn, err = srs.Vk.WriteRawTo(w)
	return pn + vn, err
}

// ReadFrom decodes ProvingKey data from reader.
func (pk *ProvingKey) ReadFrom(r io.Reader) (int64, error) {
	// decode the ProvingKey
	dec := bls12378.NewDecoder(r)
	if err := dec.Decode(&pk.G1); err != nil {
		return dec.BytesRead(), err
	}
	return dec.BytesRead(), nil
}

// UnsafeReadFrom decodes ProvingKey data from reader without checking
// that point are in the correct subgroup.
func (pk *ProvingKey) UnsafeReadFrom(r io.Reader) (int64, error) {
	// decode the ProvingKey
	dec := bls12378.NewDecoder(r, bls12378.NoSubgroupChecks())
	if err := dec.Decode(&pk.G1); err != nil {
		return dec.BytesRead(), err
	}
	return dec.BytesRead(), nil
}

// ReadFrom decodes VerifyingKey data from reader.
func (vk *VerifyingKey) ReadFrom(r io.Reader) (int64, error) {
	// decode the VerifyingKey
	dec := bls12378.NewDecoder(r)
	nLines := 63
	toDecode := make([]interface{}, 0, 4*nLines+3)
	toDecode = append(toDecode, &vk.G2[0])
	toDecode = append(toDecode, &vk.G2[1])
	toDecode = append(toDecode, &vk.G1)
	for k := 0; k < 2; k++ {
		for j := 0; j < 2; j++ {
			for i := nLines - 1; i >= 0; i-- {
				toDecode = append(toDecode, &vk.Lines[k][j][i].R0)
				toDecode = append(toDecode, &vk.Lines[k][j][i].R1)
			}
		}
	}

	for _, v := range toDecode {
		if err := dec.Decode(v); err != nil {
			return dec.BytesRead(), err
		}
	}

	return dec.BytesRead(), nil
}

// ReadFrom decodes SRS data from reader.
func (srs *SRS) ReadFrom(r io.Reader) (int64, error) {
	// decode the VerifyingKey
	var pn, vn int64
	var err error
	if pn, err = srs.Pk.ReadFrom(r); err != nil {
		return pn, err
	}
	vn, err = srs.Vk.ReadFrom(r)
	return pn + vn, err
}

// UnsafeReadFrom decodes SRS data from reader without sub group checks
func (srs *SRS) UnsafeReadFrom(r io.Reader) (int64, error) {
	// decode the VerifyingKey
	var pn, vn int64
	var err error
	if pn, err = srs.Pk.UnsafeReadFrom(r); err != nil {
		return pn, err
	}
	vn, err = srs.Vk.ReadFrom(r)
	return pn + vn, err
}

// WriteTo writes binary encoding of a OpeningProof
func (proof *OpeningProof) WriteTo(w io.Writer) (int64, error) {
	enc := bls12378.NewEncoder(w)

	toEncode := []interface{}{
		&proof.H,
		&proof.ClaimedValue,
	}

	for _, v := range toEncode {
		if err := enc.Encode(v); err != nil {
			return enc.BytesWritten(), err
		}
	}

	return enc.BytesWritten(), nil
}

// ReadFrom decodes OpeningProof data from reader.
func (proof *OpeningProof) ReadFrom(r io.Reader) (int64, error) {
	dec := bls12378.NewDecoder(r)

	toDecode := []interface{}{
		&proof.H,
		&proof.ClaimedValue,
	}

	for _, v := range toDecode {
		if err := dec.Decode(v); err != nil {
			return dec.BytesRead(), err
		}
	}

	return dec.BytesRead(), nil
}

// WriteTo writes binary encoding of a BatchOpeningProof
func (proof *BatchOpeningProof) WriteTo(w io.Writer) (int64, error) {
	enc := bls12378.NewEncoder(w)

	toEncode := []interface{}{
		&proof.H,
		proof.ClaimedValues,
	}

	for _, v := range toEncode {
		if err := enc.Encode(v); err != nil {
			return enc.BytesWritten(), err
		}
	}

	return enc.BytesWritten(), nil
}

// ReadFrom decodes BatchOpeningProof data from reader.
func (proof *BatchOpeningProof) ReadFrom(r io.Reader) (int64, error) {
	dec := bls12378.NewDecoder(r)
	toDecode := []interface{}{
		&proof.H,
		&proof.ClaimedValues,
	}

	for _, v := range toDecode {
		if err := dec.Decode(v); err != nil {
			return dec.BytesRead(), err
		}
	}

	return dec.BytesRead(), nil
}
