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

package bls12377

const MAX_BATCH_SIZE = 600

type batchOp struct {
	bucketID, pointID uint32
}

func (o batchOp) isNeg() bool {
	return o.pointID&1 == 1
}

// processChunkG1BatchAffine process a chunk of the scalars during the msm
// using affine coordinates for the buckets. To amortize the cost of the inverse in the affine addition
// we use a batch affine addition.
//
// this is derived from a PR by 0x0ece : https://github.com/ConsenSys/gnark-crypto/pull/249
// See Section 5.3: ia.cr/2022/1396
func processChunkG1BatchAffine[B ibG1Affine](chunk uint64,
	chRes chan<- g1JacExtended,
	c uint64,
	points []G1Affine,
	digits []uint32) {

	// init the buckets
	var buckets B
	for i := 0; i < len(buckets); i++ {
		buckets[i].setInfinity()
	}

	// setup for the batch affine;
	batchSize := len(buckets) / 5
	if batchSize > MAX_BATCH_SIZE {
		batchSize = MAX_BATCH_SIZE
	}
	if batchSize <= 0 {
		batchSize = 1
	}
	bucketIds := make(map[uint32]struct{}, len(buckets)/2) // TODO @gbotrel tune the capacity here
	cptP := 0                                              // count the number of point added to current batch

	var P [MAX_BATCH_SIZE]G1Affine  // allocated on the stack
	var R [MAX_BATCH_SIZE]*G1Affine // ...

	canAdd := func(bID uint32) bool {
		_, ok := bucketIds[bID]
		return !ok
	}

	isFull := func() bool {
		return cptP == batchSize
	}

	executeAndReset := func() {
		if cptP == 0 {
			return
		}
		BatchAddG1Affine(R[:cptP], P[:cptP], cptP)
		for k := range bucketIds {
			delete(bucketIds, k)
		}
		cptP = 0
	}

	add := func(op batchOp) {
		// CanAdd must be called before --> ensures bucket is not "used" in current batch

		BK := &buckets[op.bucketID]
		PP := &points[op.pointID>>1]
		if PP.IsInfinity() {
			return
		}
		// handle special cases with inf or -P / P
		if BK.IsInfinity() {
			if op.isNeg() {
				BK.Neg(PP)
			} else {
				BK.Set(PP)
			}
			return
		}
		if op.isNeg() {
			// if bucket == P --> -P == 0
			if BK.Equal(PP) {
				BK.setInfinity()
				return
			}
		} else {
			// if bucket == -P, B == 0
			if BK.X.Equal(&PP.X) && !BK.Y.Equal(&PP.Y) {
				BK.setInfinity()
				return
			}
		}

		// bucketIds[cptP] = op.bucketID
		bucketIds[op.bucketID] = struct{}{}
		R[cptP] = BK
		if op.isNeg() {
			P[cptP].Neg(PP)
		} else {
			P[cptP].Set(PP)
		}
		cptP++
	}

	queue := make([]batchOp, 0, 4096) // TODO find right capacity here.

	processQueue := func() {
		// for i := len(queue) - 1; i >= 0; i-- {
		for i := 0; i < len(queue); i++ {
			if canAdd(queue[i].bucketID) {
				add(queue[i])
				if isFull() {
					executeAndReset()
				}
				queue[i] = queue[len(queue)-1]
				queue = queue[:len(queue)-1]
				i--
			}
		}
	}

	nbBatches := 0
	for i, digit := range digits {

		if digit == 0 {
			continue
		}

		op := batchOp{pointID: uint32(i) << 1}
		// if msbWindow bit is set, we need to substract
		if digit&1 == 0 {
			// add
			op.bucketID = uint32((digit >> 1) - 1)
		} else {
			// sub
			op.bucketID = (uint32((digit >> 1)))
			op.pointID += 1
		}
		if canAdd(op.bucketID) {
			add(op)
			if isFull() {
				executeAndReset()
				nbBatches++
				if len(queue) != 0 { // TODO @gbotrel this doesn't seem to help much? should minimize queue resizing
					add(queue[len(queue)-1])
					queue = queue[:len(queue)-1]
				}
				// processQueue()
			}
		} else {
			// put it in queue.
			queue = append(queue, op)
		}
	}
	// fmt.Printf("chunk %d\nlen(queue)=%d\nnbBatches=%d\nbatchSize=%d\nnbBuckets=%d\nnbPoints=%d\n",
	// 	chunk, len(queue), nbBatches, batch.batchSize, len(buckets), len(points))
	// executeAndReset()
	for len(queue) != 0 {
		processQueue()
		executeAndReset() // execute batch even if not full.
	}

	// flush items in batch.
	executeAndReset()

	// reduce buckets into total
	// total =  bucket[0] + 2*bucket[1] + 3*bucket[2] ... + n*bucket[n-1]

	var runningSum, total g1JacExtended
	runningSum.setInfinity()
	total.setInfinity()
	for k := len(buckets) - 1; k >= 0; k-- {
		if !buckets[k].IsInfinity() {
			runningSum.addMixed(&buckets[k])
		}
		total.add(&runningSum)
	}

	chRes <- total

}

// we declare the buckets as fixed-size array types
// this allow us to allocate the buckets on the stack
type bucketG1AffineC4 [1 << (4 - 1)]G1Affine
type bucketG1AffineC5 [1 << (5 - 1)]G1Affine
type bucketG1AffineC6 [1 << (6 - 1)]G1Affine
type bucketG1AffineC7 [1 << (7 - 1)]G1Affine
type bucketG1AffineC8 [1 << (8 - 1)]G1Affine
type bucketG1AffineC9 [1 << (9 - 1)]G1Affine
type bucketG1AffineC10 [1 << (10 - 1)]G1Affine
type bucketG1AffineC11 [1 << (11 - 1)]G1Affine
type bucketG1AffineC12 [1 << (12 - 1)]G1Affine
type bucketG1AffineC13 [1 << (13 - 1)]G1Affine
type bucketG1AffineC14 [1 << (14 - 1)]G1Affine
type bucketG1AffineC15 [1 << (15 - 1)]G1Affine
type bucketG1AffineC16 [1 << (16 - 1)]G1Affine
type bucketG1AffineC20 [1 << (20 - 1)]G1Affine
type bucketG1AffineC21 [1 << (21 - 1)]G1Affine

type ibG1Affine interface {
	bucketG1AffineC4 |
		bucketG1AffineC5 |
		bucketG1AffineC6 |
		bucketG1AffineC7 |
		bucketG1AffineC8 |
		bucketG1AffineC9 |
		bucketG1AffineC10 |
		bucketG1AffineC11 |
		bucketG1AffineC12 |
		bucketG1AffineC13 |
		bucketG1AffineC14 |
		bucketG1AffineC15 |
		bucketG1AffineC16 |
		bucketG1AffineC20 |
		bucketG1AffineC21
}

// processChunkG2BatchAffine process a chunk of the scalars during the msm
// using affine coordinates for the buckets. To amortize the cost of the inverse in the affine addition
// we use a batch affine addition.
//
// this is derived from a PR by 0x0ece : https://github.com/ConsenSys/gnark-crypto/pull/249
// See Section 5.3: ia.cr/2022/1396
func processChunkG2BatchAffine[B ibG2Affine](chunk uint64,
	chRes chan<- g2JacExtended,
	c uint64,
	points []G2Affine,
	digits []uint32) {

	// init the buckets
	var buckets B
	for i := 0; i < len(buckets); i++ {
		buckets[i].setInfinity()
	}

	// setup for the batch affine;
	batchSize := len(buckets) / 5
	if batchSize > MAX_BATCH_SIZE {
		batchSize = MAX_BATCH_SIZE
	}
	if batchSize <= 0 {
		batchSize = 1
	}
	bucketIds := make(map[uint32]struct{}, len(buckets)/2) // TODO @gbotrel tune the capacity here
	cptP := 0                                              // count the number of point added to current batch

	var P [MAX_BATCH_SIZE]G2Affine  // allocated on the stack
	var R [MAX_BATCH_SIZE]*G2Affine // ...

	canAdd := func(bID uint32) bool {
		_, ok := bucketIds[bID]
		return !ok
	}

	isFull := func() bool {
		return cptP == batchSize
	}

	executeAndReset := func() {
		if cptP == 0 {
			return
		}
		BatchAddG2Affine(R[:cptP], P[:cptP], cptP)
		for k := range bucketIds {
			delete(bucketIds, k)
		}
		cptP = 0
	}

	add := func(op batchOp) {
		// CanAdd must be called before --> ensures bucket is not "used" in current batch

		BK := &buckets[op.bucketID]
		PP := &points[op.pointID>>1]
		if PP.IsInfinity() {
			return
		}
		// handle special cases with inf or -P / P
		if BK.IsInfinity() {
			if op.isNeg() {
				BK.Neg(PP)
			} else {
				BK.Set(PP)
			}
			return
		}
		if op.isNeg() {
			// if bucket == P --> -P == 0
			if BK.Equal(PP) {
				BK.setInfinity()
				return
			}
		} else {
			// if bucket == -P, B == 0
			if BK.X.Equal(&PP.X) && !BK.Y.Equal(&PP.Y) {
				BK.setInfinity()
				return
			}
		}

		// bucketIds[cptP] = op.bucketID
		bucketIds[op.bucketID] = struct{}{}
		R[cptP] = BK
		if op.isNeg() {
			P[cptP].Neg(PP)
		} else {
			P[cptP].Set(PP)
		}
		cptP++
	}

	queue := make([]batchOp, 0, 4096) // TODO find right capacity here.

	processQueue := func() {
		// for i := len(queue) - 1; i >= 0; i-- {
		for i := 0; i < len(queue); i++ {
			if canAdd(queue[i].bucketID) {
				add(queue[i])
				if isFull() {
					executeAndReset()
				}
				queue[i] = queue[len(queue)-1]
				queue = queue[:len(queue)-1]
				i--
			}
		}
	}

	nbBatches := 0
	for i, digit := range digits {

		if digit == 0 {
			continue
		}

		op := batchOp{pointID: uint32(i) << 1}
		// if msbWindow bit is set, we need to substract
		if digit&1 == 0 {
			// add
			op.bucketID = uint32((digit >> 1) - 1)
		} else {
			// sub
			op.bucketID = (uint32((digit >> 1)))
			op.pointID += 1
		}
		if canAdd(op.bucketID) {
			add(op)
			if isFull() {
				executeAndReset()
				nbBatches++
				if len(queue) != 0 { // TODO @gbotrel this doesn't seem to help much? should minimize queue resizing
					add(queue[len(queue)-1])
					queue = queue[:len(queue)-1]
				}
				// processQueue()
			}
		} else {
			// put it in queue.
			queue = append(queue, op)
		}
	}
	// fmt.Printf("chunk %d\nlen(queue)=%d\nnbBatches=%d\nbatchSize=%d\nnbBuckets=%d\nnbPoints=%d\n",
	// 	chunk, len(queue), nbBatches, batch.batchSize, len(buckets), len(points))
	// executeAndReset()
	for len(queue) != 0 {
		processQueue()
		executeAndReset() // execute batch even if not full.
	}

	// flush items in batch.
	executeAndReset()

	// reduce buckets into total
	// total =  bucket[0] + 2*bucket[1] + 3*bucket[2] ... + n*bucket[n-1]

	var runningSum, total g2JacExtended
	runningSum.setInfinity()
	total.setInfinity()
	for k := len(buckets) - 1; k >= 0; k-- {
		if !buckets[k].IsInfinity() {
			runningSum.addMixed(&buckets[k])
		}
		total.add(&runningSum)
	}

	chRes <- total

}

// we declare the buckets as fixed-size array types
// this allow us to allocate the buckets on the stack
type bucketG2AffineC4 [1 << (4 - 1)]G2Affine
type bucketG2AffineC5 [1 << (5 - 1)]G2Affine
type bucketG2AffineC6 [1 << (6 - 1)]G2Affine
type bucketG2AffineC7 [1 << (7 - 1)]G2Affine
type bucketG2AffineC8 [1 << (8 - 1)]G2Affine
type bucketG2AffineC9 [1 << (9 - 1)]G2Affine
type bucketG2AffineC10 [1 << (10 - 1)]G2Affine
type bucketG2AffineC11 [1 << (11 - 1)]G2Affine
type bucketG2AffineC12 [1 << (12 - 1)]G2Affine
type bucketG2AffineC13 [1 << (13 - 1)]G2Affine
type bucketG2AffineC14 [1 << (14 - 1)]G2Affine
type bucketG2AffineC15 [1 << (15 - 1)]G2Affine
type bucketG2AffineC16 [1 << (16 - 1)]G2Affine
type bucketG2AffineC20 [1 << (20 - 1)]G2Affine
type bucketG2AffineC21 [1 << (21 - 1)]G2Affine

type ibG2Affine interface {
	bucketG2AffineC4 |
		bucketG2AffineC5 |
		bucketG2AffineC6 |
		bucketG2AffineC7 |
		bucketG2AffineC8 |
		bucketG2AffineC9 |
		bucketG2AffineC10 |
		bucketG2AffineC11 |
		bucketG2AffineC12 |
		bucketG2AffineC13 |
		bucketG2AffineC14 |
		bucketG2AffineC15 |
		bucketG2AffineC16 |
		bucketG2AffineC20 |
		bucketG2AffineC21
}
