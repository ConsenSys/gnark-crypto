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

package bls12381

//Note: This only works for simple extensions

import (
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bls12-381/internal/fptower"
)

func g2IsogenyXNumerator(dst *fptower.E2, x *fptower.E2) {
	g2EvalPolynomial(dst,
		false,
		[]fptower.E2{
			{
				fp.Element{5185457120960601698, 494647221959407934, 8971396042087821730, 324544954362548322, 14214792730224113654, 1405280679127738945},
				fp.Element{5185457120960601698, 494647221959407934, 8971396042087821730, 324544954362548322, 14214792730224113654, 1405280679127738945},
			},
			{
				fp.Element{0},
				fp.Element{6910023028261548496, 9745789443900091043, 7668299866710145304, 2432656849393633605, 2897729527445498821, 776645607375592125},
			},
			{
				fp.Element{724047465092313539, 15783990863276714670, 12824896677063784855, 15246381572572671516, 13186611051602728692, 1485475813959743803},
				fp.Element{12678383550985550056, 4872894721950045521, 13057521970209848460, 10439700461551592610, 10672236800577525218, 388322803687796062},
			},
			{
				fp.Element{4659755689450087917, 1804066951354704782, 15570919779568036803, 15592734958806855601, 7597208057374167129, 1841438384006890194},
				fp.Element{0},
			},
		},
		x)
}

func g2IsogenyXDenominator(dst *fptower.E2, x *fptower.E2) {
	g2EvalPolynomial(dst,
		true,
		[]fptower.E2{
			{
				fp.Element{0},
				fp.Element{2250392438786206615, 17463829474098544446, 14571211649711714824, 4495761442775821336, 258811604141191305, 357646605018048850},
			},
			{
				fp.Element{4933130441833534766, 15904462746612662304, 8034115857496836953, 12755092135412849606, 7007796720291435703, 252692002104915169},
				fp.Element{8469300574244328829, 4752422838614097887, 17848302789776796362, 12930989898711414520, 16851051131888818207, 1621106615542624696},
			},
		},
		x)
}

func g2IsogenyYNumerator(dst *fptower.E2, x *fptower.E2, y *fptower.E2) {
	var _dst fptower.E2
	g2EvalPolynomial(&_dst,
		false,
		[]fptower.E2{
			{
				fp.Element{10869708750642247614, 13056187057366814946, 1750362034917495549, 6326189602300757217, 1140223926335695785, 632761649765668291},
				fp.Element{10869708750642247614, 13056187057366814946, 1750362034917495549, 6326189602300757217, 1140223926335695785, 632761649765668291},
			},
			{
				fp.Element{0},
				fp.Element{13765940311003083782, 5579209876153186557, 11349908400803699438, 11707848830955952341, 199199289641242246, 899896674917908607},
			},
			{
				fp.Element{15562563812347550836, 2436447360975022760, 6528760985104924230, 5219850230775796305, 5336118400288762609, 194161401843898031},
				fp.Element{16286611277439864375, 18220438224251737430, 906913588459157469, 2019487729638916206, 75985378181939686, 1679637215803641835},
			},
			{
				fp.Element{11849179119594500956, 13906615243538674725, 14543197362847770509, 2041759640812427310, 2879701092679313252, 1259985822978576468},
				fp.Element{0},
			},
		},
		x)

	dst.Mul(&_dst, y)
}

func g2IsogenyYDenominator(dst *fptower.E2, x *fptower.E2) {
	g2EvalPolynomial(dst,
		true,
		[]fptower.E2{
			{
				fp.Element{99923616639376095, 10339114964526300021, 6204619029868000785, 1288486622530663893, 14587509920085997152, 272081012460753233},
				fp.Element{99923616639376095, 10339114964526300021, 6204619029868000785, 1288486622530663893, 14587509920085997152, 272081012460753233},
			},
			{
				fp.Element{0},
				fp.Element{6751177316358619845, 15498000274876530106, 6820146801716041242, 13487284328327464010, 776434812423573915, 1072939815054146550},
			},
			{
				fp.Element{7399695662750302149, 14633322083064217648, 12051173786245255430, 9909266166264498601, 1288323043582377747, 379038003157372754},
				fp.Element{6002735353327561446, 6023563502162542543, 13831244861028377885, 15776815867859765525, 4123780734888324547, 1494760614490167112},
			},
		},
		x)
}

func g2Isogeny(p *G2Affine) {

	den := make([]fptower.E2, 2)

	g2IsogenyYDenominator(&den[1], &p.X)
	g2IsogenyXDenominator(&den[0], &p.X)

	g2IsogenyYNumerator(&p.Y, &p.X, &p.Y)
	g2IsogenyXNumerator(&p.X, &p.X)

	den = fptower.BatchInvert(den)

	p.X.Mul(&p.X, &den[0])
	p.Y.Mul(&p.Y, &den[1])
}

// g2SqrtRatio computes the square root of u/v and returns 0 iff u/v was indeed a quadratic residue
// if not, we get sqrt(Z * u / v). Recall that Z is non-residue
// The main idea is that since the computation of the square root involves taking large powers of u/v, the inversion of v can be avoided
func g2SqrtRatio(z *fptower.E2, u *fptower.E2, v *fptower.E2) uint64 {

	// Taken from https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/13/ F.2.1.1. for any field

	tv1 := fptower.E2{
		fp.Element{8921533702591418330, 15859389534032789116, 3389114680249073393, 15116930867080254631, 3288288975085550621, 1021049300055853010},
		fp.Element{8921533702591418330, 15859389534032789116, 3389114680249073393, 15116930867080254631, 3288288975085550621, 1021049300055853010},
	} //tv1 = c6

	var tv2, tv3, tv4, tv5 fptower.E2
	var exp big.Int
	// c4 = 7 = 2^3 - 1
	// q is odd so c1 is at least 1.
	exp.SetBytes([]byte{7})

	tv2.Exp(*v, &exp)
	tv3.Mul(&tv2, &tv2)
	tv3.Mul(&tv3, v)

	// line 5
	tv5.Mul(u, &tv3)

	// c3 = 1001205140483106588246484290269935788605945006208159541241399033561623546780709821462541004956387089373434649096260670658193992783731681621012512651314777238193313314641988297376025498093520728838658813979860931248214124593092835
	exp.SetBytes([]byte{42, 67, 122, 75, 140, 53, 252, 116, 189, 39, 142, 170, 34, 242, 94, 158, 45, 201, 14, 80, 231, 4, 107, 70, 110, 89, 228, 147, 73, 232, 189, 5, 10, 98, 207, 209, 109, 220, 166, 239, 83, 20, 147, 48, 151, 142, 240, 17, 214, 134, 25, 200, 97, 133, 199, 178, 146, 232, 90, 135, 9, 26, 4, 150, 107, 249, 30, 211, 231, 27, 116, 49, 98, 195, 56, 54, 33, 19, 207, 215, 206, 214, 177, 215, 99, 130, 234, 178, 106, 160, 0, 1, 199, 24, 227})
	tv5.Exp(tv5, &exp)
	tv5.Mul(&tv5, &tv2)
	tv2.Mul(&tv5, v)
	tv3.Mul(&tv5, u)

	// line 10
	tv4.Mul(&tv3, &tv2)

	// c5 = 4
	exp.SetBytes([]byte{4})
	tv5.Exp(tv4, &exp)
	//TODO: Optimize? Probably nah. If we did it'd look like the following comments:
	//TODO: This would however save a bunch of heap allocation
	// The optimization is useless as it does EXACTLY as the exp function
	//These branches are decided by static data. Probably optimized away by compiler.
	/*const logC5 = 3 - 1

	   // c1 is at least 3 bcz q is 1 mod 8. So logC5 is at least 2.
	   switch logC5{
	   case 0:
		   tv5.SetOne()
	   case 1:
			tv5.Set(&tv4)
	    default:
	        tv5.Double(&tv4)
	   }
	   for i := 1; i < logC5; i++ {
		   tv5.Double(&tv5)
	    }*/

	isQNr := g2NotOne(&tv5)

	tv2.Mul(&tv3, &fptower.E2{
		fp.Element{1921729236329761493, 9193968980645934504, 9862280504246317678, 6861748847800817560, 10375788487011937166, 4460107375738415},
		fp.Element{16821121318233475459, 10183025025229892778, 1779012082459463630, 3442292649700377418, 1061500799026501234, 1352426537312017168},
	})
	tv5.Mul(&tv4, &tv1)

	// line 15

	tv3.Select(int(isQNr), &tv3, &tv2)
	tv4.Select(int(isQNr), &tv4, &tv5)

	exp.Lsh(big.NewInt(1), 3-2)

	for i := 3; i >= 2; i-- {
		//line 20
		tv5.Exp(tv4, &exp)
		nE1 := g2NotOne(&tv5)

		tv2.Mul(&tv3, &tv1)
		tv1.Mul(&tv1, &tv1)
		tv5.Mul(&tv4, &tv1)

		tv3.Select(int(nE1), &tv3, &tv2)
		tv4.Select(int(nE1), &tv4, &tv5)

		exp.Rsh(&exp, 1)
	}

	*z = tv3
	return isQNr
}

// g2SetZ sets z to [-2 -1].
func g2SetZ(z *fptower.E2) {
	z.Set(&fptower.E2{
		fp.Element{9794203289623549276, 7309342082925068282, 1139538881605221074, 15659550692327388916, 16008355200866287827, 582484205531694093},
		fp.Element{4897101644811774638, 3654671041462534141, 569769440802610537, 17053147383018470266, 17227549637287919721, 291242102765847046},
	})
}

// g2MulByZ multiplies x by [-2 -1] and stores the result in z
func g2MulByZ(z *fptower.E2, x *fptower.E2) {

	panic("not implemented")
}

// From https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/13/ Pg 80
func g2SswuMap(u *fptower.E2) G2Affine {

	var tv1 fptower.E2
	tv1.Square(u)

	//mul tv1 by Z
	g2MulByZ(&tv1, &tv1)

	var tv2 fptower.E2
	tv2.Square(&tv1)
	tv2.Add(&tv2, &tv1)

	var tv3 fptower.E2
	//Standard doc line 5
	var tv4 fptower.E2
	tv4.SetOne()
	tv3.Add(&tv2, &tv4)
	tv3.Mul(&tv3, &fptower.E2{
		fp.Element{2515823342057463218, 7982686274772798116, 7934098172177393262, 8484566552980779962, 4455086327883106868, 1323173589274087377},
		fp.Element{2515823342057463218, 7982686274772798116, 7934098172177393262, 8484566552980779962, 4455086327883106868, 1323173589274087377},
	})

	tv2NZero := g2NotZero(&tv2)

	// tv4 = Z
	tv4 = fptower.E2{
		fp.Element{9794203289623549276, 7309342082925068282, 1139538881605221074, 15659550692327388916, 16008355200866287827, 582484205531694093},
		fp.Element{4897101644811774638, 3654671041462534141, 569769440802610537, 17053147383018470266, 17227549637287919721, 291242102765847046},
	}

	tv2.Neg(&tv2)
	tv4.Select(int(tv2NZero), &tv4, &tv2)
	tv2 = fptower.E2{
		fp.Element{0},
		fp.Element{16517514583386313282, 74322656156451461, 16683759486841714365, 815493829203396097, 204518332920448171, 1306242806803223655},
	}
	tv4.Mul(&tv4, &tv2)

	tv2.Square(&tv3)

	var tv6 fptower.E2
	//Standard doc line 10
	tv6.Square(&tv4)

	var tv5 fptower.E2
	tv5.Mul(&tv6, &fptower.E2{
		fp.Element{0},
		fp.Element{16517514583386313282, 74322656156451461, 16683759486841714365, 815493829203396097, 204518332920448171, 1306242806803223655},
	})

	tv2.Add(&tv2, &tv5)
	tv2.Mul(&tv2, &tv3)
	tv6.Mul(&tv6, &tv4)

	//Standards doc line 15
	tv5.Mul(&tv6, &fptower.E2{
		fp.Element{2515823342057463218, 7982686274772798116, 7934098172177393262, 8484566552980779962, 4455086327883106868, 1323173589274087377},
		fp.Element{2515823342057463218, 7982686274772798116, 7934098172177393262, 8484566552980779962, 4455086327883106868, 1323173589274087377},
	})
	tv2.Add(&tv2, &tv5)

	var x fptower.E2
	x.Mul(&tv1, &tv3)

	var y1 fptower.E2
	gx1NSquare := g2SqrtRatio(&y1, &tv2, &tv6)

	var y fptower.E2
	y.Mul(&tv1, u)

	//Standards doc line 20
	y.Mul(&y, &y1)

	x.Select(int(gx1NSquare), &tv3, &x)
	y.Select(int(gx1NSquare), &y1, &y)

	y1.Neg(&y)
	y.Select(int(g2Sgn0(u)^g2Sgn0(&y)), &y, &y1)

	//Standards doc line 25
	x.Div(&x, &tv4)

	return G2Affine{x, y}
}

// EncodeToCurveG2SSWU maps a fptower.E2 to a point on the curve using the Simplified Shallue and van de Woestijne Ulas map
//https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/13/#section-6.6.3
func EncodeToCurveG2SSWU(msg, dst []byte) (G2Affine, error) {

	var res G2Affine
	u, err := hashToFp(msg, dst, 2)
	if err != nil {
		return res, err
	}

	res = g2SswuMap(&fptower.E2{
		u[0],
		u[1],
	})

	//this is in an isogenous curve
	g2Isogeny(&res)

	res.ClearCofactor(&res)

	return res, nil
}

// HashToCurveG2SSWU hashes a byte string to the G2 curve. Usable as a random oracle.
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-3
func HashToCurveG2SSWU(msg, dst []byte) (G2Affine, error) {
	u, err := hashToFp(msg, dst, 2*2)
	if err != nil {
		return G2Affine{}, err
	}

	Q0 := g2SswuMap(&fptower.E2{
		u[0],
		u[1],
	})
	Q1 := g2SswuMap(&fptower.E2{
		u[2+0],
		u[2+1],
	})

	//TODO: Add in E' first, then apply isogeny
	g2Isogeny(&Q0)
	g2Isogeny(&Q1)

	var _Q0, _Q1 G2Jac
	_Q0.FromAffine(&Q0)
	_Q1.FromAffine(&Q1).AddAssign(&_Q0)

	_Q1.ClearCofactor(&_Q1)

	Q1.FromJacobian(&_Q1)
	return Q1, nil
}

// g2Sgn0 is an algebraic substitute for the notion of sign in ordered fields
// Namely, every non-zero quadratic residue in a finite field of characteristic =/= 2 has exactly two square roots, one of each sign
// Taken from https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/ section 4.1
// The sign of an element is not obviously related to that of its Montgomery form
func g2Sgn0(z *fptower.E2) uint64 {

	nonMont := *z
	nonMont.FromMont()

	sign := uint64(0)
	zero := uint64(1)
	var signI uint64
	var zeroI uint64

	zeroI = z.A0[0]
	signI = zeroI % 2
	zeroI = (zeroI | -zeroI) >> 63
	sign = sign | (zero & signI)
	zero = zero & zeroI

	zeroI = z.A1[0]
	signI = zeroI % 2
	zeroI = (zeroI | -zeroI) >> 63
	sign = sign | (zero & signI)
	zero = zero & zeroI

	return sign

}

func g2EvalPolynomial(z *fptower.E2, monic bool, coefficients []fptower.E2, x *fptower.E2) {
	dst := coefficients[len(coefficients)-1]

	if monic {
		dst.Add(&dst, x)
	}

	for i := len(coefficients) - 2; i >= 0; i-- {
		dst.Mul(&dst, x)
		dst.Add(&dst, &coefficients[i])
	}

	z.Set(&dst)
}

func g2NotZero(x *fptower.E2) uint64 {
	//Assuming G1 is over Fp and that if hashing is available for G2, it also is for G1
	return g1NotZero(&x.A0) | g1NotZero(&x.A1)

}

func g2NotOne(x *fptower.E2) uint64 {

	//Assuming hash is implemented for G1 and that the curve is over Fp
	return g1NotOne(&x.A0) | g1NotZero(&x.A1)

}
