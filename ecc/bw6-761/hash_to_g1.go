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

package bw6761

import (
	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-761/fp"

	"math/big"
)

//Note: This only works for simple extensions

func g1IsogenyXNumerator(dst *fp.Element, x *fp.Element) {
	g1EvalPolynomial(dst,
		false,
		[]fp.Element{
			{8909527437822417859, 3882721843441849971, 10272758549140236519, 14003949960675888070, 3082674948893537730, 15057121805406331995, 9473729116041321980, 4728985833820774787, 8655877921490901299, 4190069830971789310, 5954844214359654504, 3167048922498340},
			{18108974101656748848, 2338690158884275892, 6508602482182596491, 4145012576748966774, 3414660410706516893, 2573833883260883150, 6164918126664265221, 14548675305047459086, 16639788212446610270, 18398701356395340891, 6095608892147913870, 52477293239851844},
			{9054487050828374424, 10392717116296913754, 3254301241091298245, 11295878325229259195, 1707330205353258446, 10510288978485217383, 3082459063332132610, 7274337652523729543, 17543266143078080943, 9199350678197670445, 3047804446073956935, 26238646619925922},
		},
		x)
}

func g1IsogenyXDenominator(dst *fp.Element, x *fp.Element) {
	g1EvalPolynomial(dst,
		true,
		[]fp.Element{
			{289919226011913130, 13019990545710127566, 4409829457611675068, 13030600802816293865, 15696054586628993047, 9353078419867322391, 5664203968291172875, 5090703637405909511, 17774776443174359288, 10018561694451762270, 12632664537138156478, 46143195394855163},
		},
		x)
}

func g1IsogenyYNumerator(dst *fp.Element, x *fp.Element, y *fp.Element) {
	var _dst fp.Element
	g1EvalPolynomial(&_dst,
		false,
		[]fp.Element{
			{9682176248893090322, 2627273429413213811, 1155528216520376823, 1734722477587034670, 13988724381275734601, 17289533515091656624, 2581744904959040264, 16263110058591731584, 231510300096278344, 819211016254091825, 9584860091064199543, 19904548774929241},
			{14209419774307277534, 17047004024416446496, 12006050873920801753, 7382661640201664267, 5619017447097588016, 4097933930624713700, 13346346473479882378, 10676906847998820547, 18226515408490094624, 14642258392207702855, 11108762314101178010, 33023872084892202},
			{8716717078775571656, 12731407275181189647, 9762903723273894736, 15440890901978225969, 5121990616059775339, 13084122861746100533, 9247377189996397831, 3376268883861637013, 15736310281815139598, 9151307960883459721, 9143413338221870806, 78715939859777766},
			{4527243525414187212, 14419730595003232685, 10850522657400424930, 5647939162614629597, 10077037139531405031, 5255144489242608691, 10764601568520842113, 12860540863116640579, 17995005108393816279, 13823047375953611030, 1523902223036978467, 13119323309962961},
		},
		x)

	dst.Mul(&_dst, y)
}

func g1IsogenyYDenominator(dst *fp.Element, x *fp.Element) {
	g1EvalPolynomial(dst,
		true,
		[]fp.Element{
			{2800676018270776722, 404959871884879410, 14461481433037540995, 11679465559666498996, 9481399069190242817, 18023312492583527742, 3661347334798803492, 4152305114258814443, 3867985292375803742, 13391491194096551019, 1887398969680023676, 20806804014868441},
			{4201014027406165083, 9830811844682094923, 3245478075846759876, 8295826302644972687, 14222098603785364226, 8588224665165739997, 14715393039052981047, 6228457671388221664, 15025349975418481421, 1640492717435274912, 12054470491374811323, 31210206022302661},
			{1690257235147301491, 3999098444797791463, 11640570174130445566, 9646961545794767555, 1990010047514562840, 18364734666159086263, 16718249672545350429, 7166856194535316732, 10485397052507485351, 16714307291500037780, 4352991985123392508, 56546597402289384},
		},
		x)
}

func g1Isogeny(p *G1Affine) {

	den := make([]fp.Element, 2)

	g1IsogenyYDenominator(&den[1], &p.X)
	g1IsogenyXDenominator(&den[0], &p.X)

	g1IsogenyYNumerator(&p.Y, &p.X, &p.Y)
	g1IsogenyXNumerator(&p.X, &p.X)

	den = fp.BatchInvert(den)

	p.X.Mul(&p.X, &den[0])
	p.Y.Mul(&p.Y, &den[1])
}

// g1SqrtRatio computes the square root of u/v and returns 0 iff u/v was indeed a quadratic residue
// if not, we get sqrt(Z * u / v). Recall that Z is non-residue
// If v = 0, u/v is meaningless and the output is unspecified, without raising an error.
// The main idea is that since the computation of the square root involves taking large powers of u/v, the inversion of v can be avoided
func g1SqrtRatio(z *fp.Element, u *fp.Element, v *fp.Element) uint64 {
	// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#name-optimized-sqrt_ratio-for-q- (3 mod 4)
	var tv1 fp.Element
	tv1.Square(v) // 1. tv1 = v^2
	var tv2 fp.Element
	tv2.Mul(u, v)       // 2. tv2 = u * v
	tv1.Mul(&tv1, &tv2) // 3. tv1 = tv1 * tv2

	var y1 fp.Element
	{
		var c1 big.Int
		// c1 = 1722862596078933134849197420568914385619917228134037527378447540052405855560872934021920795822352921910216141938446653362790439780138561939837377924781325399737901274844627212593135907855899198987974925107492278210691228279767074
		c1.SetBytes([]byte{72, 186, 9, 62, 224, 243, 130, 180, 97, 242, 80, 1, 62, 191, 207, 174, 73, 134, 26, 160, 116, 81, 162, 20, 160, 157, 123, 224, 33, 239, 144, 92, 30, 233, 142, 57, 97, 58, 70, 64, 243, 174, 191, 201, 109, 8, 193, 33, 162, 114, 59, 68, 190, 127, 100, 28, 119, 52, 247, 28, 250, 255, 203, 166, 40, 69, 176, 149, 153, 234, 62, 5, 131, 62, 43, 186, 188, 41, 13, 249, 164, 79, 154, 28, 0, 0, 32, 189, 39, 64, 0, 0, 0, 0, 34}) // c1 = (q - 3) / 4     # Integer arithmetic

		y1.Exp(tv1, &c1) // 4. y1 = tv1^c1
	}

	y1.Mul(&y1, &tv2) // 5. y1 = y1 * tv2

	var y2 fp.Element
	// c2 = sqrt(-Z)
	tv3 := fp.Element{10289215067249928212, 13987875627487618797, 10154775028297877632, 5892581882377791321, 12835424790914788634, 14963278386355512102, 10283221901563449361, 9868336211540881409, 7345304935218488881, 6998778443322886180, 9453359982570584357, 56775348355244645}
	y2.Mul(&y1, &tv3)              // 6. y2 = y1 * c2
	tv3.Square(&y1)                // 7. tv3 = y1^2
	tv3.Mul(&tv3, v)               // 8. tv3 = tv3 * v
	isQNr := tv3.NotEqual(u)       // 9. isQR = tv3 == u
	z.Select(int(isQNr), &y1, &y2) // 10. y = CMOV(y2, y1, isQR)
	return isQNr
}

// g1MulByZ multiplies x by [2] and stores the result in z
func g1MulByZ(z *fp.Element, x *fp.Element) {

	res := *x

	res.Double(&res)

	*z = res
}

// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#name-simplified-swu-method
// mapToCurve1 implements the SSWU map
// No cofactor clearing or isogeny
func mapToCurve1(u *fp.Element) G1Affine {

	var sswuIsoCurveCoeffA = fp.Element{12169852093062392636, 3867460573998792965, 2540986171999662608, 3377838107874487171, 6313266756742099767, 5994530928773814047, 5007141583730923456, 2345996307867737670, 7096861766432061441, 10014420324597579745, 8416419844935780388, 63340978449966806}
	var sswuIsoCurveCoeffB = fp.Element{9514135687797572479, 9972495974968977338, 17954535578332286571, 7437044986470910914, 13903267017721129281, 1871129682978723308, 13401268269932482209, 739043012311877982, 12116264695643437343, 1632209977726909861, 3621981106970059143, 65605772132525947}

	var tv1 fp.Element
	tv1.Square(u) // 1.  tv1 = u^2

	//mul tv1 by Z
	g1MulByZ(&tv1, &tv1) // 2.  tv1 = Z * tv1

	var tv2 fp.Element
	tv2.Square(&tv1)    // 3.  tv2 = tv1^2
	tv2.Add(&tv2, &tv1) // 4.  tv2 = tv2 + tv1

	var tv3 fp.Element
	var tv4 fp.Element
	tv4.SetOne()
	tv3.Add(&tv2, &tv4)                // 5.  tv3 = tv2 + 1
	tv3.Mul(&tv3, &sswuIsoCurveCoeffB) // 6.  tv3 = B * tv3

	tv2NZero := g1NotZero(&tv2)

	// tv4 = Z
	tv4 = fp.Element{289919226011913130, 13019990545710127566, 4409829457611675068, 13030600802816293865, 15696054586628993047, 9353078419867322391, 5664203968291172875, 5090703637405909511, 17774776443174359288, 10018561694451762270, 12632664537138156478, 46143195394855163}

	tv2.Neg(&tv2)
	tv4.Select(int(tv2NZero), &tv4, &tv2) // 7.  tv4 = CMOV(Z, -tv2, tv2 != 0)
	tv4.Mul(&tv4, &sswuIsoCurveCoeffA)    // 8.  tv4 = A * tv4

	tv2.Square(&tv3) // 9.  tv2 = tv3^2

	var tv6 fp.Element
	tv6.Square(&tv4) // 10. tv6 = tv4^2

	var tv5 fp.Element
	tv5.Mul(&tv6, &sswuIsoCurveCoeffA) // 11. tv5 = A * tv6

	tv2.Add(&tv2, &tv5) // 12. tv2 = tv2 + tv5
	tv2.Mul(&tv2, &tv3) // 13. tv2 = tv2 * tv3
	tv6.Mul(&tv6, &tv4) // 14. tv6 = tv6 * tv4

	tv5.Mul(&tv6, &sswuIsoCurveCoeffB) // 15. tv5 = B * tv6
	tv2.Add(&tv2, &tv5)                // 16. tv2 = tv2 + tv5

	var x fp.Element
	x.Mul(&tv1, &tv3) // 17.   x = tv1 * tv3

	var y1 fp.Element
	gx1NSquare := g1SqrtRatio(&y1, &tv2, &tv6) // 18. (is_gx1_square, y1) = sqrt_ratio(tv2, tv6)

	var y fp.Element
	y.Mul(&tv1, u) // 19.   y = tv1 * u

	y.Mul(&y, &y1) // 20.   y = y * y1

	x.Select(int(gx1NSquare), &tv3, &x) // 21.   x = CMOV(x, tv3, is_gx1_square)
	y.Select(int(gx1NSquare), &y1, &y)  // 22.   y = CMOV(y, y1, is_gx1_square)

	y1.Neg(&y)
	y.Select(int(g1Sgn0(u)^g1Sgn0(&y)), &y, &y1)

	// 23.  e1 = sgn0(u) == sgn0(y)
	// 24.   y = CMOV(-y, y, e1)

	x.Div(&x, &tv4) // 25.   x = x / tv4

	return G1Affine{x, y}
}

func g1EvalPolynomial(z *fp.Element, monic bool, coefficients []fp.Element, x *fp.Element) {
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

// hashToFp hashes msg to count prime field elements.
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-5.2
func hashToFp(msg, dst []byte, count int) ([]fp.Element, error) {
	// 128 bits of security
	// L = ceil((ceil(log2(p)) + k) / 8), where k is the security parameter = 128
	const Bytes = 1 + (fp.Bits-1)/8
	const L = 16 + Bytes

	lenInBytes := count * L
	pseudoRandomBytes, err := ecc.ExpandMsgXmd(msg, dst, lenInBytes)
	if err != nil {
		return nil, err
	}

	res := make([]fp.Element, count)
	for i := 0; i < count; i++ {
		res[i].SetBytes(pseudoRandomBytes[i*L : (i+1)*L])
	}
	return res, nil
}

// g1Sgn0 is an algebraic substitute for the notion of sign in ordered fields
// Namely, every non-zero quadratic residue in a finite field of characteristic =/= 2 has exactly two square roots, one of each sign
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#name-the-sgn0-function
// The sign of an element is not obviously related to that of its Montgomery form
func g1Sgn0(z *fp.Element) uint64 {

	nonMont := *z
	nonMont.FromMont()
	// m == 1
	return nonMont[0] % 2

}

// MapToG1 invokes the SSWU map, and guarantees that the result is in g1
func MapToG1(u fp.Element) G1Affine {
	res := mapToCurve1(&u)
	//this is in an isogenous curve
	g1Isogeny(&res)
	res.ClearCofactor(&res)
	return res
}

// EncodeToG1 hashes a message to a point on the G1 curve using the SSWU map.
// It is faster than HashToG1, but the result is not uniformly distributed. Unsuitable as a random oracle.
// dst stands for "domain separation tag", a string unique to the construction using the hash function
//https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/13/#section-6.6.3
func EncodeToG1(msg, dst []byte) (G1Affine, error) {

	var res G1Affine
	u, err := hashToFp(msg, dst, 1)
	if err != nil {
		return res, err
	}

	res = mapToCurve1(&u[0])

	//this is in an isogenous curve
	g1Isogeny(&res)
	res.ClearCofactor(&res)
	return res, nil
}

// HashToG1 hashes a message to a point on the G1 curve using the SSWU map.
// Slower than EncodeToG1, but usable as a random oracle.
// dst stands for "domain separation tag", a string unique to the construction using the hash function
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-3
func HashToG1(msg, dst []byte) (G1Affine, error) {
	u, err := hashToFp(msg, dst, 2*1)
	if err != nil {
		return G1Affine{}, err
	}

	Q0 := mapToCurve1(&u[0])
	Q1 := mapToCurve1(&u[1])

	//TODO: Add in E' first, then apply isogeny
	g1Isogeny(&Q0)
	g1Isogeny(&Q1)

	var _Q0, _Q1 G1Jac
	_Q0.FromAffine(&Q0)
	_Q1.FromAffine(&Q1).AddAssign(&_Q0)

	_Q1.ClearCofactor(&_Q1)

	Q1.FromJacobian(&_Q1)
	return Q1, nil
}

func g1NotZero(x *fp.Element) uint64 {

	return x[0] | x[1] | x[2] | x[3] | x[4] | x[5] | x[6] | x[7] | x[8] | x[9] | x[10] | x[11]

}
