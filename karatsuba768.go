// Copyright (c) 2017 Pedro Martelletto. All rights reserved.
// Use of this source code is governed by a BSD-style license
// that can be found in the LICENSE file.

// This package implements the 768n x 768n polynomial multiplication algorithm
// presented in section 6 of https://ntruprime.cr.yp.to/ntruprime-20160511.pdf.

package karatsuba768

import "crypto/subtle"

type thinPoly []int32

// Freeze reduces x modulo 9829, for x in (-165191050,+165191050).
func Freeze(x int32) int32 {
	x -= 9829 * ((13*x) >> 17)
	x -= 9829 * ((427*x + 2097152) >> 22)
	y := x + 9829
	v := subtle.ConstantTimeLessOrEq(int(x), -1)
	return int32(subtle.ConstantTimeSelect(v, int(y), int(x)))
}

func (p thinPoly) Freeze() thinPoly {
	for i := range p {
		p[i] = Freeze(p[i])
	}
	return p
}

// Zero clears the contents of p.
func (p thinPoly) Zero() thinPoly {
	for i := range p {
		p[i] = 0
	}
	return p
}

// Set copies the contents of x to p.
func (p thinPoly) Set(x []int32) thinPoly {
	for i := range x {
		p[i] = x[i]
	}
	return p
}

// Add sets p to the addition a + b.
func (p thinPoly) Add(a, b []int32) thinPoly {
	for i := range a {
		p[i] = Freeze(a[i] + b[i])
	}
	return p
}

// Inc increments the contents of p by x.
func (p thinPoly) Inc(x []int32) thinPoly {
	for i := range x {
		p[i] += x[i]
	}
	return p
}

// Mul sets p to the multiplication of the polynomial p by the constant c.
func (p thinPoly) Mul(c int32, v []int32) thinPoly {
	for i := range v {
		p[i] = c * v[i]
	}
	return p.Freeze()
}

// x4Mul implements 4n x 4n, the lowest level of the multiplication algorithm.
func (p thinPoly) x4Mul(f, g thinPoly) thinPoly {
	p.Zero()
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			p[i+j] += Freeze(f[i] * g[j])
		}
	}
	return p
}

// Karatsuba5 uses x4Mul to implement 8n x 8xn.
func (p thinPoly) Karatsuba5(f, g thinPoly) thinPoly {
	var t = make(thinPoly, 8)
	var z = make(thinPoly, 16)
	f0, f1 := f[:4], f[4:]
	g0, g1 := g[:4], g[4:]

	t.x4Mul(f0, g0)
	z.Set(t)
	t.x4Mul(f1, g1)
	z[4:].Inc(t.Mul(-1, t))

	p.Set(z)
	p[4:].Inc(z.Mul(-1, z)[:12])
	t.x4Mul(z.Add(f0, f1), z[4:].Add(g0, g1))
	p[4:].Inc(t)

	return p
}

// Karatsuba4 uses Karatsuba5 to implement 16n x 16n.
func (p thinPoly) Karatsuba4(f, g thinPoly) thinPoly {
	var t = make(thinPoly, 16)
	var z = make(thinPoly, 32)
	f0, f1 := f[:8], f[8:]
	g0, g1 := g[:8], g[8:]

	t.Karatsuba5(f0, g0)
	z.Set(t)
	t.Karatsuba5(f1, g1)
	z[8:].Inc(t.Mul(-1, t))

	p.Set(z)
	p[8:].Inc(z.Mul(-1, z)[:24])
	t.Karatsuba5(z.Add(f0, f1), z[8:].Add(g0, g1))
	p[8:].Inc(t)

	return p
}

// Karatsuba3 uses Karatsuba4 to implement 32n x 32n.
func (p thinPoly) Karatsuba3(f, g thinPoly) thinPoly {
	var t = make(thinPoly, 32)
	var z = make(thinPoly, 64)
	f0, f1 := f[:16], f[16:]
	g0, g1 := g[:16], g[16:]

	t.Karatsuba4(f0, g0)
	z.Set(t)
	t.Karatsuba4(f1, g1)
	z[16:].Inc(t.Mul(-1, t))

	p.Set(z)
	p[16:].Inc(z.Mul(-1, z)[:48])
	t.Karatsuba4(z.Add(f0, f1), z[16:].Add(g0, g1))
	p[16:].Inc(t)

	return p
}

// Karatsuba2 uses Karatsuba3 to implement 64n x 64n.
func (p thinPoly) Karatsuba2(f, g thinPoly) thinPoly {
	var t = make(thinPoly, 64)
	var z = make(thinPoly, 128)
	f0, f1 := f[:32], f[32:]
	g0, g1 := g[:32], g[32:]

	t.Karatsuba3(f0, g0)
	z.Set(t)
	t.Karatsuba3(f1, g1)
	z[32:].Inc(t.Mul(-1, t))

	p.Set(z)
	p[32:].Inc(z.Mul(-1, z)[:96])
	t.Karatsuba3(z.Add(f0, f1), z[32:].Add(g0, g1))
	p[32:].Inc(t)

	return p
}

// Karatsuba1 uses Karatsuba2 to implement 128n x 128n.
func (p thinPoly) Karatsuba1(f, g thinPoly) thinPoly {
	var t = make(thinPoly, 128)
	var z = make(thinPoly, 256)
	f0, f1 := f[:64], f[64:]
	g0, g1 := g[:64], g[64:]

	t.Karatsuba2(f0, g0)
	z.Set(t)
	t.Karatsuba2(f1, g1)
	z[64:].Inc(t.Mul(-1, t))

	p.Set(z)
	p[64:].Inc(z.Mul(-1, z)[:192])
	t.Karatsuba2(z.Add(f0, f1), z[64:].Add(g0, g1))
	p[64:].Inc(t)

	return p.Freeze()
}

// Map with Toom6 coefficients for selected points.
var toomEvalCoeffs = map[int][]int32 {
	+1: { 1, 1, 1, 1, 1, 1 },
	-1: { 1, -1, 1, -1, 1, -1 },
	+2: { 1, 2, 4, 8, 16, 32},
	-2: { 1, -2, 4, -8, 16, -32 },
	+3: { 1, 3, 9, 27, 81, 243 },
	-3: { 1, -3, 9, -27, 81, -243 },
	+4: { 1, 4, 16, 64, 256, 1024 },
	-4: { 1, -4, 16, -64, 256, -1024 },
	+5: { 1, 5, 25, 125, 625, 3125 },
}

// toomEval evaluates the Toom6 factorization of f*g over GF(9829) at p.
func toomEval(p int, f, g *[768]int32) []int32 {
	a := make(thinPoly, 128)
	b := make(thinPoly, 128)
	t := make(thinPoly, 128)

	for i,v := range toomEvalCoeffs[p] {
		a.Inc(t.Mul(v, f[i*128:(i+1)*128]))
		b.Inc(t.Mul(v, g[i*128:(i+1)*128]))
	}

	return make(thinPoly, 256).Karatsuba1(a.Freeze(), b.Freeze())
}

// Interpolation parameters for Toom6.
var toomParam = [][]int32 {
	{ 7863, 1, 6552, 3276, 8425, 8893, 234, 5090, 4895, 3916, 6949 },
	{ 1705, 7864, 7864, 8846, 8846, 1841, 1841, 5169, 5169, 0, 576 },
	{ 9488, 9569, 7381, 7131, 33, 308, 1920, 8107, 2319, 2889, 4100 },
	{ 3328, 9228, 9228, 2041, 2041, 8027, 8027, 8527, 8527, 0, 9009 },
	{ 3266, 2727, 4935, 8102, 157, 6737, 6138, 8742, 9147, 9023, 8464 },
	{ 6655, 5993, 5993, 9515, 9515, 5365, 5365, 372, 372, 0, 273 },
	{ 8498, 2819, 5952, 901, 3916, 1018, 5776, 3309, 2826, 4301, 150 },
	{ 7969, 1488, 1488, 9085, 9085, 4425, 4425, 5590, 5590, 0, 9799 },
	{ 372, 9457, 9581, 248, 7127, 2702, 5590, 4239, 471, 9358, 9824 },
}

// toomInterpolate performs a linear interpolation of 'points' with the
// parameters passed in 'param'.
func toomInterpolate(points [][]int32, param []int32) []int32 {
	t := make(thinPoly, 256)
	u := make(thinPoly, 256)

	for i := range points {
		t.Inc(u.Mul(param[i], points[i]))
	}

	return t.Freeze()
}

// Toom6 decomposes a 768n x 768n multiplication into six instances of 128n x
// 128n. It is the highest level of the multiplication algorithm.
func (r thinPoly) Toom6(f, g *[768]int32) thinPoly {
	var e = [][]int32 {
		make(thinPoly, 256).Karatsuba1(f[0:128], g[0:128]),
		toomEval(+1, f, g),
		toomEval(-1, f, g),
		toomEval(+2, f, g),
		toomEval(-2, f, g),
		toomEval(+3, f, g),
		toomEval(-3, f, g),
		toomEval(+4, f, g),
		toomEval(-4, f, g),
		toomEval(+5, f, g),
		make(thinPoly, 256).Karatsuba1(f[640:768], g[640:768]),
	}
	var c = [][]int32 {
		e[0],
		toomInterpolate(e, toomParam[0]),
		toomInterpolate(e, toomParam[1]),
		toomInterpolate(e, toomParam[2]),
		toomInterpolate(e, toomParam[3]),
		toomInterpolate(e, toomParam[4]),
		toomInterpolate(e, toomParam[5]),
		toomInterpolate(e, toomParam[6]),
		toomInterpolate(e, toomParam[7]),
		toomInterpolate(e, toomParam[8]),
		e[10],
	}

	copy(r[:128], c[0])
	r[128:].Add(c[0][128:], c[1][:128])
	r[256:].Add(c[1][128:], c[2][:128])
	r[384:].Add(c[2][128:], c[3][:128])
	r[512:].Add(c[3][128:], c[4][:128])
	r[640:].Add(c[4][128:], c[5][:128])
	r[768:].Add(c[5][128:], c[6][:128])
	r[896:].Add(c[6][128:], c[7][:128])
	r[1024:].Add(c[7][128:], c[8][:128])
	r[1152:].Add(c[8][128:], c[9][:128])
	r[1280:].Add(c[9][128:], c[10][:128])
	copy(r[1408:], c[10][128:])

	return r
}

// Main entry point.
func Mul(h *[1536]int32, f, g *[768]int32) {
	z := thinPoly(h[:])
	z.Toom6(f, g)
}
