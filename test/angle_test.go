// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestAnp(t *testing.T) {
	vvd(t, gofa.Anp(-0.1), 6.183185307179586477, 1e-12, "Anp", "")
}

func TestAnpm(t *testing.T) {
	vvd(t, gofa.Anpm(-4.0), 2.283185307179586477, 1e-12, "Anpm", "")
}

func TestA2af(t *testing.T) {
	var idmsf [4]int
	var s byte

	gofa.A2af(4, 2.345, &s, &idmsf)

	vbv(t, s, '+', "A2af", "s")

	viv(t, idmsf[0], 134, "A2af", "0")
	viv(t, idmsf[1], 21, "A2af", "1")
	viv(t, idmsf[2], 30, "A2af", "2")
	viv(t, idmsf[3], 9706, "A2af", "3")
}

func TestA2tf(t *testing.T) {
	var ihmsf [4]int
	var s byte

	gofa.A2tf(4, -3.01234, &s, &ihmsf)

	vbv(t, s, '-', "A2tf", "s")

	viv(t, ihmsf[0], 11, "A2tf", "0")
	viv(t, ihmsf[1], 30, "A2tf", "1")
	viv(t, ihmsf[2], 22, "A2tf", "2")
	viv(t, ihmsf[3], 6484, "A2tf", "3")
}

func TestD2tf(t *testing.T) {
	var ihmsf [4]int
	var s byte

	gofa.D2tf(4, -0.987654321, &s, &ihmsf)

	vbv(t, s, '-', "D2tf", "s")

	viv(t, ihmsf[0], 23, "D2tf", "ihmsf[0]")
	viv(t, ihmsf[1], 42, "D2tf", "ihmsf[1]")
	viv(t, ihmsf[2], 13, "D2tf", "ihmsf[2]")
	viv(t, ihmsf[3], 3333, "D2tf", "ihmsf[3]")

}

func TestAf2a(t *testing.T) {
	var a float64
	var j int

	j = gofa.Af2a('-', 45, 13, 27.2, &a)

	vvd(t, a, -0.7893115794313644842, 1e-12, "Af2a", "a")
	viv(t, j, 0, "Af2a", "j")

}

func TestTf2a(t *testing.T) {
	var a float64
	var j int

	j = gofa.Tf2a('+', 4, 58, 20.2, &a)

	vvd(t, a, 1.301739278189537429, 1e-12, "Tf2a", "a")
	viv(t, j, 0, "Tf2a", "j")
}

func TestTf2d(t *testing.T) {
	var d float64
	var j int

	j = gofa.Tf2d(' ', 23, 55, 10.9, &d)

	vvd(t, d, 0.9966539351851851852, 1e-12, "Tf2d", "d")
	viv(t, j, 0, "Tf2d", "j")
}

func TestSepp(t *testing.T) {
	var a, b [3]float64
	var s float64

	a[0] = 1.0
	a[1] = 0.1
	a[2] = 0.2

	b[0] = -3.0
	b[1] = 1e-3
	b[2] = 0.2

	s = gofa.Sepp(a, b)

	vvd(t, s, 2.860391919024660768, 1e-12, "Sepp", "")
}

func TestSeps(t *testing.T) {
	var al, ap, bl, bp, s float64

	al = 1.0
	ap = 0.1

	bl = 0.2
	bp = -3.0

	s = gofa.Seps(al, ap, bl, bp)

	vvd(t, s, 2.346722016996998842, 1e-14, "Seps", "")
}

func TestPap(t *testing.T) {
	var a, b [3]float64
	var theta float64

	a[0] = 1.0
	a[1] = 0.1
	a[2] = 0.2

	b[0] = -3.0
	b[1] = 1e-3
	b[2] = 0.2

	theta = gofa.Pap(a, b)

	vvd(t, theta, 0.3671514267841113674, 1e-12, "Pap", "")

}

func TestPas(t *testing.T) {
	var al, ap, bl, bp, theta float64

	al = 1.0
	ap = 0.1
	bl = 0.2
	bp = -1.0

	theta = gofa.Pas(al, ap, bl, bp)

	vvd(t, theta, -2.724544922932270424, 1e-12, "Pas", "")
}
