// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestD2dtf(t *testing.T) {
	var j, iy, im, id int
	var ihmsf [4]int

	j = gofa.D2dtf("UTC", 5, 2400000.5, 49533.99999, &iy, &im, &id, &ihmsf)

	viv(t, iy, 1994, "D2dtf", "iy")
	viv(t, im, 6, "D2dtf", "im")
	viv(t, id, 30, "D2dtf", "id")
	viv(t, ihmsf[0], 23, "D2dtf", "ihmsf[0]")
	viv(t, ihmsf[1], 59, "D2dtf", "ihmsf[1]")
	viv(t, ihmsf[2], 60, "D2dtf", "ihmsf[2]")
	viv(t, ihmsf[3], 13599, "D2dtf", "ihmsf[3]")
	viv(t, j, 0, "D2dtf", "")
}
func TestDat(t *testing.T) {
	var j int
	var deltat float64

	j = gofa.Dat(2003, 6, 1, 0.0, &deltat)
	vvd(t, deltat, 32.0, EPS12, "Dat", "deltat")
	viv(t, j, 0, "Dat", "")

	j = gofa.Dat(2008, 1, 17, 0.0, &deltat)
	vvd(t, deltat, 33.0, EPS12, "Dat", "deltat")
	viv(t, j, 0, "Dat", "")

	j = gofa.Dat(2017, 9, 1, 0.0, &deltat)
	vvd(t, deltat, 37.0, EPS12, "Dat", "deltat")
	viv(t, j, 0, "Dat", "")
}
func TestDtdb(t *testing.T) {
	dtdb := gofa.Dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)

	vvd(t, dtdb, -0.1280368005936998991e-2, 1e-15, "Dtdb", "")
}

func TestDtf2d(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Dtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599, &u1, &u2)

	vvd(t, u1+u2, 2449534.49999, 1e-6, "Dtf2d", "u")
	viv(t, j, 0, "Dtf2d", "j")
}

func TestTaitt(t *testing.T) {
	var t1, t2 float64
	var j int

	j = gofa.Taitt(2453750.5, 0.892482639, &t1, &t2)

	vvd(t, t1, 2453750.5, 1e-6, "Taitt", "t1")
	vvd(t, t2, 0.892855139, 1e-12, "Taitt", "t2")
	viv(t, j, 0, "Taitt", "j")
}
func TestTaiut1(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Taiut1(2453750.5, 0.892482639, -32.6659, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Taiut1", "u1")
	vvd(t, u2, 0.8921045614537037037, 1e-12, "Taiut1", "u2")
	viv(t, j, 0, "Taiut1", "j")
}

func TestTaiutc(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Taiutc(2453750.5, 0.892482639, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Taiutc", "u1")
	vvd(t, u2, 0.8921006945555555556, 1e-12, "Taiutc", "u2")
	viv(t, j, 0, "Taiutc", "j")
}

func TestTcbtdb(t *testing.T) {
	var b1, b2 float64
	var j int

	j = gofa.Tcbtdb(2453750.5, 0.893019599, &b1, &b2)

	vvd(t, b1, 2453750.5, 1e-6, "Tcbtdb", "b1")
	vvd(t, b2, 0.8928551362746343397, 1e-12, "Tcbtdb", "b2")
	viv(t, j, 0, "Tcbtdb", "j")

}

func TestTcgtt(t *testing.T) {
	var t1, t2 float64
	var j int

	j = gofa.Tcgtt(2453750.5, 0.892862531, &t1, &t2)

	vvd(t, t1, 2453750.5, 1e-6, "Tcgtt", "t1")
	vvd(t, t2, 0.8928551387488816828, 1e-12, "Tcgtt", "t2")
	viv(t, j, 0, "Tcgtt", "j")
}

func TestTdbtcb(t *testing.T) {
	var b1, b2 float64
	var j int

	j = gofa.Tdbtcb(2453750.5, 0.892855137, &b1, &b2)

	vvd(t, b1, 2453750.5, 1e-6, "Tdbtcb", "b1")
	vvd(t, b2, 0.8930195997253656716, 1e-12, "Tdbtcb", "b2")
	viv(t, j, 0, "Tdbtcb", "j")

}

func TestTdbtt(t *testing.T) {
	var t1, t2 float64
	var j int

	j = gofa.Tdbtt(2453750.5, 0.892855137, -0.000201, &t1, &t2)

	vvd(t, t1, 2453750.5, 1e-6, "Tdbtt", "t1")
	vvd(t, t2, 0.8928551393263888889, 1e-12, "Tdbtt", "t2")
	viv(t, j, 0, "Tdbtt", "j")
}

func TestTttai(t *testing.T) {
	var a1, a2 float64
	var j int

	j = gofa.Tttai(2453750.5, 0.892482639, &a1, &a2)

	vvd(t, a1, 2453750.5, 1e-6, "Tttai", "a1")
	vvd(t, a2, 0.892110139, 1e-12, "Tttai", "a2")
	viv(t, j, 0, "Tttai", "j")

}

func TestTttcg(t *testing.T) {
	var g1, g2 float64
	var j int

	j = gofa.Tttcg(2453750.5, 0.892482639, &g1, &g2)

	vvd(t, g1, 2453750.5, 1e-6, "Tttcg", "g1")
	vvd(t, g2, 0.8924900312508587113, 1e-12, "Tttcg", "g2")
	viv(t, j, 0, "Tttcg", "j")
}

func TestTttdb(t *testing.T) {
	var b1, b2 float64
	var j int

	j = gofa.Tttdb(2453750.5, 0.892855139, -0.000201, &b1, &b2)

	vvd(t, b1, 2453750.5, 1e-6, "Tttdb", "b1")
	vvd(t, b2, 0.8928551366736111111, 1e-12, "Tttdb", "b2")
	viv(t, j, 0, "Tttdb", "j")
}

func TestTtut1(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Ttut1(2453750.5, 0.892855139, 64.8499, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Ttut1", "u1")
	vvd(t, u2, 0.8921045614537037037, 1e-12, "Ttut1", "u2")
	viv(t, j, 0, "Ttut1", "j")
}

func TestUt1tai(t *testing.T) {
	var a1, a2 float64
	var j int

	j = gofa.Ut1tai(2453750.5, 0.892104561, -32.6659, &a1, &a2)

	vvd(t, a1, 2453750.5, 1e-6, "Ut1tai", "a1")
	vvd(t, a2, 0.8924826385462962963, 1e-12, "Ut1tai", "a2")
	viv(t, j, 0, "Ut1tai", "j")
}

func TestUt1tt(t *testing.T) {
	var t1, t2 float64
	var j int

	j = gofa.Ut1tt(2453750.5, 0.892104561, 64.8499, &t1, &t2)

	vvd(t, t1, 2453750.5, 1e-6, "Ut1tt", "t1")
	vvd(t, t2, 0.8928551385462962963, 1e-12, "Ut1tt", "t2")
	viv(t, j, 0, "Ut1tt", "j")
}

func TestUt1utc(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Ut1utc(2453750.5, 0.892104561, 0.3341, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Ut1utc", "u1")
	vvd(t, u2, 0.8921006941018518519, 1e-12, "Ut1utc", "u2")
	viv(t, j, 0, "Ut1utc", "j")
}

func TestUtctai(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Utctai(2453750.5, 0.892100694, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Utctai", "u1")
	vvd(t, u2, 0.8924826384444444444, 1e-12, "Utctai", "u2")
	viv(t, j, 0, "Utctai", "j")

}

func TestUtcut1(t *testing.T) {
	var u1, u2 float64
	var j int

	j = gofa.Utcut1(2453750.5, 0.892100694, 0.3341, &u1, &u2)

	vvd(t, u1, 2453750.5, 1e-6, "Utcut1", "u1")
	vvd(t, u2, 0.8921045608981481481, 1e-12, "Utcut1", "u2")
	viv(t, j, 0, "Utcut1", "j")
}
