// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestIcrs2g(t *testing.T) {
	var dr, dd, dl, db float64

	dr = 5.9338074302227188048671087
	dd = -1.1784870613579944551540570
	gofa.Icrs2g(dr, dd, &dl, &db)
	vvd(t, dl, 5.5850536063818546461558, 1e-14, "Icrs2g", "L")
	vvd(t, db, -0.7853981633974483096157, 1e-14, "Icrs2g", "B")
}
func TestG2icrs(t *testing.T) {
	var dl, db, dr, dd float64

	dl = 5.5850536063818546461558105
	db = -0.7853981633974483096156608
	gofa.G2icrs(dl, db, &dr, &dd)
	vvd(t, dr, 5.9338074302227188048671, 1e-14, "G2icrs", "R")
	vvd(t, dd, -1.1784870613579944551541, 1e-14, "G2icrs", "D")
}
func TestAe2hd(t *testing.T) {
	var a, e, p, h, d float64

	a = 5.5
	e = 1.1
	p = 0.7

	gofa.Ae2hd(a, e, p, &h, &d)

	vvd(t, h, 0.5933291115507309663, 1e-14, "Ae2hd", "h")
	vvd(t, d, 0.9613934761647817620, 1e-14, "Ae2hd", "d")
}
func TestHd2ae(t *testing.T) {
	var h, d, p, a, e float64

	h = 1.1
	d = 1.2
	p = 0.3

	gofa.Hd2ae(h, d, p, &a, &e)

	vvd(t, a, 5.916889243730066194, 1e-13, "Hd2ae", "a")
	vvd(t, e, 0.4472186304990486228, 1e-14, "Hd2ae", "e")
}
func TestHd2pa(t *testing.T) {
	var h, d, p, q float64

	h = 1.1
	d = 1.2
	p = 0.3

	q = gofa.Hd2pa(h, d, p)

	vvd(t, q, 1.906227428001995580, 1e-13, "Hd2pa", "q")
}
func TestEcm06(t *testing.T) {
	var date1, date2 float64
	var rm [3][3]float64

	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Ecm06(date1, date2, &rm)

	vvd(t, rm[0][0], 0.9999952427708701137, 1e-14, "Ecm06", "rm11")
	vvd(t, rm[0][1], -0.2829062057663042347e-2, 1e-14, "Ecm06", "rm12")
	vvd(t, rm[0][2], -0.1229163741100017629e-2, 1e-14, "Ecm06", "rm13")
	vvd(t, rm[1][0], 0.3084546876908653562e-2, 1e-14, "Ecm06", "rm21")
	vvd(t, rm[1][1], 0.9174891871550392514, 1e-14, "Ecm06", "rm22")
	vvd(t, rm[1][2], 0.3977487611849338124, 1e-14, "Ecm06", "rm23")
	vvd(t, rm[2][0], 0.2488512951527405928e-5, 1e-14, "Ecm06", "rm31")
	vvd(t, rm[2][1], -0.3977506604161195467, 1e-14, "Ecm06", "rm32")
	vvd(t, rm[2][2], 0.9174935488232863071, 1e-14, "Ecm06", "rm33")
}
func TestEqec06(t *testing.T) {
	var date1, date2, dr, dd, dl, db float64

	date1 = 1234.5
	date2 = 2440000.5
	dr = 1.234
	dd = 0.987

	gofa.Eqec06(date1, date2, dr, dd, &dl, &db)

	vvd(t, dl, 1.342509918994654619, 1e-14, "Eqec06", "dl")
	vvd(t, db, 0.5926215259704608132, 1e-14, "Eqec06", "db")

}
func TestEceq06(t *testing.T) {
	var date1, date2, dl, db, dr, dd float64

	date1 = 2456165.5
	date2 = 0.401182685
	dl = 5.1
	db = -0.9

	gofa.Eceq06(date1, date2, dl, db, &dr, &dd)

	vvd(t, dr, 5.533459733613627767, 1e-14, "Eceq06", "dr")
	vvd(t, dd, -1.246542932554480576, 1e-14, "Eceq06", "dd")
}
func TestLtecm(t *testing.T) {
	var epj float64
	var rm [3][3]float64

	epj = -3000.0

	gofa.Ltecm(epj, &rm)

	vvd(t, rm[0][0], 0.3564105644859788825, 1e-14, "Ltecm", "rm11")
	vvd(t, rm[0][1], 0.8530575738617682284, 1e-14, "Ltecm", "rm12")
	vvd(t, rm[0][2], 0.3811355207795060435, 1e-14, "Ltecm", "rm13")
	vvd(t, rm[1][0], -0.9343283469640709942, 1e-14, "Ltecm", "rm21")
	vvd(t, rm[1][1], 0.3247830597681745976, 1e-14, "Ltecm", "rm22")
	vvd(t, rm[1][2], 0.1467872751535940865, 1e-14, "Ltecm", "rm23")
	vvd(t, rm[2][0], 0.1431636191201167793e-2, 1e-14, "Ltecm", "rm31")
	vvd(t, rm[2][1], -0.4084222566960599342, 1e-14, "Ltecm", "rm32")
	vvd(t, rm[2][2], 0.9127919865189030899, 1e-14, "Ltecm", "rm33")
}
func TestLteqec(t *testing.T) {
	var epj, dr, dd, dl, db float64

	epj = -1500.0
	dr = 1.234
	dd = 0.987

	gofa.Lteqec(epj, dr, dd, &dl, &db)

	vvd(t, dl, 0.5039483649047114859, 1e-14, "Lteqec", "dl")
	vvd(t, db, 0.5848534459726224882, 1e-14, "Lteqec", "db")
}
func TestLteceq(t *testing.T) {
	var epj, dl, db, dr, dd float64

	epj = 2500.0
	dl = 1.5
	db = 0.6

	gofa.Lteceq(epj, dl, db, &dr, &dd)

	vvd(t, dr, 1.275156021861921167, 1e-14, "Lteceq", "dr")
	vvd(t, dd, 0.9966573543519204791, 1e-14, "Lteceq", "dd")
}
func TestEform(t *testing.T) {
	var j int
	var a, f float64

	j = gofa.Eform(0, &a, &f)

	viv(t, j, -1, "Eform", "j0")

	j = gofa.Eform(gofa.WGS84, &a, &f)

	viv(t, j, 0, "Eform", "j1")
	vvd(t, a, 6378137.0, 1e-10, "Eform", "a1")
	vvd(t, f, 0.3352810664747480720e-2, 1e-18, "Eform", "f1")

	j = gofa.Eform(gofa.GRS80, &a, &f)

	viv(t, j, 0, "Eform", "j2")
	vvd(t, a, 6378137.0, 1e-10, "Eform", "a2")
	vvd(t, f, 0.3352810681182318935e-2, 1e-18, "Eform", "f2")

	j = gofa.Eform(gofa.WGS72, &a, &f)

	viv(t, j, 0, "Eform", "j2")
	vvd(t, a, 6378135.0, 1e-10, "Eform", "a3")
	vvd(t, f, 0.3352779454167504862e-2, 1e-18, "Eform", "f3")

	j = gofa.Eform(4, &a, &f)
	viv(t, j, -1, "Eform", "j3")
}
func TestGc2gd(t *testing.T) {
	var j int
	xyz := [3]float64{2e6, 3e6, 5.244e6}
	var e, p, h float64

	j = gofa.Gc2gd(0, xyz, &e, &p, &h)

	viv(t, j, -1, "Gc2gd", "j0")

	j = gofa.Gc2gd(gofa.WGS84, xyz, &e, &p, &h)

	viv(t, j, 0, "Gc2gd", "j1")
	vvd(t, e, 0.9827937232473290680, 1e-14, "Gc2gd", "e1")
	vvd(t, p, 0.97160184819075459, 1e-14, "Gc2gd", "p1")
	vvd(t, h, 331.4172461426059892, 1e-8, "Gc2gd", "h1")

	j = gofa.Gc2gd(gofa.GRS80, xyz, &e, &p, &h)

	viv(t, j, 0, "Gc2gd", "j2")
	vvd(t, e, 0.9827937232473290680, 1e-14, "Gc2gd", "e2")
	vvd(t, p, 0.97160184820607853, 1e-14, "Gc2gd", "p2")
	vvd(t, h, 331.41731754844348, 1e-8, "Gc2gd", "h2")

	j = gofa.Gc2gd(gofa.WGS72, xyz, &e, &p, &h)

	viv(t, j, 0, "Gc2gd", "j3")
	vvd(t, e, 0.9827937232473290680, 1e-14, "Gc2gd", "e3")
	vvd(t, p, 0.9716018181101511937, 1e-14, "Gc2gd", "p3")
	vvd(t, h, 333.2770726130318123, 1e-8, "Gc2gd", "h3")

	j = gofa.Gc2gd(4, xyz, &e, &p, &h)

	viv(t, j, -1, "Gc2gd", "j4")
}
func TestGc2gde(t *testing.T) {
	var j int
	a := 6378136.0
	f := 0.0033528
	xyz := [3]float64{2e6, 3e6, 5.244e6}
	var e, p, h float64

	j = gofa.Gc2gde(a, f, xyz, &e, &p, &h)

	viv(t, j, 0, "Gc2gde", "j")
	vvd(t, e, 0.9827937232473290680, 1e-14, "Gc2gde", "e")
	vvd(t, p, 0.9716018377570411532, 1e-14, "Gc2gde", "p")
	vvd(t, h, 332.36862495764397, 1e-8, "Gc2gde", "h")
}
func TestGd2gc(t *testing.T) {
	var j int
	e := 3.1
	p := -0.5
	h := 2500.0
	var xyz [3]float64

	j = gofa.Gd2gc(0, e, p, h, &xyz)

	viv(t, j, -1, "Gd2gc", "j0")

	j = gofa.Gd2gc(gofa.WGS84, e, p, h, &xyz)

	viv(t, j, 0, "Gd2gc", "j1")
	vvd(t, xyz[0], -5599000.5577049947, 1e-7, "Gd2gc", "1/1")
	vvd(t, xyz[1], 233011.67223479203, 1e-7, "Gd2gc", "2/1")
	vvd(t, xyz[2], -3040909.4706983363, 1e-7, "Gd2gc", "3/1")

	j = gofa.Gd2gc(gofa.GRS80, e, p, h, &xyz)

	viv(t, j, 0, "Gd2gc", "j2")
	vvd(t, xyz[0], -5599000.5577260984, 1e-7, "Gd2gc", "1/2")
	vvd(t, xyz[1], 233011.6722356702949, 1e-7, "Gd2gc", "2/2")
	vvd(t, xyz[2], -3040909.4706095476, 1e-7, "Gd2gc", "3/2")

	j = gofa.Gd2gc(gofa.WGS72, e, p, h, &xyz)

	viv(t, j, 0, "Gd2gc", "j3")
	vvd(t, xyz[0], -5598998.7626301490, 1e-7, "Gd2gc", "1/3")
	vvd(t, xyz[1], 233011.5975297822211, 1e-7, "Gd2gc", "2/3")
	vvd(t, xyz[2], -3040908.6861467111, 1e-7, "Gd2gc", "3/3")

	j = gofa.Gd2gc(4, e, p, h, &xyz)

	viv(t, j, -1, "Gd2gc", "j4")
}
func TestGd2gce(t *testing.T) {
	var j int
	a := 6378136.0
	f := 0.0033528
	e := 3.1
	p := -0.5
	h := 2500.0
	var xyz [3]float64

	j = gofa.Gd2gce(a, f, e, p, h, &xyz)

	viv(t, j, 0, "Gd2gce", "j")
	vvd(t, xyz[0], -5598999.6665116328, 1e-7, "Gd2gce", "1")
	vvd(t, xyz[1], 233011.6351463057189, 1e-7, "Gd2gce", "2")
	vvd(t, xyz[2], -3040909.0517314132, 1e-7, "Gd2gce", "3")
}
