// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestTpors(t *testing.T) {
	var xi, eta, ra, dec, az1, bz1, az2, bz2 float64
	var n int

	xi = -0.03
	eta = 0.07
	ra = 1.3
	dec = 1.5

	n = gofa.Tpors(xi, eta, ra, dec, &az1, &bz1, &az2, &bz2)

	vvd(t, az1, 1.736621577783208748, 1e-13, "Tpors", "az1")
	vvd(t, bz1, 1.436736561844090323, 1e-13, "Tpors", "bz1")

	vvd(t, az2, 4.004971075806584490, 1e-13, "Tpors", "az2")
	vvd(t, bz2, 1.565084088476417917, 1e-13, "Tpors", "bz2")

	viv(t, n, 2, "Tpors", "n")
}
func TestTporv(t *testing.T) {
	var xi, eta, ra, dec float64
	var v, vz1, vz2 [3]float64
	var n int

	xi = -0.03
	eta = 0.07
	ra = 1.3
	dec = 1.5
	gofa.S2c(ra, dec, &v)

	n = gofa.Tporv(xi, eta, v, &vz1, &vz2)

	vvd(t, vz1[0], -0.02206252822366888610, 1e-15, "Tporv", "x1")
	vvd(t, vz1[1], 0.1318251060359645016, 1e-14, "Tporv", "y1")
	vvd(t, vz1[2], 0.9910274397144543895, 1e-14, "Tporv", "z1")

	vvd(t, vz2[0], -0.003712211763801968173, 1e-16, "Tporv", "x2")
	vvd(t, vz2[1], -0.004341519956299836813, 1e-16, "Tporv", "y2")
	vvd(t, vz2[2], 0.9999836852110587012, 1e-14, "Tporv", "z2")

	viv(t, n, 2, "Tporv", "n")
}
func TestTpsts(t *testing.T) {
	var xi, eta, raz, decz, ra, dec float64

	xi = -0.03
	eta = 0.07
	raz = 2.3
	decz = 1.5

	gofa.Tpsts(xi, eta, raz, decz, &ra, &dec)

	vvd(t, ra, 0.7596127167359629775, 1e-14, "Tpsts", "ra")
	vvd(t, dec, 1.540864645109263028, 1e-13, "Tpsts", "dec")
}
func TestTpstv(t *testing.T) {
	var xi, eta, raz, decz float64
	var vz, v [3]float64

	xi = -0.03
	eta = 0.07
	raz = 2.3
	decz = 1.5
	gofa.S2c(raz, decz, &vz)

	gofa.Tpstv(xi, eta, vz, &v)

	vvd(t, v[0], 0.02170030454907376677, 1e-15, "Tpstv", "x")
	vvd(t, v[1], 0.02060909590535367447, 1e-15, "Tpstv", "y")
	vvd(t, v[2], 0.9995520806583523804, 1e-14, "Tpstv", "z")
}
func TestTpxes(t *testing.T) {
	var ra, dec, raz, decz, xi, eta float64
	var j int

	ra = 1.3
	dec = 1.55
	raz = 2.3
	decz = 1.5

	j = gofa.Tpxes(ra, dec, raz, decz, &xi, &eta)

	vvd(t, xi, -0.01753200983236980595, 1e-15, "Tpxes", "xi")
	vvd(t, eta, 0.05962940005778712891, 1e-15, "Tpxes", "eta")

	viv(t, j, 0, "Tpxes", "j")
}
func TestTpxev(t *testing.T) {
	var ra, dec, raz, decz, xi, eta float64
	var v, vz [3]float64
	var j int

	ra = 1.3
	dec = 1.55
	raz = 2.3
	decz = 1.5
	gofa.S2c(ra, dec, &v)
	gofa.S2c(raz, decz, &vz)

	j = gofa.Tpxev(v, vz, &xi, &eta)

	vvd(t, xi, -0.01753200983236980595, 1e-15, "Tpxev", "xi")
	vvd(t, eta, 0.05962940005778712891, 1e-15, "Tpxev", "eta")

	viv(t, j, 0, "Tpxev", "j")

}
