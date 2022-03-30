// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestCal2jd(t *testing.T) {
	var j int
	var djm0, djm float64

	j = gofa.Cal2jd(2003, 06, 01, &djm0, &djm)

	vvd(t, djm0, 2400000.5, EPS7, "Cal2jd", "djm0")
	vvd(t, djm, 52791.0, EPS7, "Cal2jd", "djm")

	viv(t, 0, j, "Cal2jd", "status")
}

func TestJd2cal(t *testing.T) {
	var dj1, dj2, fd float64
	var iy, im, id, j int

	dj1 = 2400000.5
	dj2 = 50123.9999

	j = gofa.Jd2cal(dj1, dj2, &iy, &im, &id, &fd)

	viv(t, iy, 1996, "Jd2cal", "y")
	viv(t, im, 2, "Jd2cal", "m")
	viv(t, id, 10, "Jd2cal", "d")
	vvd(t, fd, 0.9999, EPS7, "Jd2cal", "fd")
	viv(t, 0, j, "Cal2jd", "")

}

func TestJdcalf(t *testing.T) {
	var dj1, dj2 float64
	var iydmf [4]int
	var j int

	dj1 = 2400000.5
	dj2 = 50123.9999

	j = gofa.Jdcalf(4, dj1, dj2, &iydmf)

	viv(t, iydmf[0], 1996, "iauJdcalf", "y")
	viv(t, iydmf[1], 2, "iauJdcalf", "m")
	viv(t, iydmf[2], 10, "iauJdcalf", "d")
	viv(t, iydmf[3], 9999, "iauJdcalf", "f")

	viv(t, j, 0, "iauJdcalf", "j")
}

func TestEpb(t *testing.T) {
	var epb float64

	epb = gofa.Epb(2415019.8135, 30103.18648)

	vvd(t, epb, 1982.418424159278580, 1e-12, "iauEpb", "")
}

func TestEpj(t *testing.T) {
	var epj float64

	epj = gofa.Epj(2451545, -7392.5)

	vvd(t, epj, 1979.760438056125941, 1e-12, "iauEpj", "")
}

func TestEpb2jd(t *testing.T) {
	var epb, djm0, djm float64

	epb = 1957.3

	gofa.Epb2jd(epb, &djm0, &djm)

	vvd(t, djm0, 2400000.5, 1e-9, "iauEpb2jd", "djm0")
	vvd(t, djm, 35948.1915101513, 1e-9, "iauEpb2jd", "mjd")
}

func TestEpj2jd(t *testing.T) {
	var epj, djm0, djm float64

	epj = 1996.8

	gofa.Epj2jd(epj, &djm0, &djm)

	vvd(t, djm0, 2400000.5, 1e-9, "iauEpj2jd", "djm0")
	vvd(t, djm, 50375.7, 1e-9, "iauEpj2jd", "mjd")
}
