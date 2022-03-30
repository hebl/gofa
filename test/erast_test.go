// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestEe00(t *testing.T) {
	var epsa, dpsi, ee float64

	epsa = 0.4090789763356509900
	dpsi = -0.9630909107115582393e-5

	ee = gofa.Ee00(2400000.5, 53736.0, epsa, dpsi)

	vvd(t, ee, -0.8834193235367965479e-5, 1e-18, "Ee00", "")
}

func TestEe00a(t *testing.T) {
	ee := gofa.Ee00a(2400000.5, 53736.0)

	vvd(t, ee, -0.8834192459222588227e-5, 1e-18, "Ee00a", "")
}

func TestEe00b(t *testing.T) {
	ee := gofa.Ee00b(2400000.5, 53736.0)

	vvd(t, ee, -0.8835700060003032831e-5, 1e-18, "Ee00b", "")
}
func TestEe06a(t *testing.T) {
	ee := gofa.Ee06a(2400000.5, 53736.0)

	vvd(t, ee, -0.8834195072043790156e-5, 1e-15, "Ee06a", "")
}

func TestEect00(t *testing.T) {
	eect := gofa.Eect00(2400000.5, 53736.0)

	vvd(t, eect, 0.2046085004885125264e-8, 1e-20, "Eect00", "")
}

func TestEqeq94(t *testing.T) {
	eqeq := gofa.Eqeq94(2400000.5, 41234.0)

	vvd(t, eqeq, 0.5357758254609256894e-4, 1e-17, "Eqeq94", "")
}

func TestEra00(t *testing.T) {
	era00 := gofa.Era00(2400000.5, 54388.0)

	vvd(t, era00, 0.4022837240028158102, 1e-12, "Era00", "")
}

func TestGmst00(t *testing.T) {
	theta := gofa.Gmst00(2400000.5, 53736.0, 2400000.5, 53736.0)

	vvd(t, theta, 1.754174972210740592, 1e-12, "Gmst00", "")
}

func TestGmst06(t *testing.T) {
	theta := gofa.Gmst06(2400000.5, 53736.0, 2400000.5, 53736.0)

	vvd(t, theta, 1.754174971870091203, 1e-12, "Gmst06", "")
}

func TestGmst82(t *testing.T) {
	theta := gofa.Gmst82(2400000.5, 53736.0)

	vvd(t, theta, 1.754174981860675096, 1e-12, "Gmst82", "")
}

func TestGst00a(t *testing.T) {
	theta := gofa.Gst00a(2400000.5, 53736.0, 2400000.5, 53736.0)

	vvd(t, theta, 1.754166138018281369, 1e-12, "Gst00a", "")
}

func TestGst00b(t *testing.T) {
	theta := gofa.Gst00b(2400000.5, 53736.0)

	vvd(t, theta, 1.754166136510680589, 1e-12, "Gst00b", "")
}

func TestGst06(t *testing.T) {
	var rnpb [3][3]float64

	rnpb[0][0] = 0.9999989440476103608
	rnpb[0][1] = -0.1332881761240011518e-2
	rnpb[0][2] = -0.5790767434730085097e-3

	rnpb[1][0] = 0.1332858254308954453e-2
	rnpb[1][1] = 0.9999991109044505944
	rnpb[1][2] = -0.4097782710401555759e-4

	rnpb[2][0] = 0.5791308472168153320e-3
	rnpb[2][1] = 0.4020595661593994396e-4
	rnpb[2][2] = 0.9999998314954572365

	theta := gofa.Gst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb)

	vvd(t, theta, 1.754166138018167568, 1e-12, "Gst06", "")
}

func TestGst06a(t *testing.T) {
	theta := gofa.Gst06a(2400000.5, 53736.0, 2400000.5, 53736.0)

	vvd(t, theta, 1.754166137675019159, 1e-12, "Gst06a", "")
}

func TestGst94(t *testing.T) {
	theta := gofa.Gst94(2400000.5, 53736.0)

	vvd(t, theta, 1.754166136020645203, 1e-12, "Gst94", "")
}
