// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestIr(t *testing.T) {
	var r [3][3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Ir(&r)

	vvd(t, r[0][0], 1.0, 0.0, "Ir", "11")
	vvd(t, r[0][1], 0.0, 0.0, "Ir", "12")
	vvd(t, r[0][2], 0.0, 0.0, "Ir", "13")

	vvd(t, r[1][0], 0.0, 0.0, "Ir", "21")
	vvd(t, r[1][1], 1.0, 0.0, "Ir", "22")
	vvd(t, r[1][2], 0.0, 0.0, "Ir", "23")

	vvd(t, r[2][0], 0.0, 0.0, "Ir", "31")
	vvd(t, r[2][1], 0.0, 0.0, "Ir", "32")
	vvd(t, r[2][2], 1.0, 0.0, "Ir", "33")
}

func TestZp(t *testing.T) {
	var p [3]float64

	p[0] = 0.3
	p[1] = 1.2
	p[2] = -2.5

	gofa.Zp(&p)

	vvd(t, p[0], 0.0, 0.0, "Zp", "1")
	vvd(t, p[1], 0.0, 0.0, "Zp", "2")
	vvd(t, p[2], 0.0, 0.0, "Zp", "3")
}

func TestZr(t *testing.T) {
	var r [3][3]float64

	r[0][0] = 2.0
	r[1][0] = 3.0
	r[2][0] = 2.0

	r[0][1] = 3.0
	r[1][1] = 2.0
	r[2][1] = 3.0

	r[0][2] = 3.0
	r[1][2] = 4.0
	r[2][2] = 5.0

	gofa.Zr(&r)

	vvd(t, r[0][0], 0.0, 0.0, "Zr", "00")
	vvd(t, r[1][0], 0.0, 0.0, "Zr", "01")
	vvd(t, r[2][0], 0.0, 0.0, "Zr", "02")

	vvd(t, r[0][1], 0.0, 0.0, "Zr", "10")
	vvd(t, r[1][1], 0.0, 0.0, "Zr", "11")
	vvd(t, r[2][1], 0.0, 0.0, "Zr", "12")

	vvd(t, r[0][2], 0.0, 0.0, "Zr", "20")
	vvd(t, r[1][2], 0.0, 0.0, "Zr", "21")
	vvd(t, r[2][2], 0.0, 0.0, "Zr", "22")
}

func TestZpv(t *testing.T) {
	var pv [2][3]float64

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = -0.5
	pv[1][1] = 3.1
	pv[1][2] = 0.9

	gofa.Zpv(&pv)

	vvd(t, pv[0][0], 0.0, 0.0, "Zpv", "p1")
	vvd(t, pv[0][1], 0.0, 0.0, "Zpv", "p2")
	vvd(t, pv[0][2], 0.0, 0.0, "Zpv", "p3")

	vvd(t, pv[1][0], 0.0, 0.0, "Zpv", "v1")
	vvd(t, pv[1][1], 0.0, 0.0, "Zpv", "v2")
	vvd(t, pv[1][2], 0.0, 0.0, "Zpv", "v3")

}

func TestCp(t *testing.T) {
	var p, c [3]float64

	p[0] = 0.3
	p[1] = 1.2
	p[2] = -2.5

	gofa.Cp(p, &c)

	vvd(t, c[0], 0.3, 0.0, "Cp", "1")
	vvd(t, c[1], 1.2, 0.0, "Cp", "2")
	vvd(t, c[2], -2.5, 0.0, "Cp", "3")
}

func TestCr(t *testing.T) {
	var r, c [3][3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Cr(r, &c)

	vvd(t, c[0][0], 2.0, 0.0, "Cr", "11")
	vvd(t, c[0][1], 3.0, 0.0, "Cr", "12")
	vvd(t, c[0][2], 2.0, 0.0, "Cr", "13")

	vvd(t, c[1][0], 3.0, 0.0, "Cr", "21")
	vvd(t, c[1][1], 2.0, 0.0, "Cr", "22")
	vvd(t, c[1][2], 3.0, 0.0, "Cr", "23")

	vvd(t, c[2][0], 3.0, 0.0, "Cr", "31")
	vvd(t, c[2][1], 4.0, 0.0, "Cr", "32")
	vvd(t, c[2][2], 5.0, 0.0, "Cr", "33")
}

func TestCpv(t *testing.T) {
	var pv, c [2][3]float64

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = -0.5
	pv[1][1] = 3.1
	pv[1][2] = 0.9

	gofa.Cpv(pv, &c)

	vvd(t, c[0][0], 0.3, 0.0, "Cpv", "p1")
	vvd(t, c[0][1], 1.2, 0.0, "Cpv", "p2")
	vvd(t, c[0][2], -2.5, 0.0, "Cpv", "p3")

	vvd(t, c[1][0], -0.5, 0.0, "Cpv", "v1")
	vvd(t, c[1][1], 3.1, 0.0, "Cpv", "v2")
	vvd(t, c[1][2], 0.9, 0.0, "Cpv", "v3")

}

func TestP2pv(t *testing.T) {
	var p [3]float64
	var pv [2][3]float64

	p[0] = 0.25
	p[1] = 1.2
	p[2] = 3.0

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = -0.5
	pv[1][1] = 3.1
	pv[1][2] = 0.9

	gofa.P2pv(p, &pv)

	vvd(t, pv[0][0], 0.25, 0.0, "P2pv", "p1")
	vvd(t, pv[0][1], 1.2, 0.0, "P2pv", "p2")
	vvd(t, pv[0][2], 3.0, 0.0, "P2pv", "p3")

	vvd(t, pv[1][0], 0.0, 0.0, "P2pv", "v1")
	vvd(t, pv[1][1], 0.0, 0.0, "P2pv", "v2")
	vvd(t, pv[1][2], 0.0, 0.0, "P2pv", "v3")
}

func TestPv2p(t *testing.T) {
	var pv [2][3]float64
	var p [3]float64

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = -0.5
	pv[1][1] = 3.1
	pv[1][2] = 0.9

	gofa.Pv2p(pv, &p)

	vvd(t, p[0], 0.3, 0.0, "Pv2p", "1")
	vvd(t, p[1], 1.2, 0.0, "Pv2p", "2")
	vvd(t, p[2], -2.5, 0.0, "Pv2p", "3")

}

func TestRx(t *testing.T) {
	var phi float64
	var r [3][3]float64

	phi = 0.3456789

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Rx(phi, &r)

	vvd(t, r[0][0], 2.0, 0.0, "Rx", "11")
	vvd(t, r[0][1], 3.0, 0.0, "Rx", "12")
	vvd(t, r[0][2], 2.0, 0.0, "Rx", "13")

	vvd(t, r[1][0], 3.839043388235612460, 1e-12, "Rx", "21")
	vvd(t, r[1][1], 3.237033249594111899, 1e-12, "Rx", "22")
	vvd(t, r[1][2], 4.516714379005982719, 1e-12, "Rx", "23")

	vvd(t, r[2][0], 1.806030415924501684, 1e-12, "Rx", "31")
	vvd(t, r[2][1], 3.085711545336372503, 1e-12, "Rx", "32")
	vvd(t, r[2][2], 3.687721683977873065, 1e-12, "Rx", "33")
}

func TestRy(t *testing.T) {
	var theta float64
	var r [3][3]float64

	theta = 0.3456789

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Ry(theta, &r)

	vvd(t, r[0][0], 0.8651847818978159930, 1e-12, "Ry", "11")
	vvd(t, r[0][1], 1.467194920539316554, 1e-12, "Ry", "12")
	vvd(t, r[0][2], 0.1875137911274457342, 1e-12, "Ry", "13")

	vvd(t, r[1][0], 3, 1e-12, "Ry", "21")
	vvd(t, r[1][1], 2, 1e-12, "Ry", "22")
	vvd(t, r[1][2], 3, 1e-12, "Ry", "23")

	vvd(t, r[2][0], 3.500207892850427330, 1e-12, "Ry", "31")
	vvd(t, r[2][1], 4.779889022262298150, 1e-12, "Ry", "32")
	vvd(t, r[2][2], 5.381899160903798712, 1e-12, "Ry", "33")
}

func TestRz(t *testing.T) {
	var psi float64
	var r [3][3]float64

	psi = 0.3456789

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Rz(psi, &r)

	vvd(t, r[0][0], 2.898197754208926769, 1e-12, "Rz", "11")
	vvd(t, r[0][1], 3.500207892850427330, 1e-12, "Rz", "12")
	vvd(t, r[0][2], 2.898197754208926769, 1e-12, "Rz", "13")

	vvd(t, r[1][0], 2.144865911309686813, 1e-12, "Rz", "21")
	vvd(t, r[1][1], 0.865184781897815993, 1e-12, "Rz", "22")
	vvd(t, r[1][2], 2.144865911309686813, 1e-12, "Rz", "23")

	vvd(t, r[2][0], 3.0, 1e-12, "Rz", "31")
	vvd(t, r[2][1], 4.0, 1e-12, "Rz", "32")
	vvd(t, r[2][2], 5.0, 1e-12, "Rz", "33")
}

func TestS2c(t *testing.T) {
	var c [3]float64

	gofa.S2c(3.0123, -0.999, &c)

	vvd(t, c[0], -0.5366267667260523906, 1e-12, "S2c", "1")
	vvd(t, c[1], 0.0697711109765145365, 1e-12, "S2c", "2")
	vvd(t, c[2], -0.8409302618566214041, 1e-12, "S2c", "3")

}

func TestC2s(t *testing.T) {
	var p [3]float64
	var theta, phi float64

	p[0] = 100.0
	p[1] = -50.0
	p[2] = 25.0

	gofa.C2s(p, &theta, &phi)

	vvd(t, theta, -0.4636476090008061162, 1e-14, "C2s", "theta")
	vvd(t, phi, 0.2199879773954594463, 1e-14, "C2s", "phi")
}

func TestP2s(t *testing.T) {
	var p [3]float64
	var theta, phi, r float64

	p[0] = 100.0
	p[1] = -50.0
	p[2] = 25.0

	gofa.P2s(p, &theta, &phi, &r)

	vvd(t, theta, -0.4636476090008061162, 1e-12, "P2s", "theta")
	vvd(t, phi, 0.2199879773954594463, 1e-12, "P2s", "phi")
	vvd(t, r, 114.5643923738960002, 1e-9, "P2s", "r")
}

func TestS2p(t *testing.T) {
	var p [3]float64

	gofa.S2p(-3.21, 0.123, 0.456, &p)

	vvd(t, p[0], -0.4514964673880165228, 1e-12, "S2p", "x")
	vvd(t, p[1], 0.0309339427734258688, 1e-12, "S2p", "y")
	vvd(t, p[2], 0.0559466810510877933, 1e-12, "S2p", "z")

}

func TestS2pv(t *testing.T) {
	var pv [2][3]float64

	gofa.S2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5, &pv)

	vvd(t, pv[0][0], -0.4514964673880165228, 1e-12, "S2pv", "x")
	vvd(t, pv[0][1], 0.0309339427734258688, 1e-12, "S2pv", "y")
	vvd(t, pv[0][2], 0.0559466810510877933, 1e-12, "S2pv", "z")

	vvd(t, pv[1][0], 0.1292270850663260170e-4, 1e-16, "S2pv", "vx")
	vvd(t, pv[1][1], 0.2652814182060691422e-5, 1e-16, "S2pv", "vy")
	vvd(t, pv[1][2], 0.2568431853930292259e-5, 1e-16, "S2pv", "vz")
}

func TestPv2s(t *testing.T) {
	var pv [2][3]float64
	var theta, phi, r, td, pd, rd float64

	pv[0][0] = -0.4514964673880165
	pv[0][1] = 0.03093394277342585
	pv[0][2] = 0.05594668105108779

	pv[1][0] = 1.292270850663260e-5
	pv[1][1] = 2.652814182060692e-6
	pv[1][2] = 2.568431853930293e-6

	gofa.Pv2s(pv, &theta, &phi, &r, &td, &pd, &rd)

	vvd(t, theta, 3.073185307179586515, 1e-12, "Pv2s", "theta")
	vvd(t, phi, 0.1229999999999999992, 1e-12, "Pv2s", "phi")
	vvd(t, r, 0.4559999999999999757, 1e-12, "Pv2s", "r")
	vvd(t, td, -0.7800000000000000364e-5, 1e-16, "Pv2s", "td")
	vvd(t, pd, 0.9010000000000001639e-5, 1e-16, "Pv2s", "pd")
	vvd(t, rd, -0.1229999999999999832e-4, 1e-16, "Pv2s", "rd")
}

func TestPpp(t *testing.T) {
	var a, b, apb [3]float64

	a[0] = 2.0
	a[1] = 2.0
	a[2] = 3.0

	b[0] = 1.0
	b[1] = 3.0
	b[2] = 4.0

	gofa.Ppp(a, b, &apb)

	vvd(t, apb[0], 3.0, 1e-12, "Ppp", "0")
	vvd(t, apb[1], 5.0, 1e-12, "Ppp", "1")
	vvd(t, apb[2], 7.0, 1e-12, "Ppp", "2")
}

func TestPmp(t *testing.T) {
	var a, b, amb [3]float64

	a[0] = 2.0
	a[1] = 2.0
	a[2] = 3.0

	b[0] = 1.0
	b[1] = 3.0
	b[2] = 4.0

	gofa.Pmp(a, b, &amb)

	vvd(t, amb[0], 1.0, 1e-12, "Pmp", "0")
	vvd(t, amb[1], -1.0, 1e-12, "Pmp", "1")
	vvd(t, amb[2], -1.0, 1e-12, "Pmp", "2")

}

func TestPpsp(t *testing.T) {
	var a, b, apsb [3]float64
	var s float64

	a[0] = 2.0
	a[1] = 2.0
	a[2] = 3.0

	s = 5.0

	b[0] = 1.0
	b[1] = 3.0
	b[2] = 4.0

	gofa.Ppsp(a, s, b, &apsb)

	vvd(t, apsb[0], 7.0, 1e-12, "Ppsp", "0")
	vvd(t, apsb[1], 17.0, 1e-12, "Ppsp", "1")
	vvd(t, apsb[2], 23.0, 1e-12, "Ppsp", "2")
}

func TestPdp(t *testing.T) {
	var a, b [3]float64
	var adb float64

	a[0] = 2.0
	a[1] = 2.0
	a[2] = 3.0

	b[0] = 1.0
	b[1] = 3.0
	b[2] = 4.0

	adb = gofa.Pdp(a, b)

	vvd(t, adb, 20, 1e-12, "Pdp", "")

}

func TestPxp(t *testing.T) {
	var a, b, axb [3]float64

	a[0] = 2.0
	a[1] = 2.0
	a[2] = 3.0

	b[0] = 1.0
	b[1] = 3.0
	b[2] = 4.0

	gofa.Pxp(a, b, &axb)

	vvd(t, axb[0], -1.0, 1e-12, "Pxp", "1")
	vvd(t, axb[1], -5.0, 1e-12, "Pxp", "2")
	vvd(t, axb[2], 4.0, 1e-12, "Pxp", "3")
}

func TestPm(t *testing.T) {
	var p [3]float64
	var r float64

	p[0] = 0.3
	p[1] = 1.2
	p[2] = -2.5

	r = gofa.Pm(p)

	vvd(t, r, 2.789265136196270604, 1e-12, "Pm", "")
}

func TestPn(t *testing.T) {
	var p, u [3]float64
	var r float64

	p[0] = 0.3
	p[1] = 1.2
	p[2] = -2.5

	gofa.Pn(p, &r, &u)

	vvd(t, r, 2.789265136196270604, 1e-12, "Pn", "r")

	vvd(t, u[0], 0.1075552109073112058, 1e-12, "Pn", "u1")
	vvd(t, u[1], 0.4302208436292448232, 1e-12, "Pn", "u2")
	vvd(t, u[2], -0.8962934242275933816, 1e-12, "Pn", "u3")
}

func TestSxp(t *testing.T) {
	var s float64
	var p, sp [3]float64

	s = 2.0

	p[0] = 0.3
	p[1] = 1.2
	p[2] = -2.5

	gofa.Sxp(s, p, &sp)

	vvd(t, sp[0], 0.6, 0.0, "Sxp", "1")
	vvd(t, sp[1], 2.4, 0.0, "Sxp", "2")
	vvd(t, sp[2], -5.0, 0.0, "Sxp", "3")
}

func TestPvppv(t *testing.T) {
	var a, b, apb [2][3]float64

	a[0][0] = 2.0
	a[0][1] = 2.0
	a[0][2] = 3.0

	a[1][0] = 5.0
	a[1][1] = 6.0
	a[1][2] = 3.0

	b[0][0] = 1.0
	b[0][1] = 3.0
	b[0][2] = 4.0

	b[1][0] = 3.0
	b[1][1] = 2.0
	b[1][2] = 1.0

	gofa.Pvppv(a, b, &apb)

	vvd(t, apb[0][0], 3.0, 1e-12, "Pvppv", "p1")
	vvd(t, apb[0][1], 5.0, 1e-12, "Pvppv", "p2")
	vvd(t, apb[0][2], 7.0, 1e-12, "Pvppv", "p3")

	vvd(t, apb[1][0], 8.0, 1e-12, "Pvppv", "v1")
	vvd(t, apb[1][1], 8.0, 1e-12, "Pvppv", "v2")
	vvd(t, apb[1][2], 4.0, 1e-12, "Pvppv", "v3")
}

func TestPvmpv(t *testing.T) {
	var a, b, amb [2][3]float64

	a[0][0] = 2.0
	a[0][1] = 2.0
	a[0][2] = 3.0

	a[1][0] = 5.0
	a[1][1] = 6.0
	a[1][2] = 3.0

	b[0][0] = 1.0
	b[0][1] = 3.0
	b[0][2] = 4.0

	b[1][0] = 3.0
	b[1][1] = 2.0
	b[1][2] = 1.0

	gofa.Pvmpv(a, b, &amb)

	vvd(t, amb[0][0], 1.0, 1e-12, "Pvmpv", "11")
	vvd(t, amb[0][1], -1.0, 1e-12, "Pvmpv", "21")
	vvd(t, amb[0][2], -1.0, 1e-12, "Pvmpv", "31")

	vvd(t, amb[1][0], 2.0, 1e-12, "Pvmpv", "12")
	vvd(t, amb[1][1], 4.0, 1e-12, "Pvmpv", "22")
	vvd(t, amb[1][2], 2.0, 1e-12, "Pvmpv", "32")
}

func TestPvdpv(t *testing.T) {
	var a, b [2][3]float64
	var adb [2]float64

	a[0][0] = 2.0
	a[0][1] = 2.0
	a[0][2] = 3.0

	a[1][0] = 6.0
	a[1][1] = 0.0
	a[1][2] = 4.0

	b[0][0] = 1.0
	b[0][1] = 3.0
	b[0][2] = 4.0

	b[1][0] = 0.0
	b[1][1] = 2.0
	b[1][2] = 8.0

	gofa.Pvdpv(a, b, &adb)

	vvd(t, adb[0], 20.0, 1e-12, "Pvdpv", "1")
	vvd(t, adb[1], 50.0, 1e-12, "Pvdpv", "2")
}

func TestPvxpv(t *testing.T) {
	var a, b, axb [2][3]float64

	a[0][0] = 2.0
	a[0][1] = 2.0
	a[0][2] = 3.0

	a[1][0] = 6.0
	a[1][1] = 0.0
	a[1][2] = 4.0

	b[0][0] = 1.0
	b[0][1] = 3.0
	b[0][2] = 4.0

	b[1][0] = 0.0
	b[1][1] = 2.0
	b[1][2] = 8.0

	gofa.Pvxpv(a, b, &axb)

	vvd(t, axb[0][0], -1.0, 1e-12, "Pvxpv", "p1")
	vvd(t, axb[0][1], -5.0, 1e-12, "Pvxpv", "p2")
	vvd(t, axb[0][2], 4.0, 1e-12, "Pvxpv", "p3")

	vvd(t, axb[1][0], -2.0, 1e-12, "Pvxpv", "v1")
	vvd(t, axb[1][1], -36.0, 1e-12, "Pvxpv", "v2")
	vvd(t, axb[1][2], 22.0, 1e-12, "Pvxpv", "v3")
}

func TestPvm(t *testing.T) {
	var pv [2][3]float64
	var r, s float64

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = 0.45
	pv[1][1] = -0.25
	pv[1][2] = 1.1

	gofa.Pvm(pv, &r, &s)

	vvd(t, r, 2.789265136196270604, 1e-12, "Pvm", "r")
	vvd(t, s, 1.214495780149111922, 1e-12, "Pvm", "s")
}

func TestSxpv(t *testing.T) {
	var s float64
	var pv, spv [2][3]float64

	s = 2.0

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = 0.5
	pv[1][1] = 3.2
	pv[1][2] = -0.7

	gofa.Sxpv(s, pv, &spv)

	vvd(t, spv[0][0], 0.6, 0.0, "Sxpv", "p1")
	vvd(t, spv[0][1], 2.4, 0.0, "Sxpv", "p2")
	vvd(t, spv[0][2], -5.0, 0.0, "Sxpv", "p3")

	vvd(t, spv[1][0], 1.0, 0.0, "Sxpv", "v1")
	vvd(t, spv[1][1], 6.4, 0.0, "Sxpv", "v2")
	vvd(t, spv[1][2], -1.4, 0.0, "Sxpv", "v3")
}

func TestS2xpv(t *testing.T) {
	var s1, s2 float64
	var pv, spv [2][3]float64

	s1 = 2.0
	s2 = 3.0

	pv[0][0] = 0.3
	pv[0][1] = 1.2
	pv[0][2] = -2.5

	pv[1][0] = 0.5
	pv[1][1] = 2.3
	pv[1][2] = -0.4

	gofa.S2xpv(s1, s2, pv, &spv)

	vvd(t, spv[0][0], 0.6, 1e-12, "S2xpv", "p1")
	vvd(t, spv[0][1], 2.4, 1e-12, "S2xpv", "p2")
	vvd(t, spv[0][2], -5.0, 1e-12, "S2xpv", "p3")

	vvd(t, spv[1][0], 1.5, 1e-12, "S2xpv", "v1")
	vvd(t, spv[1][1], 6.9, 1e-12, "S2xpv", "v2")
	vvd(t, spv[1][2], -1.2, 1e-12, "S2xpv", "v3")
}

func TestPvu(t *testing.T) {
	var pv, upv [2][3]float64

	pv[0][0] = 126668.5912743160734
	pv[0][1] = 2136.792716839935565
	pv[0][2] = -245251.2339876830229

	pv[1][0] = -0.4051854035740713039e-2
	pv[1][1] = -0.6253919754866175788e-2
	pv[1][2] = 0.1189353719774107615e-1

	gofa.Pvu(2920.0, pv, &upv)

	vvd(t, upv[0][0], 126656.7598605317105, 1e-6, "Pvu", "p1")
	vvd(t, upv[0][1], 2118.531271155726332, 1e-8, "Pvu", "p2")
	vvd(t, upv[0][2], -245216.5048590656190, 1e-6, "Pvu", "p3")

	vvd(t, upv[1][0], -0.4051854035740713039e-2, 1e-12, "Pvu", "v1")
	vvd(t, upv[1][1], -0.6253919754866175788e-2, 1e-12, "Pvu", "v2")
	vvd(t, upv[1][2], 0.1189353719774107615e-1, 1e-12, "Pvu", "v3")
}

func TestPvup(t *testing.T) {
	var pv [2][3]float64
	var p [3]float64

	pv[0][0] = 126668.5912743160734
	pv[0][1] = 2136.792716839935565
	pv[0][2] = -245251.2339876830229

	pv[1][0] = -0.4051854035740713039e-2
	pv[1][1] = -0.6253919754866175788e-2
	pv[1][2] = 0.1189353719774107615e-1

	gofa.Pvup(2920.0, pv, &p)

	vvd(t, p[0], 126656.7598605317105, 1e-6, "Pvup", "1")
	vvd(t, p[1], 2118.531271155726332, 1e-8, "Pvup", "2")
	vvd(t, p[2], -245216.5048590656190, 1e-6, "Pvup", "3")
}

func TestRxr(t *testing.T) {
	var a, b, atb [3][3]float64

	a[0][0] = 2.0
	a[0][1] = 3.0
	a[0][2] = 2.0

	a[1][0] = 3.0
	a[1][1] = 2.0
	a[1][2] = 3.0

	a[2][0] = 3.0
	a[2][1] = 4.0
	a[2][2] = 5.0

	b[0][0] = 1.0
	b[0][1] = 2.0
	b[0][2] = 2.0

	b[1][0] = 4.0
	b[1][1] = 1.0
	b[1][2] = 1.0

	b[2][0] = 3.0
	b[2][1] = 0.0
	b[2][2] = 1.0

	gofa.Rxr(a, b, &atb)

	vvd(t, atb[0][0], 20.0, 1e-12, "Rxr", "11")
	vvd(t, atb[0][1], 7.0, 1e-12, "Rxr", "12")
	vvd(t, atb[0][2], 9.0, 1e-12, "Rxr", "13")

	vvd(t, atb[1][0], 20.0, 1e-12, "Rxr", "21")
	vvd(t, atb[1][1], 8.0, 1e-12, "Rxr", "22")
	vvd(t, atb[1][2], 11.0, 1e-12, "Rxr", "23")

	vvd(t, atb[2][0], 34.0, 1e-12, "Rxr", "31")
	vvd(t, atb[2][1], 10.0, 1e-12, "Rxr", "32")
	vvd(t, atb[2][2], 15.0, 1e-12, "Rxr", "33")
}

func TestTr(t *testing.T) {
	var r, rt [3][3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	gofa.Tr(r, &rt)

	vvd(t, rt[0][0], 2.0, 0.0, "Tr", "11")
	vvd(t, rt[0][1], 3.0, 0.0, "Tr", "12")
	vvd(t, rt[0][2], 3.0, 0.0, "Tr", "13")

	vvd(t, rt[1][0], 3.0, 0.0, "Tr", "21")
	vvd(t, rt[1][1], 2.0, 0.0, "Tr", "22")
	vvd(t, rt[1][2], 4.0, 0.0, "Tr", "23")

	vvd(t, rt[2][0], 2.0, 0.0, "Tr", "31")
	vvd(t, rt[2][1], 3.0, 0.0, "Tr", "32")
	vvd(t, rt[2][2], 5.0, 0.0, "Tr", "33")
}

func TestRxp(t *testing.T) {
	var r [3][3]float64
	var p, rp [3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	p[0] = 0.2
	p[1] = 1.5
	p[2] = 0.1

	gofa.Rxp(r, p, &rp)

	vvd(t, rp[0], 5.1, 1e-12, "iauRxp", "1")
	vvd(t, rp[1], 3.9, 1e-12, "iauRxp", "2")
	vvd(t, rp[2], 7.1, 1e-12, "iauRxp", "3")
}

func TestTrxp(t *testing.T) {
	var r [3][3]float64
	var p, trp [3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	p[0] = 0.2
	p[1] = 1.5
	p[2] = 0.1

	gofa.Trxp(r, p, &trp)

	vvd(t, trp[0], 5.2, 1e-12, "iauTrxp", "1")
	vvd(t, trp[1], 4.0, 1e-12, "iauTrxp", "2")
	vvd(t, trp[2], 5.4, 1e-12, "iauTrxp", "3")
}

func TestRxpv(t *testing.T) {
	var r [3][3]float64
	var pv, rpv [2][3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	pv[0][0] = 0.2
	pv[0][1] = 1.5
	pv[0][2] = 0.1

	pv[1][0] = 1.5
	pv[1][1] = 0.2
	pv[1][2] = 0.1

	gofa.Rxpv(r, pv, &rpv)

	vvd(t, rpv[0][0], 5.1, 1e-12, "iauRxpv", "11")
	vvd(t, rpv[1][0], 3.8, 1e-12, "iauRxpv", "12")

	vvd(t, rpv[0][1], 3.9, 1e-12, "iauRxpv", "21")
	vvd(t, rpv[1][1], 5.2, 1e-12, "iauRxpv", "22")

	vvd(t, rpv[0][2], 7.1, 1e-12, "iauRxpv", "31")
	vvd(t, rpv[1][2], 5.8, 1e-12, "iauRxpv", "32")
}

func TestTrxpv(t *testing.T) {
	var r [3][3]float64
	var pv, trpv [2][3]float64

	r[0][0] = 2.0
	r[0][1] = 3.0
	r[0][2] = 2.0

	r[1][0] = 3.0
	r[1][1] = 2.0
	r[1][2] = 3.0

	r[2][0] = 3.0
	r[2][1] = 4.0
	r[2][2] = 5.0

	pv[0][0] = 0.2
	pv[0][1] = 1.5
	pv[0][2] = 0.1

	pv[1][0] = 1.5
	pv[1][1] = 0.2
	pv[1][2] = 0.1

	gofa.Trxpv(r, pv, &trpv)

	vvd(t, trpv[0][0], 5.2, 1e-12, "iauTrxpv", "p1")
	vvd(t, trpv[0][1], 4.0, 1e-12, "iauTrxpv", "p1")
	vvd(t, trpv[0][2], 5.4, 1e-12, "iauTrxpv", "p1")

	vvd(t, trpv[1][0], 3.9, 1e-12, "iauTrxpv", "v1")
	vvd(t, trpv[1][1], 5.3, 1e-12, "iauTrxpv", "v2")
	vvd(t, trpv[1][2], 4.1, 1e-12, "iauTrxpv", "v3")
}

func TestRv2m(t *testing.T) {
	var w [3]float64
	var r [3][3]float64

	w[0] = 0.0
	w[1] = 1.41371669
	w[2] = -1.88495559

	gofa.Rv2m(w, &r)

	vvd(t, r[0][0], -0.7071067782221119905, 1e-14, "iauRv2m", "11")
	vvd(t, r[0][1], -0.5656854276809129651, 1e-14, "iauRv2m", "12")
	vvd(t, r[0][2], -0.4242640700104211225, 1e-14, "iauRv2m", "13")

	vvd(t, r[1][0], 0.5656854276809129651, 1e-14, "iauRv2m", "21")
	vvd(t, r[1][1], -0.0925483394532274246, 1e-14, "iauRv2m", "22")
	vvd(t, r[1][2], -0.8194112531408833269, 1e-14, "iauRv2m", "23")

	vvd(t, r[2][0], 0.4242640700104211225, 1e-14, "iauRv2m", "31")
	vvd(t, r[2][1], -0.8194112531408833269, 1e-14, "iauRv2m", "32")
	vvd(t, r[2][2], 0.3854415612311154341, 1e-14, "iauRv2m", "33")
}

func TestRm2v(t *testing.T) {
	var r [3][3]float64
	var w [3]float64

	r[0][0] = 0.00
	r[0][1] = -0.80
	r[0][2] = -0.60

	r[1][0] = 0.80
	r[1][1] = -0.36
	r[1][2] = 0.48

	r[2][0] = 0.60
	r[2][1] = 0.48
	r[2][2] = -0.64

	gofa.Rm2v(r, &w)

	vvd(t, w[0], 0.0, 1e-12, "iauRm2v", "1")
	vvd(t, w[1], 1.413716694115406957, 1e-12, "iauRm2v", "2")
	vvd(t, w[2], -1.884955592153875943, 1e-12, "iauRm2v", "3")
}
