// Copyright 2022 HE Boliang
// All rights reserved.

package gofa_test

import (
	"testing"

	"github.com/hebl/gofa"
)

func TestAb(t *testing.T) {
	var pnat, v, ppr [3]float64
	var s, bm1 float64

	pnat[0] = -0.76321968546737951
	pnat[1] = -0.60869453983060384
	pnat[2] = -0.21676408580639883
	v[0] = 2.1044018893653786e-5
	v[1] = -8.9108923304429319e-5
	v[2] = -3.8633714797716569e-5
	s = 0.99980921395708788
	bm1 = 0.99999999506209258

	gofa.Ab(pnat, v, s, bm1, &ppr)

	vvd(t, ppr[0], -0.7631631094219556269, 1e-12, "Ab", "1")
	vvd(t, ppr[1], -0.6087553082505590832, 1e-12, "Ab", "2")
	vvd(t, ppr[2], -0.2167926269368471279, 1e-12, "Ab", "3")
}
func TestApcg(t *testing.T) {
	var date1, date2 float64
	var ebpv [2][3]float64
	var ehp [3]float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	ebpv[0][0] = 0.901310875
	ebpv[0][1] = -0.417402664
	ebpv[0][2] = -0.180982288
	ebpv[1][0] = 0.00742727954
	ebpv[1][1] = 0.0140507459
	ebpv[1][2] = 0.00609045792
	ehp[0] = 0.903358544
	ehp[1] = -0.415395237
	ehp[2] = -0.180084014

	gofa.Apcg(date1, date2, ebpv, ehp, &astrom)

	vvd(t, astrom.Pmt, 12.65133794027378508, 1e-11, "Apcg", "pmt")
	vvd(t, astrom.Eb[0], 0.901310875, 1e-12, "Apcg", "eb(1)")
	vvd(t, astrom.Eb[1], -0.417402664, 1e-12, "Apcg", "eb(2)")
	vvd(t, astrom.Eb[2], -0.180982288, 1e-12, "Apcg", "eb(3)")
	vvd(t, astrom.Eh[0], 0.8940025429324143045, 1e-12, "Apcg", "eh(1)")
	vvd(t, astrom.Eh[1], -0.4110930268679817955, 1e-12, "Apcg", "eh(2)")
	vvd(t, astrom.Eh[2], -0.1782189004872870264, 1e-12, "Apcg", "eh(3)")
	vvd(t, astrom.Em, 1.010465295811013146, 1e-12, "Apcg", "em")
	vvd(t, astrom.V[0], 0.4289638913597693554e-4, 1e-16, "Apcg", "v(1)")
	vvd(t, astrom.V[1], 0.8115034051581320575e-4, 1e-16, "Apcg", "v(2)")
	vvd(t, astrom.V[2], 0.3517555136380563427e-4, 1e-16, "Apcg", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999951686012981, 1e-12, "Apcg", "bm1")
	vvd(t, astrom.Bpn[0][0], 1.0, 0.0, "Apcg", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.0, 0.0, "Apcg", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.0, 0.0, "Apcg", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], 0.0, 0.0, "Apcg", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 1.0, 0.0, "Apcg", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], 0.0, 0.0, "Apcg", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], 0.0, 0.0, "Apcg", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.0, 0.0, "Apcg", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 1.0, 0.0, "Apcg", "bpn(3,3)")
}

func TestApcg13(t *testing.T) {
	var date1, date2 float64

	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Apcg13(date1, date2, &astrom)

	vvd(t, astrom.Pmt, 12.65133794027378508, 1e-11, "Apcg13", "pmt")
	vvd(t, astrom.Eb[0], 0.9013108747340644755, 1e-12, "Apcg13", "eb(1)")
	vvd(t, astrom.Eb[1], -0.4174026640406119957, 1e-12, "Apcg13", "eb(2)")
	vvd(t, astrom.Eb[2], -0.1809822877867817771, 1e-12, "Apcg13", "eb(3)")
	vvd(t, astrom.Eh[0], 0.8940025429255499549, 1e-12, "Apcg13", "eh(1)")
	vvd(t, astrom.Eh[1], -0.4110930268331896318, 1e-12, "Apcg13", "eh(2)")
	vvd(t, astrom.Eh[2], -0.1782189006019749850, 1e-12, "Apcg13", "eh(3)")
	vvd(t, astrom.Em, 1.010465295964664178, 1e-12, "Apcg13", "em")
	vvd(t, astrom.V[0], 0.4289638912941341125e-4, 1e-16, "Apcg13", "v(1)")
	vvd(t, astrom.V[1], 0.8115034032405042132e-4, 1e-16, "Apcg13", "v(2)")
	vvd(t, astrom.V[2], 0.3517555135536470279e-4, 1e-16, "Apcg13", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999951686013142, 1e-12, "Apcg13", "bm1")
	vvd(t, astrom.Bpn[0][0], 1.0, 0.0, "Apcg13", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.0, 0.0, "Apcg13", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.0, 0.0, "Apcg13", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], 0.0, 0.0, "Apcg13", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 1.0, 0.0, "Apcg13", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], 0.0, 0.0, "Apcg13", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], 0.0, 0.0, "Apcg13", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.0, 0.0, "Apcg13", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 1.0, 0.0, "Apcg13", "bpn(3,3)")
}

func TestApci(t *testing.T) {
	var date1, date2 float64
	var ebpv [2][3]float64
	var ehp [3]float64
	var x, y, s float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	ebpv[0][0] = 0.901310875
	ebpv[0][1] = -0.417402664
	ebpv[0][2] = -0.180982288
	ebpv[1][0] = 0.00742727954
	ebpv[1][1] = 0.0140507459
	ebpv[1][2] = 0.00609045792
	ehp[0] = 0.903358544
	ehp[1] = -0.415395237
	ehp[2] = -0.180084014
	x = 0.0013122272
	y = -2.92808623e-5
	s = 3.05749468e-8

	gofa.Apci(date1, date2, ebpv, ehp, x, y, s, &astrom)

	vvd(t, astrom.Pmt, 12.65133794027378508, 1e-11, "Apci", "pmt")
	vvd(t, astrom.Eb[0], 0.901310875, 1e-12, "Apci", "eb(1)")
	vvd(t, astrom.Eb[1], -0.417402664, 1e-12, "Apci", "eb(2)")
	vvd(t, astrom.Eb[2], -0.180982288, 1e-12, "Apci", "eb(3)")
	vvd(t, astrom.Eh[0], 0.8940025429324143045, 1e-12, "Apci", "eh(1)")
	vvd(t, astrom.Eh[1], -0.4110930268679817955, 1e-12, "Apci", "eh(2)")
	vvd(t, astrom.Eh[2], -0.1782189004872870264, 1e-12, "Apci", "eh(3)")
	vvd(t, astrom.Em, 1.010465295811013146, 1e-12, "Apci", "em")
	vvd(t, astrom.V[0], 0.4289638913597693554e-4, 1e-16, "Apci", "v(1)")
	vvd(t, astrom.V[1], 0.8115034051581320575e-4, 1e-16, "Apci", "v(2)")
	vvd(t, astrom.V[2], 0.3517555136380563427e-4, 1e-16, "Apci", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999951686012981, 1e-12, "Apci", "bm1")
	vvd(t, astrom.Bpn[0][0], 0.9999991390295159156, 1e-12, "Apci", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.4978650072505016932e-7, 1e-12, "Apci", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.1312227200000000000e-2, 1e-12, "Apci", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], -0.1136336653771609630e-7, 1e-12, "Apci", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 0.9999999995713154868, 1e-12, "Apci", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], -0.2928086230000000000e-4, 1e-12, "Apci", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], -0.1312227200895260194e-2, 1e-12, "Apci", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.2928082217872315680e-4, 1e-12, "Apci", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 0.9999991386008323373, 1e-12, "Apci", "bpn(3,3)")
}

func TestApci13(t *testing.T) {
	var date1, date2 float64
	var eo float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Apci13(date1, date2, &astrom, &eo)

	vvd(t, astrom.Pmt, 12.65133794027378508, 1e-11, "Apci13", "pmt")
	vvd(t, astrom.Eb[0], 0.9013108747340644755, 1e-12, "Apci13", "eb(1)")
	vvd(t, astrom.Eb[1], -0.4174026640406119957, 1e-12, "Apci13", "eb(2)")
	vvd(t, astrom.Eb[2], -0.1809822877867817771, 1e-12, "Apci13", "eb(3)")
	vvd(t, astrom.Eh[0], 0.8940025429255499549, 1e-12, "Apci13", "eh(1)")
	vvd(t, astrom.Eh[1], -0.4110930268331896318, 1e-12, "Apci13", "eh(2)")
	vvd(t, astrom.Eh[2], -0.1782189006019749850, 1e-12, "Apci13", "eh(3)")
	vvd(t, astrom.Em, 1.010465295964664178, 1e-12, "Apci13", "em")
	vvd(t, astrom.V[0], 0.4289638912941341125e-4, 1e-16, "Apci13", "v(1)")
	vvd(t, astrom.V[1], 0.8115034032405042132e-4, 1e-16, "Apci13", "v(2)")
	vvd(t, astrom.V[2], 0.3517555135536470279e-4, 1e-16, "Apci13", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999951686013142, 1e-12, "Apci13", "bm1")
	vvd(t, astrom.Bpn[0][0], 0.9999992060376761710, 1e-12, "Apci13", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.4124244860106037157e-7, 1e-12, "Apci13", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.1260128571051709670e-2, 1e-12, "Apci13", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], -0.1282291987222130690e-7, 1e-12, "Apci13", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 0.9999999997456835325, 1e-12, "Apci13", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], -0.2255288829420524935e-4, 1e-12, "Apci13", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], -0.1260128571661374559e-2, 1e-12, "Apci13", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.2255285422953395494e-4, 1e-12, "Apci13", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 0.9999992057833604343, 1e-12, "Apci13", "bpn(3,3)")
	vvd(t, eo, -0.2900618712657375647e-2, 1e-12, "Apci13", "eo")

}

func TestApco(t *testing.T) {
	var date1, date2 float64
	var ebpv [2][3]float64
	var ehp [3]float64
	var x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb float64
	var astrom gofa.ASTROM

	date1 = 2456384.5
	date2 = 0.970031644
	ebpv[0][0] = -0.974170438
	ebpv[0][1] = -0.211520082
	ebpv[0][2] = -0.0917583024
	ebpv[1][0] = 0.00364365824
	ebpv[1][1] = -0.0154287319
	ebpv[1][2] = -0.00668922024
	ehp[0] = -0.973458265
	ehp[1] = -0.209215307
	ehp[2] = -0.0906996477
	x = 0.0013122272
	y = -2.92808623e-5
	s = 3.05749468e-8
	theta = 3.14540971
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	sp = -3.01974337e-11
	refa = 0.000201418779
	refb = -2.36140831e-7

	gofa.Apco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb, &astrom)

	vvd(t, astrom.Pmt, 13.25248468622587269, 1e-11, "Apco", "pmt")
	vvd(t, astrom.Eb[0], -0.9741827110630322720, 1e-12, "Apco", "eb(1)")
	vvd(t, astrom.Eb[1], -0.2115130190135344832, 1e-12, "Apco", "eb(2)")
	vvd(t, astrom.Eb[2], -0.09179840186949532298, 1e-12, "Apco", "eb(3)")
	vvd(t, astrom.Eh[0], -0.9736425571689739035, 1e-12, "Apco", "eh(1)")
	vvd(t, astrom.Eh[1], -0.2092452125849330936, 1e-12, "Apco", "eh(2)")
	vvd(t, astrom.Eh[2], -0.09075578152243272599, 1e-12, "Apco", "eh(3)")
	vvd(t, astrom.Em, 0.9998233241709957653, 1e-12, "Apco", "em")
	vvd(t, astrom.V[0], 0.2078704992916728762e-4, 1e-16, "Apco", "v(1)")
	vvd(t, astrom.V[1], -0.8955360107151952319e-4, 1e-16, "Apco", "v(2)")
	vvd(t, astrom.V[2], -0.3863338994288951082e-4, 1e-16, "Apco", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999950277561236, 1e-12, "Apco", "bm1")
	vvd(t, astrom.Bpn[0][0], 0.9999991390295159156, 1e-12, "Apco", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.4978650072505016932e-7, 1e-12, "Apco", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.1312227200000000000e-2, 1e-12, "Apco", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], -0.1136336653771609630e-7, 1e-12, "Apco", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 0.9999999995713154868, 1e-12, "Apco", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], -0.2928086230000000000e-4, 1e-12, "Apco", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], -0.1312227200895260194e-2, 1e-12, "Apco", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.2928082217872315680e-4, 1e-12, "Apco", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 0.9999991386008323373, 1e-12, "Apco", "bpn(3,3)")
	vvd(t, astrom.Along, -0.5278008060295995734, 1e-12, "Apco", "along")
	vvd(t, astrom.Xpl, 0.1133427418130752958e-5, 1e-17, "Apco", "xpl")
	vvd(t, astrom.Ypl, 0.1453347595780646207e-5, 1e-17, "Apco", "ypl")
	vvd(t, astrom.Sphi, -0.9440115679003211329, 1e-12, "Apco", "sphi")
	vvd(t, astrom.Cphi, 0.3299123514971474711, 1e-12, "Apco", "cphi")
	vvd(t, astrom.Diurab, 0, 0, "Apco", "diurab")
	vvd(t, astrom.Eral, 2.617608903970400427, 1e-12, "Apco", "eral")
	vvd(t, astrom.Refa, 0.2014187790000000000e-3, 1e-15, "Apco", "refa")
	vvd(t, astrom.Refb, -0.2361408310000000000e-6, 1e-18, "Apco", "refb")
}

func TestApco13(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, eo float64
	var astrom gofa.ASTROM
	var j int

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	j = gofa.Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom, &eo)

	vvd(t, astrom.Pmt, 13.25248468622475727, 1e-11, "Apco13", "pmt")
	vvd(t, astrom.Eb[0], -0.9741827107320875162, 1e-12, "Apco13", "eb(1)")
	vvd(t, astrom.Eb[1], -0.2115130190489716682, 1e-12, "Apco13", "eb(2)")
	vvd(t, astrom.Eb[2], -0.09179840189496755339, 1e-12, "Apco13", "eb(3)")
	vvd(t, astrom.Eh[0], -0.9736425572586935247, 1e-12, "Apco13", "eh(1)")
	vvd(t, astrom.Eh[1], -0.2092452121603336166, 1e-12, "Apco13", "eh(2)")
	vvd(t, astrom.Eh[2], -0.09075578153885665295, 1e-12, "Apco13", "eh(3)")
	vvd(t, astrom.Em, 0.9998233240913898141, 1e-12, "Apco13", "em")
	vvd(t, astrom.V[0], 0.2078704994520489246e-4, 1e-16, "Apco13", "v(1)")
	vvd(t, astrom.V[1], -0.8955360133238868938e-4, 1e-16, "Apco13", "v(2)")
	vvd(t, astrom.V[2], -0.3863338993055887398e-4, 1e-16, "Apco13", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999950277561004, 1e-12, "Apco13", "bm1")
	vvd(t, astrom.Bpn[0][0], 0.9999991390295147999, 1e-12, "Apco13", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0.4978650075315529277e-7, 1e-12, "Apco13", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0.001312227200850293372, 1e-12, "Apco13", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], -0.1136336652812486604e-7, 1e-12, "Apco13", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 0.9999999995713154865, 1e-12, "Apco13", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], -0.2928086230975367296e-4, 1e-12, "Apco13", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], -0.001312227201745553566, 1e-12, "Apco13", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0.2928082218847679162e-4, 1e-12, "Apco13", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 0.9999991386008312212, 1e-12, "Apco13", "bpn(3,3)")
	vvd(t, astrom.Along, -0.5278008060295995733, 1e-12, "Apco13", "along")
	vvd(t, astrom.Xpl, 0.1133427418130752958e-5, 1e-17, "Apco13", "xpl")
	vvd(t, astrom.Ypl, 0.1453347595780646207e-5, 1e-17, "Apco13", "ypl")
	vvd(t, astrom.Sphi, -0.9440115679003211329, 1e-12, "Apco13", "sphi")
	vvd(t, astrom.Cphi, 0.3299123514971474711, 1e-12, "Apco13", "cphi")
	vvd(t, astrom.Diurab, 0, 0, "Apco13", "diurab")
	vvd(t, astrom.Eral, 2.617608909189664000, 1e-12, "Apco13", "eral")
	vvd(t, astrom.Refa, 0.2014187785940396921e-3, 1e-15, "Apco13", "refa")
	vvd(t, astrom.Refb, -0.2361408314943696227e-6, 1e-18, "Apco13", "refb")
	vvd(t, eo, -0.003020548354802412839, 1e-14, "Apco13", "eo")
	viv(t, j, 0, "Apco13", "j")

}

func TestApcs(t *testing.T) {
	var date1, date2 float64
	var pv, ebpv [2][3]float64
	var ehp [3]float64
	var astrom gofa.ASTROM

	date1 = 2456384.5
	date2 = 0.970031644
	pv[0][0] = -1836024.09
	pv[0][1] = 1056607.72
	pv[0][2] = -5998795.26
	pv[1][0] = -77.0361767
	pv[1][1] = -133.310856
	pv[1][2] = 0.0971855934
	ebpv[0][0] = -0.974170438
	ebpv[0][1] = -0.211520082
	ebpv[0][2] = -0.0917583024
	ebpv[1][0] = 0.00364365824
	ebpv[1][1] = -0.0154287319
	ebpv[1][2] = -0.00668922024
	ehp[0] = -0.973458265
	ehp[1] = -0.209215307
	ehp[2] = -0.0906996477

	gofa.Apcs(date1, date2, pv, ebpv, ehp, &astrom)

	vvd(t, astrom.Pmt, 13.25248468622587269, 1e-11, "Apcs", "pmt")
	vvd(t, astrom.Eb[0], -0.9741827110629881886, 1e-12, "Apcs", "eb(1)")
	vvd(t, astrom.Eb[1], -0.2115130190136415986, 1e-12, "Apcs", "eb(2)")
	vvd(t, astrom.Eb[2], -0.09179840186954412099, 1e-12, "Apcs", "eb(3)")
	vvd(t, astrom.Eh[0], -0.9736425571689454706, 1e-12, "Apcs", "eh(1)")
	vvd(t, astrom.Eh[1], -0.2092452125850435930, 1e-12, "Apcs", "eh(2)")
	vvd(t, astrom.Eh[2], -0.09075578152248299218, 1e-12, "Apcs", "eh(3)")
	vvd(t, astrom.Em, 0.9998233241709796859, 1e-12, "Apcs", "em")
	vvd(t, astrom.V[0], 0.2078704993282685510e-4, 1e-16, "Apcs", "v(1)")
	vvd(t, astrom.V[1], -0.8955360106989405683e-4, 1e-16, "Apcs", "v(2)")
	vvd(t, astrom.V[2], -0.3863338994289409097e-4, 1e-16, "Apcs", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999950277561237, 1e-12, "Apcs", "bm1")
	vvd(t, astrom.Bpn[0][0], 1, 0, "Apcs", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0, 0, "Apcs", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0, 0, "Apcs", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], 0, 0, "Apcs", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 1, 0, "Apcs", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], 0, 0, "Apcs", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], 0, 0, "Apcs", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0, 0, "Apcs", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 1, 0, "Apcs", "bpn(3,3)")
}

func TestApcs13(t *testing.T) {
	var date1, date2 float64
	var pv [2][3]float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	pv[0][0] = -6241497.16
	pv[0][1] = 401346.896
	pv[0][2] = -1251136.04
	pv[1][0] = -29.264597
	pv[1][1] = -455.021831
	pv[1][2] = 0.0266151194

	gofa.Apcs13(date1, date2, pv, &astrom)

	vvd(t, astrom.Pmt, 12.65133794027378508, 1e-11, "Apcs13", "pmt")
	vvd(t, astrom.Eb[0], 0.9012691529025250644, 1e-12, "Apcs13", "eb(1)")
	vvd(t, astrom.Eb[1], -0.4173999812023194317, 1e-12, "Apcs13", "eb(2)")
	vvd(t, astrom.Eb[2], -0.1809906511146429670, 1e-12, "Apcs13", "eb(3)")
	vvd(t, astrom.Eh[0], 0.8939939101760130792, 1e-12, "Apcs13", "eh(1)")
	vvd(t, astrom.Eh[1], -0.4111053891734021478, 1e-12, "Apcs13", "eh(2)")
	vvd(t, astrom.Eh[2], -0.1782336880636997374, 1e-12, "Apcs13", "eh(3)")
	vvd(t, astrom.Em, 1.010428384373491095, 1e-12, "Apcs13", "em")
	vvd(t, astrom.V[0], 0.4279877294121697570e-4, 1e-16, "Apcs13", "v(1)")
	vvd(t, astrom.V[1], 0.7963255087052120678e-4, 1e-16, "Apcs13", "v(2)")
	vvd(t, astrom.V[2], 0.3517564013384691531e-4, 1e-16, "Apcs13", "v(3)")
	vvd(t, astrom.Bm1, 0.9999999952947980978, 1e-12, "Apcs13", "bm1")
	vvd(t, astrom.Bpn[0][0], 1, 0, "Apcs13", "bpn(1,1)")
	vvd(t, astrom.Bpn[1][0], 0, 0, "Apcs13", "bpn(2,1)")
	vvd(t, astrom.Bpn[2][0], 0, 0, "Apcs13", "bpn(3,1)")
	vvd(t, astrom.Bpn[0][1], 0, 0, "Apcs13", "bpn(1,2)")
	vvd(t, astrom.Bpn[1][1], 1, 0, "Apcs13", "bpn(2,2)")
	vvd(t, astrom.Bpn[2][1], 0, 0, "Apcs13", "bpn(3,2)")
	vvd(t, astrom.Bpn[0][2], 0, 0, "Apcs13", "bpn(1,3)")
	vvd(t, astrom.Bpn[1][2], 0, 0, "Apcs13", "bpn(2,3)")
	vvd(t, astrom.Bpn[2][2], 1, 0, "Apcs13", "bpn(3,3)")
}

func TestAper(t *testing.T) {
	var theta float64
	var astrom gofa.ASTROM

	astrom.Along = 1.234
	theta = 5.678

	gofa.Aper(theta, &astrom)

	vvd(t, astrom.Eral, 6.912000000000000000, 1e-12, "Aper", "pmt")
}

func TestAper13(t *testing.T) {
	var ut11, ut12 float64
	var astrom gofa.ASTROM

	astrom.Along = 1.234
	ut11 = 2456165.5
	ut12 = 0.401182685

	gofa.Aper13(ut11, ut12, &astrom)

	vvd(t, astrom.Eral, 3.316236661789694933, 1e-12, "Aper13", "pmt")
}

func TestApio(t *testing.T) {
	var sp, theta, elong, phi, hm, xp, yp, refa, refb float64
	var astrom gofa.ASTROM

	sp = -3.01974337e-11
	theta = 3.14540971
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	refa = 0.000201418779
	refb = -2.36140831e-7

	gofa.Apio(sp, theta, elong, phi, hm, xp, yp, refa, refb, &astrom)

	vvd(t, astrom.Along, -0.5278008060295995734, 1e-12, "Apio", "along")
	vvd(t, astrom.Xpl, 0.1133427418130752958e-5, 1e-17, "Apio", "xpl")
	vvd(t, astrom.Ypl, 0.1453347595780646207e-5, 1e-17, "Apio", "ypl")
	vvd(t, astrom.Sphi, -0.9440115679003211329, 1e-12, "Apio", "sphi")
	vvd(t, astrom.Cphi, 0.3299123514971474711, 1e-12, "Apio", "cphi")
	vvd(t, astrom.Diurab, 0.5135843661699913529e-6, 1e-12, "Apio", "diurab")
	vvd(t, astrom.Eral, 2.617608903970400427, 1e-12, "Apio", "eral")
	vvd(t, astrom.Refa, 0.2014187790000000000e-3, 1e-15, "Apio", "refa")
	vvd(t, astrom.Refb, -0.2361408310000000000e-6, 1e-18, "Apio", "refb")
}

func TestApio13(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl float64
	var j int
	var astrom gofa.ASTROM

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	j = gofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom)

	vvd(t, astrom.Along, -0.5278008060295995733, 1e-12, "Apio13", "along")
	vvd(t, astrom.Xpl, 0.1133427418130752958e-5, 1e-17, "Apio13", "xpl")
	vvd(t, astrom.Ypl, 0.1453347595780646207e-5, 1e-17, "Apio13", "ypl")
	vvd(t, astrom.Sphi, -0.9440115679003211329, 1e-12, "Apio13", "sphi")
	vvd(t, astrom.Cphi, 0.3299123514971474711, 1e-12, "Apio13", "cphi")
	vvd(t, astrom.Diurab, 0.5135843661699913529e-6, 1e-12, "Apio13", "diurab")
	vvd(t, astrom.Eral, 2.617608909189664000, 1e-12, "Apio13", "eral")
	vvd(t, astrom.Refa, 0.2014187785940396921e-3, 1e-15, "Apio13", "refa")
	vvd(t, astrom.Refb, -0.2361408314943696227e-6, 1e-18, "Apio13", "refb")
	viv(t, j, 0, "Apio13", "j")
}

func TestAtcc13(t *testing.T) {
	var rc, dc, pr, pd, px, rv, date1, date2, ra, da float64

	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0
	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Atcc13(rc, dc, pr, pd, px, rv, date1, date2, &ra, &da)

	vvd(t, ra, 2.710126504531372384, 1e-12, "Atcc13", "ra")
	vvd(t, da, 0.1740632537628350152, 1e-12, "Atcc13", "da")
}

func TestAtccq(t *testing.T) {
	var date1, date2 float64
	var eo, rc, dc, pr, pd, px, rv, ra, da float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0

	gofa.Atccq(rc, dc, pr, pd, px, rv, &astrom, &ra, &da)

	vvd(t, ra, 2.710126504531372384, 1e-12, "Atccq", "ra")
	vvd(t, da, 0.1740632537628350152, 1e-12, "Atccq", "da")
}

func TestAtci13(t *testing.T) {
	var rc, dc, pr, pd, px, rv, date1, date2, ri, di, eo float64

	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0
	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Atci13(rc, dc, pr, pd, px, rv, date1, date2, &ri, &di, &eo)

	vvd(t, ri, 2.710121572968696744, 1e-12, "Atci13", "ri")
	vvd(t, di, 0.1729371367219539137, 1e-12, "Atci13", "di")
	vvd(t, eo, -0.002900618712657375647, 1e-14, "Atci13", "eo")
}

func TestAtciq(t *testing.T) {
	var date1, date2 float64
	var eo, rc, dc, pr, pd, px, rv, ri, di float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0

	gofa.Atciq(rc, dc, pr, pd, px, rv, &astrom, &ri, &di)

	vvd(t, ri, 2.710121572968696744, 1e-12, "Atciq", "ri")
	vvd(t, di, 0.1729371367219539137, 1e-12, "Atciq", "di")
}

func TestAtciqn(t *testing.T) {
	// var b [3]gofa.LDBODY
	b := make([]gofa.LDBODY, 3)
	var date1, date2 float64
	var eo, rc, dc, pr, pd, px, rv, ri, di float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0
	b[0].Bm = 0.00028574
	b[0].Dl = 3e-10
	b[0].Pv[0][0] = -7.81014427
	b[0].Pv[0][1] = -5.60956681
	b[0].Pv[0][2] = -1.98079819
	b[0].Pv[1][0] = 0.0030723249
	b[0].Pv[1][1] = -0.00406995477
	b[0].Pv[1][2] = -0.00181335842
	b[1].Bm = 0.00095435
	b[1].Dl = 3e-9
	b[1].Pv[0][0] = 0.738098796
	b[1].Pv[0][1] = 4.63658692
	b[1].Pv[0][2] = 1.9693136
	b[1].Pv[1][0] = -0.00755816922
	b[1].Pv[1][1] = 0.00126913722
	b[1].Pv[1][2] = 0.000727999001
	b[2].Bm = 1.0
	b[2].Dl = 6e-6
	b[2].Pv[0][0] = -0.000712174377
	b[2].Pv[0][1] = -0.00230478303
	b[2].Pv[0][2] = -0.00105865966
	b[2].Pv[1][0] = 6.29235213e-6
	b[2].Pv[1][1] = -3.30888387e-7
	b[2].Pv[1][2] = -2.96486623e-7

	gofa.Atciqn(rc, dc, pr, pd, px, rv, &astrom, 3, b, &ri, &di)

	vvd(t, ri, 2.710122008104983335, 1e-12, "Atciqn", "ri")
	vvd(t, di, 0.1729371916492767821, 1e-12, "Atciqn", "di")

}

func TestAtciqz(t *testing.T) {
	var date1, date2 float64
	var eo, rc, dc, ri, di float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	rc = 2.71
	dc = 0.174

	gofa.Atciqz(rc, dc, &astrom, &ri, &di)

	vvd(t, ri, 2.709994899247256984, 1e-12, "Atciqz", "ri")
	vvd(t, di, 0.1728740720984931891, 1e-12, "Atciqz", "di")
}

func TestAtco13(t *testing.T) {
	var rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, aob, zob, hob, dob, rob, eo float64
	var j int

	rc = 2.71
	dc = 0.174
	pr = 1e-5
	pd = 5e-6
	px = 0.1
	rv = 55.0
	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	j = gofa.Atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &aob, &zob, &hob, &dob, &rob, &eo)

	vvd(t, aob, 0.9251774485485515207e-1, 1e-12, "Atco13", "aob")
	vvd(t, zob, 1.407661405256499357, 1e-12, "Atco13", "zob")
	vvd(t, hob, -0.9265154431529724692e-1, 1e-12, "Atco13", "hob")
	vvd(t, dob, 0.1716626560072526200, 1e-12, "Atco13", "dob")
	vvd(t, rob, 2.710260453504961012, 1e-12, "Atco13", "rob")
	vvd(t, eo, -0.003020548354802412839, 1e-14, "Atco13", "eo")
	viv(t, j, 0, "Atco13", "j")
}

func TestAtic13(t *testing.T) {
	var ri, di, date1, date2, rc, dc, eo float64

	ri = 2.710121572969038991
	di = 0.1729371367218230438
	date1 = 2456165.5
	date2 = 0.401182685

	gofa.Atic13(ri, di, date1, date2, &rc, &dc, &eo)

	vvd(t, rc, 2.710126504531716819, 1e-12, "Atic13", "rc")
	vvd(t, dc, 0.1740632537627034482, 1e-12, "Atic13", "dc")
	vvd(t, eo, -0.002900618712657375647, 1e-14, "Atic13", "eo")
}

func TestAticq(t *testing.T) {
	var date1, date2 float64
	var eo, ri, di, rc, dc float64
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	ri = 2.710121572969038991
	di = 0.1729371367218230438

	gofa.Aticq(ri, di, &astrom, &rc, &dc)

	vvd(t, rc, 2.710126504531716819, 1e-12, "Aticq", "rc")
	vvd(t, dc, 0.1740632537627034482, 1e-12, "Aticq", "dc")
}

func TestAticqn(t *testing.T) {
	var date1, date2 float64
	var eo, ri, di, rc, dc float64
	b := make([]gofa.LDBODY, 3)
	var astrom gofa.ASTROM

	date1 = 2456165.5
	date2 = 0.401182685
	gofa.Apci13(date1, date2, &astrom, &eo)
	ri = 2.709994899247599271
	di = 0.1728740720983623469
	b[0].Bm = 0.00028574
	b[0].Dl = 3e-10
	b[0].Pv[0][0] = -7.81014427
	b[0].Pv[0][1] = -5.60956681
	b[0].Pv[0][2] = -1.98079819
	b[0].Pv[1][0] = 0.0030723249
	b[0].Pv[1][1] = -0.00406995477
	b[0].Pv[1][2] = -0.00181335842
	b[1].Bm = 0.00095435
	b[1].Dl = 3e-9
	b[1].Pv[0][0] = 0.738098796
	b[1].Pv[0][1] = 4.63658692
	b[1].Pv[0][2] = 1.9693136
	b[1].Pv[1][0] = -0.00755816922
	b[1].Pv[1][1] = 0.00126913722
	b[1].Pv[1][2] = 0.000727999001
	b[2].Bm = 1.0
	b[2].Dl = 6e-6
	b[2].Pv[0][0] = -0.000712174377
	b[2].Pv[0][1] = -0.00230478303
	b[2].Pv[0][2] = -0.00105865966
	b[2].Pv[1][0] = 6.29235213e-6
	b[2].Pv[1][1] = -3.30888387e-7
	b[2].Pv[1][2] = -2.96486623e-7

	gofa.Aticqn(ri, di, &astrom, 3, b, &rc, &dc)

	vvd(t, rc, 2.709999575033027333, 1e-12, "Atciqn", "rc")
	vvd(t, dc, 0.1739999656316469990, 1e-12, "Atciqn", "dc")
}

func TestAtio13(t *testing.T) {
	var ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, aob, zob, hob, dob, rob float64
	var j int

	ri = 2.710121572969038991
	di = 0.1729371367218230438
	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	j = gofa.Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &aob, &zob, &hob, &dob, &rob)

	vvd(t, aob, 0.9233952224895122499e-1, 1e-12, "Atio13", "aob")
	vvd(t, zob, 1.407758704513549991, 1e-12, "Atio13", "zob")
	vvd(t, hob, -0.9247619879881698140e-1, 1e-12, "Atio13", "hob")
	vvd(t, dob, 0.1717653435756234676, 1e-12, "Atio13", "dob")
	vvd(t, rob, 2.710085107988480746, 1e-12, "Atio13", "rob")
	viv(t, j, 0, "Atio13", "j")
}

func TestAtioq(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ri, di, aob, zob, hob, dob, rob float64
	var astrom gofa.ASTROM

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55
	gofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom)
	ri = 2.710121572969038991
	di = 0.1729371367218230438

	gofa.Atioq(ri, di, &astrom, &aob, &zob, &hob, &dob, &rob)

	vvd(t, aob, 0.9233952224895122499e-1, 1e-12, "Atioq", "aob")
	vvd(t, zob, 1.407758704513549991, 1e-12, "Atioq", "zob")
	vvd(t, hob, -0.9247619879881698140e-1, 1e-12, "Atioq", "hob")
	vvd(t, dob, 0.1717653435756234676, 1e-12, "Atioq", "dob")
	vvd(t, rob, 2.710085107988480746, 1e-12, "Atioq", "rob")
}

func TestAtoc13(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ob1, ob2, rc, dc float64
	var j int

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	ob1 = 2.710085107986886201
	ob2 = 0.1717653435758265198
	j = gofa.Atoc13("R", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &rc, &dc)
	vvd(t, rc, 2.709956744659136129, 1e-12, "Atoc13", "R/rc")
	vvd(t, dc, 0.1741696500898471362, 1e-12, "Atoc13", "R/dc")
	viv(t, j, 0, "Atoc13", "R/j")

	ob1 = -0.09247619879782006106
	ob2 = 0.1717653435758265198
	j = gofa.Atoc13("H", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &rc, &dc)
	vvd(t, rc, 2.709956744659734086, 1e-12, "Atoc13", "H/rc")
	vvd(t, dc, 0.1741696500898471362, 1e-12, "Atoc13", "H/dc")
	viv(t, j, 0, "Atoc13", "H/j")

	ob1 = 0.09233952224794989993
	ob2 = 1.407758704513722461
	j = gofa.Atoc13("A", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &rc, &dc)
	vvd(t, rc, 2.709956744659734086, 1e-12, "Atoc13", "A/rc")
	vvd(t, dc, 0.1741696500898471366, 1e-12, "Atoc13", "A/dc")
	viv(t, j, 0, "Atoc13", "A/j")

}

func TestAtoi13(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ob1, ob2, ri, di float64
	var j int

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55

	ob1 = 2.710085107986886201
	ob2 = 0.1717653435758265198
	j = gofa.Atoi13("R", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &ri, &di)
	vvd(t, ri, 2.710121574447540810, 1e-12, "Atoi13", "R/ri")
	vvd(t, di, 0.1729371839116608778, 1e-12, "Atoi13", "R/di")
	viv(t, j, 0, "Atoi13", "R/J")

	ob1 = -0.09247619879782006106
	ob2 = 0.1717653435758265198
	j = gofa.Atoi13("H", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &ri, &di)
	vvd(t, ri, 2.710121574448138676, 1e-12, "Atoi13", "H/ri")
	vvd(t, di, 0.1729371839116608778, 1e-12, "Atoi13", "H/di")
	viv(t, j, 0, "Atoi13", "H/J")

	ob1 = 0.09233952224794989993
	ob2 = 1.407758704513722461
	j = gofa.Atoi13("A", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &ri, &di)
	vvd(t, ri, 2.710121574448138676, 1e-12, "Atoi13", "A/ri")
	vvd(t, di, 0.1729371839116608781, 1e-12, "Atoi13", "A/di")
	viv(t, j, 0, "Atoi13", "A/J")
}

func TestAtoiq(t *testing.T) {
	var utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ob1, ob2, ri, di float64
	var astrom gofa.ASTROM

	utc1 = 2456384.5
	utc2 = 0.969254051
	dut1 = 0.1550675
	elong = -0.527800806
	phi = -1.2345856
	hm = 2738.0
	xp = 2.47230737e-7
	yp = 1.82640464e-6
	phpa = 731.0
	tc = 12.8
	rh = 0.59
	wl = 0.55
	gofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom)

	ob1 = 2.710085107986886201
	ob2 = 0.1717653435758265198
	gofa.Atoiq("R", ob1, ob2, &astrom, &ri, &di)
	vvd(t, ri, 2.710121574447540810, 1e-12, "Atoiq", "R/ri")
	vvd(t, di, 0.17293718391166087785, 1e-12, "Atoiq", "R/di")

	ob1 = -0.09247619879782006106
	ob2 = 0.1717653435758265198
	gofa.Atoiq("H", ob1, ob2, &astrom, &ri, &di)
	vvd(t, ri, 2.710121574448138676, 1e-12, "Atoiq", "H/ri")
	vvd(t, di, 0.1729371839116608778, 1e-12, "Atoiq", "H/di")

	ob1 = 0.09233952224794989993
	ob2 = 1.407758704513722461
	gofa.Atoiq("A", ob1, ob2, &astrom, &ri, &di)
	vvd(t, ri, 2.710121574448138676, 1e-12, "Atoiq", "A/ri")
	vvd(t, di, 0.1729371839116608781, 1e-12, "Atoiq", "A/di")
}

func TestLd(t *testing.T) {
	var bm, em, dlim float64
	var p, q, e, p1 [3]float64

	bm = 0.00028574
	p[0] = -0.763276255
	p[1] = -0.608633767
	p[2] = -0.216735543
	q[0] = -0.763276255
	q[1] = -0.608633767
	q[2] = -0.216735543
	e[0] = 0.76700421
	e[1] = 0.605629598
	e[2] = 0.211937094
	em = 8.91276983
	dlim = 3e-10

	gofa.Ld(bm, p, q, e, em, dlim, &p1)

	vvd(t, p1[0], -0.7632762548968159627, 1e-12, "Ld", "1")
	vvd(t, p1[1], -0.6086337670823762701, 1e-12, "Ld", "2")
	vvd(t, p1[2], -0.2167355431320546947, 1e-12, "Ld", "3")
}

func TestLdn(t *testing.T) {
	var n int
	b := make([]gofa.LDBODY, 3)
	var ob, sc, sn [3]float64

	n = 3
	b[0].Bm = 0.00028574
	b[0].Dl = 3e-10
	b[0].Pv[0][0] = -7.81014427
	b[0].Pv[0][1] = -5.60956681
	b[0].Pv[0][2] = -1.98079819
	b[0].Pv[1][0] = 0.0030723249
	b[0].Pv[1][1] = -0.00406995477
	b[0].Pv[1][2] = -0.00181335842
	b[1].Bm = 0.00095435
	b[1].Dl = 3e-9
	b[1].Pv[0][0] = 0.738098796
	b[1].Pv[0][1] = 4.63658692
	b[1].Pv[0][2] = 1.9693136
	b[1].Pv[1][0] = -0.00755816922
	b[1].Pv[1][1] = 0.00126913722
	b[1].Pv[1][2] = 0.000727999001
	b[2].Bm = 1.0
	b[2].Dl = 6e-6
	b[2].Pv[0][0] = -0.000712174377
	b[2].Pv[0][1] = -0.00230478303
	b[2].Pv[0][2] = -0.00105865966
	b[2].Pv[1][0] = 6.29235213e-6
	b[2].Pv[1][1] = -3.30888387e-7
	b[2].Pv[1][2] = -2.96486623e-7
	ob[0] = -0.974170437
	ob[1] = -0.2115201
	ob[2] = -0.0917583114
	sc[0] = -0.763276255
	sc[1] = -0.608633767
	sc[2] = -0.216735543

	gofa.Ldn(n, b, ob, sc, &sn)

	vvd(t, sn[0], -0.7632762579693333866, 1e-12, "Ldn", "1")
	vvd(t, sn[1], -0.6086337636093002660, 1e-12, "Ldn", "2")
	vvd(t, sn[2], -0.2167355420646328159, 1e-12, "Ldn", "3")
}

func TestLdsun(t *testing.T) {
	var p, e, p1 [3]float64
	var em float64

	p[0] = -0.763276255
	p[1] = -0.608633767
	p[2] = -0.216735543
	e[0] = -0.973644023
	e[1] = -0.20925523
	e[2] = -0.0907169552
	em = 0.999809214

	gofa.Ldsun(p, e, em, &p1)

	vvd(t, p1[0], -0.7632762580731413169, 1e-12, "Ldsun", "1")
	vvd(t, p1[1], -0.6086337635262647900, 1e-12, "Ldsun", "2")
	vvd(t, p1[2], -0.2167355419322321302, 1e-12, "Ldsun", "3")
}

func TestPmpx(t *testing.T) {
	var rc, dc, pr, pd, px, rv, pmt float64
	var pob, pco [3]float64

	rc = 1.234
	dc = 0.789
	pr = 1e-5
	pd = -2e-5
	px = 1e-2
	rv = 10.0
	pmt = 8.75
	pob[0] = 0.9
	pob[1] = 0.4
	pob[2] = 0.1

	gofa.Pmpx(rc, dc, pr, pd, px, rv, pmt, pob, &pco)

	vvd(t, pco[0], 0.2328137623960308438, 1e-12, "Pmpx", "1")
	vvd(t, pco[1], 0.6651097085397855328, 1e-12, "Pmpx", "2")
	vvd(t, pco[2], 0.7095257765896359837, 1e-12, "Pmpx", "3")
}

func TestPmsafe(t *testing.T) {
	var j int
	var ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b, ra2, dec2, pmr2, pmd2, px2, rv2 float64

	ra1 = 1.234
	dec1 = 0.789
	pmr1 = 1e-5
	pmd1 = -2e-5
	px1 = 1e-2
	rv1 = 10.0
	ep1a = 2400000.5
	ep1b = 48348.5625
	ep2a = 2400000.5
	ep2b = 51544.5

	j = gofa.Pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b, &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2)

	vvd(t, ra2, 1.234087484501017061, 1e-12, "Pmsafe", "ra2")
	vvd(t, dec2, 0.7888249982450468567, 1e-12, "Pmsafe", "dec2")
	vvd(t, pmr2, 0.9996457663586073988e-5, 1e-12, "Pmsafe", "pmr2")
	vvd(t, pmd2, -0.2000040085106754565e-4, 1e-16, "Pmsafe", "pmd2")
	vvd(t, px2, 0.9999997295356830666e-2, 1e-12, "Pmsafe", "px2")
	vvd(t, rv2, 10.38468380293920069, 1e-10, "Pmsafe", "rv2")
	viv(t, j, 0, "Pmsafe", "j")
}

func TestPvstar(t *testing.T) {
	var pv [2][3]float64
	var ra, dec, pmr, pmd, px, rv float64
	var j int

	pv[0][0] = 126668.5912743160601
	pv[0][1] = 2136.792716839935195
	pv[0][2] = -245251.2339876830091

	pv[1][0] = -0.4051854035740712739e-2
	pv[1][1] = -0.6253919754866173866e-2
	pv[1][2] = 0.1189353719774107189e-1

	j = gofa.Pvstar(pv, &ra, &dec, &pmr, &pmd, &px, &rv)

	vvd(t, ra, 0.1686756e-1, 1e-12, "Pvstar", "ra")
	vvd(t, dec, -1.093989828, 1e-12, "Pvstar", "dec")
	vvd(t, pmr, -0.1783235160000472788e-4, 1e-16, "Pvstar", "pmr")
	vvd(t, pmd, 0.2336024047000619347e-5, 1e-16, "Pvstar", "pmd")
	vvd(t, px, 0.74723, 1e-12, "Pvstar", "px")
	vvd(t, rv, -21.60000010107306010, 1e-11, "Pvstar", "rv")

	viv(t, j, 0, "Pvstar", "j")
}

func TestPvtob(t *testing.T) {
	var elong, phi, hm, xp, yp, sp, theta float64
	var pv [2][3]float64

	elong = 2.0
	phi = 0.5
	hm = 3000.0
	xp = 1e-6
	yp = -0.5e-6
	sp = 1e-8
	theta = 5.0

	gofa.Pvtob(elong, phi, hm, xp, yp, sp, theta, &pv)

	vvd(t, pv[0][0], 4225081.367071159207, 1e-5, "Pvtob", "p(1)")
	vvd(t, pv[0][1], 3681943.215856198144, 1e-5, "Pvtob", "p(2)")
	vvd(t, pv[0][2], 3041149.399241260785, 1e-5, "Pvtob", "p(3)")
	vvd(t, pv[1][0], -268.4915389365998787, 1e-9, "Pvtob", "v(1)")
	vvd(t, pv[1][1], 308.0977983288903123, 1e-9, "Pvtob", "v(2)")
	vvd(t, pv[1][2], 0, 0, "Pvtob", "v(3)")
}

func TestRefco(t *testing.T) {
	var phpa, tc, rh, wl, refa, refb float64

	phpa = 800.0
	tc = 10.0
	rh = 0.9
	wl = 0.4

	gofa.Refco(phpa, tc, rh, wl, &refa, &refb)

	vvd(t, refa, 0.2264949956241415009e-3, 1e-15, "Refco", "refa")
	vvd(t, refb, -0.2598658261729343970e-6, 1e-18, "Refco", "refb")

}
func TestStarpm(t *testing.T) {
	var ra1, dec1, pmr1, pmd1, px1, rv1 float64
	var ra2, dec2, pmr2, pmd2, px2, rv2 float64
	var j int

	ra1 = 0.01686756
	dec1 = -1.093989828
	pmr1 = -1.78323516e-5
	pmd1 = 2.336024047e-6
	px1 = 0.74723
	rv1 = -21.6

	j = gofa.Starpm(ra1, dec1, pmr1, pmd1, px1, rv1, 2400000.5, 50083.0, 2400000.5, 53736.0, &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2)

	vvd(t, ra2, 0.01668919069414256149, 1e-13, "Starpm", "ra")
	vvd(t, dec2, -1.093966454217127897, 1e-13, "Starpm", "dec")
	vvd(t, pmr2, -0.1783662682153176524e-4, 1e-17, "Starpm", "pmr")
	vvd(t, pmd2, 0.2338092915983989595e-5, 1e-17, "Starpm", "pmd")
	vvd(t, px2, 0.7473533835317719243, 1e-13, "Starpm", "px")
	vvd(t, rv2, -21.59905170476417175, 1e-11, "Starpm", "rv")

	viv(t, j, 0, "Starpm", "j")
}

func TestStarpv(t *testing.T) {
	var ra, dec, pmr, pmd, px, rv float64
	var pv [2][3]float64
	var j int

	ra = 0.01686756
	dec = -1.093989828
	pmr = -1.78323516e-5
	pmd = 2.336024047e-6
	px = 0.74723
	rv = -21.6

	j = gofa.Starpv(ra, dec, pmr, pmd, px, rv, &pv)

	vvd(t, pv[0][0], 126668.5912743160601, 1e-10, "Starpv", "11")
	vvd(t, pv[0][1], 2136.792716839935195, 1e-12, "Starpv", "12")
	vvd(t, pv[0][2], -245251.2339876830091, 1e-10, "Starpv", "13")

	vvd(t, pv[1][0], -0.4051854008955659551e-2, 1e-13, "Starpv", "21")
	vvd(t, pv[1][1], -0.6253919754414777970e-2, 1e-15, "Starpv", "22")
	vvd(t, pv[1][2], 0.1189353714588109341e-1, 1e-13, "Starpv", "23")

	viv(t, j, 0, "Starpv", "j")
}
