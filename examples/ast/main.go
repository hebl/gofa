// Copyright 2022 HE Boliang
// All rights reserved.

package main

//astrometry example sofa_ast_c.pdf P17.

import (
	"fmt"
	"math"

	"github.com/hebl/gofa"
)

/**
C out:

     ICRS, epoch J2000.0: 14 34 16.8118300 -12 31 10.396500
  catalog -> astrometric: 14 34 16.4960283 -12 31 02.523786
     astrometric -> CIRS: 14 34 20.2370587 -12 34 36.381654
         catalog -> CIRS: 14 34 20.2370587 -12 34 36.381654
     geocentric apparent: 14 35 01.7725802 -12 34 36.381654
     CIRS -> topocentric: 14 34 20.2570287 -12 34 36.141207
        CIRS -> observed: 14 34 16.9649101 -12 34 44.643091
        ICRS -> observed: 14 34 16.9649106 -12 34 44.643094
ICRS -> CIRS (JPL, IERS): 14 34 20.2370639 -12 34 36.381756
ICRS -> CIRS (+ planets): 14 34 20.2370658 -12 34 36.381784
     CIRS -> astrometric: 14 34 16.4960283 -12 31 02.523786

This out:

     ICRS, epoch J2000.0: 14 34 16.8118300 -12 31 10.396500
  catalog -> astrometric: 14 34 16.4960283 -12 31 02.523786
     astrometric -> CIRS: 14 34 20.2370587 -12 34 36.381654
         catalog -> CIRS: 14 34 20.2370587 -12 34 36.381654
     geocentric apparent: 14 35 01.7725802 -12 34 36.381654
     CIRS -> topocentric: 14 34 20.2570287 -12 34 36.141207
        CIRS -> observed: 14 34 16.9649101 -12 34 44.643091
        ICRS -> observed: 14 34 16.9649106 -12 34 44.643094
ICRS -> CIRS (JPL, IERS): 14 34 20.2370639 -12 34 36.381756
ICRS -> CIRS (+ planets): 14 34 20.2370658 -12 34 36.381784
     CIRS -> astrometric: 14 34 16.4960283 -12 31 02.523786
*/

const (
	DAS2R = 4.848136811095359935899141e-6 /* arcsec to radians */
)

func main() {
	var astrom gofa.ASTROM
	b := make([]gofa.LDBODY, 3)
	var phi, elong, hm, phpa, tc, rh, wl, utc1, utc2, tai1, tai2 float64
	var tt1, tt2, xp, yp, dut1, dx, dy, rc, dc, pr, pd, px, rv, rca, dca float64
	var ri, di, eo, ra, da, aot, zot, hot, dot, rot, aob, zob, hob, dob, rob float64
	var pvh, pvb [2][3]float64
	var r [3][3]float64
	var x, y, s float64

	/* Site longitude, latitude (radians) and height above the geoid (m). */
	gofa.Af2a('-', 5, 41, 54.2, &elong)
	gofa.Af2a('-', 15, 57, 42.8, &phi)
	hm = 625.0

	/* Ambient pressure (HPa), temperature (C) and rel. humidity (frac). */
	phpa = 952.0
	tc = 18.5
	rh = 0.83

	/* Effective color (microns). */
	wl = 0.55

	/* UTC date. */
	if gofa.Dtf2d("UTC", 2013, 4, 2, 23, 15, 43.55, &utc1, &utc2) != 0 {
		return
	}

	/* TT date. */
	if gofa.Utctai(utc1, utc2, &tai1, &tai2) != 0 {
		return
	}
	if gofa.Taitt(tai1, tai2, &tt1, &tt2) != 0 {
		return
	}

	/* EOPs: polar motion in radians, UT1-UTC in seconds. */
	xp = 50.995e-3 * DAS2R
	yp = 376.723e-3 * DAS2R
	dut1 = 155.0675e-3

	/* Corrections to IAU 2000A CIP (radians). */
	dx = 0.269e-3 * DAS2R
	dy = -0.274e-3 * DAS2R

	/* Star ICRS RA,Dec (radians). */
	if gofa.Tf2a(' ', 14, 34, 16.81183, &rc) != 0 {
		return
	}
	if gofa.Af2a('-', 12, 31, 10.3965, &dc) != 0 {
		return
	}
	reprd("ICRS, epoch J2000.0:", rc, dc)

	/* Annual proper motion: RA/Dec derivatives, epoch J2000.0. */
	pr = math.Atan2(-354.45e-3*DAS2R, math.Cos(dc))
	pd = 595.35e-3 * DAS2R

	/* Parallax (arcsec) and recession speed (km/s). */
	px = 164.99e-3
	rv = 0.0

	/* ICRS catalog to astrometric place... */
	gofa.Atcc13(rc, dc, pr, pd, px, rv, tt1, tt2, &rca, &dca)
	reprd("catalog -> astrometric:", rca, dca)

	/* ...then to CIRS (geocentric observer). */
	gofa.Atci13(rca, dca, 0.0, 0.0, 0.0, 0.0, tt1, tt2, &ri, &di, &eo)
	reprd("astrometric -> CIRS:", ri, di)

	/* ICRS catalog directly to CIRS (geocentric observer). */
	gofa.Atci13(rc, dc, pr, pd, px, rv, tt1, tt2, &ri, &di, &eo)
	reprd("catalog -> CIRS:", ri, di)

	/* Apparent place. */
	ra = gofa.Anp(ri - eo)
	da = di
	reprd("geocentric apparent:", ra, da)

	/* CIRS to topocentric. */
	if gofa.Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp,
		0.0, 0.0, 0.0, 0.0,
		&aot, &zot, &hot, &dot, &rot) != 0 {
		return
	}
	reprd("CIRS -> topocentric:", rot, dot)

	/* CIRS to observed. */
	if gofa.Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp,
		phpa, tc, rh, wl,
		&aob, &zob, &hob, &dob, &rob) != 0 {
		return
	}
	reprd("CIRS -> observed:", rob, dob)

	/* ICRS to observed. */
	if gofa.Atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1,
		elong, phi, hm, xp, yp, phpa, tc, rh, wl,
		&aob, &zob, &hob, &dob, &rob, &eo) != 0 {
		return
	}
	reprd("ICRS -> observed:", rob, dob)

	/* ICRS to CIRS using some user-supplied parameters. */

	/* SOFA heliocentric Earth ephemeris. */
	if gofa.Epv00(tt1, tt2, &pvh, &pvb) != 0 {
		return
	}

	/* JPL DE405 barycentric Earth ephemeris. */
	pvb[0][0] = -0.9741704366519668
	pvb[0][1] = -0.2115201000882231
	pvb[0][2] = -0.0917583114068277
	pvb[1][0] = 0.0036436589347388
	pvb[1][1] = -0.0154287318503146
	pvb[1][2] = -0.0066892203821059

	/* IAU 2000A CIP. */
	gofa.Pnm00a(tt1, tt2, &r)
	gofa.Bpn2xy(r, &x, &y)

	/* Apply IERS corrections. */
	x += dx
	y += dy

	/* SOFA CIO locator. */
	s = gofa.S06(tt1, tt2, x, y)

	/* Populate the context. */
	gofa.Apci(tt1, tt2, pvb, pvh[0], x, y, s, &astrom)

	/* Carry out the transformation and report the results. */
	gofa.Atciq(rc, dc, pr, pd, px, rv, &astrom, &ri, &di)
	reprd("ICRS -> CIRS (JPL, IERS):", ri, di)

	/* The same but with Saturn then Jupiter then Sun light deflection. */
	b[0].Bm = 0.00028574
	b[0].Dl = 3e-10
	b[0].Pv[0][0] = -7.8101442680818964
	b[0].Pv[0][1] = -5.6095668114887358
	b[0].Pv[0][2] = -1.9807981923749924
	b[0].Pv[1][0] = 0.0030723248971152
	b[0].Pv[1][1] = -0.0040699547707598
	b[0].Pv[1][2] = -0.0018133584165345

	b[1].Bm = 0.00095435
	b[1].Dl = 3e-9
	b[1].Pv[0][0] = 0.7380987962351833
	b[1].Pv[0][1] = 4.6365869247538951
	b[1].Pv[0][2] = 1.9693136030111202
	b[1].Pv[1][0] = -0.0075581692172088
	b[1].Pv[1][1] = 0.0012691372216750
	b[1].Pv[1][2] = 0.0007279990012801

	b[2].Bm = 1.0
	b[2].Dl = 6e-6
	b[2].Pv[0][0] = -0.0007121743770509
	b[2].Pv[0][1] = -0.0023047830339257
	b[2].Pv[0][2] = -0.0010586596574639
	b[2].Pv[1][0] = 0.0000062923521264
	b[2].Pv[1][1] = -0.0000003308883872
	b[2].Pv[1][2] = -0.0000002964866231

	gofa.Atciqn(rc, dc, pr, pd, px, rv, &astrom, 3, b, &ri, &di)
	reprd("ICRS -> CIRS (+ planets):", ri, di)
	/* CIRS to ICRS (astrometric). */
	gofa.Aticqn(ri, di, &astrom, 3, b, &rca, &dca)
	reprd("CIRS -> astrometric:", rca, dca)

}

func reprd(s string, ra, dc float64) {
	var pm byte
	var i [4]int
	fmt.Printf("%25s", s)
	gofa.A2tf(7, ra, &pm, &i)
	fmt.Printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3])
	gofa.A2af(6, dc, &pm, &i)
	fmt.Printf(" %c%2.2d %2.2d %2.2d.%6.6d\n", pm, i[0], i[1], i[2], i[3])
}
