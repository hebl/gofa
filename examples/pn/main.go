// Copyright 2022 HE Boliang
// All rights reserved.

package main

import (
	"fmt"
	"math"

	"github.com/hebl/gofa"
)

// gofa_pn_c.pdf examples

func printMatrix(m [3][3]float64) {
	for i := 0; i < 3; i++ {
		fmt.Printf("%19.15f  %19.15f  %19.15f\n", m[i][0], m[i][1], m[i][2])
	}
	fmt.Println(" ")
}

func fmtDeg(v float64) string {
	return fmt.Sprintf("%.15f", v*gofa.DR2D)
}

func fmtHms(v float64) string {
	var sign byte
	var ihmsf [4]int
	gofa.A2tf(9, v, &sign, &ihmsf)
	return fmt.Sprintf("%s%02dh %02dm %02d.%09ds", string(sign), ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])
}

func main() {
	var iy, im, id, ih, min int
	var sec, xp, yp, dut1,
		ddp80, dde80, dx00, dy00, dx06, dy06,
		utc1, utc2, tai1, tai2, tt1, tt2, ut11, ut12, dp80, de80,
		dpsi, deps, epsa,
		ee, gst, x, y, s,
		era, sp, dp00, de00, ddp00, dde00 float64
	var v1, v2 [3]float64
	var rp, rn, rnpb, rc2ti, rpom, rc2it, rc2i, rb, rpb [3][3]float64

	/* UTC. */
	iy = 2007
	im = 4
	id = 5
	ih = 12
	min = 0
	sec = 0.0

	/* Polar motion (arcsec->radians). */
	xp = 0.0349282 * gofa.DAS2R
	yp = 0.4833163 * gofa.DAS2R

	/* UT1-UTC (s). */
	dut1 = -0.072073685
	/* Nutation corrections wrt IAU 1976/1980 (mas->radians). */
	ddp80 = -55.0655 * gofa.DMAS2R
	dde80 = -6.3580 * gofa.DMAS2R

	/* CIP offsets wrt IAU 2000A (mas->radians). */
	dx00 = 0.1725 * gofa.DMAS2R
	dy00 = -0.2650 * gofa.DMAS2R

	/* CIP offsets wrt IAU 2006/2000A (mas->radians). */
	dx06 = 0.1750 * gofa.DMAS2R
	dy06 = -0.2259 * gofa.DMAS2R

	/* From UTC get TT and UT1. */
	gofa.Dtf2d("utc", iy, im, id, ih, min, sec, &utc1, &utc2)
	gofa.Utctai(utc1, utc2, &tai1, &tai2)
	gofa.Taitt(tai1, tai2, &tt1, &tt2)
	gofa.Utcut1(utc1, utc2, dut1, &ut11, &ut12)

	fmt.Printf("TT  = %.1f + %.15f\n", tt1, tt2)
	fmt.Printf("UT1 = %.1f + %.15f\n\n", ut11, ut12)

	/* ============= */
	/* IAU 1976/1980 */
	/* ============= */
	/* IAU 1976 precession matrix, J2000.0 to date. */
	gofa.Pmat76(tt1, tt2, &rp)

	/* IAU 1980 nutation. */
	gofa.Nut80(tt1, tt2, &dp80, &de80)

	/* Add adjustments: frame bias, precession-rates, geophysical. */
	dpsi = dp80 + ddp80
	deps = de80 + dde80
	/* Mean obliquity. */
	epsa = gofa.Obl80(tt1, tt2)

	/* Build the rotation matrix. */
	gofa.Numat(epsa, dpsi, deps, &rn)

	/* Combine the matrices: PN = N x P. */
	gofa.Rxr(rn, rp, &rnpb)

	/* Equation of the equinoxes, including nutation correction. */
	ee = gofa.Eqeq94(tt1, tt2) + ddp80*math.Cos(epsa)

	/* Greenwich apparent sidereal time (IAU 1982/1994). */
	gst = gofa.Anp(gofa.Gmst82(ut11, ut12) + ee)

	/* Form celestial-terrestrial matrix (no polar motion yet). */
	gofa.Cr(rnpb, &rc2ti)
	gofa.Rz(gst, &rc2ti)

	/* Polar motion matrix (TIRS->ITRS, IERS 1996). */
	gofa.Ir(&rpom)
	gofa.Rx(-yp, &rpom)
	gofa.Ry(-xp, &rpom)

	/* Form celestial-terrestrial matrix (including polar motion). */
	gofa.Rxr(rpom, rc2ti, &rc2it)

	//
	fmt.Printf("IAU 1976/1980\n\n")
	//
	fmt.Printf("NPB matrix, equinox based:\n")
	printMatrix(rnpb)

	fmt.Printf("GST = %s°\n    = %s\n\n", fmtDeg(gst), fmtHms(gst))

	//
	fmt.Printf("celestial to terrestrial matrix (no polar motion)\n")
	printMatrix(rc2ti)

	//
	fmt.Printf("celestial to terrestrial matrix\n")
	printMatrix(rc2it) //

	//

	/* ==================== */
	/* IAU 2000A, CIO based */
	/* ==================== */

	/* CIP and CIO, IAU 2000A. */
	gofa.Xys00a(tt1, tt2, &x, &y, &s)

	/* Add CIP corrections. */
	x = x + dx00
	y = y + dy00

	/* GCRS to CIRS matrix. */
	gofa.C2ixys(x, y, s, &rc2i)

	/* Earth rotation angle. */
	era = gofa.Era00(ut11, ut12)

	/* Form celestial-terrestrial matrix (no polar motion yet). */
	gofa.Cr(rc2i, &rc2ti)
	gofa.Rz(era, &rc2ti)

	/* Polar motion matrix (TIRS->ITRS, IERS 2003). */
	sp = gofa.Sp00(tt1, tt2)
	gofa.Pom00(xp, yp, sp, &rpom)

	/* Form celestial-terrestrial matrix (including polar motion). */
	gofa.Rxr(rpom, rc2ti, &rc2it)

	//
	fmt.Printf("IAU 2000A, CIO based\n\n")
	fmt.Printf("X = %.15f\nY = %.15f\ns = %.9f\"\n\n", x, y, s*gofa.DR2AS)

	//
	fmt.Printf("NPB matrix, CIO based\n")
	printMatrix(rc2i)

	//
	fmt.Printf("ERA = %s°\n    = %s\n\n", fmtDeg(era), fmtHms(era))

	fmt.Printf("celestial to terrestrial matrix (no polar motion)\n")
	printMatrix(rc2ti)

	fmt.Printf("celestial to terrestrial matrix\n")
	printMatrix(rc2it) //

	//

	/* ======================== */
	/* IAU 2000A, equinox based */
	/* ======================== */

	/* Nutation, IAU 2000A. */
	gofa.Nut00a(tt1, tt2, &dp00, &de00)

	/* Precession-nutation quantities, IAU 2000. */
	gofa.Pn00(tt1, tt2, dp00, de00, &epsa, &rb, &rp, &rpb, &rn, &rnpb)

	/* Transform dX,dY corrections from GCRS to mean of date. */
	v1[0] = dx00
	v1[1] = dy00
	v1[2] = 0.0
	gofa.Rxp(rnpb, v1, &v2)
	ddp00 = v2[0] / math.Sin(epsa)
	dde00 = v2[1]

	/* Corrected nutation. */
	dpsi = dp00 + ddp00
	deps = de00 + dde00

	/* Build the rotation matrix. */
	gofa.Numat(epsa, dpsi, deps, &rn)

	/* Combine the matrices: PN = N x P. */
	gofa.Rxr(rn, rpb, &rnpb)

	/* Greenwich apparent sidereal time (IAU 2000). */
	gst = gofa.Anp(gofa.Gmst00(ut11, ut12, tt1, tt2) + gofa.Ee00(tt1, tt2, epsa, dpsi))

	/* Form celestial-terrestrial matrix (no polar motion yet). */
	gofa.Cr(rnpb, &rc2ti)
	gofa.Rz(gst, &rc2ti)

	/* Polar motion matrix (TIRS->ITRS, IERS 2003). */
	sp = gofa.Sp00(tt1, tt2)
	gofa.Pom00(xp, yp, sp, &rpom)

	/* Form celestial-terrestrial matrix (including polar motion). */
	gofa.Rxr(rpom, rc2ti, &rc2it)

	//
	fmt.Printf("IAU 2000A, equinox based\n\n")

	//
	fmt.Printf("NPB matrix, CIO based\n")
	printMatrix(rnpb)

	//
	fmt.Printf("GST = %s°\n    = %s\n\n", fmtDeg(gst), fmtHms(gst))

	fmt.Printf("celestial to terrestrial matrix (no polar motion)\n")
	printMatrix(rc2ti)

	fmt.Printf("celestial to terrestrial matrix\n")
	printMatrix(rc2it) //

	//
	/* ========================= */
	/* IAU 2006/2000A, CIO based */
	/* ========================= */

	/* CIP and CIO, IAU 2006/2000A. */
	gofa.Xys06a(tt1, tt2, &x, &y, &s)

	/* Add CIP corrections. */
	x += dx06
	y += dy06

	/* GCRS to CIRS matrix. */
	gofa.C2ixys(x, y, s, &rc2i)

	/* Earth rotation angle. */
	era = gofa.Era00(ut11, ut12)

	/* Form celestial-terrestrial matrix (no polar motion yet). */
	gofa.Cr(rc2i, &rc2ti)
	gofa.Rz(era, &rc2ti)

	/* Polar motion matrix (TIRS->ITRS, IERS 2003). */
	sp = gofa.Sp00(tt1, tt2)
	gofa.Pom00(xp, yp, sp, &rpom)

	/* Form celestial-terrestrial matrix (including polar motion). */
	gofa.Rxr(rpom, rc2ti, &rc2it)

	/* =========================================== */
	/* IAU 2006/2000A, CIO based, using X,Y series */
	/* =========================================== */

	/* CIP and CIO, IAU 2006/2000A. */
	gofa.Xy06(tt1, tt2, &x, &y)
	s = gofa.S06(tt1, tt2, x, y)

	/* Add CIP corrections. */
	x += dx06
	y += dy06

	/* GCRS to CIRS matrix. */
	gofa.C2ixys(x, y, s, &rc2i)

	/* Earth rotation angle. */
	era = gofa.Era00(ut11, ut12)

	/* Form celestial-terrestrial matrix (no polar motion yet). */
	gofa.Cr(rc2i, &rc2ti)
	gofa.Rz(era, &rc2ti)

	/* Polar motion matrix (TIRS->ITRS, IERS 2003). */
	sp = gofa.Sp00(tt1, tt2)
	gofa.Pom00(xp, yp, sp, &rpom)

	/* Form celestial-terrestrial matrix (including polar motion). */
	gofa.Rxr(rpom, rc2ti, &rc2it)

	//
	fmt.Printf("IAU 2006/2000A, CIO based\n\n")
	fmt.Printf("X = %.15f\nY = %.15f\ns = %.9f\"\n\n", x, y, s*gofa.DR2AS)

	//
	fmt.Printf("NPB matrix, CIO based\n")
	printMatrix(rc2i)

	//
	fmt.Printf("ERA = %s°\n    = %s\n\n", fmtDeg(era), fmtHms(era))

	fmt.Printf("celestial to terrestrial matrix (no polar motion)\n")
	printMatrix(rc2ti)

	fmt.Printf("celestial to terrestrial matrix\n")
	printMatrix(rc2it) //
}
