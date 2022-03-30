// Copyright 2022 HE Boliang
// All rights reserved.

package main

import (
	"fmt"
	"math"

	"github.com/hebl/gofa"
)

func main() {
	p03()
	fmt.Println()

	p05()
	fmt.Println()

	p09()
	fmt.Println()

	p11()
}

func p03() {
	// C out:
	// 17011 km,  61 deg

	// This out:
	// 17011 km,  61 deg

	const R2D = 57.2957795 /* radians to degrees */
	const RKM = 6375.0     /* Earth radius in km */

	/* Longitudes and latitudes (radians) for London and Sydney. */
	al := -0.2 / R2D
	bl := 51.5 / R2D
	as := 151.2 / R2D
	bs := -33.9 / R2D

	/*Great-circle distance and initial heading. */
	fmt.Printf("%5.0f km,%4.0f deg\n",
		gofa.Seps(al, bl, as, bs)*RKM,
		gofa.Pas(al, bl, as, bs)*R2D)

}

func p05() {
	// C out:
	// -00 59 59.999
	// -01 00 00.00

	// This out:
	// 	-00 59 59.999
	// -01 00 00.00

	var ha float64
	var sign byte
	var ihmsf [4]int

	ha = -0.261799315

	gofa.A2tf(3, ha, &sign, &ihmsf)
	fmt.Printf("%c%2.2d %2.2d %2.2d.%3.3d\n", sign, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	gofa.A2tf(2, ha, &sign, &ihmsf)
	fmt.Printf("%c%2.2d %2.2d %2.2d.%2.2d\n", sign, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])
}

func p09() {
	// C out:
	// The two positions are  0.000000 arcsec apart.
	// The two positions are  0.000612 arcsec apart.

	// This out:
	// The two positions are  0.000000 arcsec apart.
	// The two positions are  0.000612 arcsec apart.

	var ra, da, rb, db, theta, s, c float64
	var a, b, axb [3]float64

	/* RA,Dec of source ICRF J044238.6-001743 in the ICRF2 catalog. */
	gofa.Tf2a('+', 04, 42, 38.66073910, &ra)
	gofa.Af2a('-', 00, 17, 43.4203921, &da)

	/* RA,Dec of the same source in the ICRF3 (S/X) catalog. */
	gofa.Tf2a('+', 04, 42, 38.66072366, &rb)
	gofa.Af2a('-', 00, 17, 43.4209582, &db)

	/* Method 1: spherical trigonometry (cosine rule). */
	theta = math.Acos(math.Sin(da)*math.Sin(db) + math.Cos(da)*math.Cos(db)*math.Cos(rb-ra))
	fmt.Printf("The two positions are %9.6f arcsec apart.\n", theta*206264.80624709635515647335733)

	/* Method 2: vectors (sine and cosine from cross and dot product). */
	gofa.S2c(ra, da, &a)
	gofa.S2c(rb, db, &b)
	gofa.Pxp(a, b, &axb)
	s = gofa.Pm(axb)
	c = gofa.Pdp(a, b)
	theta = math.Atan2(s, c)
	fmt.Printf("The two positions are %9.6f arcsec apart.\n", theta*206264.80624709635515647335733)
}

func p11() {
	// sofa_vm_c.pdf p11
	// C out:
	// Frame bias is 23.1 mas around 19 29 14.8 -39 06 19 GCRS.

	// This out:
	// Frame bias is 23.1 mas around 19 29 14.8 -39 06 19 GCRS.

	var sign byte
	var i4 [4]int
	var rb [3][3]float64
	var vb [3]float64
	var angle, ra, dec float64

	/* Generate frame bias matrix. */
	gofa.Ir(&rb)
	gofa.Rz(-0.0146*gofa.DAS2R, &rb)
	gofa.Ry(-0.041775*gofa.DAS2R*math.Sin(84381.448*gofa.DAS2R), &rb)
	gofa.Rx(0.0068192*gofa.DAS2R, &rb)

	/* Convert into r-vector form. */
	gofa.Rm2v(rb, &vb)

	/* Report. */
	gofa.P2s(vb, &ra, &dec, &angle)
	fmt.Printf("Frame bias is %4.1f mas around", angle*1e3/gofa.DAS2R)

	gofa.A2tf(1, gofa.Anp(ra), &sign, &i4)
	fmt.Printf("%3.2d %2.2d %2.2d.%d", i4[0], i4[1], i4[2], i4[3])

	gofa.A2af(0, dec, &sign, &i4)
	fmt.Printf(" %c%2.2d %2.2d %2.2d GCRS. \n", sign, i4[0], i4[1], i4[2])
}
