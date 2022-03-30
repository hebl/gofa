// Copyright 2022 HE Boliang
// All rights reserved.

package main

import (
	"fmt"

	"github.com/hebl/gofa"
)

// Gnomonic projection example sofa_ast_c.pdf P26.

func main() {
	const fl = 3430.0 /* Focal length (mm) */
	var sign byte
	var j, n int
	var i [4]int
	var raz, decz, ra, dec, xi, eta, x, y, az1, bz1, az2, bz2 float64

	/* Observed [RA,Dec] of tangent point (radians). */
	raz = 9.0 * 15.0 * gofa.DD2R
	decz = 89.5 * gofa.DD2R

	//
	/* Observed [RA,Dec] of a star (radians). */
	ra = 5.0 * 15.0 * gofa.DD2R
	dec = 89.75 * gofa.DD2R
	/* Plate coordinates of star image. */
	j = gofa.Tpxes(ra, dec, raz, decz, &xi, &eta)
	if j != 0 {
		return
	}

	fmt.Printf("image [x,y]/mm = %+8.3f, %+8.3f\n", xi*fl, eta*fl)
	// C out: 		image [x,y]/mm =  -12.961,  +22.450
	// This out:	image [x,y]/mm =  -12.961,  +22.450

	/* Plate coordinates of a different star image. */
	x = -20.0
	y = 60.0

	/* Celestial coordinates of that star. */
	xi = x / fl
	eta = y / fl
	gofa.Tpsts(xi, eta, raz, decz, &ra, &dec)
	gofa.A2tf(2, ra, &sign, &i)
	fmt.Printf("star [RA,Dec] = %2.2d %2.2d %2.2d.%2.2d",
		i[0], i[1], i[2], i[3])
	gofa.A2af(1, dec, &sign, &i)
	fmt.Printf(" %c%2.2d %2.2d %2.2d.%1d\n", sign, i[0], i[1], i[2], i[3])
	// C out: 		star [RA,Dec] = 23 14 31.74 +89 23 48.8
	// This out:	star [RA,Dec] = 23 14 31.74 +89 23 48.8

	//
	/* Given that star's [RA,Dec] and [x,y], solve for tangent point. */
	n = gofa.Tpors(xi, eta, ra, dec, &az1, &bz1, &az2, &bz2)
	if n != 0 {
		gofa.A2tf(2, az1, &sign, &i)
		fmt.Printf("TP1 = %2.2d %2.2d %2.2d.%2.2d",
			i[0], i[1], i[2], i[3])
		gofa.A2af(1, bz1, &sign, &i)
		fmt.Printf(" %c%2.2d %2.2d %2.2d.%1d\n", sign, i[0], i[1], i[2], i[3])
		if n == 2 {
			gofa.A2tf(2, az2, &sign, &i)
			fmt.Printf("TP2 = %2.2d %2.2d %2.2d.%2.2d",
				i[0], i[1], i[2], i[3])
			gofa.A2af(1, bz2, &sign, &i)
			fmt.Printf(" %c%2.2d %2.2d %2.2d.%1d\n", sign, i[0], i[1], i[2], i[3])
		}
	}
	// C out:
	//   TP1 = 01 29 03.48 +88 29 44.5
	//   TP2 = 09 00 00.00 +89 30 00.0
	// This out:
	//   TP1 = 01 29 03.48 +88 29 44.5
	//   TP2 = 09 00 00.00 +89 30 00.0

}
