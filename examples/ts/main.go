// Copyright 2022 HE Boliang
// All rights reserved.

package main

import (
	"fmt"
	"math"

	"github.com/hebl/gofa"
)

func main() {
	ts()
}

/**
C version out:

UTC 2006/01/15 21:24:37.500000
UT1 2006/01/15 21:24:37.834100
TAI 2006/01/15 21:25:10.500000
TT  2006/01/15 21:25:42.684000
TCG 2006/01/15 21:25:43.322690
TDB 2006/01/15 21:25:42.684373
TCB 2006/01/15 21:25:56.893952

This version out:
UTC 2006/01/15 21:24:37.500000
UT1 2006/01/15 21:24:37.834100
TAI 2006/01/15 21:25:10.500000
TT  2006/01/15 21:25:42.684000
TCG 2006/01/15 21:25:43.322690
TDB 2006/01/15 21:25:42.684373
TCB 2006/01/15 21:25:56.893952

It's same!!!
*/

func ts() int {
	var latnd, latnm, lonwd, lonwm, j, iy, mo, id, ih, im int
	var ihmsf [4]int
	var slatn, slonw, hm, elon, phi, u, v, sec float64
	var utc1, utc2, dut, ut11, ut12, ut, tai1, tai2, tt1, tt2, tcg1, tcg2, dtr, tdb1, tdb2, tcb1, tcb2 float64
	var xyz [3]float64

	/* Site terrestrial coordinates (WGS84). */
	latnd = 19
	latnm = 28
	slatn = 52.5
	lonwd = 155
	lonwm = 55
	slonw = 59.6
	hm = 0.0

	/* Transform to geocentric. */
	j = gofa.Af2a('+', latnd, latnm, slatn, &phi)
	if j != 0 {
		return 1
	}
	j = gofa.Af2a('-', lonwd, lonwm, slonw, &elon)
	if j != 0 {
		return 1
	}
	j = gofa.Gd2gc(1, elon, phi, hm, &xyz)
	if j != 0 {
		return 1
	}

	u = math.Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1])
	v = xyz[2]

	/* UTC date and time. */
	iy = 2006
	mo = 1
	id = 15
	ih = 21
	im = 24
	sec = 37.5

	/* Transform into internal format. */
	j = gofa.Dtf2d("UTC", iy, mo, id, ih, im, sec, &utc1, &utc2)
	if j != 0 {
		return 1
	}

	/* UT1-UTC (s, from IERS). */
	dut = 0.3341

	/* UTC -> UT1. */
	j = gofa.Utcut1(utc1, utc2, dut, &ut11, &ut12)
	if j != 0 {
		return 1
	}

	/* Extract fraction for TDB-TT calculation, later. */
	ut = math.Mod(math.Mod(ut11, 1.0)+math.Mod(ut12, 1.0), 1.0) + 0.5

	/* UTC -> TAI -> TT -> TCG. */
	j = gofa.Utctai(utc1, utc2, &tai1, &tai2)
	if j != 0 {
		return 1
	}
	j = gofa.Taitt(tai1, tai2, &tt1, &tt2)
	if j != 0 {
		return 1
	}
	j = gofa.Tttcg(tt1, tt2, &tcg1, &tcg2)
	if j != 0 {
		return 1
	}

	/*TDB-TT (using TT as a substitute for TDB). */
	dtr = gofa.Dtdb(tt1, tt2, ut, elon, u/1e3, v/1e3)

	/* TT -> TDB -> TCB. */
	j = gofa.Tttdb(tt1, tt2, dtr, &tdb1, &tdb2)
	if j != 0 {
		return 1
	}
	j = gofa.Tdbtcb(tdb1, tdb2, &tcb1, &tcb2)
	if j != 0 {
		return 1
	}

	/* Report. */
	j = gofa.D2dtf("UTC", 6, utc1, utc2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("UTC%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("ut1", 6, ut11, ut12, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("UT1%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("tai", 6, tai1, tai2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("TAI%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("tt", 6, tt1, tt2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("TT %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("tcg", 6, tcg1, tcg2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("TCG%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("tdb", 6, tdb1, tdb2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("TDB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	j = gofa.D2dtf("tcb", 6, tcb1, tcb2, &iy, &mo, &id, &ihmsf)
	if j != 0 {
		return 1
	}
	fmt.Printf("TCB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
		iy, mo, id, ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3])

	return 0
}
