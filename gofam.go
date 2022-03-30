// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Constant value and math functions

import "math"

// Const
const (
	/* Pi */
	DPI = 3.141592653589793238462643

	/* 2Pi */
	D2PI = 6.283185307179586476925287

	/* Pi/2 */
	DPI2 = 1.570796326794896619231322

	/* Radians to degrees */
	DR2D = 57.29577951308232087679815

	/* Degrees to radians */
	DD2R = 1.745329251994329576923691e-2

	/* Radians to arcseconds */
	DR2AS = 206264.8062470963551564734

	/* Arcseconds to radians */
	DAS2R = 4.848136811095359935899141e-6

	/* Seconds of time to radians */
	DS2R = 7.272205216643039903848712e-5

	/* Arcseconds in a full circle */
	TURNAS = 1296000.0

	/* Milliarcseconds to radians */
	DMAS2R = (DAS2R / 1e3)

	/* Length of tropical year B1900 (days) */
	DTY = 365.242198781

	/* Seconds per day. */
	DAYSEC = 86400.0

	/* Days per Julian year */
	DJY = 365.25

	/* Days per Julian century */
	DJC = 36525.0

	/* Days per Julian millennium */
	DJM = 365250.0

	/* Reference epoch (J2000.0), Julian Date */
	DJ00 = 2451545.0

	/* Julian Date of Modified Julian Date zero */
	DJM0 = 2400000.5

	/* Reference epoch (J2000.0), Modified Julian Date */
	DJM00 = 51544.5

	/* 1977 Jan 1.0 as MJD */
	DJM77 = 43144.0

	/* TT minus TAI (s) */
	TTMTAI = 32.184

	/* Astronomical unit (m, IAU 2012) */
	DAU = 149597870.7e3

	/* Speed of light (m/s) */
	CMPS = 299792458.0

	/* Light time for 1 au (s) */
	AULT = (DAU / CMPS)

	/* Speed of light (au per day) */
	DC = (DAYSEC / AULT)

	/* L_G = 1 - d(TT)/d(TCG) */
	ELG = 6.969290134e-10

	/* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
	ELB  = 1.550519768e-8
	TDB0 = -6.55e-5

	/* Schwarzschild radius of the Sun (au) */
	/* = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 */
	SRS = 1.97412574336e-8

	/* Reference ellipsoids */
	WGS84 = 1
	GRS80 = 2
	WGS72 = 3

	//
	DBL_EPSILON = 2.220446049250313080847263336181640625e-16
)

/* dint(A) - truncate to nearest whole number towards zero (double) */
func dint(A float64) float64 {
	if A < 0.0 {
		return math.Ceil(A)
	} else {
		return math.Floor(A)
	}
}

/* dnint(A) - round to nearest whole number (double) */
func dnint(A float64) float64 {
	if math.Abs(A) < 0.5 {
		return 0.0
	} else {
		if A < 0.0 {
			return math.Ceil(A - 0.5)
		} else {
			return math.Floor(A + 0.5)
		}
	}
}

/* dsign(A,B) - magnitude of A with sign of B (double) */
func dsign(A, B float64) float64 {
	if B < 0.0 {
		return -math.Abs(A)
	} else {
		return math.Abs(A)
	}
}

/* max(A,B) - larger (most +ve) of two numbers (generic) */
func gmax(A, B float64) float64 {
	if A > B {
		return A
	}
	return B
}

/* min(A,B) - smaller (least +ve) of two numbers (generic) */
func gmin(A, B float64) float64 {
	if A < B {
		return A
	}
	return B
}

func abs(A int) int {
	if A < 0 {
		return -A
	} else {
		return A
	}
}

//math function
func fabs(x float64) float64     { return math.Abs(x) }
func sin(x float64) float64      { return math.Sin(x) }
func cos(x float64) float64      { return math.Cos(x) }
func tan(x float64) float64      { return math.Tan(x) }
func asin(x float64) float64     { return math.Asin(x) }
func atan(x float64) float64     { return math.Atan(x) }
func atan2(x, y float64) float64 { return math.Atan2(x, y) }
func sqrt(x float64) float64     { return math.Sqrt(x) }
func fmod(x, y float64) float64  { return math.Mod(x, y) }
func pow10(n int) float64        { return math.Pow10(n) }
func pow(x, y float64) float64   { return math.Pow(x, y) }

const (
	sqrt2 = math.Sqrt2
)
