// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

const (
	VERSION = 1.18
)

/*
ASTROM Star-independent astrometry parameters

(Vectors Eb, Eh, Em and V are all with respect to BCRS axes.)
*/
type ASTROM struct {
	Pmt    float64       /* PM time interval (SSB, Julian years) */
	Eb     [3]float64    /* SSB to observer (vector, au) */
	Eh     [3]float64    /* Sun to observer (unit vector) */
	Em     float64       /* distance from Sun to observer (au) */
	V      [3]float64    /* barycentric observer velocity (vector, c) */
	Bm1    float64       /* sqrt(1-|v|^2): reciprocal of Lorenz factor */
	Bpn    [3][3]float64 /* bias-precession-nutation matrix */
	Along  float64       /* longitude + s' + dERA(DUT) (radians) */
	Phi    float64       /* geodetic latitude (radians) */
	Xpl    float64       /* polar motion xp wrt local meridian (radians) */
	Ypl    float64       /* polar motion yp wrt local meridian (radians) */
	Sphi   float64       /* sine of geodetic latitude */
	Cphi   float64       /* cosine of geodetic latitude */
	Diurab float64       /* magnitude of diurnal aberration vector */
	Eral   float64       /* "local" Earth rotation angle (radians) */
	Refa   float64       /* refraction constant A (radians) */
	Refb   float64       /* refraction constant B (radians) */
}

/*
LDBODY Body parameters for light deflection
*/
type LDBODY struct {
	Bm float64       /* mass of the body (solar masses) */
	Dl float64       /* deflection limiter (radians^2/2) */
	Pv [2][3]float64 /* barycentric PV of the body (au, au/day) */
}
