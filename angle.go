// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Angle

// Operations on Angles

/*
Anp   Normalize angle into the range 0 <= a < 2pi.

Given:
    a  float64     angle (radians)

Returned (function value):
       float64     angle in range 0-2pi
*/
func Anp(a float64) float64 {
	var w float64

	w = fmod(a, D2PI)
	if w < 0 {
		w += D2PI
	}

	return w

}

/*
Anpm  Normalize angle into the range -pi <= a < +pi.

Given:
    a   float64     angle (radians)

Returned (function value):
        float64     angle in range +/-pi
*/
func Anpm(a float64) float64 {
	var w float64

	w = fmod(a, D2PI)
	if fabs(w) >= DPI {
		w -= dsign(D2PI, a)
	}

	return w
}

/*
A2af Decompose radians into degrees, arcminutes, arcseconds, fraction.

Given:
    ndp     int      resolution (Note 1)
    angle   float64  angle in radians

Returned:
    sign    byte    '+' or '-'
    idmsf   [4]int  degrees, arcminutes, arcseconds, fraction

Notes:

 1) The argument ndp is interpreted as follows:

    ndp         resolution
     :      ...0000 00 00
    -7         1000 00 00
    -6          100 00 00
    -5           10 00 00
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
     0            0 00 01
     1            0 00 00.1
     2            0 00 00.01
     3            0 00 00.001
     :            0 00 00.000...

 2) The largest positive useful value for ndp is determined by the
    size of angle, the format of float64s on the target platform, and
    the risk of overflowing idmsf[3].  On a typical platform, for
    angle up to 2pi, the available floating-point precision might
    correspond to ndp=12.  However, the practical limit is typically
    ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    only 16 bits.

 3) The absolute value of angle may exceed 2pi.  In cases where it
    does not, it is up to the caller to test for and handle the
    case where angle is very nearly 2pi and rounds up to 360 degrees,
    by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.

Called:
    D2tf      decompose days to hms

*/
func A2af(ndp int, angle float64, sign *byte, idmsf *[4]int) {
	/* Hours to degrees * radians to turns */
	F := 15.0 / D2PI

	/* Scale then use days to h,m,s function. */
	D2tf(ndp, angle*F, sign, idmsf)

}

/*
A2tf Decompose radians into hours, minutes, seconds, fraction.

Given:
    ndp     int      resolution (Note 1)
    angle   float64  angle in radians

Returned:
    sign    byte    '+' or '-'
    ihmsf   [4]int  hours, minutes, seconds, fraction

Notes:

 1) The argument ndp is interpreted as follows:

    ndp         resolution
     :      ...0000 00 00
    -7         1000 00 00
    -6          100 00 00
    -5           10 00 00
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
     0            0 00 01
     1            0 00 00.1
     2            0 00 00.01
     3            0 00 00.001
     :            0 00 00.000...

 2) The largest positive useful value for ndp is determined by the
    size of angle, the format of float64s on the target platform, and
    the risk of overflowing ihmsf[3].  On a typical platform, for
    angle up to 2pi, the available floating-point precision might
    correspond to ndp=12.  However, the practical limit is typically
    ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    only 16 bits.

 3) The absolute value of angle may exceed 2pi.  In cases where it
    does not, it is up to the caller to test for and handle the
    case where angle is very nearly 2pi and rounds up to 24 hours,
    by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.

Called:
    D2tf      decompose days to hms
*/
func A2tf(ndp int, angle float64, sign *byte, ihmsf *[4]int) {
	D2tf(ndp, angle/D2PI, sign, ihmsf)
}

/*
D2tf Decompose days to hours, minutes, seconds, fraction.

Given:
    ndp     int      resolution (Note 1)
    days    float64  interval in days

Returned:
    sign    byte    '+' or '-'
    ihmsf   [4]int  hours, minutes, seconds, fraction

Notes:

 1) The argument ndp is interpreted as follows:

    ndp         resolution
     :      ...0000 00 00
    -7         1000 00 00
    -6          100 00 00
    -5           10 00 00
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
     0            0 00 01
     1            0 00 00.1
     2            0 00 00.01
     3            0 00 00.001
     :            0 00 00.000...

 2) The largest positive useful value for ndp is determined by the
    size of days, the format of float64 on the target platform, and
    the risk of overflowing ihmsf[3].  On a typical platform, for
    days up to 1.0, the available floating-point precision might
    correspond to ndp=12.  However, the practical limit is typically
    ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
    only 16 bits.

 3) The absolute value of days may exceed 1.0.  In cases where it
    does not, it is up to the caller to test for and handle the
    case where days is very nearly 1.0 and rounds up to 24 hours,
    by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
*/
func D2tf(ndp int, days float64, sign *byte, ihmsf *[4]int) {
	var nrs, n int
	var rs, rm, rh, a, w, ah, am, as, af float64

	/* Handle sign. */
	if days >= 0.0 {
		*sign = '+'

	} else {
		*sign = '-'
	}

	/* Interval in seconds. */
	a = DAYSEC * fabs(days)

	/* Pre-round if resolution coarser than 1s (then pretend ndp=1). */
	if ndp < 0 {
		nrs = 1
		for n = 1; n <= -ndp; n++ {
			if n == 2 || n == 4 {
				nrs *= 6
			} else {
				nrs *= 4
			}

		}
		rs = float64(nrs)
		w = a / rs
		a = rs * dnint(w)
	}

	/* Express the unit of each field in resolution units. */
	nrs = 1
	for n = 1; n <= ndp; n++ {
		nrs *= 10
	}
	rs = float64(nrs)
	rm = rs * 60.0
	rh = rm * 60.0

	/* Round the interval and express in resolution units. */
	a = dnint(rs * a)

	/* Break into fields. */
	ah = a / rh
	ah = dint(ah)
	a -= ah * rh
	am = a / rm
	am = dint(am)
	a -= am * rm
	as = a / rs
	as = dint(as)
	af = a - as*rs

	/* Return results. */
	ihmsf[0] = int(ah)
	ihmsf[1] = int(am)
	ihmsf[2] = int(as)
	ihmsf[3] = int(af)
}

/*
Af2a    Convert degrees, arcminutes, arcseconds to radians.

Given:
    s         byte     sign:  '-' = negative, otherwise positive
    ideg      int      degrees
    iamin     int      arcminutes
    asec      float64  arcseconds

Returned:
    rad       float64  angle in radians

Returned (function value):
    int    status:  0 = OK
                    1 = ideg outside range 0-359
                    2 = iamin outside range 0-59
                    3 = asec outside range 0-59.999...

Notes:

 1)  The result is computed even if any of the range checks fail.

 2)  Negative ideg, iamin and/or asec produce a warning status, but
     the absolute value is used in the conversion.

 3)  If there are multiple errors, the status value reflects only the
     first, the smallest taking precedence.
*/
func Af2a(s byte, ideg, iamin int, asec float64, rad *float64) int {
	/* Compute the interval. */
	if s == '-' {
		*rad = -1.0
	} else {
		*rad = 1.0
	}
	*rad *= (60.0*(60.0*(fabs(float64(ideg)))+(fabs(float64(iamin)))) + fabs(asec)) * DAS2R

	/* Validate arguments and return status. */
	if ideg < 0 || ideg > 359 {
		return 1
	}
	if iamin < 0 || iamin > 59 {
		return 2
	}
	if asec < 0.0 || asec >= 60.0 {
		return 3
	}
	return 0
}

/*
Tf2a  Convert hours, minutes, seconds to radians.

Given:
    s         byte     sign:  '-' = negative, otherwise positive
    ihour     int      hours
    imin      int      minutes
    sec       float64  seconds

Returned:
    rad       float64  angle in radians

Returned (function value):
    int     status:  0 = OK
                     1 = ihour outside range 0-23
                     2 = imin outside range 0-59
                     3 = sec outside range 0-59.999...

Notes:

 1) The result is computed even if any of the range checks fail.

 2) Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

 3) If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.
*/
func Tf2a(s byte, ihour, imin int, sec float64, rad *float64) int {
	/* Compute the interval. */
	if s == '-' {
		*rad = -1.0
	} else {
		*rad = 1.0
	}

	*rad *= (60.0*(60.0*(fabs(float64(ihour)))+(fabs(float64(imin)))) + fabs(sec)) * DS2R

	/* Validate arguments and return status. */
	if ihour < 0 || ihour > 23 {
		return 1
	}
	if imin < 0 || imin > 59 {
		return 2
	}
	if sec < 0.0 || sec >= 60.0 {
		return 3
	}
	return 0
}

/*
Tf2d  Convert hours, minutes, seconds to days.

Given:
    s         byte     sign:  '-' = negative, otherwise positive
    ihour     int      hours
    imin      int      minutes
    sec       float64  seconds

Returned:
    days      float64  interval in days

Returned (function value):
    int    status:  0 = OK
                    1 = ihour outside range 0-23
                    2 = imin outside range 0-59
                    3 = sec outside range 0-59.999...

Notes:

 1) The result is computed even if any of the range checks fail.

 2) Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

 3) If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.
*/
func Tf2d(s byte, ihour, imin int, sec float64, days *float64) int {
	/* Compute the interval. */
	// *days  = ( s == '-' ? -1.0 : 1.0 ) *
	if s == '-' {
		*days = -1.0
	} else {
		*days = 1.0
	}
	*days *= (60.0*(60.0*(fabs(float64(ihour)))+(fabs(float64(imin)))) + fabs(sec)) / DAYSEC

	/* Validate arguments and return status. */

	if ihour < 0 || ihour > 23 {
		return 1
	}
	if imin < 0 || imin > 59 {
		return 2
	}
	if sec < 0.0 || sec >= 60.0 {
		return 3
	}
	return 0

}

// Separation and position-angle

/*
Sepp  Angular separation between two p-vectors.

Given:
    a      [3]float64    first p-vector (not necessarily unit length)
    b      [3]float64    second p-vector (not necessarily unit length)

Returned (function value):
           float64       angular separation (radians, always positive)

Notes:

 1) If either vector is null, a zero result is returned.

 2) The angular separation is most simply formulated in terms of
    scalar product.  However, this gives poor accuracy for angles
    near zero and pi.  The present algorithm uses both cross product
    and dot product, to deliver full accuracy whatever the size of
    the angle.

Called:
    Pxp       vector product of two p-vectors
    Pm        modulus of p-vector
    Pdp       scalar product of two p-vectors
*/
func Sepp(a, b [3]float64) float64 {
	var axb [3]float64
	var ss, cs, s float64

	/* Sine of angle between the vectors, multiplied by the two moduli. */
	Pxp(a, b, &axb)
	ss = Pm(axb)

	/* Cosine of the angle, multiplied by the two moduli. */
	cs = Pdp(a, b)

	/* The angle. */
	if (ss != 0.0) || (cs != 0.0) {
		s = atan2(ss, cs)
	} else {
		s = 0.0
	}

	return s
}

/*
Seps Angular separation between two sets of spherical coordinates.

Given:
    al     float64       first longitude (radians)
    ap     float64       first latitude (radians)
    bl     float64       second longitude (radians)
    bp     float64       second latitude (radians)

Returned (function value):
           float64       angular separation (radians)

Called:
    S2c       spherical coordinates to unit vector
    Sepp      angular separation between two p-vectors
*/
func Seps(al, ap, bl, bp float64) float64 {
	var ac, bc [3]float64
	var s float64

	/* Spherical to Cartesian. */
	S2c(al, ap, &ac)
	S2c(bl, bp, &bc)

	/* Angle between the vectors. */
	s = Sepp(ac, bc)

	return s
}

/*
Pap Position-angle from two p-vectors.

Given:
    a      [3]float64  direction of reference point
    b      [3]float64  direction of point whose PA is required

Returned (function value):
           float64     position angle of b with respect to a (radians)

Notes:

 1) The result is the position angle, in radians, of direction b with
    respect to direction a.  It is in the range -pi to +pi.  The
    sense is such that if b is a small distance "north" of a the
    position angle is approximately zero, and if b is a small
    distance "east" of a the position angle is approximately +pi/2.

 2) The vectors a and b need not be of unit length.

 3) Zero is returned if the two directions are the same or if either
    vector is null.

 4) If vector a is at a pole, the result is ill-defined.

Called:
    Pn        decompose p-vector into modulus and direction
    Pm        modulus of p-vector
    Pxp       vector product of two p-vectors
    Pmp       p-vector minus p-vector
    Pdp       scalar product of two p-vectors
*/
func Pap(a, b [3]float64) float64 {
	var am, bm, st, ct, xa, ya, za, pa float64
	var au, eta, xi, a2b [3]float64

	/* Modulus and direction of the a vector. */
	Pn(a, &am, &au)

	/* Modulus of the b vector. */
	bm = Pm(b)

	/* Deal with the case of a null vector. */
	if (am == 0.0) || (bm == 0.0) {
		st = 0.0
		ct = 1.0
	} else {

		/* The "north" axis tangential from a (arbitrary length). */
		xa = a[0]
		ya = a[1]
		za = a[2]
		eta[0] = -xa * za
		eta[1] = -ya * za
		eta[2] = xa*xa + ya*ya

		/* The "east" axis tangential from a (same length). */
		Pxp(eta, au, &xi)

		/* The vector from a to b. */
		Pmp(b, a, &a2b)

		/* Resolve into components along the north and east axes. */
		st = Pdp(a2b, xi)
		ct = Pdp(a2b, eta)

		/* Deal with degenerate cases. */
		if (st == 0.0) && (ct == 0.0) {
			ct = 1.0
		}
	}

	/* Position angle. */
	pa = atan2(st, ct)

	return pa
}

/*
Pas Position-angle from spherical coordinates.

Given:
    al     float64     longitude of point A (e.g. RA) in radians
    ap     float64     latitude of point A (e.g. Dec) in radians
    bl     float64     longitude of point B
    bp     float64     latitude of point B

Returned (function value):
          float64     position angle of B with respect to A

Notes:

 1) The result is the bearing (position angle), in radians, of point
    B with respect to point A.  It is in the range -pi to +pi.  The
    sense is such that if B is a small distance "east" of point A,
    the bearing is approximately +pi/2.

 2) Zero is returned if the two points are coincident.
*/
func Pas(al, ap, bl, bp float64) float64 {
	var dl, x, y, pa float64

	dl = bl - al
	y = sin(dl) * cos(bp)
	x = sin(bp)*cos(ap) - cos(bp)*sin(ap)*cos(dl)
	// pa = ((x != 0.0) || (y != 0.0)) ? atan2(y, x) : 0.0;
	if (x != 0.0) || (y != 0.0) {
		pa = atan2(y, x)
	} else {
		pa = 0.0
	}

	return pa
}
