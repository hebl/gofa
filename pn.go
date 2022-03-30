// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// SOFA Precession / Nutation / Polar Motion

/*
Bi00    Frame bias components, IAU 2000 Frame bias components of IAU 2000
precession-nutation models;  part of the Mathews-Herring-Buffett (MHB2000)
nutation series, with additions.

Returned:
    dpsibi,depsbi  float64  longitude and obliquity corrections
    dra            float64  the ICRS RA of the J2000.0 mean equinox

Notes:

 1) The frame bias corrections in longitude and obliquity (radians)
    are required in order to correct for the offset between the GCRS
    pole and the mean J2000.0 pole.  They define, with respect to the
    GCRS frame, a J2000.0 mean pole that is consistent with the rest
    of the IAU 2000A precession-nutation model.

 2) In addition to the displacement of the pole, the complete
    description of the frame bias requires also an offset in right
    ascension.  This is not part of the IAU 2000A model, and is from
    Chapront et al. (2002).  It is returned in radians.

 3) This is a supplemented implementation of one aspect of the IAU
    2000A nutation model, formally adopted by the IAU General
    Assembly in 2000, namely MHB2000 (Mathews et al. 2002).

References:

    Chapront, J., Chapront-Touze, M. & Francou, G., Astron.
    Astrophys., 387, 700, 2002.

    Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
    and precession:  New nutation series for nonrigid Earth and
    insights into the Earth's interior", J.Geophys.Res., 107, B4,
    2002.  The MHB2000 code itself was obtained on 2002 September 9
    from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
*/
func Bi00(dpsibi, depsbi, dra *float64) {
	/* The frame bias corrections in longitude and obliquity */
	const DPBIAS = -0.041775 * DAS2R
	const DEBIAS = -0.0068192 * DAS2R

	/* The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) */
	const DRA0 = -0.0146 * DAS2R

	/* Return the results (which are fixed). */
	*dpsibi = DPBIAS
	*depsbi = DEBIAS
	*dra = DRA0
}

/*
Bp00 Frame bias and precession matrices, IAU 2000 Frame bias and precession, IAU 2000.

Given:
    date1,date2  float64         TT as a 2-part Julian Date (Note 1)

Returned:
    rb           [3][3]float64   frame bias matrix (Note 2)
    rp           [3][3]float64   precession matrix (Note 3)
    rbp          [3][3]float64   bias-precession matrix (Note 4)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1         date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
    applying frame bias.

 3) The matrix rp transforms vectors from J2000.0 mean equator and
    equinox to mean equator and equinox of date by applying
    precession.

 4) The matrix rbp transforms vectors from GCRS to mean equator and
    equinox of date by applying frame bias then precession.  It is
    the product rp x rb.

 5) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the order given.

Called:
    Bi00      frame bias components, IAU 2000
    Pr00      IAU 2000 precession adjustments
    Ir        initialize r-matrix to identity
    Rx        rotate around X-axis
    Ry        rotate around Y-axis
    Rz        rotate around Z-axis
    Cr        copy r-matrix
    Rxr       product of two r-matrices

Reference:
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
      intermediate origin" (CIO) by IAU 2006 Resolution 2.
*/
func Bp00(date1, date2 float64, rb, rp, rbp *[3][3]float64) {
	/* J2000.0 obliquity (Lieske et al. 1977) */
	const EPS0 = 84381.448 * DAS2R

	var t, dpsibi, depsbi, dra0, psia77, oma77, chia,
		dpsipr, depspr, psia, oma float64
	var rbw [3][3]float64

	/* Interval between fundamental epoch J2000.0 and current date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Frame bias. */
	Bi00(&dpsibi, &depsbi, &dra0)

	/* Precession angles (Lieske et al. 1977) */
	psia77 = (5038.7784 + (-1.07259+(-0.001147)*t)*t) * t * DAS2R
	oma77 = EPS0 + ((0.05127+(-0.007726)*t)*t)*t*DAS2R
	chia = (10.5526 + (-2.38064+(-0.001125)*t)*t) * t * DAS2R

	/* Apply IAU 2000 precession corrections. */
	Pr00(date1, date2, &dpsipr, &depspr)
	psia = psia77 + dpsipr
	oma = oma77 + depspr

	/* Frame bias matrix: GCRS to J2000.0. */
	Ir(&rbw)
	Rz(dra0, &rbw)
	Ry(dpsibi*sin(EPS0), &rbw)
	Rx(-depsbi, &rbw)
	Cr(rbw, rb)

	/* Precession matrix: J2000.0 to mean of date. */
	Ir(rp)
	Rx(EPS0, rp)
	Rz(-psia, rp)
	Rx(-oma, rp)
	Rz(chia, rp)

	/* Bias-precession matrix: GCRS to mean of date. */
	Rxr(*rp, rbw, rbp)
}

/*
Bp06  Frame bias and precession matrices, IAU 2006 Frame bias and precession, IAU 2006.

Given:
    date1,date2  float64         TT as a 2-part Julian Date (Note 1)

Returned:
    rb           [3][3]float64   frame bias matrix (Note 2)
    rp           [3][3]float64   precession matrix (Note 3)
    rbp          [3][3]float64   bias-precession matrix (Note 4)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1         date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
    applying frame bias.

 3) The matrix rp transforms vectors from mean J2000.0 to mean of
    date by applying precession.

 4) The matrix rbp transforms vectors from GCRS to mean of date by
    applying frame bias then precession.  It is the product rp x rb.

 5) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the order given.

Called:
    Pfw06     bias-precession F-W angles, IAU 2006
    Fw2m      F-W angles to r-matrix
    Pmat06    PB matrix, IAU 2006
    Tr        transpose r-matrix
    Rxr       product of two r-matrices
    Cr        copy r-matrix

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Bp06(date1, date2 float64, rb, rp, rbp *[3][3]float64) {
	var gamb, phib, psib, epsa float64
	var rbpw, rbt [3][3]float64

	/* B matrix. */
	Pfw06(DJM0, DJM00, &gamb, &phib, &psib, &epsa)
	Fw2m(gamb, phib, psib, epsa, rb)

	/* PxB matrix (temporary). */
	Pmat06(date1, date2, &rbpw)

	/* P matrix. */
	Tr(*rb, &rbt)
	Rxr(rbpw, rbt, rp)

	/* PxB matrix. */
	Cr(rbpw, rbp)
}

/*
Bpn2xy CIP XY given Bias-precession-nutation matrix Extract from the
bias-precession-nutation matrix the X,Y coordinates of the Celestial Intermediate Pole.

Given:
    rbpn      [3][3]float64  celestial-to-true matrix (Note 1)

Returned:
    x,y       float64        Celestial Intermediate Pole (Note 2)

Notes:

 1) The matrix rbpn transforms vectors from GCRS to true equator (and
    CIO or equinox) of date, and therefore the Celestial Intermediate
    Pole unit vector is the bottom row of the matrix.

 2) The arguments x,y are components of the Celestial Intermediate
    Pole unit vector in the Geocentric Celestial Reference System.

Reference:

    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154
    (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
      intermediate origin" (CIO) by IAU 2006 Resolution 2.
*/
func Bpn2xy(rbpn [3][3]float64, x, y *float64) {
	/* Extract the X,Y coordinates. */
	*x = rbpn[2][0]
	*y = rbpn[2][1]
}

/*
C2i00a Celestial-to-intermediate matrix, IAU 2000A Form the
celestial-to-intermediate matrix for a given date using the
IAU 2000A precession-nutation model.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rc2i        [3][3]float64 celestial-to-intermediate matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

      date1             date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

           =  rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

 3) A faster, but slightly less accurate, result (about 1 mas) can be
    obtained by using instead the iauC2i00b function.

Called:
    Pnm00a    classical NPB matrix, IAU 2000A
    C2ibpn    celestial-to-intermediate matrix, given NPB matrix

References:

    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154
    (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
      intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2i00a(date1, date2 float64, rc2i *[3][3]float64) {
	var rbpn [3][3]float64

	/* Obtain the celestial-to-true matrix (IAU 2000A). */
	Pnm00a(date1, date2, &rbpn)

	/* Form the celestial-to-intermediate matrix. */
	C2ibpn(date1, date2, rbpn, rc2i)
}

/*
C2i00b Celestial-to-intermediate matrix, IAU 2000B
Form the celestial-to-intermediate matrix for a given date using the
IAU 2000B precession-nutation model.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rc2i        [3][3]float64 celestial-to-intermediate matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

      date1             date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

           =  rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

 3) The present function is faster, but slightly less accurate (about
    1 mas), than the iauC2i00a function.

Called:
    Pnm00b    classical NPB matrix, IAU 2000B
    C2ibpn    celestial-to-intermediate matrix, given NPB matrix

References:

    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154
    (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
      intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2i00b(date1, date2 float64, rc2i *[3][3]float64) {
	var rbpn [3][3]float64

	/* Obtain the celestial-to-true matrix (IAU 2000B). */
	Pnm00b(date1, date2, &rbpn)

	/* Form the celestial-to-intermediate matrix. */
	C2ibpn(date1, date2, rbpn, rc2i)
}

/*
C2i06a Celestial-to-intermediate matrix, IAU 2006/2000A
Form the celestial-to-intermediate matrix for a given date using the
IAU 2006 precession and IAU 2000A nutation models.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rc2i        [3][3]float64 celestial-to-intermediate matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

      date1             date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

           =  RC2T * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

Called:
    Pnm06a    classical NPB matrix, IAU 2006/2000A
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S06       the CIO locator s, given X,Y, IAU 2006
    C2ixys    celestial-to-intermediate matrix, given X,Y and s

References:

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG
*/
func C2i06a(date1, date2 float64, rc2i *[3][3]float64) {
	var rbpn [3][3]float64
	var x, y, s float64

	/* Obtain the celestial-to-true matrix (IAU 2006/2000A). */
	Pnm06a(date1, date2, &rbpn)

	/* Extract the X,Y coordinates. */
	Bpn2xy(rbpn, &x, &y)

	/* Obtain the CIO locator. */
	s = S06(date1, date2, x, y)

	/* Form the celestial-to-intermediate matrix. */
	C2ixys(x, y, s, rc2i)
}

/*
C2ibpn  Celestial-to-intermediate matrix given B-P-N
Form the celestial-to-intermediate matrix for a given date given
the bias-precession-nutation matrix.  IAU 2000.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)
    rbpn        [3][3]float64 celestial-to-true matrix (Note 2)

Returned:
    rc2i        [3][3]float64 celestial-to-intermediate matrix (Note 3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

      date1             date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix rbpn transforms vectors from GCRS to true equator (and
    CIO or equinox) of date.  Only the CIP (bottom row) is used.

 3) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

          = RC2T * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

 4) Although its name does not include "00", This function is in fact
    specific to the IAU 2000 models.

Called:
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    C2ixy     celestial-to-intermediate matrix, given X,Y

References:
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
      intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2ibpn(date1, date2 float64, rbpn [3][3]float64, rc2i *[3][3]float64) {
	var x, y float64

	/* Extract the X,Y coordinates. */
	Bpn2xy(rbpn, &x, &y)

	/* Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific). */
	C2ixy(date1, date2, x, y, rc2i)
}

/*
C2ixy Celestial-to-intermediate matrix given CIP
Form the celestial to intermediate-frame-of-date matrix for a given
date when the CIP X,Y coordinates are known. IAU 2000.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)
    x,y         float64       Celestial Intermediate Pole (Note 2)

Returned:
    rc2i        [3][3]float64 celestial-to-intermediate matrix (Note 3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

      date1             date2

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The Celestial Intermediate Pole coordinates are the x,y components
    of the unit vector in the Geocentric Celestial Reference System.

 3) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

          = RC2T * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

 4) Although its name does not include "00", This function is in fact
    specific to the IAU 2000 models.

Called:
    C2ixys    celestial-to-intermediate matrix, given X,Y and s
    S00       the CIO locator s, given X,Y, IAU 2000A

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2ixy(date1, date2 float64, x, y float64, rc2i *[3][3]float64) {
	/* Compute s and then the matrix. */
	C2ixys(x, y, S00(date1, date2, x, y), rc2i)
}

/*
C2ixys Celestial-to-intermediate matrix given CIP and s
Form the celestial to intermediate-frame-of-date matrix given the CIP
X,Y and the CIO locator s.

Given:
    x,y      float64         Celestial Intermediate Pole (Note 1)
    s        float64         the CIO locator s (Note 2)

Returned:
    rc2i     [3][3]float64   celestial-to-intermediate matrix (Note 3)

Notes:

 1) The Celestial Intermediate Pole coordinates are the x,y
    components of the unit vector in the Geocentric Celestial
    Reference System.

 2) The CIO locator s (in radians) positions the Celestial
    Intermediate Origin on the equator of the CIP.

 3) The matrix rc2i is the first stage in the transformation from
    celestial to terrestrial coordinates:

    [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

          = RC2T * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

Called:
    Ir        initialize r-matrix to identity
    Rz        rotate around Z-axis
    Ry        rotate around Y-axis

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2ixys(x, y, s float64, rc2i *[3][3]float64) {
	var r2, e, d float64

	/* Obtain the spherical angles E and d. */
	r2 = x*x + y*y
	// e = (r2 > 0.0) ? atan2(y, x) : 0.0;
	if r2 > 0.0 {
		e = atan2(y, x)
	} else {
		e = 0.0
	}
	d = atan(sqrt(r2 / (1.0 - r2)))

	/* Form the matrix. */
	Ir(rc2i)
	Rz(e, rc2i)
	Ry(d, rc2i)
	Rz(-(e + s), rc2i)
}

/*
C2t00a Celestial-to-terrestrial matrix, IAU 2000A
Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000A precession-nutation model.

Given:
    tta,ttb  float64         TT as a 2-part Julian Date (Note 1)
    uta,utb  float64         UT1 as a 2-part Julian Date (Note 1)
    xp,yp    float64         CIP coordinates (radians, Note 2)

Returned:
    rc2t     [3][3]float64   celestial-to-terrestrial matrix (Note 3)

Notes:

 1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
    apportioned in any convenient way between the arguments uta and
    utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
    these ways, among others:

         uta            utb

     2450123.7           0.0       (JD method)
     2451545.0       -1421.3       (J2000 method)
     2400000.5       50123.2       (MJD method)
     2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  In the case of uta,utb, the
    date & time method is best matched to the Earth rotation angle
    algorithm used:  maximum precision is delivered when the uta
    argument is for 0hrs UT1 on the day in question and the utb
    argument lies in the range 0 to 1, or vice versa.

 2) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 3) The matrix rc2t transforms from celestial to terrestrial
    coordinates:

    [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

          = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), RC2I is the
    celestial-to-intermediate matrix, ERA is the Earth rotation
    angle and RPOM is the polar motion matrix.

 4) A faster, but slightly less accurate, result (about 1 mas) can
    be obtained by using instead the iauC2t00b function.

Called:
    C2i00a    celestial-to-intermediate matrix, IAU 2000A
    Era00     Earth rotation angle, IAU 2000
    Sp00      the TIO locator s', IERS 2000
    Pom00     polar motion matrix
    C2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2t00a(tta, ttb float64, uta, utb float64, xp, yp float64, rc2t *[3][3]float64) {
	var rc2i, rpom [3][3]float64
	var era, sp float64

	/* Form the celestial-to-intermediate matrix for this TT (IAU 2000A). */
	C2i00a(tta, ttb, &rc2i)

	/* Predict the Earth rotation angle for this UT1. */
	era = Era00(uta, utb)

	/* Estimate s'. */
	sp = Sp00(tta, ttb)

	/* Form the polar motion matrix. */
	Pom00(xp, yp, sp, &rpom)

	/* Combine to form the celestial-to-terrestrial matrix. */
	C2tcio(rc2i, era, rpom, rc2t)
}

/*
C2t00b Celestial-to-terrestrial matrix, IAU 2000B
Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000B precession-nutation model.

Given:
    tta,ttb  float64         TT as a 2-part Julian Date (Note 1)
    uta,utb  float64         UT1 as a 2-part Julian Date (Note 1)
    xp,yp    float64         coordinates of the pole (radians, Note 2)

Returned:
    rc2t     [3][3]float64   celestial-to-terrestrial matrix (Note 3)

Notes:

 1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
    apportioned in any convenient way between the arguments uta and
    utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
    these ways, among others:

          uta            utb

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  In the case of uta,utb, the
    date & time method is best matched to the Earth rotation angle
    algorithm used:  maximum precision is delivered when the uta
    argument is for 0hrs UT1 on the day in question and the utb
    argument lies in the range 0 to 1, or vice versa.

 2) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 3) The matrix rc2t transforms from celestial to terrestrial
    coordinates:

    [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

          = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), RC2I is the
    celestial-to-intermediate matrix, ERA is the Earth rotation
    angle and RPOM is the polar motion matrix.

 4) The present function is faster, but slightly less accurate (about
    1 mas), than the iauC2t00a function.

Called:
    C2i00b    celestial-to-intermediate matrix, IAU 2000B
    Era00     Earth rotation angle, IAU 2000
    Pom00     polar motion matrix
    C2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2t00b(tta, ttb float64, uta, utb float64, xp, yp float64, rc2t *[3][3]float64) {
	var rc2i, rpom [3][3]float64
	var era float64

	/* Form the celestial-to-intermediate matrix for this TT (IAU 2000B). */
	C2i00b(tta, ttb, &rc2i)

	/* Predict the Earth rotation angle for this UT1. */
	era = Era00(uta, utb)

	/* Form the polar motion matrix (neglecting s'). */
	Pom00(xp, yp, 0.0, &rpom)

	/* Combine to form the celestial-to-terrestrial matrix. */
	C2tcio(rc2i, era, rpom, rc2t)
}

/*
C2t06a Celestial-to-terrestrial matrix, IAU 2006/2000A
Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2006/2000A precession-nutation
model.

Given:
    tta,ttb  float64         TT as a 2-part Julian Date (Note 1)
    uta,utb  float64         UT1 as a 2-part Julian Date (Note 1)
    xp,yp    float64         coordinates of the pole (radians, Note 2)

Returned:
    rc2t     [3][3]float64   celestial-to-terrestrial matrix (Note 3)

Notes:

 1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
    apportioned in any convenient way between the arguments uta and
    utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
    these ways, among others:

          uta            utb

      2450123.7           0.0       (JD method)
      2451545.0       -1421.3       (J2000 method)
      2400000.5       50123.2       (MJD method)
      2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  In the case of uta,utb, the
    date & time method is best matched to the Earth rotation angle
    algorithm used:  maximum precision is delivered when the uta
    argument is for 0hrs UT1 on the day in question and the utb
    argument lies in the range 0 to 1, or vice versa.

 2) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 3) The matrix rc2t transforms from celestial to terrestrial
    coordinates:

    [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

          = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), RC2I is the
    celestial-to-intermediate matrix, ERA is the Earth rotation
    angle and RPOM is the polar motion matrix.

Called:
    C2i06a    celestial-to-intermediate matrix, IAU 2006/2000A
    Era00     Earth rotation angle, IAU 2000
    Sp00      the TIO locator s', IERS 2000
    Pom00     polar motion matrix
    C2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG
*/
func C2t06a(tta, ttb float64, uta, utb float64, xp, yp float64, rc2t *[3][3]float64) {
	var rc2i, rpom [3][3]float64
	var era, sp float64

	/* Form the celestial-to-intermediate matrix for this TT. */
	C2i06a(tta, ttb, &rc2i)

	/* Predict the Earth rotation angle for this UT1. */
	era = Era00(uta, utb)

	/* Estimate s'. */
	sp = Sp00(tta, ttb)

	/* Form the polar motion matrix. */
	Pom00(xp, yp, sp, &rpom)

	/* Combine to form the celestial-to-terrestrial matrix. */
	C2tcio(rc2i, era, rpom, rc2t)
}

/*
C2tcio Form CIO-based Celestial-to-terrestrial matrix
Assemble the celestial to terrestrial matrix from CIO-based
components (the celestial-to-intermediate matrix, the Earth Rotation
Angle and the polar motion matrix).

Given:
    rc2i     [3][3]float64    celestial-to-intermediate matrix
    era      float64          Earth rotation angle (radians)
    rpom     [3][3]float64    polar-motion matrix

Returned:
    rc2t     [3][3]float64    celestial-to-terrestrial matrix

Notes:

 1) This function constructs the rotation matrix that transforms
    vectors in the celestial system into vectors in the terrestrial
    system.  It does so starting from precomputed components, namely
    the matrix which rotates from celestial coordinates to the
    intermediate frame, the Earth rotation angle and the polar motion
    matrix.  One use of the present function is when generating a
    series of celestial-to-terrestrial matrices where only the Earth
    Rotation Angle changes, avoiding the considerable overhead of
    recomputing the precession-nutation more often than necessary to
    achieve given accuracy objectives.

 2) The relationship between the arguments is as follows:

      [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

            = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003).

Called:
    Cr        copy r-matrix
    Rz        rotate around Z-axis
    Rxr       product of two r-matrices

Reference:

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG
*/
func C2tcio(rc2i [3][3]float64, era float64, rpom [3][3]float64, rc2t *[3][3]float64) {
	var r [3][3]float64

	/* Construct the matrix. */
	Cr(rc2i, &r)
	Rz(era, &r)
	Rxr(rpom, r, rc2t)
}

/*
C2teqx Celestial-to-terrestrial matrix, classical
Assemble the celestial to terrestrial matrix from equinox-based
components (the celestial-to-true matrix, the Greenwich Apparent
Sidereal Time and the polar motion matrix).

Given:
       rbpn   [3][3]float64  celestial-to-true matrix
       gst    float64        Greenwich (apparent) Sidereal Time (radians)
       rpom   [3][3]float64  polar-motion matrix

Returned:
    rc2t   [3][3]float64  celestial-to-terrestrial matrix (Note 2)

Notes:

 1) This function constructs the rotation matrix that transforms
    vectors in the celestial system into vectors in the terrestrial
    system.  It does so starting from precomputed components, namely
    the matrix which rotates from celestial coordinates to the
    true equator and equinox of date, the Greenwich Apparent Sidereal
    Time and the polar motion matrix.  One use of the present function
    is when generating a series of celestial-to-terrestrial matrices
    where only the Sidereal Time changes, avoiding the considerable
    overhead of recomputing the precession-nutation more often than
    necessary to achieve given accuracy objectives.

 2) The relationship between the arguments is as follows:

     [TRS] = rpom * R_3(gst) * rbpn * [CRS]

           = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003).

Called:
    Cr        copy r-matrix
    Rz        rotate around Z-axis
    Rxr       product of two r-matrices

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2teqx(rbpn [3][3]float64, gst float64, rpom [3][3]float64, rc2t *[3][3]float64) {
	var r [3][3]float64

	/* Construct the matrix. */
	Cr(rbpn, &r)
	Rz(gst, &r)
	Rxr(rpom, r, rc2t)
}

/*
C2tpe Celestial-to-terrestrial matrix given nutation
Form the celestial to terrestrial matrix given the date, the UT1,
the nutation and the polar motion.  IAU 2000.

Given:
    tta,ttb    float64        TT as a 2-part Julian Date (Note 1)
    uta,utb    float64        UT1 as a 2-part Julian Date (Note 1)
    dpsi,deps  float64        nutation (Note 2)
    xp,yp      float64        coordinates of the pole (radians, Note 3)

Returned:
    rc2t       [3][3]float64  celestial-to-terrestrial matrix (Note 4)

Notes:

 1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
    apportioned in any convenient way between the arguments uta and
    utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
    these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  In the case of uta,utb, the
    date & time method is best matched to the Earth rotation angle
    algorithm used:  maximum precision is delivered when the uta
    argument is for 0hrs UT1 on the day in question and the utb
    argument lies in the range 0 to 1, or vice versa.

 2) The caller is responsible for providing the nutation components;
    they are in longitude and obliquity, in radians and are with
    respect to the equinox and ecliptic of date.  For high-accuracy
    applications, free core nutation should be included as well as
    any other relevant corrections to the position of the CIP.

 3) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 4) The matrix rc2t transforms from celestial to terrestrial
    coordinates:

     [TRS] = RPOM * R_3(GST) * RBPN * [CRS]

           = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), RBPN is the
    bias-precession-nutation matrix, GST is the Greenwich (apparent)
    Sidereal Time and RPOM is the polar motion matrix.

 5) Although its name does not include "00", This function is in fact
    specific to the IAU 2000 models.

Called:
    Pn00      bias/precession/nutation results, IAU 2000
    Gmst00    Greenwich mean sidereal time, IAU 2000
    Sp00      the TIO locator s', IERS 2000
    Ee00      equation of the equinoxes, IAU 2000
    Pom00     polar motion matrix
    C2teqx    form equinox-based celestial-to-terrestrial matrix

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2tpe(tta, ttb float64, uta, utb float64, dpsi, deps float64, xp, yp float64, rc2t *[3][3]float64) {
	var epsa, gmst, ee, sp float64
	var rb, rp, rbp, rn, rbpn, rpom [3][3]float64

	/* Form the celestial-to-true matrix for this TT. */
	Pn00(tta, ttb, dpsi, deps, &epsa, &rb, &rp, &rbp, &rn, &rbpn)

	/* Predict the Greenwich Mean Sidereal Time for this UT1 and TT. */
	gmst = Gmst00(uta, utb, tta, ttb)

	/* Predict the equation of the equinoxes given TT and nutation. */
	ee = Ee00(tta, ttb, epsa, dpsi)

	/* Estimate s'. */
	sp = Sp00(tta, ttb)

	/* Form the polar motion matrix. */
	Pom00(xp, yp, sp, &rpom)

	/* Combine to form the celestial-to-terrestrial matrix. */
	C2teqx(rbpn, gmst+ee, rpom, rc2t)
}

/*
C2txy Celestial-to-terrestrial matrix given CIP
Form the celestial to terrestrial matrix given the date, the UT1,
the CIP coordinates and the polar motion.  IAU 2000.

Given:
    tta,ttb  float64         TT as a 2-part Julian Date (Note 1)
    uta,utb  float64         UT1 as a 2-part Julian Date (Note 1)
    x,y      float64         Celestial Intermediate Pole (Note 2)
    xp,yp    float64         coordinates of the pole (radians, Note 3)

Returned:
    rc2t     [3][3]float64   celestial-to-terrestrial matrix (Note 4)

Notes:

 1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
    apportioned in any convenient way between the arguments uta and
    utb.  For example, JD(UT1)=2450123.7 could be expressed in any o
    these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  In the case of uta,utb, the
    date & time method is best matched to the Earth rotation angle
    algorithm used:  maximum precision is delivered when the uta
    argument is for 0hrs UT1 on the day in question and the utb
    argument lies in the range 0 to 1, or vice versa.

 2) The Celestial Intermediate Pole coordinates are the x,y
    components of the unit vector in the Geocentric Celestial
    Reference System.

 3) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 4) The matrix rc2t transforms from celestial to terrestrial
    coordinates:

      [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

            = rc2t * [CRS]

    where [CRS] is a vector in the Geocentric Celestial Reference
    System and [TRS] is a vector in the International Terrestrial
    Reference System (see IERS Conventions 2003), ERA is the Earth
    Rotation Angle and RPOM is the polar motion matrix.

 5) Although its name does not include "00", This function is in fact
    specific to the IAU 2000 models.

Called:
    C2ixy     celestial-to-intermediate matrix, given X,Y
    Era00     Earth rotation angle, IAU 2000
    Sp00      the TIO locator s', IERS 2000
    Pom00     polar motion matrix
    C2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func C2txy(tta, ttb float64, uta, utb float64, x, y float64, xp, yp float64, rc2t *[3][3]float64) {
	var rc2i, rpom [3][3]float64
	var era, sp float64

	/* Form the celestial-to-intermediate matrix for this TT. */
	C2ixy(tta, ttb, x, y, &rc2i)

	/* Predict the Earth rotation angle for this UT1. */
	era = Era00(uta, utb)

	/* Estimate s'. */
	sp = Sp00(tta, ttb)

	/* Form the polar motion matrix. */
	Pom00(xp, yp, sp, &rpom)

	/* Combine to form the celestial-to-terrestrial matrix. */
	C2tcio(rc2i, era, rpom, rc2t)
}

/*
Eo06a Equation of the origins, IAU 2006/2000A
Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    the equation of the origins in radians

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The equation of the origins is the distance between the true
    equinox and the celestial intermediate origin and, equivalently,
    the difference between Earth rotation angle and Greenwich
    apparent sidereal time (ERA-GST).  It comprises the precession
    (since J2000.0) in right ascension plus the equation of the
    equinoxes (including the small correction terms).

Called:
    Pnm06a    classical NPB matrix, IAU 2006/2000A
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S06       the CIO locator s, given X,Y, IAU 2006
    Eors      equation of the origins, given NPB matrix and s

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Eo06a(date1, date2 float64) float64 {
	var r [3][3]float64
	var x, y, s, eo float64

	/* Classical nutation x precession x bias matrix. */
	Pnm06a(date1, date2, &r)

	/* Extract CIP coordinates. */
	Bpn2xy(r, &x, &y)

	/* The CIO locator, s. */
	s = S06(date1, date2, x, y)

	/* Solve for the EO. */
	eo = Eors(r, s)

	return eo

}

/*
Eors Equation of the origins, given NPB matrix and s
Equation of the origins, given the classical NPB matrix and the
quantity s.

Given:
    rnpb  [3][3]float64  classical nutation x precession x bias matrix
    s     float64        the quantity s (the CIO locator) in radians

Returned (function value):
    float64        the equation of the origins in radians

Notes:

 1) The equation of the origins is the distance between the true
    equinox and the celestial intermediate origin and, equivalently,
    the difference between Earth rotation angle and Greenwich
    apparent sidereal time (ERA-GST).  It comprises the precession
    (since J2000.0) in right ascension plus the equation of the
    equinoxes (including the small correction terms).

 2) The algorithm is from Wallace & Capitaine (2006).

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Eors(rnpb [3][3]float64, s float64) float64 {
	var x, ax, xs, ys, zs, p, q, eo float64

	/* Evaluate Wallace & Capitaine (2006) expression (16). */
	x = rnpb[2][0]
	ax = x / (1.0 + rnpb[2][2])
	xs = 1.0 - ax*x
	ys = -ax * rnpb[2][1]
	zs = -x
	p = rnpb[0][0]*xs + rnpb[0][1]*ys + rnpb[0][2]*zs
	q = rnpb[1][0]*xs + rnpb[1][1]*ys + rnpb[1][2]*zs
	// eo = ((p != 0) || (q != 0)) ? s - atan2(q, p) : s;
	if (p != 0) || (q != 0) {
		eo = s - atan2(q, p)
	} else {
		eo = s
	}

	return eo
}

/*
Fw2m Fukushima-Williams angles to r-matrix
Form rotation matrix given the Fukushima-Williams angles.

Given:
    gamb     float64         F-W angle gamma_bar (radians)
    phib     float64         F-W angle phi_bar (radians)
    psi      float64         F-W angle psi (radians)
    eps      float64         F-W angle epsilon (radians)

Returned:
    r        [3][3]float64   rotation matrix

Notes:

 1) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole,
         E = ecliptic pole of date,
    and   P = CIP,

    the four Fukushima-Williams angles are as follows:

      gamb = gamma = epE
      phib = phi = pE
      psi = psi = pEP
      eps = epsilon = EP

 2) The matrix representing the combined effects of frame bias,
    precession and nutation is:

      NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)

 3) The present function can construct three different matrices,
    depending on which angles are supplied as the arguments gamb,
    phib, psi and eps:

    o  To obtain the nutation x precession x frame bias matrix,
      first generate the four precession angles known conventionally
      as gamma_bar, phi_bar, psi_bar and epsilon_A, then generate
      the nutation components Dpsi and Depsilon and add them to
      psi_bar and epsilon_A, and finally call the present function
      using those four angles as arguments.

    o  To obtain the precession x frame bias matrix, generate the
      four precession angles and call the present function.

    o  To obtain the frame bias matrix, generate the four precession
      angles for date J2000.0 and call the present function.

    The nutation-only and precession-only matrices can if necessary
    be obtained by combining these three appropriately.

Called:
    Ir        initialize r-matrix to identity
    Rz        rotate around Z-axis
    Rx        rotate around X-axis

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
*/
func Fw2m(gamb, phib, psi, eps float64, r *[3][3]float64) {
	Ir(r)
	Rz(gamb, r)
	Rx(phib, r)
	Rz(-psi, r)
	Rx(-eps, r)
}

/*
Fw2xy Fukushima-Williams angles to X,Y
CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

Given:
    gamb     float64    F-W angle gamma_bar (radians)
    phib     float64    F-W angle phi_bar (radians)
    psi      float64    F-W angle psi (radians)
    eps      float64    F-W angle epsilon (radians)

Returned:
    x,y      float64    CIP unit vector X,Y

Notes:

 1) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole
         E = ecliptic pole of date,
    and   P = CIP,

    the four Fukushima-Williams angles are as follows:

      gamb = gamma = epE
      phib = phi = pE
      psi = psi = pEP
      eps = epsilon = EP

 2) The matrix representing the combined effects of frame bias,
    precession and nutation is:

      NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)

    The returned values x,y are elements [2][0] and [2][1] of the
    matrix.  Near J2000.0, they are essentially angles in radians.

Called:
    Fw2m      F-W angles to r-matrix
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix

Reference:

    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
*/
func Fw2xy(gamb, phib float64, psi, eps float64, x, y *float64) {
	var r [3][3]float64

	/* Form NxPxB matrix. */
	Fw2m(gamb, phib, psi, eps, &r)

	/* Extract CIP X,Y. */
	Bpn2xy(r, x, y)
}

/*
Ltp Long-term precession matrix.

Given:
    epj     float64         Julian epoch (TT)

Returned:
    rp      [3][3]float64   precession matrix, J2000.0 to date

Notes:

 1) The matrix is in the sense

      P_date = rp x P_J2000,

    where P_J2000 is a vector with respect to the J2000.0 mean
    equator and equinox and P_date is the same vector with respect to
    the equator and equinox of epoch epj.

 2) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

Called:
    Ltpequ    equator pole, long term
    Ltpecl    ecliptic pole, long term
    Pxp       vector product
    Pn        normalize vector

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Ltp(epj float64, rp *[3][3]float64) {

	var peqr, pecl, v, eqx [3]float64
	var w float64

	/* Equator pole (bottom row of matrix). */
	Ltpequ(epj, &peqr)

	/* Ecliptic pole. */
	Ltpecl(epj, &pecl)

	/* Equinox (top row of matrix). */
	Pxp(peqr, pecl, &v)
	Pn(v, &w, &eqx)

	/* Middle row of matrix. */
	Pxp(peqr, eqx, &v)

	/* Assemble the matrix. */
	for i := 0; i < 3; i++ {
		rp[0][i] = eqx[i]
		rp[1][i] = v[i]
		rp[2][i] = peqr[i]
	}

}

/*
Ltpb Long-term precession matrix, including ICRS frame bias.

Given:
    epj     float64         Julian epoch (TT)

Returned:
    rpb     [3][3]float64   precession-bias matrix, J2000.0 to date

Notes:

 1) The matrix is in the sense

      P_date = rpb x P_ICRS,

    where P_ICRS is a vector in the Geocentric Celestial Reference
    System, and P_date is the vector with respect to the Celestial
    Intermediate Reference System at that date but with nutation
    neglected.

 2) A first order frame bias formulation is used, of sub-
    microarcsecond accuracy compared with a full 3D rotation.

 3) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Ltpb(epj float64, rpb *[3][3]float64) {
	/* Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33) */
	dx := -0.016617 * DAS2R
	de := -0.0068192 * DAS2R
	dr := -0.0146 * DAS2R

	var rp [3][3]float64

	/* Precession matrix. */
	Ltp(epj, &rp)

	/* Apply the bias. */
	for i := 0; i < 3; i++ {
		rpb[i][0] = rp[i][0] - rp[i][1]*dr + rp[i][2]*dx
		rpb[i][1] = rp[i][0]*dr + rp[i][1] + rp[i][2]*de
		rpb[i][2] = -rp[i][0]*dx - rp[i][1]*de + rp[i][2]
	}
}

/*
Ltpecl Long-term precession of the ecliptic.

Given:
    epj     float64         Julian epoch (TT)

Returned:
    vec     [3]float64      ecliptic pole unit vector

Notes:

 1) The returned vector is with respect to the J2000.0 mean equator
    and equinox.

 2) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Ltpecl(epj float64, vec *[3]float64) {

	/* Obliquity at J2000.0 (radians). */
	const eps0 = 84381.406 * DAS2R

	/* Polynomial coefficients */
	const NPOL = 4
	pqpol := [2][4]float64{
		{5851.607687, -0.1189000, -0.00028913, 0.000000101},
		{-1600.886300, 1.1689818, -0.00000020, -0.000000437},
	}

	/* Periodic coefficients */
	pqper := [][5]float64{
		{708.15, -5486.751211, -684.661560, 667.666730, -5523.863691},
		{2309.00, -17.127623, 2446.283880, -2354.886252, -549.747450},
		{1620.00, -617.517403, 399.671049, -428.152441, -310.998056},
		{492.20, 413.442940, -356.652376, 376.202861, 421.535876},
		{1183.00, 78.614193, -186.387003, 184.778874, -36.776172},
		{622.00, -180.732815, -316.800070, 335.321713, -145.278396},
		{882.00, -87.676083, 198.296701, -185.138669, -34.744450},
		{547.00, 46.140315, 101.135679, -120.972830, 22.885731},
	}
	NPER := len(pqper)

	/* Miscellaneous */

	var t, p, q, w, a, s, c float64

	/* Centuries since J2000. */
	t = (epj - 2000.0) / 100.0

	/* Initialize P_A and Q_A accumulators. */
	p = 0.0
	q = 0.0

	/* Periodic terms. */
	w = D2PI * t
	for i := 0; i < NPER; i++ {
		a = w / pqper[i][0]
		s = sin(a)
		c = cos(a)
		p += c*pqper[i][1] + s*pqper[i][3]
		q += c*pqper[i][2] + s*pqper[i][4]
	}

	/* Polynomial terms. */
	w = 1.0
	for i := 0; i < NPOL; i++ {
		p += pqpol[0][i] * w
		q += pqpol[1][i] * w
		w *= t
	}

	/* P_A and Q_A (radians). */
	p *= DAS2R
	q *= DAS2R

	/* Form the ecliptic pole vector. */
	w = 1.0 - p*p - q*q
	if w < 0.0 {
		w = 0.0
	} else {
		w = sqrt(w)
	}

	s = sin(eps0)
	c = cos(eps0)
	vec[0] = p
	vec[1] = -q*c - w*s
	vec[2] = -q*s + w*c
}

/*
Ltpequ Long-term precession of the equator.

Given:
    epj     float64         Julian epoch (TT)

Returned:
    veq     [3]float64      equator pole unit vector

Notes:

 1) The returned vector is with respect to the J2000.0 mean equator
    and equinox.

 2) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Ltpequ(epj float64, veq *[3]float64) {
	/* Polynomial coefficients */
	const NPOL = 4
	xypol := [2][4]float64{
		{5453.282155, 0.4252841, -0.00037173, -0.000000152},
		{-73750.930350, -0.7675452, -0.00018725, 0.000000231},
	}

	/* Periodic coefficients */
	xyper := [][5]float64{
		{256.75, -819.940624, 75004.344875, 81491.287984, 1558.515853},
		{708.15, -8444.676815, 624.033993, 787.163481, 7774.939698},
		{274.20, 2600.009459, 1251.136893, 1251.296102, -2219.534038},
		{241.45, 2755.175630, -1102.212834, -1257.950837, -2523.969396},
		{2309.00, -167.659835, -2660.664980, -2966.799730, 247.850422},
		{492.20, 871.855056, 699.291817, 639.744522, -846.485643},
		{396.10, 44.769698, 153.167220, 131.600209, -1393.124055},
		{288.90, -512.313065, -950.865637, -445.040117, 368.526116},
		{231.10, -819.415595, 499.754645, 584.522874, 749.045012},
		{1610.00, -538.071099, -145.188210, -89.756563, 444.704518},
		{620.00, -189.793622, 558.116553, 524.429630, 235.934465},
		{157.87, -402.922932, -23.923029, -13.549067, 374.049623},
		{220.30, 179.516345, -165.405086, -210.157124, -171.330180},
		{1200.00, -9.814756, 9.344131, -44.919798, -22.899655},
	}
	NPER := len(xyper)

	/* Miscellaneous */

	var t, x, y, w, a, s, c float64

	/* Centuries since J2000. */
	t = (epj - 2000.0) / 100.0

	/* Initialize X and Y accumulators. */
	x = 0.0
	y = 0.0

	/* Periodic terms. */
	w = D2PI * t
	for i := 0; i < NPER; i++ {
		a = w / xyper[i][0]
		s = sin(a)
		c = cos(a)
		x += c*xyper[i][1] + s*xyper[i][3]
		y += c*xyper[i][2] + s*xyper[i][4]
	}

	/* Polynomial terms. */
	w = 1.0
	for i := 0; i < NPOL; i++ {
		x += xypol[0][i] * w
		y += xypol[1][i] * w
		w *= t
	}

	/* X and Y (direction cosines). */
	x *= DAS2R
	y *= DAS2R

	/* Form the equator pole vector. */
	veq[0] = x
	veq[1] = y
	w = 1.0 - x*x - y*y
	if w < 0.0 {
		veq[2] = 0.0
	} else {
		veq[2] = sqrt(w)
	}

}

/*
Num00a Nutation matrix, IAU 2000A
Form the matrix of nutation for a given date, IAU 2000A model.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    rmatn        [3][3]float64    nutation matrix

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(true) = rmatn * V(mean), where
    the p-vector V(true) is with respect to the true equatorial triad
    of date and the p-vector V(mean) is with respect to the mean
    equatorial triad of date.

 3) A faster, but slightly less accurate, result (about 1 mas) can be
    obtained by using instead the iauNum00b function.

Called:
    Pn00a     bias/precession/nutation, IAU 2000A

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.222-3 (p114).
*/
func Num00a(date1, date2 float64, rmatn *[3][3]float64) {
	var dpsi, deps, epsa float64
	var rb, rp, rbp, rbpn [3][3]float64

	/* Obtain the required matrix (discarding other results). */
	Pn00a(date1, date2,
		&dpsi, &deps, &epsa, &rb, &rp, &rbp, rmatn, &rbpn)
}

/*
Num00b Nutation matrix, IAU 2000B
Form the matrix of nutation for a given date, IAU 2000B model.

Given:
    date1,date2  float64         TT as a 2-part Julian Date (Note 1)

Returned:
    rmatn        [3][3]float64   nutation matrix

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(true) = rmatn * V(mean), where
    the p-vector V(true) is with respect to the true equatorial triad
    of date and the p-vector V(mean) is with respect to the mean
    equatorial triad of date.

 3) The present function is faster, but slightly less accurate (about
    1 mas), than the iauNum00a function.

Called:
    Pn00b     bias/precession/nutation, IAU 2000B

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.222-3 (p114).
*/
func Num00b(date1, date2 float64, rmatn *[3][3]float64) {
	var dpsi, deps, epsa float64
	var rb, rp, rbp, rbpn [3][3]float64

	/* Obtain the required matrix (discarding other results). */
	Pn00b(date1, date2,
		&dpsi, &deps, &epsa, &rb, &rp, &rbp, rmatn, &rbpn)
}

/*
Num06a Nutation matrix, IAU 2006/2000A
Form the matrix of nutation for a given date, IAU 2006/2000A model.

Given:
    date1,date2   float64          TT as a 2-part Julian Date (Note 1)

Returned:
    rmatn         [3][3]float64    nutation matrix

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(true) = rmatn * V(mean), where
    the p-vector V(true) is with respect to the true equatorial triad
    of date and the p-vector V(mean) is with respect to the mean
    equatorial triad of date.

Called:
    Obl06     mean obliquity, IAU 2006
    Nut06a    nutation, IAU 2006/2000A
    Numat     form nutation matrix

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.222-3 (p114).
*/
func Num06a(date1, date2 float64, rmatn *[3][3]float64) {
	var eps, dp, de float64

	/* Mean obliquity. */
	eps = Obl06(date1, date2)

	/* Nutation components. */
	Nut06a(date1, date2, &dp, &de)

	/* Nutation matrix. */
	Numat(eps, dp, de, rmatn)

}

/*
Numat Nutation matrix, generic Form the matrix of nutation.

Given:
    epsa        float64         mean obliquity of date (Note 1)
    dpsi,deps   float64         nutation (Note 2)

Returned:
    rmatn       [3][3]float64   nutation matrix (Note 3)

Notes:


 1) The supplied mean obliquity epsa, must be consistent with the
    precession-nutation models from which dpsi and deps were obtained.

 2) The caller is responsible for providing the nutation components;
    they are in longitude and obliquity, in radians and are with
    respect to the equinox and ecliptic of date.

 3) The matrix operates in the sense V(true) = rmatn * V(mean),
    where the p-vector V(true) is with respect to the true
    equatorial triad of date and the p-vector V(mean) is with
    respect to the mean equatorial triad of date.

Called:
    Ir        initialize r-matrix to identity
    Rx        rotate around X-axis
    Rz        rotate around Z-axis

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.222-3 (p114).
*/
func Numat(epsa, dpsi, deps float64, rmatn *[3][3]float64) {
	/* Build the rotation matrix. */
	Ir(rmatn)
	Rx(epsa, rmatn)
	Rz(-dpsi, rmatn)
	Rx(-(epsa + deps), rmatn)
}

/*
Nut00a Nutation, IAU 2000A Nutation, IAU 2000A model
(MHB2000 luni-solar and planetary nutation with free core nutation omitted).

Given:
    date1,date2   float64   TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps     float64   nutation, luni-solar + planetary (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components in longitude and obliquity are in radians
    and with respect to the equinox and ecliptic of date.  The
    obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
    value of 84381.448 arcsec.

    Both the luni-solar and planetary nutations are included.  The
    latter are due to direct planetary nutations and the
    perturbations of the lunar and terrestrial orbits.

 3) The function computes the MHB2000 nutation series with the
    associated corrections for planetary nutations.  It is an
    implementation of the nutation part of the IAU 2000A precession-
    nutation model, formally adopted by the IAU General Assembly in
    2000, namely MHB2000 (Mathews et al. 2002), but with the free
    core nutation (FCN - see Note 4) omitted.

 4) The full MHB2000 model also contains contributions to the
    nutations in longitude and obliquity due to the free-excitation
    of the free-core-nutation during the period 1979-2000.  These FCN
    terms, which are time-dependent and unpredictable, are NOT
    included in the present function and, if required, must be
    independently computed.  With the FCN corrections included, the
    present function delivers a pole which is at current epochs
    accurate to a few hundred microarcseconds.  The omission of FCN
    introduces further errors of about that size.

 5) The present function provides classical nutation.  The MHB2000
    algorithm, from which it is adapted, deals also with (i) the
    offsets between the GCRS and mean poles and (ii) the adjustments
    in longitude and obliquity due to the changed precession rates.
    These additional functions, namely frame bias and precession
    adjustments, are supported by the SOFA functions iauBi00  and
    Pr00.

 6) The MHB2000 algorithm also provides "total" nutations, comprising
    the arithmetic sum of the frame bias, precession adjustments,
    luni-solar nutation and planetary nutation.  These total
    nutations can be used in combination with an existing IAU 1976
    precession implementation, such as iauPmat76,  to deliver GCRS-
    to-true predictions of sub-mas accuracy at current dates.
    However, there are three shortcomings in the MHB2000 model that
    must be taken into account if more accurate or definitive results
    are required (see Wallace 2002):

     (i) The MHB2000 total nutations are simply arithmetic sums,
         yet in reality the various components are successive Euler
         rotations.  This slight lack of rigor leads to cross terms
         that exceed 1 mas after a century.  The rigorous procedure
         is to form the GCRS-to-true rotation matrix by applying the
         bias, precession and nutation in that order.

    (ii) Although the precession adjustments are stated to be with
         respect to Lieske et al. (1977), the MHB2000 model does
         not specify which set of Euler angles are to be used and
         how the adjustments are to be applied.  The most literal
         and straightforward procedure is to adopt the 4-rotation
         epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
         to psi_A and DEPSPR to both omega_A and eps_A.

    (iii) The MHB2000 model predates the determination by Chapront
         et al. (2002) of a 14.6 mas displacement between the
         J2000.0 mean equinox and the origin of the ICRS frame.  It
         should, however, be noted that neglecting this displacement
         when calculating star coordinates does not lead to a
         14.6 mas change in right ascension, only a small second-
         order distortion in the pattern of the precession-nutation
         effect.

    For these reasons, the SOFA functions do not generate the "total
    nutations" directly, though they can of course easily be
    generated by calling iauBi00, iauPr00 and the present function
    and adding the results.

 7) The MHB2000 model contains 41 instances where the same frequency
    appears multiple times, of which 38 are duplicates and three are
    triplicates.  To keep the present code close to the original MHB
    algorithm, this small inefficiency has not been corrected.

Called:
    Fal03     mean anomaly of the Moon
    Faf03     mean argument of the latitude of the Moon
    Faom03    mean longitude of the Moon's ascending node
    Fame03    mean longitude of Mercury
    Fave03    mean longitude of Venus
    Fae03     mean longitude of Earth
    Fama03    mean longitude of Mars
    Faju03    mean longitude of Jupiter
    Fasa03    mean longitude of Saturn
    Faur03    mean longitude of Uranus
    Fapa03    general accumulated precession in longitude

References:

    Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
    Astron.Astrophys. 387, 700

    Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    Astron.Astrophys. 58, 1-16

    Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
    107, B4.  The MHB_2000 code itself was obtained on 9th September
    2002 from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

    Wallace, P.T., "Software for Implementing the IAU 2000
    Resolutions", in IERS Workshop 5.1 (2002)
*/
func Nut00a(date1, date2 float64, dpsi, deps *float64) {
	var i int
	var t, el, elp, f, d, om, arg, dp, de, sarg, carg,
		al, af, ad, aom, alme, alve, alea, alma,
		alju, alsa, alur, alne, apa, dpsils, depsls,
		dpsipl, depspl float64

	/* Units of 0.1 microarcsecond to radians */
	const U2R = DAS2R / 1e7

	/* ------------------------- */
	/* Luni-Solar nutation model */
	/* ------------------------- */

	/* The units for the sine and cosine coefficients are */
	/* 0.1 microarcsecond and the same per Julian century */

	xls := []struct {
		nl, nlp, nf, nd, nom int     /* coefficients of l,l',F,D,Om */
		sp, spt, cp          float64 /* longitude sin, t*sin, cos coefficients */
		ce, cet, se          float64 /* obliquity cos, t*cos, sin coefficients */
	}{

		/* 1- 10 */
		{0, 0, 0, 0, 1,
			-172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
		{0, 0, 2, -2, 2,
			-13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
		{0, 0, 2, 0, 2, -2276413.0, -234.0, 2796.0, 978459.0, -485.0, 1374.0},
		{0, 0, 0, 0, 2, 2074554.0, 207.0, -698.0, -897492.0, 470.0, -291.0},
		{0, 1, 0, 0, 0, 1475877.0, -3633.0, 11817.0, 73871.0, -184.0, -1924.0},
		{0, 1, 2, -2, 2, -516821.0, 1226.0, -524.0, 224386.0, -677.0, -174.0},
		{1, 0, 0, 0, 0, 711159.0, 73.0, -872.0, -6750.0, 0.0, 358.0},
		{0, 0, 2, 0, 1, -387298.0, -367.0, 380.0, 200728.0, 18.0, 318.0},
		{1, 0, 2, 0, 2, -301461.0, -36.0, 816.0, 129025.0, -63.0, 367.0},
		{0, -1, 2, -2, 2, 215829.0, -494.0, 111.0, -95929.0, 299.0, 132.0},

		/* 11-20 */
		{0, 0, 2, -2, 1, 128227.0, 137.0, 181.0, -68982.0, -9.0, 39.0},
		{-1, 0, 2, 0, 2, 123457.0, 11.0, 19.0, -53311.0, 32.0, -4.0},
		{-1, 0, 0, 2, 0, 156994.0, 10.0, -168.0, -1235.0, 0.0, 82.0},
		{1, 0, 0, 0, 1, 63110.0, 63.0, 27.0, -33228.0, 0.0, -9.0},
		{-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0, 31429.0, 0.0, -75.0},
		{-1, 0, 2, 2, 2, -59641.0, -11.0, 149.0, 25543.0, -11.0, 66.0},
		{1, 0, 2, 0, 1, -51613.0, -42.0, 129.0, 26366.0, 0.0, 78.0},
		{-2, 0, 2, 0, 1, 45893.0, 50.0, 31.0, -24236.0, -10.0, 20.0},
		{0, 0, 0, 2, 0, 63384.0, 11.0, -150.0, -1220.0, 0.0, 29.0},
		{0, 0, 2, 2, 2, -38571.0, -1.0, 158.0, 16452.0, -11.0, 68.0},

		/* 21-30 */
		{0, -2, 2, -2, 2, 32481.0, 0.0, 0.0, -13870.0, 0.0, 0.0},
		{-2, 0, 0, 2, 0, -47722.0, 0.0, -18.0, 477.0, 0.0, -25.0},
		{2, 0, 2, 0, 2, -31046.0, -1.0, 131.0, 13238.0, -11.0, 59.0},
		{1, 0, 2, -2, 2, 28593.0, 0.0, -1.0, -12338.0, 10.0, -3.0},
		{-1, 0, 2, 0, 1, 20441.0, 21.0, 10.0, -10758.0, 0.0, -3.0},
		{2, 0, 0, 0, 0, 29243.0, 0.0, -74.0, -609.0, 0.0, 13.0},
		{0, 0, 2, 0, 0, 25887.0, 0.0, -66.0, -550.0, 0.0, 11.0},
		{0, 1, 0, 0, 1, -14053.0, -25.0, 79.0, 8551.0, -2.0, -45.0},
		{-1, 0, 0, 2, 1, 15164.0, 10.0, 11.0, -8001.0, 0.0, -1.0},
		{0, 2, 2, -2, 2, -15794.0, 72.0, -16.0, 6850.0, -42.0, -5.0},

		/* 31-40 */
		{0, 0, -2, 2, 0, 21783.0, 0.0, 13.0, -167.0, 0.0, 13.0},
		{1, 0, 0, -2, 1, -12873.0, -10.0, -37.0, 6953.0, 0.0, -14.0},
		{0, -1, 0, 0, 1, -12654.0, 11.0, 63.0, 6415.0, 0.0, 26.0},
		{-1, 0, 2, 2, 1, -10204.0, 0.0, 25.0, 5222.0, 0.0, 15.0},
		{0, 2, 0, 0, 0, 16707.0, -85.0, -10.0, 168.0, -1.0, 10.0},
		{1, 0, 2, 2, 2, -7691.0, 0.0, 44.0, 3268.0, 0.0, 19.0},
		{-2, 0, 2, 0, 0, -11024.0, 0.0, -14.0, 104.0, 0.0, 2.0},
		{0, 1, 2, 0, 2, 7566.0, -21.0, -11.0, -3250.0, 0.0, -5.0},
		{0, 0, 2, 2, 1, -6637.0, -11.0, 25.0, 3353.0, 0.0, 14.0},
		{0, -1, 2, 0, 2, -7141.0, 21.0, 8.0, 3070.0, 0.0, 4.0},

		/* 41-50 */
		{0, 0, 0, 2, 1, -6302.0, -11.0, 2.0, 3272.0, 0.0, 4.0},
		{1, 0, 2, -2, 1, 5800.0, 10.0, 2.0, -3045.0, 0.0, -1.0},
		{2, 0, 2, -2, 2, 6443.0, 0.0, -7.0, -2768.0, 0.0, -4.0},
		{-2, 0, 0, 2, 1, -5774.0, -11.0, -15.0, 3041.0, 0.0, -5.0},
		{2, 0, 2, 0, 1, -5350.0, 0.0, 21.0, 2695.0, 0.0, 12.0},
		{0, -1, 2, -2, 1, -4752.0, -11.0, -3.0, 2719.0, 0.0, -3.0},
		{0, 0, 0, -2, 1, -4940.0, -11.0, -21.0, 2720.0, 0.0, -9.0},
		{-1, -1, 0, 2, 0, 7350.0, 0.0, -8.0, -51.0, 0.0, 4.0},
		{2, 0, 0, -2, 1, 4065.0, 0.0, 6.0, -2206.0, 0.0, 1.0},
		{1, 0, 0, 2, 0, 6579.0, 0.0, -24.0, -199.0, 0.0, 2.0},

		/* 51-60 */
		{0, 1, 2, -2, 1, 3579.0, 0.0, 5.0, -1900.0, 0.0, 1.0},
		{1, -1, 0, 0, 0, 4725.0, 0.0, -6.0, -41.0, 0.0, 3.0},
		{-2, 0, 2, 0, 2, -3075.0, 0.0, -2.0, 1313.0, 0.0, -1.0},
		{3, 0, 2, 0, 2, -2904.0, 0.0, 15.0, 1233.0, 0.0, 7.0},
		{0, -1, 0, 2, 0, 4348.0, 0.0, -10.0, -81.0, 0.0, 2.0},
		{1, -1, 2, 0, 2, -2878.0, 0.0, 8.0, 1232.0, 0.0, 4.0},
		{0, 0, 0, 1, 0, -4230.0, 0.0, 5.0, -20.0, 0.0, -2.0},
		{-1, -1, 2, 2, 2, -2819.0, 0.0, 7.0, 1207.0, 0.0, 3.0},
		{-1, 0, 2, 0, 0, -4056.0, 0.0, 5.0, 40.0, 0.0, -2.0},
		{0, -1, 2, 2, 2, -2647.0, 0.0, 11.0, 1129.0, 0.0, 5.0},

		/* 61-70 */
		{-2, 0, 0, 0, 1, -2294.0, 0.0, -10.0, 1266.0, 0.0, -4.0},
		{1, 1, 2, 0, 2, 2481.0, 0.0, -7.0, -1062.0, 0.0, -3.0},
		{2, 0, 0, 0, 1, 2179.0, 0.0, -2.0, -1129.0, 0.0, -2.0},
		{-1, 1, 0, 1, 0, 3276.0, 0.0, 1.0, -9.0, 0.0, 0.0},
		{1, 1, 0, 0, 0, -3389.0, 0.0, 5.0, 35.0, 0.0, -2.0},
		{1, 0, 2, 0, 0, 3339.0, 0.0, -13.0, -107.0, 0.0, 1.0},
		{-1, 0, 2, -2, 1, -1987.0, 0.0, -6.0, 1073.0, 0.0, -2.0},
		{1, 0, 0, 0, 2, -1981.0, 0.0, 0.0, 854.0, 0.0, 0.0},
		{-1, 0, 0, 1, 0, 4026.0, 0.0, -353.0, -553.0, 0.0, -139.0},
		{0, 0, 2, 1, 2, 1660.0, 0.0, -5.0, -710.0, 0.0, -2.0},

		/* 71-80 */
		{-1, 0, 2, 4, 2, -1521.0, 0.0, 9.0, 647.0, 0.0, 4.0},
		{-1, 1, 0, 1, 1, 1314.0, 0.0, 0.0, -700.0, 0.0, 0.0},
		{0, -2, 2, -2, 1, -1283.0, 0.0, 0.0, 672.0, 0.0, 0.0},
		{1, 0, 2, 2, 1, -1331.0, 0.0, 8.0, 663.0, 0.0, 4.0},
		{-2, 0, 2, 2, 2, 1383.0, 0.0, -2.0, -594.0, 0.0, -2.0},
		{-1, 0, 0, 0, 2, 1405.0, 0.0, 4.0, -610.0, 0.0, 2.0},
		{1, 1, 2, -2, 2, 1290.0, 0.0, 0.0, -556.0, 0.0, 0.0},
		{-2, 0, 2, 4, 2, -1214.0, 0.0, 5.0, 518.0, 0.0, 2.0},
		{-1, 0, 4, 0, 2, 1146.0, 0.0, -3.0, -490.0, 0.0, -1.0},
		{2, 0, 2, -2, 1, 1019.0, 0.0, -1.0, -527.0, 0.0, -1.0},

		/* 81-90 */
		{2, 0, 2, 2, 2, -1100.0, 0.0, 9.0, 465.0, 0.0, 4.0},
		{1, 0, 0, 2, 1, -970.0, 0.0, 2.0, 496.0, 0.0, 1.0},
		{3, 0, 0, 0, 0, 1575.0, 0.0, -6.0, -50.0, 0.0, 0.0},
		{3, 0, 2, -2, 2, 934.0, 0.0, -3.0, -399.0, 0.0, -1.0},
		{0, 0, 4, -2, 2, 922.0, 0.0, -1.0, -395.0, 0.0, -1.0},
		{0, 1, 2, 0, 1, 815.0, 0.0, -1.0, -422.0, 0.0, -1.0},
		{0, 0, -2, 2, 1, 834.0, 0.0, 2.0, -440.0, 0.0, 1.0},
		{0, 0, 2, -2, 3, 1248.0, 0.0, 0.0, -170.0, 0.0, 1.0},
		{-1, 0, 0, 4, 0, 1338.0, 0.0, -5.0, -39.0, 0.0, 0.0},
		{2, 0, -2, 0, 1, 716.0, 0.0, -2.0, -389.0, 0.0, -1.0},

		/* 91-100 */
		{-2, 0, 0, 4, 0, 1282.0, 0.0, -3.0, -23.0, 0.0, 1.0},
		{-1, -1, 0, 2, 1, 742.0, 0.0, 1.0, -391.0, 0.0, 0.0},
		{-1, 0, 0, 1, 1, 1020.0, 0.0, -25.0, -495.0, 0.0, -10.0},
		{0, 1, 0, 0, 2, 715.0, 0.0, -4.0, -326.0, 0.0, 2.0},
		{0, 0, -2, 0, 1, -666.0, 0.0, -3.0, 369.0, 0.0, -1.0},
		{0, -1, 2, 0, 1, -667.0, 0.0, 1.0, 346.0, 0.0, 1.0},
		{0, 0, 2, -1, 2, -704.0, 0.0, 0.0, 304.0, 0.0, 0.0},
		{0, 0, 2, 4, 2, -694.0, 0.0, 5.0, 294.0, 0.0, 2.0},
		{-2, -1, 0, 2, 0, -1014.0, 0.0, -1.0, 4.0, 0.0, -1.0},
		{1, 1, 0, -2, 1, -585.0, 0.0, -2.0, 316.0, 0.0, -1.0},

		/* 101-110 */
		{-1, 1, 0, 2, 0, -949.0, 0.0, 1.0, 8.0, 0.0, -1.0},
		{-1, 1, 0, 1, 2, -595.0, 0.0, 0.0, 258.0, 0.0, 0.0},
		{1, -1, 0, 0, 1, 528.0, 0.0, 0.0, -279.0, 0.0, 0.0},
		{1, -1, 2, 2, 2, -590.0, 0.0, 4.0, 252.0, 0.0, 2.0},
		{-1, 1, 2, 2, 2, 570.0, 0.0, -2.0, -244.0, 0.0, -1.0},
		{3, 0, 2, 0, 1, -502.0, 0.0, 3.0, 250.0, 0.0, 2.0},
		{0, 1, -2, 2, 0, -875.0, 0.0, 1.0, 29.0, 0.0, 0.0},
		{-1, 0, 0, -2, 1, -492.0, 0.0, -3.0, 275.0, 0.0, -1.0},
		{0, 1, 2, 2, 2, 535.0, 0.0, -2.0, -228.0, 0.0, -1.0},
		{-1, -1, 2, 2, 1, -467.0, 0.0, 1.0, 240.0, 0.0, 1.0},

		/* 111-120 */
		{0, -1, 0, 0, 2, 591.0, 0.0, 0.0, -253.0, 0.0, 0.0},
		{1, 0, 2, -4, 1, -453.0, 0.0, -1.0, 244.0, 0.0, -1.0},
		{-1, 0, -2, 2, 0, 766.0, 0.0, 1.0, 9.0, 0.0, 0.0},
		{0, -1, 2, 2, 1, -446.0, 0.0, 2.0, 225.0, 0.0, 1.0},
		{2, -1, 2, 0, 2, -488.0, 0.0, 2.0, 207.0, 0.0, 1.0},
		{0, 0, 0, 2, 2, -468.0, 0.0, 0.0, 201.0, 0.0, 0.0},
		{1, -1, 2, 0, 1, -421.0, 0.0, 1.0, 216.0, 0.0, 1.0},
		{-1, 1, 2, 0, 2, 463.0, 0.0, 0.0, -200.0, 0.0, 0.0},
		{0, 1, 0, 2, 0, -673.0, 0.0, 2.0, 14.0, 0.0, 0.0},
		{0, -1, -2, 2, 0, 658.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 121-130 */
		{0, 3, 2, -2, 2, -438.0, 0.0, 0.0, 188.0, 0.0, 0.0},
		{0, 0, 0, 1, 1, -390.0, 0.0, 0.0, 205.0, 0.0, 0.0},
		{-1, 0, 2, 2, 0, 639.0, -11.0, -2.0, -19.0, 0.0, 0.0},
		{2, 1, 2, 0, 2, 412.0, 0.0, -2.0, -176.0, 0.0, -1.0},
		{1, 1, 0, 0, 1, -361.0, 0.0, 0.0, 189.0, 0.0, 0.0},
		{1, 1, 2, 0, 1, 360.0, 0.0, -1.0, -185.0, 0.0, -1.0},
		{2, 0, 0, 2, 0, 588.0, 0.0, -3.0, -24.0, 0.0, 0.0},
		{1, 0, -2, 2, 0, -578.0, 0.0, 1.0, 5.0, 0.0, 0.0},
		{-1, 0, 0, 2, 2, -396.0, 0.0, 0.0, 171.0, 0.0, 0.0},
		{0, 1, 0, 1, 0, 565.0, 0.0, -1.0, -6.0, 0.0, 0.0},

		/* 131-140 */
		{0, 1, 0, -2, 1, -335.0, 0.0, -1.0, 184.0, 0.0, -1.0},
		{-1, 0, 2, -2, 2, 357.0, 0.0, 1.0, -154.0, 0.0, 0.0},
		{0, 0, 0, -1, 1, 321.0, 0.0, 1.0, -174.0, 0.0, 0.0},
		{-1, 1, 0, 0, 1, -301.0, 0.0, -1.0, 162.0, 0.0, 0.0},
		{1, 0, 2, -1, 2, -334.0, 0.0, 0.0, 144.0, 0.0, 0.0},
		{1, -1, 0, 2, 0, 493.0, 0.0, -2.0, -15.0, 0.0, 0.0},
		{0, 0, 0, 4, 0, 494.0, 0.0, -2.0, -19.0, 0.0, 0.0},
		{1, 0, 2, 1, 2, 337.0, 0.0, -1.0, -143.0, 0.0, -1.0},
		{0, 0, 2, 1, 1, 280.0, 0.0, -1.0, -144.0, 0.0, 0.0},
		{1, 0, 0, -2, 2, 309.0, 0.0, 1.0, -134.0, 0.0, 0.0},

		/* 141-150 */
		{-1, 0, 2, 4, 1, -263.0, 0.0, 2.0, 131.0, 0.0, 1.0},
		{1, 0, -2, 0, 1, 253.0, 0.0, 1.0, -138.0, 0.0, 0.0},
		{1, 1, 2, -2, 1, 245.0, 0.0, 0.0, -128.0, 0.0, 0.0},
		{0, 0, 2, 2, 0, 416.0, 0.0, -2.0, -17.0, 0.0, 0.0},
		{-1, 0, 2, -1, 1, -229.0, 0.0, 0.0, 128.0, 0.0, 0.0},
		{-2, 0, 2, 2, 1, 231.0, 0.0, 0.0, -120.0, 0.0, 0.0},
		{4, 0, 2, 0, 2, -259.0, 0.0, 2.0, 109.0, 0.0, 1.0},
		{2, -1, 0, 0, 0, 375.0, 0.0, -1.0, -8.0, 0.0, 0.0},
		{2, 1, 2, -2, 2, 252.0, 0.0, 0.0, -108.0, 0.0, 0.0},
		{0, 1, 2, 1, 2, -245.0, 0.0, 1.0, 104.0, 0.0, 0.0},

		/* 151-160 */
		{1, 0, 4, -2, 2, 243.0, 0.0, -1.0, -104.0, 0.0, 0.0},
		{-1, -1, 0, 0, 1, 208.0, 0.0, 1.0, -112.0, 0.0, 0.0},
		{0, 1, 0, 2, 1, 199.0, 0.0, 0.0, -102.0, 0.0, 0.0},
		{-2, 0, 2, 4, 1, -208.0, 0.0, 1.0, 105.0, 0.0, 0.0},
		{2, 0, 2, 0, 0, 335.0, 0.0, -2.0, -14.0, 0.0, 0.0},
		{1, 0, 0, 1, 0, -325.0, 0.0, 1.0, 7.0, 0.0, 0.0},
		{-1, 0, 0, 4, 1, -187.0, 0.0, 0.0, 96.0, 0.0, 0.0},
		{-1, 0, 4, 0, 1, 197.0, 0.0, -1.0, -100.0, 0.0, 0.0},
		{2, 0, 2, 2, 1, -192.0, 0.0, 2.0, 94.0, 0.0, 1.0},
		{0, 0, 2, -3, 2, -188.0, 0.0, 0.0, 83.0, 0.0, 0.0},

		/* 161-170 */
		{-1, -2, 0, 2, 0, 276.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 1, 0, 0, 0, -286.0, 0.0, 1.0, 6.0, 0.0, 0.0},
		{0, 0, 4, 0, 2, 186.0, 0.0, -1.0, -79.0, 0.0, 0.0},
		{0, 0, 0, 0, 3, -219.0, 0.0, 0.0, 43.0, 0.0, 0.0},
		{0, 3, 0, 0, 0, 276.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 2, -4, 1, -153.0, 0.0, -1.0, 84.0, 0.0, 0.0},
		{0, -1, 0, 2, 1, -156.0, 0.0, 0.0, 81.0, 0.0, 0.0},
		{0, 0, 0, 4, 1, -154.0, 0.0, 1.0, 78.0, 0.0, 0.0},
		{-1, -1, 2, 4, 2, -174.0, 0.0, 1.0, 75.0, 0.0, 0.0},
		{1, 0, 2, 4, 2, -163.0, 0.0, 2.0, 69.0, 0.0, 1.0},

		/* 171-180 */
		{-2, 2, 0, 2, 0, -228.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-2, -1, 2, 0, 1, 91.0, 0.0, -4.0, -54.0, 0.0, -2.0},
		{-2, 0, 0, 2, 2, 175.0, 0.0, 0.0, -75.0, 0.0, 0.0},
		{-1, -1, 2, 0, 2, -159.0, 0.0, 0.0, 69.0, 0.0, 0.0},
		{0, 0, 4, -2, 1, 141.0, 0.0, 0.0, -72.0, 0.0, 0.0},
		{3, 0, 2, -2, 1, 147.0, 0.0, 0.0, -75.0, 0.0, 0.0},
		{-2, -1, 0, 2, 1, -132.0, 0.0, 0.0, 69.0, 0.0, 0.0},
		{1, 0, 0, -1, 1, 159.0, 0.0, -28.0, -54.0, 0.0, 11.0},
		{0, -2, 0, 2, 0, 213.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{-2, 0, 0, 4, 1, 123.0, 0.0, 0.0, -64.0, 0.0, 0.0},

		/* 181-190 */
		{-3, 0, 0, 0, 1, -118.0, 0.0, -1.0, 66.0, 0.0, 0.0},
		{1, 1, 2, 2, 2, 144.0, 0.0, -1.0, -61.0, 0.0, 0.0},
		{0, 0, 2, 4, 1, -121.0, 0.0, 1.0, 60.0, 0.0, 0.0},
		{3, 0, 2, 2, 2, -134.0, 0.0, 1.0, 56.0, 0.0, 1.0},
		{-1, 1, 2, -2, 1, -105.0, 0.0, 0.0, 57.0, 0.0, 0.0},
		{2, 0, 0, -4, 1, -102.0, 0.0, 0.0, 56.0, 0.0, 0.0},
		{0, 0, 0, -2, 2, 120.0, 0.0, 0.0, -52.0, 0.0, 0.0},
		{2, 0, 2, -4, 1, 101.0, 0.0, 0.0, -54.0, 0.0, 0.0},
		{-1, 1, 0, 2, 1, -113.0, 0.0, 0.0, 59.0, 0.0, 0.0},
		{0, 0, 2, -1, 1, -106.0, 0.0, 0.0, 61.0, 0.0, 0.0},

		/* 191-200 */
		{0, -2, 2, 2, 2, -129.0, 0.0, 1.0, 55.0, 0.0, 0.0},
		{2, 0, 0, 2, 1, -114.0, 0.0, 0.0, 57.0, 0.0, 0.0},
		{4, 0, 2, -2, 2, 113.0, 0.0, -1.0, -49.0, 0.0, 0.0},
		{2, 0, 0, -2, 2, -102.0, 0.0, 0.0, 44.0, 0.0, 0.0},
		{0, 2, 0, 0, 1, -94.0, 0.0, 0.0, 51.0, 0.0, 0.0},
		{1, 0, 0, -4, 1, -100.0, 0.0, -1.0, 56.0, 0.0, 0.0},
		{0, 2, 2, -2, 1, 87.0, 0.0, 0.0, -47.0, 0.0, 0.0},
		{-3, 0, 0, 4, 0, 161.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-1, 1, 2, 0, 1, 96.0, 0.0, 0.0, -50.0, 0.0, 0.0},
		{-1, -1, 0, 4, 0, 151.0, 0.0, -1.0, -5.0, 0.0, 0.0},

		/* 201-210 */
		{-1, -2, 2, 2, 2, -104.0, 0.0, 0.0, 44.0, 0.0, 0.0},
		{-2, -1, 2, 4, 2, -110.0, 0.0, 0.0, 48.0, 0.0, 0.0},
		{1, -1, 2, 2, 1, -100.0, 0.0, 1.0, 50.0, 0.0, 0.0},
		{-2, 1, 0, 2, 0, 92.0, 0.0, -5.0, 12.0, 0.0, -2.0},
		{-2, 1, 2, 0, 1, 82.0, 0.0, 0.0, -45.0, 0.0, 0.0},
		{2, 1, 0, -2, 1, 82.0, 0.0, 0.0, -45.0, 0.0, 0.0},
		{-3, 0, 2, 0, 1, -78.0, 0.0, 0.0, 41.0, 0.0, 0.0},
		{-2, 0, 2, -2, 1, -77.0, 0.0, 0.0, 43.0, 0.0, 0.0},
		{-1, 1, 0, 2, 2, 2.0, 0.0, 0.0, 54.0, 0.0, 0.0},
		{0, -1, 2, -1, 2, 94.0, 0.0, 0.0, -40.0, 0.0, 0.0},

		/* 211-220 */
		{-1, 0, 4, -2, 2, -93.0, 0.0, 0.0, 40.0, 0.0, 0.0},
		{0, -2, 2, 0, 2, -83.0, 0.0, 10.0, 40.0, 0.0, -2.0},
		{-1, 0, 2, 1, 2, 83.0, 0.0, 0.0, -36.0, 0.0, 0.0},
		{2, 0, 0, 0, 2, -91.0, 0.0, 0.0, 39.0, 0.0, 0.0},
		{0, 0, 2, 0, 3, 128.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-2, 0, 4, 0, 2, -79.0, 0.0, 0.0, 34.0, 0.0, 0.0},
		{-1, 0, -2, 0, 1, -83.0, 0.0, 0.0, 47.0, 0.0, 0.0},
		{-1, 1, 2, 2, 1, 84.0, 0.0, 0.0, -44.0, 0.0, 0.0},
		{3, 0, 0, 0, 1, 83.0, 0.0, 0.0, -43.0, 0.0, 0.0},
		{-1, 0, 2, 3, 2, 91.0, 0.0, 0.0, -39.0, 0.0, 0.0},

		/* 221-230 */
		{2, -1, 2, 0, 1, -77.0, 0.0, 0.0, 39.0, 0.0, 0.0},
		{0, 1, 2, 2, 1, 84.0, 0.0, 0.0, -43.0, 0.0, 0.0},
		{0, -1, 2, 4, 2, -92.0, 0.0, 1.0, 39.0, 0.0, 0.0},
		{2, -1, 2, 2, 2, -92.0, 0.0, 1.0, 39.0, 0.0, 0.0},
		{0, 2, -2, 2, 0, -94.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 2, -1, 1, 68.0, 0.0, 0.0, -36.0, 0.0, 0.0},
		{0, -2, 0, 0, 1, -61.0, 0.0, 0.0, 32.0, 0.0, 0.0},
		{1, 0, 2, -4, 2, 71.0, 0.0, 0.0, -31.0, 0.0, 0.0},
		{1, -1, 0, -2, 1, 62.0, 0.0, 0.0, -34.0, 0.0, 0.0},
		{-1, -1, 2, 0, 1, -63.0, 0.0, 0.0, 33.0, 0.0, 0.0},

		/* 231-240 */
		{1, -1, 2, -2, 2, -73.0, 0.0, 0.0, 32.0, 0.0, 0.0},
		{-2, -1, 0, 4, 0, 115.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 0, 0, 3, 0, -103.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, -1, 2, 2, 2, 63.0, 0.0, 0.0, -28.0, 0.0, 0.0},
		{0, 2, 2, 0, 2, 74.0, 0.0, 0.0, -32.0, 0.0, 0.0},
		{1, 1, 0, 2, 0, -103.0, 0.0, -3.0, 3.0, 0.0, -1.0},
		{2, 0, 2, -1, 2, -69.0, 0.0, 0.0, 30.0, 0.0, 0.0},
		{1, 0, 2, 1, 1, 57.0, 0.0, 0.0, -29.0, 0.0, 0.0},
		{4, 0, 0, 0, 0, 94.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{2, 1, 2, 0, 1, 64.0, 0.0, 0.0, -33.0, 0.0, 0.0},

		/* 241-250 */
		{3, -1, 2, 0, 2, -63.0, 0.0, 0.0, 26.0, 0.0, 0.0},
		{-2, 2, 0, 2, 1, -38.0, 0.0, 0.0, 20.0, 0.0, 0.0},
		{1, 0, 2, -3, 1, -43.0, 0.0, 0.0, 24.0, 0.0, 0.0},
		{1, 1, 2, -4, 1, -45.0, 0.0, 0.0, 23.0, 0.0, 0.0},
		{-1, -1, 2, -2, 1, 47.0, 0.0, 0.0, -24.0, 0.0, 0.0},
		{0, -1, 0, -1, 1, -48.0, 0.0, 0.0, 25.0, 0.0, 0.0},
		{0, -1, 0, -2, 1, 45.0, 0.0, 0.0, -26.0, 0.0, 0.0},
		{-2, 0, 0, 0, 2, 56.0, 0.0, 0.0, -25.0, 0.0, 0.0},
		{-2, 0, -2, 2, 0, 88.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 0, -2, 4, 0, -75.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 251-260 */
		{1, -2, 0, 0, 0, 85.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 0, 1, 1, 49.0, 0.0, 0.0, -26.0, 0.0, 0.0},
		{-1, 2, 0, 2, 0, -74.0, 0.0, -3.0, -1.0, 0.0, -1.0},
		{1, -1, 2, -2, 1, -39.0, 0.0, 0.0, 21.0, 0.0, 0.0},
		{1, 2, 2, -2, 2, 45.0, 0.0, 0.0, -20.0, 0.0, 0.0},
		{2, -1, 2, -2, 2, 51.0, 0.0, 0.0, -22.0, 0.0, 0.0},
		{1, 0, 2, -1, 1, -40.0, 0.0, 0.0, 21.0, 0.0, 0.0},
		{2, 1, 2, -2, 1, 41.0, 0.0, 0.0, -21.0, 0.0, 0.0},
		{-2, 0, 0, -2, 1, -42.0, 0.0, 0.0, 24.0, 0.0, 0.0},
		{1, -2, 2, 0, 2, -51.0, 0.0, 0.0, 22.0, 0.0, 0.0},

		/* 261-270 */
		{0, 1, 2, 1, 1, -42.0, 0.0, 0.0, 22.0, 0.0, 0.0},
		{1, 0, 4, -2, 1, 39.0, 0.0, 0.0, -21.0, 0.0, 0.0},
		{-2, 0, 4, 2, 2, 46.0, 0.0, 0.0, -18.0, 0.0, 0.0},
		{1, 1, 2, 1, 2, -53.0, 0.0, 0.0, 22.0, 0.0, 0.0},
		{1, 0, 0, 4, 0, 82.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{1, 0, 2, 2, 0, 81.0, 0.0, -1.0, -4.0, 0.0, 0.0},
		{2, 0, 2, 1, 2, 47.0, 0.0, 0.0, -19.0, 0.0, 0.0},
		{3, 1, 2, 0, 2, 53.0, 0.0, 0.0, -23.0, 0.0, 0.0},
		{4, 0, 2, 0, 1, -45.0, 0.0, 0.0, 22.0, 0.0, 0.0},
		{-2, -1, 2, 0, 0, -44.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 271-280 */
		{0, 1, -2, 2, 1, -33.0, 0.0, 0.0, 16.0, 0.0, 0.0},
		{1, 0, -2, 1, 0, -61.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{0, -1, -2, 2, 1, 28.0, 0.0, 0.0, -15.0, 0.0, 0.0},
		{2, -1, 0, -2, 1, -38.0, 0.0, 0.0, 19.0, 0.0, 0.0},
		{-1, 0, 2, -1, 2, -33.0, 0.0, 0.0, 21.0, 0.0, 0.0},
		{1, 0, 2, -3, 2, -60.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 2, -2, 3, 48.0, 0.0, 0.0, -10.0, 0.0, 0.0},
		{0, 0, 2, -3, 1, 27.0, 0.0, 0.0, -14.0, 0.0, 0.0},
		{-1, 0, -2, 2, 1, 38.0, 0.0, 0.0, -20.0, 0.0, 0.0},
		{0, 0, 2, -4, 2, 31.0, 0.0, 0.0, -13.0, 0.0, 0.0},

		/* 281-290 */
		{-2, 1, 0, 0, 1, -29.0, 0.0, 0.0, 15.0, 0.0, 0.0},
		{-1, 0, 0, -1, 1, 28.0, 0.0, 0.0, -15.0, 0.0, 0.0},
		{2, 0, 2, -4, 2, -32.0, 0.0, 0.0, 15.0, 0.0, 0.0},
		{0, 0, 4, -4, 4, 45.0, 0.0, 0.0, -8.0, 0.0, 0.0},
		{0, 0, 4, -4, 2, -44.0, 0.0, 0.0, 19.0, 0.0, 0.0},
		{-1, -2, 0, 2, 1, 28.0, 0.0, 0.0, -15.0, 0.0, 0.0},
		{-2, 0, 0, 3, 0, -51.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, -2, 2, 1, -36.0, 0.0, 0.0, 20.0, 0.0, 0.0},
		{-3, 0, 2, 2, 2, 44.0, 0.0, 0.0, -19.0, 0.0, 0.0},
		{-3, 0, 2, 2, 1, 26.0, 0.0, 0.0, -14.0, 0.0, 0.0},

		/* 291-300 */
		{-2, 0, 2, 2, 0, -60.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, -1, 0, 0, 1, 35.0, 0.0, 0.0, -18.0, 0.0, 0.0},
		{-2, 1, 2, 2, 2, -27.0, 0.0, 0.0, 11.0, 0.0, 0.0},
		{1, 1, 0, 1, 0, 47.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, 1, 4, -2, 2, 36.0, 0.0, 0.0, -15.0, 0.0, 0.0},
		{-1, 1, 0, -2, 1, -36.0, 0.0, 0.0, 20.0, 0.0, 0.0},
		{0, 0, 0, -4, 1, -35.0, 0.0, 0.0, 19.0, 0.0, 0.0},
		{1, -1, 0, 2, 1, -37.0, 0.0, 0.0, 19.0, 0.0, 0.0},
		{1, 1, 0, 2, 1, 32.0, 0.0, 0.0, -16.0, 0.0, 0.0},
		{-1, 2, 2, 2, 2, 35.0, 0.0, 0.0, -14.0, 0.0, 0.0},

		/* 301-310 */
		{3, 1, 2, -2, 2, 32.0, 0.0, 0.0, -13.0, 0.0, 0.0},
		{0, -1, 0, 4, 0, 65.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, -1, 0, 2, 0, 47.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, 0, 4, 0, 1, 32.0, 0.0, 0.0, -16.0, 0.0, 0.0},
		{2, 0, 4, -2, 2, 37.0, 0.0, 0.0, -16.0, 0.0, 0.0},
		{-1, -1, 2, 4, 1, -30.0, 0.0, 0.0, 15.0, 0.0, 0.0},
		{1, 0, 0, 4, 1, -32.0, 0.0, 0.0, 16.0, 0.0, 0.0},
		{1, -2, 2, 2, 2, -31.0, 0.0, 0.0, 13.0, 0.0, 0.0},
		{0, 0, 2, 3, 2, 37.0, 0.0, 0.0, -16.0, 0.0, 0.0},
		{-1, 1, 2, 4, 2, 31.0, 0.0, 0.0, -13.0, 0.0, 0.0},

		/* 311-320 */
		{3, 0, 0, 2, 0, 49.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 0, 4, 2, 2, 32.0, 0.0, 0.0, -13.0, 0.0, 0.0},
		{1, 1, 2, 2, 1, 23.0, 0.0, 0.0, -12.0, 0.0, 0.0},
		{-2, 0, 2, 6, 2, -43.0, 0.0, 0.0, 18.0, 0.0, 0.0},
		{2, 1, 2, 2, 2, 26.0, 0.0, 0.0, -11.0, 0.0, 0.0},
		{-1, 0, 2, 6, 2, -32.0, 0.0, 0.0, 14.0, 0.0, 0.0},
		{1, 0, 2, 4, 1, -29.0, 0.0, 0.0, 14.0, 0.0, 0.0},
		{2, 0, 2, 4, 2, -27.0, 0.0, 0.0, 12.0, 0.0, 0.0},
		{1, 1, -2, 1, 0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, 1, 2, 1, 2, -11.0, 0.0, 0.0, 5.0, 0.0, 0.0},

		/* 321-330 */
		{2, 0, -2, 0, 2, -21.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{-1, 0, 0, 1, 2, -34.0, 0.0, 0.0, 15.0, 0.0, 0.0},
		{-4, 0, 2, 2, 1, -10.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{-1, -1, 0, 1, 0, -36.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, -2, 2, 2, -9.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{1, 0, 0, -1, 2, -12.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{0, -1, 2, -2, 3, -21.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{-2, 1, 2, 0, 0, -29.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, 0, 2, -2, 4, -15.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-2, -2, 0, 2, 0, -20.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 331-340 */
		{-2, 0, -2, 4, 0, 28.0, 0.0, 0.0, 0.0, 0.0, -2.0},
		{0, -2, -2, 2, 0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 2, 0, -2, 1, -22.0, 0.0, 0.0, 12.0, 0.0, 0.0},
		{3, 0, 0, -4, 1, -14.0, 0.0, 0.0, 7.0, 0.0, 0.0},
		{-1, 1, 2, -2, 2, 24.0, 0.0, 0.0, -11.0, 0.0, 0.0},
		{1, -1, 2, -4, 1, 11.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{1, 1, 0, -2, 2, 14.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{-3, 0, 2, 0, 0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, 0, 2, 0, 2, 18.0, 0.0, 0.0, -8.0, 0.0, 0.0},
		{-2, 0, 0, 1, 0, -38.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 341-350 */
		{0, 0, -2, 1, 0, -31.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, 0, 0, 2, 1, -16.0, 0.0, 0.0, 8.0, 0.0, 0.0},
		{-1, -1, -2, 2, 0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 2, -4, 1, -18.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{2, 1, 0, -4, 1, -10.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{0, 2, 0, -2, 1, -17.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{1, 0, 0, -3, 1, 9.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{-2, 0, 2, -2, 2, 16.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{-2, -1, 0, 0, 1, 22.0, 0.0, 0.0, -12.0, 0.0, 0.0},
		{-4, 0, 0, 2, 0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 351-360 */
		{1, 1, 0, -4, 1, -13.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{-1, 0, 2, -4, 1, -17.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{0, 0, 4, -4, 1, -14.0, 0.0, 0.0, 8.0, 0.0, 0.0},
		{0, 3, 2, -2, 2, 0.0, 0.0, 0.0, -7.0, 0.0, 0.0},
		{-3, -1, 0, 4, 0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, 0, 0, 4, 1, 19.0, 0.0, 0.0, -10.0, 0.0, 0.0},
		{1, -1, -2, 2, 0, -34.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 0, 2, 2, -20.0, 0.0, 0.0, 8.0, 0.0, 0.0},
		{1, -2, 0, 0, 1, 9.0, 0.0, 0.0, -5.0, 0.0, 0.0},
		{1, -1, 0, 0, 2, -18.0, 0.0, 0.0, 7.0, 0.0, 0.0},

		/* 361-370 */
		{0, 0, 0, 1, 2, 13.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{-1, -1, 2, 0, 0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, -2, 2, -2, 2, -12.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{0, -1, 2, -1, 1, 15.0, 0.0, 0.0, -8.0, 0.0, 0.0},
		{-1, 0, 2, 0, 3, -11.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{1, 1, 0, 0, 2, 13.0, 0.0, 0.0, -5.0, 0.0, 0.0},
		{-1, 1, 2, 0, 0, -18.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 2, 0, 0, 0, -35.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 2, 2, 0, 2, 9.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{-1, 0, 4, -2, 1, -19.0, 0.0, 0.0, 10.0, 0.0, 0.0},

		/* 371-380 */
		{3, 0, 2, -4, 2, -26.0, 0.0, 0.0, 11.0, 0.0, 0.0},
		{1, 2, 2, -2, 1, 8.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{1, 0, 4, -4, 2, -10.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{-2, -1, 0, 4, 1, 10.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{0, -1, 0, 2, 2, -21.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{-2, 1, 0, 4, 0, -15.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, -1, 2, 2, 1, 9.0, 0.0, 0.0, -5.0, 0.0, 0.0},
		{2, 0, -2, 2, 0, -29.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 0, 1, 1, -19.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{0, 1, 0, 2, 2, 12.0, 0.0, 0.0, -5.0, 0.0, 0.0},

		/* 381-390 */
		{1, -1, 2, -1, 2, 22.0, 0.0, 0.0, -9.0, 0.0, 0.0},
		{-2, 0, 4, 0, 1, -10.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{2, 1, 0, 0, 1, -20.0, 0.0, 0.0, 11.0, 0.0, 0.0},
		{0, 1, 2, 0, 0, -20.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, -1, 4, -2, 2, -17.0, 0.0, 0.0, 7.0, 0.0, 0.0},
		{0, 0, 4, -2, 4, 15.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{0, 2, 2, 0, 1, 8.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{-3, 0, 0, 6, 0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 0, 4, 1, -12.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{1, -2, 0, 2, 0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 391-400 */
		{-1, 0, 0, 4, 2, -13.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{-1, -2, 2, 2, 1, -14.0, 0.0, 0.0, 8.0, 0.0, 0.0},
		{-1, 0, 0, -2, 2, 13.0, 0.0, 0.0, -5.0, 0.0, 0.0},
		{1, 0, -2, -2, 1, -17.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{0, 0, -2, -2, 1, -12.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{-2, 0, -2, 0, 1, -10.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{0, 0, 0, 3, 1, 10.0, 0.0, 0.0, -6.0, 0.0, 0.0},
		{0, 0, 0, 3, 0, -15.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 1, 0, 4, 0, -22.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 2, 2, 0, 28.0, 0.0, 0.0, -1.0, 0.0, 0.0},

		/* 401-410 */
		{-2, 0, 2, 3, 2, 15.0, 0.0, 0.0, -7.0, 0.0, 0.0},
		{1, 0, 0, 2, 2, 23.0, 0.0, 0.0, -10.0, 0.0, 0.0},
		{0, -1, 2, 1, 2, 12.0, 0.0, 0.0, -5.0, 0.0, 0.0},
		{3, -1, 0, 0, 0, 29.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{2, 0, 0, 1, 0, -25.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{1, -1, 2, 0, 0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, 2, 1, 0, -18.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 2, 0, 3, 15.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{3, 1, 0, 0, 0, -23.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{3, -1, 2, -2, 2, 12.0, 0.0, 0.0, -5.0, 0.0, 0.0},

		/* 411-420 */
		{2, 0, 2, -1, 1, -8.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{1, 1, 2, 0, 0, -19.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, 4, -1, 2, -10.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{1, 2, 2, 0, 2, 21.0, 0.0, 0.0, -9.0, 0.0, 0.0},
		{-2, 0, 0, 6, 0, 23.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, -1, 0, 4, 1, -16.0, 0.0, 0.0, 8.0, 0.0, 0.0},
		{-2, -1, 2, 4, 1, -19.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{0, -2, 2, 2, 1, -22.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{0, -1, 2, 2, 0, 27.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-1, 0, 2, 3, 1, 16.0, 0.0, 0.0, -8.0, 0.0, 0.0},

		/* 421-430 */
		{-2, 1, 2, 4, 2, 19.0, 0.0, 0.0, -8.0, 0.0, 0.0},
		{2, 0, 0, 2, 2, 9.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{2, -2, 2, 0, 2, -9.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{-1, 1, 2, 3, 2, -9.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{3, 0, 2, -1, 2, -8.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{4, 0, 2, -2, 1, 18.0, 0.0, 0.0, -9.0, 0.0, 0.0},
		{-1, 0, 0, 6, 0, 16.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-1, -2, 2, 4, 2, -10.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{-3, 0, 2, 6, 2, -23.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{-1, 0, 2, 4, 0, 16.0, 0.0, 0.0, -1.0, 0.0, 0.0},

		/* 431-440 */
		{3, 0, 0, 2, 1, -12.0, 0.0, 0.0, 6.0, 0.0, 0.0},
		{3, -1, 2, 0, 1, -8.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{3, 0, 2, 0, 0, 30.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{1, 0, 4, 0, 2, 24.0, 0.0, 0.0, -10.0, 0.0, 0.0},
		{5, 0, 2, -2, 2, 10.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{0, -1, 2, 4, 1, -16.0, 0.0, 0.0, 7.0, 0.0, 0.0},
		{2, -1, 2, 2, 1, -16.0, 0.0, 0.0, 7.0, 0.0, 0.0},
		{0, 1, 2, 4, 2, 17.0, 0.0, 0.0, -7.0, 0.0, 0.0},
		{1, -1, 2, 4, 2, -24.0, 0.0, 0.0, 10.0, 0.0, 0.0},
		{3, -1, 2, 2, 2, -12.0, 0.0, 0.0, 5.0, 0.0, 0.0},

		/* 441-450 */
		{3, 0, 2, 2, 1, -24.0, 0.0, 0.0, 11.0, 0.0, 0.0},
		{5, 0, 2, 0, 2, -23.0, 0.0, 0.0, 9.0, 0.0, 0.0},
		{0, 0, 2, 6, 2, -13.0, 0.0, 0.0, 5.0, 0.0, 0.0},
		{4, 0, 2, 2, 2, -15.0, 0.0, 0.0, 7.0, 0.0, 0.0},
		{0, -1, 1, -1, 1, 0.0, 0.0, -1988.0, 0.0, 0.0, -1679.0},
		{-1, 0, 1, 0, 3, 0.0, 0.0, -63.0, 0.0, 0.0, -27.0},
		{0, -2, 2, -2, 3, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, -1, 0, 1, 0.0, 0.0, 5.0, 0.0, 0.0, 4.0},
		{2, -2, 0, -2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-1, 0, 1, 0, 2, 0.0, 0.0, 364.0, 0.0, 0.0, 176.0},

		/* 451-460 */
		{-1, 0, 1, 0, 1, 0.0, 0.0, -1044.0, 0.0, 0.0, -891.0},
		{-1, -1, 2, -1, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-2, 2, 0, 2, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 0, 1, 0, 0, 0.0, 0.0, 330.0, 0.0, 0.0, 0.0},
		{-4, 1, 2, 2, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-3, 0, 2, 1, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-2, -1, 2, 0, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{1, 0, -2, 1, 1, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, -1, -2, 0, 1, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-4, 0, 2, 2, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 461-470 */
		{-3, 1, 0, 3, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 0, -1, 2, 0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0},
		{0, -2, 0, 0, 2, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{0, -2, 0, 0, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-3, 0, 0, 3, 0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, -1, 0, 2, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 0, -2, 3, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-4, 0, 0, 4, 0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 1, -2, 0, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{2, -1, 0, -2, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},

		/* 471-480 */
		{0, 0, 1, -1, 0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 2, 0, 1, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 1, 2, 0, 2, -7.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{1, 1, 0, -1, 1, 7.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{1, 0, 1, -2, 1, 0.0, 0.0, -12.0, 0.0, 0.0, -10.0},
		{0, 2, 0, 0, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{1, -1, 2, -3, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 1, 2, -1, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, 0, 4, -2, 2, -7.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-2, 0, 4, -2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},

		/* 481-490 */
		{-2, -2, 0, 2, 1, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-2, 0, -2, 4, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 2, 2, -4, 1, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{1, 1, 2, -4, 2, 7.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-1, 2, 2, -2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, 0, 0, -3, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 2, 0, 0, 1, -5.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{0, 0, 0, -2, 0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 2, -2, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 1, 0, 0, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 491-500 */
		{0, 0, 0, -1, 2, -8.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-2, 1, 0, 1, 0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, -2, 0, -2, 1, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{1, 0, -2, 0, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-3, 1, 0, 2, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 1, -2, 2, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 0, 0, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-3, 0, 0, 2, 0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, -1, 0, 2, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 0, 2, -6, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},

		/* 501-510 */
		{0, 1, 2, -4, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 0, 0, -4, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-2, 1, 2, -2, 1, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, -1, 2, -4, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{0, 1, 0, -2, 2, 9.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-1, 0, 0, -2, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 0, -2, -2, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-4, 0, 2, 0, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, -1, 0, -1, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, -2, 0, 2, 9.0, 0.0, 0.0, -3.0, 0.0, 0.0},

		/* 511-520 */
		{-3, 0, 0, 1, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 0, -2, 1, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 0, -2, 2, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{0, 0, -4, 2, 0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, -1, -2, 2, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 2, -6, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 0, 2, -4, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{1, 0, 0, -4, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{2, 1, 2, -4, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{2, 1, 2, -4, 1, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},

		/* 521-530 */
		{0, 1, 4, -4, 4, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 4, -4, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-1, -1, -2, 4, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -3, 0, 2, 0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 0, -2, 4, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, -1, 0, 3, 0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, -2, 3, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 0, 0, 3, 1, -5.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{0, -1, 0, 1, 0, -13.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, 0, 2, 2, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 531-540 */
		{1, 1, -2, 2, 0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 1, 0, 2, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{1, -2, 2, -2, 1, 10.0, 0.0, 13.0, 6.0, 0.0, -5.0},
		{0, 0, 1, 0, 2, 0.0, 0.0, 30.0, 0.0, 0.0, 14.0},
		{0, 0, 1, 0, 1, 0.0, 0.0, -162.0, 0.0, 0.0, -138.0},
		{0, 0, 1, 0, 0, 0.0, 0.0, 75.0, 0.0, 0.0, 0.0},
		{-1, 2, 0, 2, 1, -7.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{0, 0, 2, 0, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, 0, 2, 0, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 0, 0, -1, 1, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 541-550 */
		{3, 0, 0, -2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{1, 0, 2, -2, 3, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 2, 0, 0, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, 0, 2, -3, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 1, 4, -2, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, -2, 0, 4, 0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, -3, 0, 2, 0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, -2, 4, 0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 0, 3, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 0, 0, 4, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},

		/* 551-560 */
		{-1, 0, 0, 3, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, -2, 0, 0, 0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, -1, 0, 1, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 0, 0, 2, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, -2, 2, 0, 1, -6.0, 0.0, -3.0, 3.0, 0.0, 1.0},
		{-1, 0, 1, 2, 1, 0.0, 0.0, -3.0, 0.0, 0.0, -2.0},
		{-1, 1, 0, 3, 0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -1, 2, 1, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, -1, 2, 0, 0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 1, 2, 2, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},

		/* 561-570 */
		{2, -2, 2, -2, 2, -1.0, 0.0, 3.0, 3.0, 0.0, -1.0},
		{1, 1, 0, 1, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{1, 0, 1, 0, 1, 0.0, 0.0, -13.0, 0.0, 0.0, -11.0},
		{1, 0, 1, 0, 0, 3.0, 0.0, 6.0, 0.0, 0.0, 0.0},
		{0, 2, 0, 2, 0, -7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, -1, 2, -2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{0, -1, 4, -2, 1, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{0, 0, 4, -2, 3, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 4, -2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{4, 0, 2, -4, 2, -7.0, 0.0, 0.0, 3.0, 0.0, 0.0},

		/* 571-580 */
		{2, 2, 2, -2, 2, 8.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{2, 0, 4, -4, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, -2, 0, 4, 0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, -3, 2, 2, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-3, 0, 2, 4, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-3, 0, 2, -2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, -1, 0, -2, 1, 8.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{-3, 0, 0, 0, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-3, 0, -2, 2, 0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 0, -4, 1, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},

		/* 581-590 */
		{-2, 1, 0, -2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-4, 0, 0, 0, 1, -8.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{-1, 0, 0, -4, 1, -7.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-3, 0, 0, -2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 0, 3, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-1, 1, 0, 4, 1, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{1, -2, 2, 0, 1, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{0, 1, 0, 3, 0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1, 0, 2, 2, 3, 6.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0, 0, 2, 2, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 591-600 */
		{-2, 0, 2, 2, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 1, 2, 2, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{3, 0, 0, 0, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, 1, 0, 1, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, -1, 2, -1, 2, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{0, 0, 2, 0, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 3, 0, 3, 0.0, 0.0, -26.0, 0.0, 0.0, -11.0},
		{0, 0, 3, 0, 2, 0.0, 0.0, -10.0, 0.0, 0.0, -5.0},
		{-1, 2, 2, 2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-1, 0, 4, 0, 0, -13.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 601-610 */
		{1, 2, 2, 0, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{3, 1, 2, -2, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{1, 1, 4, -2, 2, 7.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-2, -1, 0, 6, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, -2, 0, 4, 0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-2, 0, 0, 6, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, -2, 2, 4, 2, -6.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, -3, 2, 2, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 0, 4, 2, -7.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-1, -1, 2, 3, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 611-620 */
		{-2, 0, 2, 4, 0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, -1, 0, 2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{1, 0, 0, 3, 0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 0, 4, 1, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{0, 1, 0, 4, 0, -11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, -1, 2, 1, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{0, 0, 2, 2, 3, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 2, 2, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-1, 0, 2, 2, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-2, 0, 4, 2, 1, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},

		/* 621-630 */
		{2, 1, 0, 2, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 1, 0, 2, 0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, -1, 2, 0, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 2, 1, 0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 1, 2, 2, 0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 0, 2, 0, 3, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{3, 0, 2, 0, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{1, 0, 2, 0, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{1, 0, 3, 0, 3, 0.0, 0.0, -5.0, 0.0, 0.0, -2.0},
		{1, 1, 2, 1, 1, -7.0, 0.0, 0.0, 4.0, 0.0, 0.0},

		/* 631-640 */
		{0, 2, 2, 2, 2, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{2, 1, 2, 0, 0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 0, 4, -2, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{4, 1, 2, -2, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{-1, -1, 0, 6, 0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-3, -1, 2, 6, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{-1, 0, 0, 6, 1, -5.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-3, 0, 2, 6, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{1, -1, 0, 4, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{1, -1, 0, 4, 0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 641-650 */
		{-2, 0, 2, 5, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{1, -2, 2, 2, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{3, -1, 0, 2, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, -1, 2, 2, 0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, 2, 3, 1, 5.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{-1, 1, 2, 4, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{0, 1, 2, 3, 2, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-1, 0, 4, 2, 1, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 0, 2, 1, 1, 6.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{5, 0, 0, 0, 0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 651-660 */
		{2, 1, 2, 1, 2, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{1, 0, 4, 0, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{3, 1, 2, 0, 1, 7.0, 0.0, 0.0, -4.0, 0.0, 0.0},
		{3, 0, 4, -2, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{-2, -1, 2, 6, 2, -5.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 0, 6, 0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, -2, 2, 4, 2, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{-2, 0, 2, 6, 1, -6.0, 0.0, 0.0, 3.0, 0.0, 0.0},
		{2, 0, 0, 4, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{2, 0, 0, 4, 0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0},

		/* 661-670 */
		{2, -2, 2, 2, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0, 0, 2, 4, 0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1, 0, 2, 3, 2, 7.0, 0.0, 0.0, -3.0, 0.0, 0.0},
		{4, 0, 0, 2, 0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2, 0, 2, 2, 0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0, 0, 4, 2, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{4, -1, 2, 0, 2, -6.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{3, 0, 2, 1, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{2, 1, 2, 2, 1, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{4, 1, 2, 0, 2, 5.0, 0.0, 0.0, -2.0, 0.0, 0.0},

		/* 671-678 */
		{-1, -1, 2, 6, 2, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{-1, 0, 2, 6, 1, -4.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{1, -1, 2, 4, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{1, 1, 2, 4, 2, 4.0, 0.0, 0.0, -2.0, 0.0, 0.0},
		{3, 1, 2, 2, 2, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{5, 0, 2, 0, 1, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{2, -1, 2, 4, 2, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0},
		{2, 0, 2, 4, 1, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0},
	}

	/* Number of terms in the luni-solar nutation model */
	NLS := len(xls)

	/* ------------------------ */
	/* Planetary nutation model */
	/* ------------------------ */

	/* The units for the sine and cosine coefficients are */
	/* 0.1 microarcsecond                                 */

	xpl := []struct {
		nl, /* coefficients of l, F, D and Omega */
		nf,
		nd,
		nom,
		nme, /* coefficients of planetary longitudes */
		nve,
		nea,
		nma,
		nju,
		nsa,
		nur,
		nne,
		npa int /* coefficient of general precession */
		sp, cp int /* longitude sin, cos coefficients */
		se, ce int /* obliquity sin, cos coefficients */
	}{

		/* 1-10 */
		{0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 0, 1440, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, -8, 16, -4, -5, 0, 0, 2, 56, -117, -42, -40},
		{0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 2, 125, -43, 0, -54},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, 2, 0, 5, 0, 0},
		{0, 0, 0, 0, 0, 0, -4, 8, -1, -5, 0, 0, 2, 3, -7, -3, 0},
		{0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 1, 3, 0, 0, -2},
		{0, 1, -1, 1, 0, 0, 3, -8, 3, 0, 0, 0, 0, -114, 0, 0, 61},
		{-1, 0, 0, 0, 0, 10, -3, 0, 0, 0, 0, 0, 0, -219, 89, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, -2, 6, -3, 0, 2, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, -462, 1604, 0, 0},

		/* 11-20 */
		{0, 1, -1, 1, 0, 0, -5, 8, -3, 0, 0, 0, 0, 99, 0, 0, -53},
		{0, 0, 0, 0, 0, 0, -4, 8, -3, 0, 0, 0, 1, -3, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 4, -8, 1, 5, 0, 0, 2, 0, 6, 2, 0},
		{0, 0, 0, 0, 0, -5, 6, 4, 0, 0, 0, 0, 2, 3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 2, -12, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 1, 14, -218, 117, 8},
		{0, 1, -1, 1, 0, 0, -1, 0, 2, -5, 0, 0, 0, 31, -481, -257, -17},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0, -491, 128, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, -2, 5, 0, 0, 0, -3084, 5123, 2735, 1647},
		{0, 0, 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 1, -1444, 2409, -1286, -771},

		/* 21-30 */
		{0, 0, 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 2, 11, -24, -11, -9},
		{2, -1, -1, 0, 0, 0, 3, -7, 0, 0, 0, 0, 0, 26, -9, 0, 0},
		{1, 0, -2, 0, 0, 19, -21, 3, 0, 0, 0, 0, 0, 103, -60, 0, 0},
		{0, 1, -1, 1, 0, 2, -4, 0, -3, 0, 0, 0, 0, 0, -13, -7, 0},
		{1, 0, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0, -26, -29, -16, 14},
		{0, 1, -1, 1, 0, 0, -1, 0, -4, 10, 0, 0, 0, 9, -27, -14, -5},
		{-2, 0, 2, 1, 0, 0, 2, 0, 0, -5, 0, 0, 0, 12, 0, 0, -6},
		{0, 0, 0, 0, 0, 3, -7, 4, 0, 0, 0, 0, 0, -7, 0, 0, 0},
		{0, -1, 1, 0, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0, 24, 0, 0},
		{-2, 0, 2, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0, 284, 0, 0, -151},

		/* 31-40 */
		{-1, 0, 0, 0, 0, 18, -16, 0, 0, 0, 0, 0, 0, 226, 101, 0, 0},
		{-2, 1, 1, 2, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, -8, -2, 0},
		{-1, 1, -1, 1, 0, 18, -17, 0, 0, 0, 0, 0, 0, 0, -6, -3, 0},
		{-1, 0, 1, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0, 5, 0, 0, -3},
		{0, 0, 0, 0, 0, -8, 13, 0, 0, 0, 0, 0, 2, -41, 175, 76, 17},
		{0, 2, -2, 2, 0, -8, 11, 0, 0, 0, 0, 0, 0, 0, 15, 6, 0},
		{0, 0, 0, 0, 0, -8, 13, 0, 0, 0, 0, 0, 1, 425, 212, -133, 269},
		{0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0, 1200, 598, 319, -641},
		{0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 0, 235, 334, 0, 0},
		{0, 1, -1, 1, 0, 8, -14, 0, 0, 0, 0, 0, 0, 11, -12, -7, -6},

		/* 41-50 */
		{0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 1, 5, -6, 3, 3},
		{-2, 0, 2, 1, 0, 0, 2, 0, -4, 5, 0, 0, 0, -5, 0, 0, 3},
		{-2, 0, 2, 2, 0, 3, -3, 0, 0, 0, 0, 0, 0, 6, 0, 0, -3},
		{-2, 0, 2, 0, 0, 0, 2, 0, -3, 1, 0, 0, 0, 15, 0, 0, 0},
		{0, 0, 0, 1, 0, 3, -5, 0, 2, 0, 0, 0, 0, 13, 0, 0, -7},
		{-2, 0, 2, 0, 0, 0, 2, 0, -4, 3, 0, 0, 0, -6, -9, 0, 0},
		{0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 266, -78, 0, 0},
		{0, 0, 0, 1, 0, 0, -1, 2, 0, 0, 0, 0, 0, -460, -435, -232, 246},
		{0, 1, -1, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0, 15, 7, 0},
		{-1, 1, 0, 1, 0, 3, -5, 0, 0, 0, 0, 0, 0, -3, 0, 0, 2},

		/* 51-60 */
		{-1, 0, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 131, 0, 0},
		{-2, 0, 2, 0, 0, 0, 2, 0, -2, -2, 0, 0, 0, 4, 0, 0, 0},
		{-2, 2, 0, 2, 0, 0, -5, 9, 0, 0, 0, 0, 0, 0, 3, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 2, 0, -17, -19, -10, 9},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -9, -11, 6, -5},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -6, 0, 0, 3},
		{-1, 0, 1, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, -16, 8, 0, 0},
		{0, -1, 1, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0},

		/* 61-70 */
		{0, 1, -1, 2, 0, 0, -1, 0, 0, 2, 0, 0, 0, 11, 24, 11, -5},
		{0, 0, 0, 1, 0, 0, -9, 17, 0, 0, 0, 0, 0, -3, -4, -2, 1},
		{0, 0, 0, 2, 0, -3, 5, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 1, -1, 1, 0, 0, -1, 0, -1, 2, 0, 0, 0, 0, -8, -4, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 3, 0, 0},
		{1, 0, -2, 0, 0, 17, -16, 0, -2, 0, 0, 0, 0, 0, 5, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, 1, -3, 0, 0, 0, 0, 3, 2, 0},
		{-2, 0, 2, 1, 0, 0, 5, -6, 0, 0, 0, 0, 0, -6, 4, 2, 3},
		{0, -2, 2, 0, 0, 0, 9, -13, 0, 0, 0, 0, 0, -3, -5, 0, 0},
		{0, 1, -1, 2, 0, 0, -1, 0, 0, 1, 0, 0, 0, -5, 0, 0, 2},

		/* 71-80 */
		{0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 4, 24, 13, -2},
		{0, -1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -42, 20, 0, 0},
		{0, -2, 2, 0, 0, 5, -6, 0, 0, 0, 0, 0, 0, -10, 233, 0, 0},
		{0, -1, 1, 1, 0, 5, -7, 0, 0, 0, 0, 0, 0, -3, 0, 0, 1},
		{-2, 0, 2, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0, 78, -18, 0, 0},
		{2, 1, -3, 1, 0, -6, 7, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -3, -1, 0},
		{0, -1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -4, -2, 1},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 2, 0, 0, 0, -8, -4, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, -5, 3, 0},

		/* 81-90 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, -7, 0, 0, 3},
		{0, 0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 2, -14, 8, 3, 6},
		{0, 0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 1, 0, 8, -4, 0},
		{0, 1, -1, 1, 0, 0, -9, 15, 0, 0, 0, 0, 0, 0, 19, 10, 0},
		{0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0, 45, -22, 0, 0},
		{1, -1, -1, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0, -3, 0, 0, 0},
		{2, 0, -2, 0, 0, 2, -5, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0},
		{-2, 0, 2, 0, 0, 0, 2, 0, -5, 5, 0, 0, 0, 0, 3, 0, 0},
		{2, 0, -2, 1, 0, 0, -6, 8, 0, 0, 0, 0, 0, 3, 5, 3, -2},
		{2, 0, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0, 89, -16, -9, -48},

		/* 91-100 */
		{-2, 1, 1, 0, 0, 0, 1, 0, -3, 0, 0, 0, 0, 0, 3, 0, 0},
		{-2, 1, 1, 1, 0, 0, 1, 0, -3, 0, 0, 0, 0, -3, 7, 4, 2},
		{-2, 0, 2, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0, -349, -62, 0, 0},
		{-2, 0, 2, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, -15, 22, 0, 0},
		{-2, 0, 2, 0, 0, 0, 2, 0, -1, -5, 0, 0, 0, -3, 0, 0, 0},
		{-1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -53, 0, 0, 0},
		{-1, 1, 1, 1, 0, -20, 20, 0, 0, 0, 0, 0, 0, 5, 0, 0, -3},
		{1, 0, -2, 0, 0, 20, -21, 0, 0, 0, 0, 0, 0, 0, -8, 0, 0},
		{0, 0, 0, 1, 0, 0, 8, -15, 0, 0, 0, 0, 0, 15, -7, -4, -8},
		{0, 2, -2, 1, 0, 0, -10, 15, 0, 0, 0, 0, 0, -3, 0, 0, 1},

		/* 101-110 */
		{0, -1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, -21, -78, 0, 0},
		{0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 20, -70, -37, -11},
		{0, 1, -1, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 6, 3, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, -2, 4, 0, 0, 0, 5, 3, 2, -2},
		{2, 0, -2, 1, 0, -6, 8, 0, 0, 0, 0, 0, 0, -17, -4, -2, 9},
		{0, -2, 2, 1, 0, 5, -6, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 32, 15, -8, 17},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 174, 84, 45, -93},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 11, 56, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, -66, -12, -6, 35},

		/* 111-120 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 47, 8, 4, -25},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 8, 4, 0},
		{0, 2, -2, 1, 0, 0, -9, 13, 0, 0, 0, 0, 0, 10, -22, -12, -5},
		{0, 0, 0, 1, 0, 0, 7, -13, 0, 0, 0, 0, 0, -3, 0, 0, 2},
		{-2, 0, 2, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, -24, 12, 0, 0},
		{0, 0, 0, 0, 0, 0, 9, -17, 0, 0, 0, 0, 0, 5, -6, 0, 0},
		{0, 0, 0, 0, 0, 0, -9, 17, 0, 0, 0, 0, 2, 3, 0, 0, -2},
		{1, 0, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0, 4, 3, 1, -2},
		{1, 0, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0, 0, 29, 15, 0},
		{0, 0, 0, 2, 0, 0, -1, 2, 0, 0, 0, 0, 0, -5, -4, -2, 2},

		/* 121-130 */
		{0, -1, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 8, -3, -1, -5},
		{0, -2, 2, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0},
		{0, 0, 0, 0, 0, 3, -5, 0, 2, 0, 0, 0, 0, 10, 0, 0, 0},
		{-2, 0, 2, 1, 0, 0, 2, 0, -3, 1, 0, 0, 0, 3, 0, 0, -2},
		{-2, 0, 2, 1, 0, 3, -3, 0, 0, 0, 0, 0, 0, -5, 0, 0, 3},
		{0, 0, 0, 1, 0, 8, -13, 0, 0, 0, 0, 0, 0, 46, 66, 35, -25},
		{0, -1, 1, 0, 0, 8, -12, 0, 0, 0, 0, 0, 0, -14, 7, 0, 0},
		{0, 2, -2, 1, 0, -8, 11, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0},
		{-1, 0, 1, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, -5, 0, 0, 0},
		{-1, 0, 0, 1, 0, 18, -16, 0, 0, 0, 0, 0, 0, -68, -34, -18, 36},

		/* 131-140 */
		{0, 1, -1, 1, 0, 0, -1, 0, -1, 1, 0, 0, 0, 0, 14, 7, 0},
		{0, 0, 0, 1, 0, 3, -7, 4, 0, 0, 0, 0, 0, 10, -6, -3, -5},
		{-2, 1, 1, 1, 0, 0, -3, 7, 0, 0, 0, 0, 0, -5, -4, -2, 3},
		{0, 1, -1, 2, 0, 0, -1, 0, -2, 5, 0, 0, 0, -3, 5, 2, 1},
		{0, 0, 0, 1, 0, 0, 0, 0, -2, 5, 0, 0, 0, 76, 17, 9, -41},
		{0, 0, 0, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0, 84, 298, 159, -45},
		{1, 0, 0, 1, 0, -10, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 2, -2, 1, 0, 0, -2, 0, 0, 0, 0, 0, 0, -3, 0, 0, 2},
		{-1, 0, 0, 1, 0, 10, -3, 0, 0, 0, 0, 0, 0, -3, 0, 0, 1},
		{0, 0, 0, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0, -82, 292, 156, 44},

		/* 141-150 */
		{0, 0, 0, 1, 0, 0, 0, 0, 2, -5, 0, 0, 0, -73, 17, 9, 39},
		{0, -1, 1, 0, 0, 0, 1, 0, 2, -5, 0, 0, 0, -9, -16, 0, 0},
		{2, -1, -1, 1, 0, 0, 3, -7, 0, 0, 0, 0, 0, 3, 0, -1, -2},
		{-2, 0, 2, 0, 0, 0, 2, 0, 0, -5, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 1, 0, -3, 7, -4, 0, 0, 0, 0, 0, -9, -5, -3, 5},
		{-2, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, -439, 0, 0, 0},
		{1, 0, 0, 1, 0, -18, 16, 0, 0, 0, 0, 0, 0, 57, -28, -15, -30},
		{-2, 1, 1, 1, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, -6, -3, 0},
		{0, 1, -1, 2, 0, -8, 12, 0, 0, 0, 0, 0, 0, -4, 0, 0, 2},
		{0, 0, 0, 1, 0, -8, 13, 0, 0, 0, 0, 0, 0, -40, 57, 30, 21},

		/* 151-160 */
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 1, 23, 7, 3, -13},
		{0, 1, -1, 1, 0, 0, 0, -2, 0, 0, 0, 0, 0, 273, 80, 43, -146},
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -449, 430, 0, 0},
		{0, 1, -1, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0, -8, -47, -25, 4},
		{0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 1, 6, 47, 25, -3},
		{-1, 0, 1, 1, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 23, 13, 0},
		{-1, 0, 1, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0, -3, 0, 0, 2},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, -2, 0, 0, 0, 3, -4, -2, -2},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 2, 0, 0, 0, -48, -110, -59, 26},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 51, 114, 61, -27},

		/* 161-170 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, -133, 0, 0, 57},
		{0, 1, -1, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0},
		{0, 0, 0, 1, 0, -3, 5, 0, 0, 0, 0, 0, 0, -21, -6, -3, 11},
		{0, 1, -1, 2, 0, -3, 4, 0, 0, 0, 0, 0, 0, 0, -3, -1, 0},
		{0, 0, 0, 1, 0, 0, -2, 4, 0, 0, 0, 0, 0, -11, -21, -11, 6},
		{0, 2, -2, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0, -18, -436, -233, 9},
		{0, -1, 1, 0, 0, 5, -7, 0, 0, 0, 0, 0, 0, 35, -7, 0, 0},
		{0, 0, 0, 1, 0, 5, -8, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0},
		{-2, 0, 2, 1, 0, 6, -8, 0, 0, 0, 0, 0, 0, 11, -3, -1, -6},
		{0, 0, 0, 1, 0, 0, -8, 15, 0, 0, 0, 0, 0, -5, -3, -1, 3},

		/* 171-180 */
		{-2, 0, 2, 1, 0, 0, 2, 0, -3, 0, 0, 0, 0, -53, -9, -5, 28},
		{-2, 0, 2, 1, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0, 3, 2, 1},
		{1, 0, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 4, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, -4, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0, -50, 194, 103, 27},
		{0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, -13, 52, 28, 7},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -91, 248, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 6, 49, 26, -3},
		{0, 1, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, -6, -47, -25, 3},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 5, 3, 0},

		/* 181-190 */
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 52, 23, 10, -23},
		{0, 1, -1, 2, 0, 0, -1, 0, 0, -1, 0, 0, 0, -3, 0, 0, 1},
		{0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 5, 3, 0},
		{0, -1, 1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, -4, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, -7, 13, 0, 0, 0, 0, 2, -4, 8, 3, 2},
		{0, 0, 0, 0, 0, 0, 7, -13, 0, 0, 0, 0, 0, 10, 0, 0, 0},
		{2, 0, -2, 1, 0, 0, -5, 6, 0, 0, 0, 0, 0, 3, 0, 0, -2},
		{0, 2, -2, 1, 0, 0, -8, 11, 0, 0, 0, 0, 0, 0, 8, 4, 0},
		{0, 2, -2, 1, -1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 8, 4, 1},
		{-2, 0, 2, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, -4, 0, 0, 0},

		/* 191-200 */
		{0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, -4, 0, 0, 0},
		{0, 1, -1, 1, 0, 0, -1, 0, 0, 3, 0, 0, 0, -8, 4, 2, 4},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 8, -4, -2, -4},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0, 15, 7, 0},
		{-2, 0, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, -138, 0, 0, 0},
		{0, 0, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0, 0, -7, -3, 0},
		{0, 0, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, -7, -3, 0},
		{2, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0, 54, 0, 0, -29},
		{0, 1, -1, 2, 0, 0, -1, 0, 2, 0, 0, 0, 0, 0, 10, 4, 0},
		{0, 1, -1, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, -7, 0, 0, 3},

		/* 201-210 */
		{0, 0, 0, 1, 0, 0, 1, -2, 0, 0, 0, 0, 0, -37, 35, 19, 20},
		{0, -1, 1, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 4, 0, 0},
		{0, -1, 1, 0, 0, 0, 1, 0, 0, -2, 0, 0, 0, -4, 9, 0, 0},
		{0, 2, -2, 1, 0, 0, -2, 0, 0, 2, 0, 0, 0, 8, 0, 0, -4},
		{0, 1, -1, 1, 0, 3, -6, 0, 0, 0, 0, 0, 0, -9, -14, -8, 5},
		{0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 1, -3, -9, -5, 3},
		{0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 0, -145, 47, 0, 0},
		{0, 1, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0, -10, 40, 21, 5},
		{0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 1, 11, -49, -26, -7},
		{0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 2, -2150, 0, 0, 932},

		/* 211-220 */
		{0, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0, -12, 0, 0, 5},
		{0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 2, 85, 0, 0, -37},
		{0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 1, 4, 0, 0, -2},
		{0, 1, -1, 1, 0, 0, 1, -4, 0, 0, 0, 0, 0, 3, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, -86, 153, 0, 0},
		{0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 1, -6, 9, 5, 3},
		{0, 1, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0, 9, -13, -7, -5},
		{0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 1, -8, 12, 6, 4},
		{0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 2, -51, 0, 0, 22},
		{0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 2, -11, -268, -116, 5},

		/* 221-230 */
		{0, 2, -2, 2, 0, -5, 6, 0, 0, 0, 0, 0, 0, 0, 12, 5, 0},
		{0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 2, 0, 7, 3, 0},
		{0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 1, 31, 6, 3, -17},
		{0, 1, -1, 1, 0, -5, 7, 0, 0, 0, 0, 0, 0, 140, 27, 14, -75},
		{0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 1, 57, 11, 6, -30},
		{0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, 0, -14, -39, 0, 0},
		{0, 1, -1, 2, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -6, -2, 0},
		{0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 4, 15, 8, -2},
		{0, -1, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 4, 0, 0},
		{0, 2, -2, 1, 0, 0, -2, 0, 1, 0, 0, 0, 0, -3, 0, 0, 1},

		/* 231-240 */
		{0, 0, 0, 0, 0, 0, -6, 11, 0, 0, 0, 0, 2, 0, 11, 5, 0},
		{0, 0, 0, 0, 0, 0, 6, -11, 0, 0, 0, 0, 0, 9, 6, 0, 0},
		{0, 0, 0, 0, -1, 0, 4, 0, 0, 0, 0, 0, 2, -4, 10, 4, 2},
		{0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0, 0, 0, 5, 3, 0, 0},
		{2, 0, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0, 16, 0, 0, -9},
		{-2, 0, 2, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, -3, 0, 0, 0},
		{0, 2, -2, 1, 0, 0, -7, 9, 0, 0, 0, 0, 0, 0, 3, 2, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 2, 7, 0, 0, -3},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -25, 22, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 42, 223, 119, -22},

		/* 241-250 */
		{0, 1, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0, -27, -143, -77, 14},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 9, 49, 26, -5},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, -1166, 0, 0, 505},
		{0, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 2, -6, 0, 0, 3},
		{0, 0, 0, 1, 0, 3, -5, 0, 0, 0, 0, 0, 0, -8, 0, 1, 4},
		{0, -1, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0},
		{0, 2, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0, 117, 0, 0, -63},
		{0, 0, 0, 1, 0, 0, 2, -4, 0, 0, 0, 0, 0, -4, 8, 4, 2},
		{0, 2, -2, 1, 0, 0, -4, 4, 0, 0, 0, 0, 0, 3, 0, 0, -2},

		/* 251-260 */
		{0, 1, -1, 2, 0, -5, 7, 0, 0, 0, 0, 0, 0, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0, 31, 0, 0},
		{0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 1, -5, 0, 1, 3},
		{0, 1, -1, 1, 0, 0, -4, 6, 0, 0, 0, 0, 0, 4, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 1, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 2, -24, -13, -6, 10},
		{0, -1, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0},
		{0, 0, 0, 1, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, -32, -17, 0},
		{0, 0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 2, 8, 12, 5, -3},
		{0, 0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 1, 3, 0, 0, -1},

		/* 261-270 */
		{0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, 0, 7, 13, 0, 0},
		{0, -1, 1, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0, -3, 16, 0, 0},
		{0, 2, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0, 50, 0, 0, -27},
		{-2, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -5, -3, 0},
		{0, -2, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0},
		{0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0, 1, 0, 5, 3, 1},
		{0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0, 2, 24, 5, 2, -11},
		{0, 0, 0, 0, 0, -2, 3, 0, 0, 0, 0, 0, 2, 5, -11, -5, -2},
		{0, 0, 0, 0, 0, -2, 3, 0, 0, 0, 0, 0, 1, 30, -3, -2, -16},
		{0, 1, -1, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0, 18, 0, 0, -9},

		/* 271-280 */
		{0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 0, 8, 614, 0, 0},
		{0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 1, 3, -3, -1, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 6, 17, 9, -3},
		{0, 1, -1, 1, 0, 0, -1, 0, 3, 0, 0, 0, 0, -3, -9, -5, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 6, 3, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, -127, 21, 9, 55},
		{0, 0, 0, 0, 0, 0, 4, -8, 0, 0, 0, 0, 0, 3, 5, 0, 0},
		{0, 0, 0, 0, 0, 0, -4, 8, 0, 0, 0, 0, 2, -6, -10, -4, 3},
		{0, -2, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 2, 16, 9, 4, -7},

		/* 281-290 */
		{0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 1, 3, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, 0, 0, 22, 0, 0},
		{0, 0, 0, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0, 0, 19, 10, 0},
		{0, 2, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0, 7, 0, 0, -4},
		{0, 0, 0, 0, 0, 0, -5, 10, 0, 0, 0, 0, 2, 0, -5, -2, 0},
		{0, 0, 0, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2, -9, 3, 1, 4},
		{0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 2, 17, 0, 0, -7},
		{0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 1, 0, -3, -2, -1},
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, -20, 34, 0, 0},

		/* 291-300 */
		{0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 1, -10, 0, 1, 5},
		{0, 1, -1, 1, 0, 1, -3, 0, 0, 0, 0, 0, 0, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 22, -87, 0, 0},
		{0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 1, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 2, -3, -6, -2, 1},
		{0, 0, 0, 0, 0, -7, 11, 0, 0, 0, 0, 0, 2, -16, -3, -1, 7},
		{0, 0, 0, 0, 0, -7, 11, 0, 0, 0, 0, 0, 1, 0, -3, -2, 0},
		{0, -2, 2, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, -68, 39, 0, 0},
		{0, 2, -2, 1, 0, -4, 4, 0, 0, 0, 0, 0, 0, 27, 0, 0, -14},

		/* 301-310 */
		{0, -1, 1, 0, 0, 4, -5, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -25, 0, 0, 0},
		{0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 0, 1, -12, -3, -2, 6},
		{0, 1, -1, 1, 0, -4, 6, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 0, 2, 3, 66, 29, -1},
		{0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 2, 490, 0, 0, -213},
		{0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 1, -22, 93, 49, 12},
		{0, 1, -1, 1, 0, -4, 5, 0, 0, 0, 0, 0, 0, -7, 28, 15, 4},
		{0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 1, -3, 13, 7, 2},
		{0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, 0, -46, 14, 0, 0},

		/* 311-320 */
		{-2, 0, 2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 1, 0, 0},
		{0, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0},
		{0, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0, -28, 0, 0, 15},
		{0, 0, 0, 0, 0, 0, -1, 0, 5, 0, 0, 0, 2, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, 0, 3, 0, 0},
		{0, 0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 2, -11, 0, 0, 5},
		{0, 0, 0, 0, 0, 0, -7, 12, 0, 0, 0, 0, 2, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 25, 106, 57, -13},

		/* 321-330 */
		{0, 1, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 5, 21, 11, -3},
		{0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1485, 0, 0, 0},
		{0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1, -7, -32, -17, 4},
		{0, 1, -1, 1, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0},
		{0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 2, -6, -3, -2, 3},
		{0, 0, 0, 0, 0, 0, -1, 0, 4, 0, 0, 0, 2, 30, -6, -2, -13},
		{0, 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0, -4, 4, 0, 0},
		{0, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0, -19, 0, 0, 10},
		{0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 2, 0, 4, 2, -1},
		{0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0, 0, 3, 0, 0},

		/* 331-340 */
		{0, 2, -2, 1, 0, 0, -3, 0, 3, 0, 0, 0, 0, 4, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, -3, 7, 0, 0, 0, 0, 2, 0, -3, -1, 0},
		{-2, 0, 2, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 2, 5, 3, 1, -2},
		{0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, 0, 11, 0, 0},
		{0, 0, 0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 2, 118, 0, 0, -52},
		{0, 0, 0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 1, 0, -5, -3, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, -3, 0, 0, 0, 0, -28, 36, 0, 0},
		{0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0, 5, -5, 0, 0},
		{0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 0, 1, 14, -59, -31, -8},

		/* 341-350 */
		{0, 1, -1, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0, 0, 9, 5, 1},
		{0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 0, 2, -458, 0, 0, 198},
		{0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 0, 2, 0, -45, -20, 0},
		{0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 0, 1, 9, 0, 0, -5},
		{0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0},
		{0, 0, 0, 1, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, -4, -2, -1},
		{0, 2, -2, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0, 11, 0, 0, -6},
		{0, 0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 2, 6, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, -16, 23, 0, 0},
		{0, 0, 0, 1, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, -4, -2, 0},

		/* 351-360 */
		{0, 0, 0, 0, 0, 0, -1, 0, 2, 0, 0, 0, 2, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0, -166, 269, 0, 0},
		{0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 15, 0, 0, -8},
		{0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 0, 2, 10, 0, 0, -4},
		{0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, -78, 45, 0, 0},
		{0, 0, 0, 0, 0, -3, 4, 0, 0, 0, 0, 0, 2, 0, -5, -2, 0},
		{0, 0, 0, 0, 0, -3, 4, 0, 0, 0, 0, 0, 1, 7, 0, 0, -4},
		{0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0, -5, 328, 0, 0},
		{0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 1, 3, 0, 0, -2},
		{0, 0, 0, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0, 5, 0, 0, -2},

		/* 361-370 */
		{0, 0, 0, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, -3, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 1, -5, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, -4, -2, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1223, -26, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 1, 0, 7, 3, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, -3, 5, 0, 0, 0, 3, 0, 0, 0},
		{0, 0, 0, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 0, 0, 0, -6, 20, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, -368, 0, 0, 0},

		/* 371-380 */
		{0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, -75, 0, 0, 0},
		{0, 0, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 11, 0, 0, -6},
		{0, 0, 0, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0, 3, 0, 0, -2},
		{0, 0, 0, 0, 0, -8, 14, 0, 0, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 1, 0, 2, -5, 0, 0, 0, -13, -30, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 0, 21, 3, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 8, -27, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -8, 3, 0, 0, 0, 0, -19, -11, 0, 0},

		/* 381-390 */
		{0, 0, 0, 0, 0, 0, -3, 8, -3, 0, 0, 0, 2, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, 0, -2, 5, 0, 0, 2, 0, 5, 2, 0},
		{0, 0, 0, 0, 0, -8, 12, 0, 0, 0, 0, 0, 2, -6, 0, 0, 2},
		{0, 0, 0, 0, 0, -8, 12, 0, 0, 0, 0, 0, 0, -8, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 1, -2, 0, 0, 0, -1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, -14, 0, 0, 6},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 6, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, -74, 0, 0, 32},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, -3, -1, 0},
		{0, 2, -2, 1, 0, -5, 5, 0, 0, 0, 0, 0, 0, 4, 0, 0, -2},

		/* 391-400 */
		{0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 8, 11, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 3, 2, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2, -262, 0, 0, 114},
		{0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0},
		{0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 0, 1, -7, 0, 0, 4},
		{0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 0, 2, 0, -27, -12, 0},
		{0, 0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 2, -19, -8, -4, 8},
		{0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 0, 2, 202, 0, 0, -87},
		{0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 0, 1, -8, 35, 19, 5},
		{0, 1, -1, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0},

		/* 401-410 */
		{0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, 0, 16, -5, 0, 0},
		{0, 2, -2, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 5, 0, 0, -3},
		{0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -3, 0, 0},
		{0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, -35, -48, -21, 15},
		{0, 0, 0, 0, 0, 0, -2, 6, 0, 0, 0, 0, 2, -3, -5, -2, 1},
		{0, 0, 0, 1, 0, 2, -2, 0, 0, 0, 0, 0, 0, 6, 0, 0, -3},
		{0, 0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 2, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, 0, -5, 0, 0},
		{0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0, 0, 1, 12, 55, 29, -6},

		/* 411-420 */
		{0, 1, -1, 1, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0},
		{0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, -598, 0, 0, 0},
		{0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 1, -3, -13, -7, 1},
		{0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2, -5, -7, -3, 2},
		{0, 0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 2, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, 5, -7, 0, 0},
		{0, 0, 0, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0, 4, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, 16, -6, 0, 0},
		{0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, 0, 8, -3, 0, 0},
		{0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 0, 1, 8, -31, -16, -4},

		/* 421-430 */
		{0, 1, -1, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 0, 2, 113, 0, 0, -49},
		{0, 0, 0, 0, 0, -7, 10, 0, 0, 0, 0, 0, 2, 0, -24, -10, 0},
		{0, 0, 0, 0, 0, -7, 10, 0, 0, 0, 0, 0, 1, 4, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 27, 0, 0, 0},
		{0, 0, 0, 0, 0, -4, 8, 0, 0, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, -4, 5, 0, 0, 0, 0, 0, 2, 0, -4, -2, 0},
		{0, 0, 0, 0, 0, -4, 5, 0, 0, 0, 0, 0, 1, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, -13, 0, 0, 6},

		/* 431-440 */
		{0, 0, 0, 0, 0, 0, -2, 0, 5, 0, 0, 0, 2, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2, -18, -10, -4, 8},
		{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -4, -28, 0, 0},
		{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, -5, 6, 3, 2},
		{0, 0, 0, 0, 0, -9, 13, 0, 0, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, -1, 5, 0, 0, 0, 0, 2, -5, -9, -4, 2},
		{0, 0, 0, 0, 0, 0, -2, 0, 4, 0, 0, 0, 2, 17, 0, 0, -7},
		{0, 0, 0, 0, 0, 0, 2, 0, -4, 0, 0, 0, 0, 11, 4, 0, 0},
		{0, 0, 0, 0, 0, 0, -2, 7, 0, 0, 0, 0, 2, 0, -6, -2, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0, 83, 15, 0, 0},

		/* 441-450 */
		{0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 0, 1, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 0, 2, 0, -114, -49, 0},
		{0, 0, 0, 0, 0, -6, 8, 0, 0, 0, 0, 0, 2, 117, 0, 0, -51},
		{0, 0, 0, 0, 0, -6, 8, 0, 0, 0, 0, 0, 1, -5, 19, 10, 2},
		{0, 0, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0, -3, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, -3, 9, 0, 0, 0, 0, 2, 0, -3, -1, 0},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 2, 0, -6, -2, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 393, 3, 0, 0},

		/* 451-460 */
		{0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 1, -4, 21, 11, 2},
		{0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 2, -6, 0, -1, 3},
		{0, 0, 0, 0, 0, -5, 10, 0, 0, 0, 0, 0, 2, -3, 8, 4, 1},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 8, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 2, 18, -29, -13, -8},
		{0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 1, 8, 34, 18, -4},
		{0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 89, 0, 0, 0},
		{0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 1, 3, 12, 6, -1},
		{0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 2, 54, -15, -7, -24},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, -3, 0, 0, 0, 0, 3, 0, 0},

		/* 461-470 */
		{0, 0, 0, 0, 0, 0, -5, 13, 0, 0, 0, 0, 2, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 0, 0, 35, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 2, -154, -30, -13, 67},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 15, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 1, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 0, 9, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 2, 80, -71, -31, -35},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, 2, 0, -20, -9, 0},
		{0, 0, 0, 0, 0, 0, -6, 15, 0, 0, 0, 0, 2, 11, 5, 2, -5},
		{0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 0, 2, 61, -96, -42, -27},

		/* 471-480 */
		{0, 0, 0, 0, 0, -3, 9, -4, 0, 0, 0, 0, 2, 14, 9, 4, -6},
		{0, 0, 0, 0, 0, 0, 2, 0, 2, -5, 0, 0, 2, -11, -6, -3, 5},
		{0, 0, 0, 0, 0, 0, -2, 8, -1, -5, 0, 0, 2, 0, -3, -1, 0},
		{0, 0, 0, 0, 0, 0, 6, -8, 3, 0, 0, 0, 2, 123, -415, -180, -53},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -35},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 7, -32, -17, -4},
		{0, 1, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -9, -5, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, -4, 2, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, -89, 0, 0, 38},

		/* 481-490 */
		{0, 0, 0, 0, 0, 0, -6, 16, -4, -5, 0, 0, 2, 0, -86, -19, -6},
		{0, 0, 0, 0, 0, 0, -2, 8, -3, 0, 0, 0, 2, 0, 0, -19, 6},
		{0, 0, 0, 0, 0, 0, -2, 8, -3, 0, 0, 0, 2, -123, -416, -180, 53},
		{0, 0, 0, 0, 0, 0, 6, -8, 1, 5, 0, 0, 2, 0, -3, -1, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, -2, 5, 0, 0, 2, 12, -6, -3, -5},
		{0, 0, 0, 0, 0, 3, -5, 4, 0, 0, 0, 0, 2, -13, 9, 4, 6},
		{0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 2, 0, -15, -7, 0},
		{0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 1, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 2, -62, -97, -42, 27},
		{0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 2, -11, 5, 2, 5},

		/* 491-500 */
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, -19, -8, 0},
		{0, 0, 0, 0, 0, 3, -3, 0, 2, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 2, -2, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 4, 2, 0},
		{0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0},
		{0, 2, -2, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2, -85, -70, -31, 37},
		{0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2, 163, -12, -5, -72},
		{0, 0, 0, 0, 0, -3, 7, 0, 0, 0, 0, 0, 2, -63, -16, -7, 28},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 2, -21, -32, -14, 9},
		{0, 0, 0, 0, 0, -5, 6, 0, 0, 0, 0, 0, 2, 0, -3, -1, 0},

		/* 501-510 */
		{0, 0, 0, 0, 0, -5, 6, 0, 0, 0, 0, 0, 1, 3, 0, 0, -2},
		{0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0},
		{0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2, 3, 10, 4, -1},
		{0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, -1, 6, 0, 0, 0, 0, 2, 0, -7, -3, 0},
		{0, 0, 0, 0, 0, 0, 7, -9, 0, 0, 0, 0, 2, 0, -4, -2, 0},
		{0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 6, 19, 0, 0},
		{0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 2, 5, -173, -75, -2},
		{0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 2, 0, -7, -3, 0},
		{0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 2, 7, -12, -5, -3},

		/* 511-520 */
		{0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 0, 1, -3, 0, 0, 2},
		{0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 0, 2, 3, -4, -2, -1},
		{0, 0, 0, 0, 0, -7, 9, 0, 0, 0, 0, 0, 2, 74, 0, 0, -32},
		{0, 0, 0, 0, 0, -7, 9, 0, 0, 0, 0, 0, 1, -3, 12, 6, 2},
		{0, 0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 2, 26, -14, -6, -11},
		{0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 2, 19, 0, 0, -8},
		{0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 0, 1, 6, 24, 13, -3},
		{0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0, 83, 0, 0, 0},
		{0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 1, 0, -10, -5, 0},
		{0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 2, 11, -3, -1, -5},

		/* 521-530 */
		{0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2, 3, 0, 1, -1},
		{0, 0, 0, 0, 0, 0, -3, 0, 5, 0, 0, 0, 2, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, -4, 0, 0, 0},
		{0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 5, -23, -12, -3},
		{0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, -339, 0, 0, 147},
		{0, 0, 0, 0, 0, -9, 12, 0, 0, 0, 0, 0, 2, 0, -10, -5, 0},
		{0, 0, 0, 0, 0, 0, 3, 0, -4, 0, 0, 0, 0, 5, 0, 0, 0},
		{0, 2, -2, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 2, 0, -4, -2, 0},
		{0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 0, 18, -3, 0, 0},

		/* 531-540 */
		{0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 2, 9, -11, -5, -4},
		{0, 0, 0, 0, 0, -2, 6, 0, 0, 0, 0, 0, 2, -8, 0, 0, 4},
		{0, 0, 0, 0, 0, -6, 7, 0, 0, 0, 0, 0, 1, 3, 0, 0, -1},
		{0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 2, 6, -9, -4, -2},
		{0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 0, -4, -12, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 2, 67, -91, -39, -29},
		{0, 0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 2, 30, -18, -8, -13},
		{0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 2, 0, -114, -50, 0},

		/* 541-550 */
		{0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 0, 0, 2, 0, 0, 0, 23},
		{0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 0, 0, 2, 517, 16, 7, -224},
		{0, 0, 0, 0, 0, 0, 3, 0, 0, -2, 0, 0, 2, 0, -7, -3, 0},
		{0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 2, 143, -3, -1, -62},
		{0, 0, 0, 0, 0, 0, 3, 0, 0, -1, 0, 0, 2, 29, 0, 0, -13},
		{0, 2, -2, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0, -4, 0, 0, 2},
		{0, 0, 0, 0, 0, -8, 16, 0, 0, 0, 0, 0, 2, -6, 0, 0, 3},
		{0, 0, 0, 0, 0, 0, 3, 0, 2, -5, 0, 0, 2, 5, 12, 5, -2},
		{0, 0, 0, 0, 0, 0, 7, -8, 3, 0, 0, 0, 2, -25, 0, 0, 11},
		{0, 0, 0, 0, 0, 0, -5, 16, -4, -5, 0, 0, 2, -3, 0, 0, 1},

		/* 551-560 */
		{0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, 0, -1, 8, -3, 0, 0, 0, 2, -22, 12, 5, 10},
		{0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 2, 50, 0, 0, -22},
		{0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 1, 0, 7, 4, 0},
		{0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 2, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 2, -4, 4, 2, 2},
		{0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 2, -5, -11, -5, 2},
		{0, 0, 0, 0, 0, -3, 8, 0, 0, 0, 0, 0, 2, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, -5, 5, 0, 0, 0, 0, 0, 1, 4, 17, 9, -2},
		{0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 0, 59, 0, 0, 0},

		/* 561-570 */
		{0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 1, 0, -4, -2, 0},
		{0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 2, -8, 0, 0, 4},
		{0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 4, -15, -8, -2},
		{0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 370, -8, 0, -160},
		{0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 2, 0, 0, -3, 0},
		{0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 2, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, 0, 6, -5, 0, 0, 0, 0, 2, -6, 3, 1, 3},
		{0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 2, -10, 0, 0, 4},

		/* 571-580 */
		{0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 0, 2, 0, 9, 4, 0},
		{0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 2, 4, 17, 7, -2},
		{0, 0, 0, 0, 0, -9, 11, 0, 0, 0, 0, 0, 2, 34, 0, 0, -15},
		{0, 0, 0, 0, 0, -9, 11, 0, 0, 0, 0, 0, 1, 0, 5, 3, 0},
		{0, 0, 0, 0, 0, 0, 4, 0, -4, 0, 0, 0, 2, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 4, 0, -3, 0, 0, 0, 2, -37, -7, -3, 16},
		{0, 0, 0, 0, 0, -6, 6, 0, 0, 0, 0, 0, 1, 3, 13, 7, -2},
		{0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0},
		{0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 1, 0, -3, -2, 0},
		{0, 0, 0, 0, 0, 0, 4, 0, -2, 0, 0, 0, 2, -184, -3, -1, 80},

		/* 581-590 */
		{0, 0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 2, -3, 0, 0, 1},
		{0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 1, 0, -10, -6, -1},
		{0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 2, 31, -6, 0, -13},
		{0, 0, 0, 0, 0, 0, 4, 0, -1, 0, 0, 0, 2, -3, -32, -14, 1},
		{0, 0, 0, 0, 0, 0, 4, 0, 0, -2, 0, 0, 2, -7, 0, 0, 3},
		{0, 0, 0, 0, 0, 0, 5, -2, 0, 0, 0, 0, 2, 0, -8, -4, 0},
		{0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 3, -4, 0, 0},
		{0, 0, 0, 0, 0, 8, -9, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0},
		{0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 0, 2, 0, 3, 1, 0},

		/* 591-600 */
		{0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 2, 19, -23, -10, 2},
		{0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -10},
		{0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1, 0, 3, 2, 0},
		{0, 0, 0, 0, 0, -7, 7, 0, 0, 0, 0, 0, 1, 0, 9, 5, -1},
		{0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 0, 0, 28, 0, 0, 0},
		{0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 1, 0, -7, -4, 0},
		{0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 2, 8, -4, 0, -4},
		{0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0},
		{0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, 0, -4, 0, 0, 0, 2, -3, 0, 0, 1},

		/* 601-610 */
		{0, 0, 0, 0, 0, 0, 5, 0, -3, 0, 0, 0, 2, -9, 0, 1, 4},
		{0, 0, 0, 0, 0, 0, 5, 0, -2, 0, 0, 0, 2, 3, 12, 5, -1},
		{0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 17, -3, -1, 0},
		{0, 0, 0, 0, 0, -8, 8, 0, 0, 0, 0, 0, 1, 0, 7, 4, 0},
		{0, 0, 0, 0, 0, 8, -8, 0, 0, 0, 0, 0, 0, 19, 0, 0, 0},
		{0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 1, 0, -5, -3, 0},
		{0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 2, 14, -3, 0, -1},
		{0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0},
		{0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1, 0, 0, 0, -5},
		{0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1, 0, 5, 3, 0},

		/* 611-620 */
		{0, 0, 0, 0, 0, 9, -9, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0},
		{0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 0, 1, 0, -3, -2, 0},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 2, 9, 4, 3},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0, 4, 2, 0},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 6, 0, 0, -3},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0, 3, 1, 0},
		{0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 5, 0, 0, -2},

		/* 621-630 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, -1},
		{1, 0, -2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, -3, 0, 0, 0},
		{1, 0, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0},
		{1, 0, -2, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 7, 0, 0, 0},
		{1, 0, -2, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -4, 0, 0, 0},
		{-1, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0},
		{-1, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 6, 0, 0, 0},
		{-1, 0, 2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, -4, 0, 0},
		{1, 0, -2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, -4, 0, 0},
		{-2, 0, 2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, 5, 0, 0, 0},

		/* 631-640 */
		{-1, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0, -3, 0, 0, 0},
		{-1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 4, 0, 0, 0},
		{-1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0},
		{-1, 0, 2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0},
		{1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0},
		{-1, 0, 2, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0, 13, 0, 0, 0},
		{-2, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0, 21, 11, 0, 0},
		{1, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, -5, 0, 0},
		{-1, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -5, -2, 0},
		{1, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0},

		/* 641-650 */
		{-1, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, -5, 0, 0},
		{-1, 0, 2, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0, -3, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 20, 10, 0, 0},
		{-1, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, -34, 0, 0, 0},
		{-1, 0, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, -19, 0, 0, 0},
		{1, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0, 3, 0, 0, -2},
		{1, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0, -3, 0, 0, 1},
		{1, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0, -6, 0, 0, 3},
		{1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -4, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 3, 0, 0, 0},

		/* 651-660 */
		{0, 0, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0},
		{0, 0, -2, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 4, 0, 0, 0},
		{0, 2, 0, 2, 0, -2, 2, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 2, 0, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0, 6, 0, 0, -3},
		{0, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0, -8, 0, 0, 3},
		{0, 2, 0, 2, 0, -2, 3, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, -3, 0, 0, 0},
		{0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -3, -2, 0},
		{1, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 126, -63, -27, -55},
		{-1, 2, 0, 2, 0, 10, -3, 0, 0, 0, 0, 0, 0, -5, 0, 1, 2},

		/* 661-670 */
		{0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3, 28, 15, 2},
		{1, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5, 0, 1, -2},
		{0, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 9, 4, 1},
		{0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0, 0, 9, 4, -1},
		{-1, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0, -126, -63, -27, 55},
		{2, 2, -2, 2, 0, 0, -2, 0, 3, 0, 0, 0, 0, 3, 0, 0, -1},
		{1, 2, 0, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0, 21, -11, -6, -11},
		{0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0},
		{-1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -21, -11, -6, 11},
		{-2, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0, -3, 0, 0, 1},

		/* 671-680 */
		{0, 2, 0, 2, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
		{0, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0, 8, 0, 0, -4},
		{0, 2, 0, 2, 0, 0, 1, 0, -1, 0, 0, 0, 0, -6, 0, 0, 3},
		{0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 0, 0, 0, -3, 0, 0, 1},
		{-1, 2, 2, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0, 3, 0, 0, -1},
		{1, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0, -3, 0, 0, 1},
		{-1, 2, 2, 2, 0, 0, 2, 0, -3, 0, 0, 0, 0, -5, 0, 0, 2},
		{2, 2, 0, 2, 0, 0, 2, 0, -3, 0, 0, 0, 0, 24, -12, -5, -11},
		{1, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0, 0, 3, 1, 0},
		{1, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 3, 1, 0},

		/* 681-687 */
		{1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0},
		{0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, -24, -12, -5, 10},
		{2, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 0, -1, -2},
		{-1, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0, 13, 0, 0, -6},
		{-1, 2, 2, 2, 0, 3, -3, 0, 0, 0, 0, 0, 0, 7, 0, 0, -3},
		{1, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1},
		{0, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0, 3, 0, 0, -1},
	}

	/* Number of terms in the planetary nutation model */
	NPL := len(xpl)

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* ------------------- */
	/* LUNI-SOLAR NUTATION */
	/* ------------------- */

	/* Fundamental (Delaunay) arguments */

	/* Mean anomaly of the Moon (IERS 2003). */
	el = Fal03(t)

	/* Mean anomaly of the Sun (MHB2000). */
	elp = fmod(1287104.79305+
		t*(129596581.0481+
			t*(-0.5532+
				t*(0.000136+
					t*(-0.00001149)))), TURNAS) * DAS2R

	/* Mean longitude of the Moon minus that of the ascending node */
	/* (IERS 2003. */
	f = Faf03(t)

	/* Mean elongation of the Moon from the Sun (MHB2000). */
	d = fmod(1072260.70369+
		t*(1602961601.2090+
			t*(-6.3706+
				t*(0.006593+
					t*(-0.00003169)))), TURNAS) * DAS2R

	/* Mean longitude of the ascending node of the Moon (IERS 2003). */
	om = Faom03(t)

	/* Initialize the nutation values. */
	dp = 0.0
	de = 0.0

	/* Summation of luni-solar nutation series (in reverse order). */
	for i = NLS - 1; i >= 0; i-- {

		/* Argument and functions. */
		arg = fmod(float64(xls[i].nl)*el+
			float64(xls[i].nlp)*elp+
			float64(xls[i].nf)*f+
			float64(xls[i].nd)*d+
			float64(xls[i].nom)*om, D2PI)
		sarg = sin(arg)
		carg = cos(arg)

		/* Term. */
		dp += (xls[i].sp+xls[i].spt*t)*sarg + xls[i].cp*carg
		de += (xls[i].ce+xls[i].cet*t)*carg + xls[i].se*sarg
	}

	/* Convert from 0.1 microarcsec units to radians. */
	dpsils = dp * U2R
	depsls = de * U2R

	/* ------------------ */
	/* PLANETARY NUTATION */
	/* ------------------ */

	/* n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
	/* in different functions, using slightly different Delaunay */
	/* arguments in the two cases.  This behaviour is faithfully */
	/* reproduced here.  Use of the IERS 2003 expressions for both */
	/* cases leads to negligible changes, well below */
	/* 0.1 microarcsecond. */

	/* Mean anomaly of the Moon (MHB2000). */
	al = fmod(2.35555598+8328.6914269554*t, D2PI)

	/* Mean longitude of the Moon minus that of the ascending node */
	/*(MHB2000). */
	af = fmod(1.627905234+8433.466158131*t, D2PI)

	/* Mean elongation of the Moon from the Sun (MHB2000). */
	ad = fmod(5.198466741+7771.3771468121*t, D2PI)

	/* Mean longitude of the ascending node of the Moon (MHB2000). */
	aom = fmod(2.18243920-33.757045*t, D2PI)

	/* General accumulated precession in longitude (IERS 2003). */
	apa = Fapa03(t)

	/* Planetary longitudes, Mercury through Uranus (IERS 2003). */
	alme = Fame03(t)
	alve = Fave03(t)
	alea = Fae03(t)
	alma = Fama03(t)
	alju = Faju03(t)
	alsa = Fasa03(t)
	alur = Faur03(t)

	/* Neptune longitude (MHB2000). */
	alne = fmod(5.321159000+3.8127774000*t, D2PI)

	/* Initialize the nutation values. */
	dp = 0.0
	de = 0.0

	/* Summation of planetary nutation series (in reverse order). */
	for i = NPL - 1; i >= 0; i-- {

		/* Argument and functions. */
		arg = fmod(float64(xpl[i].nl)*al+
			float64(xpl[i].nf)*af+
			float64(xpl[i].nd)*ad+
			float64(xpl[i].nom)*aom+
			float64(xpl[i].nme)*alme+
			float64(xpl[i].nve)*alve+
			float64(xpl[i].nea)*alea+
			float64(xpl[i].nma)*alma+
			float64(xpl[i].nju)*alju+
			float64(xpl[i].nsa)*alsa+
			float64(xpl[i].nur)*alur+
			float64(xpl[i].nne)*alne+
			float64(xpl[i].npa)*apa, D2PI)
		sarg = sin(arg)
		carg = cos(arg)

		/* Term. */
		dp += float64(xpl[i].sp)*sarg + float64(xpl[i].cp)*carg
		de += float64(xpl[i].se)*sarg + float64(xpl[i].ce)*carg

	}

	/* Convert from 0.1 microarcsec units to radians. */
	dpsipl = dp * U2R
	depspl = de * U2R

	/* ------- */
	/* RESULTS */
	/* ------- */

	/* Add luni-solar and planetary components. */
	*dpsi = dpsils + dpsipl
	*deps = depsls + depspl
}

/*
Nut00b Nutation, IAU 2000B Nutation, IAU 2000B model.

Given:
    date1,date2   float64    TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps     float64    nutation, luni-solar + planetary (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components in longitude and obliquity are in radians
    and with respect to the equinox and ecliptic of date.  The
    obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
    value of 84381.448 arcsec.  (The errors that result from using
    this function with the IAU 2006 value of 84381.406 arcsec can be
    neglected.)

    The nutation model consists only of luni-solar terms, but
    includes also a fixed offset which compensates for certain long-
    period planetary terms (Note 7).

 3) This function is an implementation of the IAU 2000B abridged
    nutation model formally adopted by the IAU General Assembly in
    2000.  The function computes the MHB_2000_SHORT luni-solar
    nutation series (Luzum 2001), but without the associated
    corrections for the precession rate adjustments and the offset
    between the GCRS and J2000.0 mean poles.

 4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
    terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
    77 terms, plus additional simplifications, yet still delivers
    results of 1 mas accuracy at present epochs.  This combination of
    accuracy and size makes the IAU 2000B abridged nutation model
    suitable for most practical applications.

    The function delivers a pole accurate to 1 mas from 1900 to 2100
    (usually better than 1 mas, very occasionally just outside
    1 mas).  The full IAU 2000A model, which is implemented in the
    function iauNut00a (q.v.), delivers considerably greater accuracy
    at current dates;  however, to realize this improved accuracy,
    corrections for the essentially unpredictable free-core-nutation
    (FCN) must also be included.

 5) The present function provides classical nutation.  The
    MHB_2000_SHORT algorithm, from which it is adapted, deals also
    with (i) the offsets between the GCRS and mean poles and (ii) the
    adjustments in longitude and obliquity due to the changed
    precession rates.  These additional functions, namely frame bias
    and precession adjustments, are supported by the SOFA functions
    Bi00  and iauPr00.

 6) The MHB_2000_SHORT algorithm also provides "total" nutations,
    comprising the arithmetic sum of the frame bias, precession
    adjustments, and nutation (luni-solar + planetary).  These total
    nutations can be used in combination with an existing IAU 1976
    precession implementation, such as iauPmat76,  to deliver GCRS-
    to-true predictions of mas accuracy at current epochs.  However,
    for symmetry with the iauNut00a  function (q.v. for the reasons),
    the SOFA functions do not generate the "total nutations"
    directly.  Should they be required, they could of course easily
    be generated by calling iauBi00, iauPr00 and the present function
    and adding the results.

 7) The IAU 2000B model includes "planetary bias" terms that are
    fixed in size but compensate for long-period nutations.  The
    amplitudes quoted in McCarthy & Luzum (2003), namely
    Dpsi = -1.5835 mas and Depsilon = +1.6339 mas, are optimized for
    the "total nutations" method described in Note 6.  The Luzum
    (2001) values used in this SOFA implementation, namely -0.135 mas
    and +0.388 mas, are optimized for the "rigorous" method, where
    frame bias, precession and nutation are applied separately and in
    that order.  During the interval 1995-2050, the SOFA
    implementation delivers a maximum error of 1.001 mas (not
    including FCN).

References:

    Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
    for the precession quantities based upon the IAU /1976/ system of
    astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)

    Luzum, B., private communication, 2001 (Fortran code
    MHB_2000_SHORT)

    McCarthy, D.D. & Luzum, B.J., "An abridged model of the
    precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
    85, 37-49 (2003)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
*/
func Nut00b(date1, date2 float64, dpsi, deps *float64) {
	var t, el, elp, f, d, om, arg, dp, de, sarg, carg,
		dpsils, depsls, dpsipl, depspl float64
	var i int

	/* Units of 0.1 microarcsecond to radians */
	const U2R = DAS2R / 1e7

	/* ---------------------------------------- */
	/* Fixed offsets in lieu of planetary terms */
	/* ---------------------------------------- */

	const DPPLAN = -0.135 * DMAS2R
	const DEPLAN = 0.388 * DMAS2R

	/* --------------------------------------------------- */
	/* Luni-solar nutation: argument and term coefficients */
	/* --------------------------------------------------- */

	/* The units for the sine and cosine coefficients are */
	/* 0.1 microarcsec and the same per Julian century    */

	x := []struct {
		nl, nlp, nf, nd, nom int     /* coefficients of l,l',F,D,Om */
		ps, pst, pc          float64 /* longitude sin, t*sin, cos coefficients */
		ec, ect, es          float64 /* obliquity cos, t*cos, sin coefficients */

	}{

		/* 1-10 */
		{0, 0, 0, 0, 1,
			-172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
		{0, 0, 2, -2, 2,
			-13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
		{0, 0, 2, 0, 2, -2276413.0, -234.0, 2796.0, 978459.0, -485.0, 1374.0},
		{0, 0, 0, 0, 2, 2074554.0, 207.0, -698.0, -897492.0, 470.0, -291.0},
		{0, 1, 0, 0, 0, 1475877.0, -3633.0, 11817.0, 73871.0, -184.0, -1924.0},
		{0, 1, 2, -2, 2, -516821.0, 1226.0, -524.0, 224386.0, -677.0, -174.0},
		{1, 0, 0, 0, 0, 711159.0, 73.0, -872.0, -6750.0, 0.0, 358.0},
		{0, 0, 2, 0, 1, -387298.0, -367.0, 380.0, 200728.0, 18.0, 318.0},
		{1, 0, 2, 0, 2, -301461.0, -36.0, 816.0, 129025.0, -63.0, 367.0},
		{0, -1, 2, -2, 2, 215829.0, -494.0, 111.0, -95929.0, 299.0, 132.0},

		/* 11-20 */
		{0, 0, 2, -2, 1, 128227.0, 137.0, 181.0, -68982.0, -9.0, 39.0},
		{-1, 0, 2, 0, 2, 123457.0, 11.0, 19.0, -53311.0, 32.0, -4.0},
		{-1, 0, 0, 2, 0, 156994.0, 10.0, -168.0, -1235.0, 0.0, 82.0},
		{1, 0, 0, 0, 1, 63110.0, 63.0, 27.0, -33228.0, 0.0, -9.0},
		{-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0, 31429.0, 0.0, -75.0},
		{-1, 0, 2, 2, 2, -59641.0, -11.0, 149.0, 25543.0, -11.0, 66.0},
		{1, 0, 2, 0, 1, -51613.0, -42.0, 129.0, 26366.0, 0.0, 78.0},
		{-2, 0, 2, 0, 1, 45893.0, 50.0, 31.0, -24236.0, -10.0, 20.0},
		{0, 0, 0, 2, 0, 63384.0, 11.0, -150.0, -1220.0, 0.0, 29.0},
		{0, 0, 2, 2, 2, -38571.0, -1.0, 158.0, 16452.0, -11.0, 68.0},

		/* 21-30 */
		{0, -2, 2, -2, 2, 32481.0, 0.0, 0.0, -13870.0, 0.0, 0.0},
		{-2, 0, 0, 2, 0, -47722.0, 0.0, -18.0, 477.0, 0.0, -25.0},
		{2, 0, 2, 0, 2, -31046.0, -1.0, 131.0, 13238.0, -11.0, 59.0},
		{1, 0, 2, -2, 2, 28593.0, 0.0, -1.0, -12338.0, 10.0, -3.0},
		{-1, 0, 2, 0, 1, 20441.0, 21.0, 10.0, -10758.0, 0.0, -3.0},
		{2, 0, 0, 0, 0, 29243.0, 0.0, -74.0, -609.0, 0.0, 13.0},
		{0, 0, 2, 0, 0, 25887.0, 0.0, -66.0, -550.0, 0.0, 11.0},
		{0, 1, 0, 0, 1, -14053.0, -25.0, 79.0, 8551.0, -2.0, -45.0},
		{-1, 0, 0, 2, 1, 15164.0, 10.0, 11.0, -8001.0, 0.0, -1.0},
		{0, 2, 2, -2, 2, -15794.0, 72.0, -16.0, 6850.0, -42.0, -5.0},

		/* 31-40 */
		{0, 0, -2, 2, 0, 21783.0, 0.0, 13.0, -167.0, 0.0, 13.0},
		{1, 0, 0, -2, 1, -12873.0, -10.0, -37.0, 6953.0, 0.0, -14.0},
		{0, -1, 0, 0, 1, -12654.0, 11.0, 63.0, 6415.0, 0.0, 26.0},
		{-1, 0, 2, 2, 1, -10204.0, 0.0, 25.0, 5222.0, 0.0, 15.0},
		{0, 2, 0, 0, 0, 16707.0, -85.0, -10.0, 168.0, -1.0, 10.0},
		{1, 0, 2, 2, 2, -7691.0, 0.0, 44.0, 3268.0, 0.0, 19.0},
		{-2, 0, 2, 0, 0, -11024.0, 0.0, -14.0, 104.0, 0.0, 2.0},
		{0, 1, 2, 0, 2, 7566.0, -21.0, -11.0, -3250.0, 0.0, -5.0},
		{0, 0, 2, 2, 1, -6637.0, -11.0, 25.0, 3353.0, 0.0, 14.0},
		{0, -1, 2, 0, 2, -7141.0, 21.0, 8.0, 3070.0, 0.0, 4.0},

		/* 41-50 */
		{0, 0, 0, 2, 1, -6302.0, -11.0, 2.0, 3272.0, 0.0, 4.0},
		{1, 0, 2, -2, 1, 5800.0, 10.0, 2.0, -3045.0, 0.0, -1.0},
		{2, 0, 2, -2, 2, 6443.0, 0.0, -7.0, -2768.0, 0.0, -4.0},
		{-2, 0, 0, 2, 1, -5774.0, -11.0, -15.0, 3041.0, 0.0, -5.0},
		{2, 0, 2, 0, 1, -5350.0, 0.0, 21.0, 2695.0, 0.0, 12.0},
		{0, -1, 2, -2, 1, -4752.0, -11.0, -3.0, 2719.0, 0.0, -3.0},
		{0, 0, 0, -2, 1, -4940.0, -11.0, -21.0, 2720.0, 0.0, -9.0},
		{-1, -1, 0, 2, 0, 7350.0, 0.0, -8.0, -51.0, 0.0, 4.0},
		{2, 0, 0, -2, 1, 4065.0, 0.0, 6.0, -2206.0, 0.0, 1.0},
		{1, 0, 0, 2, 0, 6579.0, 0.0, -24.0, -199.0, 0.0, 2.0},

		/* 51-60 */
		{0, 1, 2, -2, 1, 3579.0, 0.0, 5.0, -1900.0, 0.0, 1.0},
		{1, -1, 0, 0, 0, 4725.0, 0.0, -6.0, -41.0, 0.0, 3.0},
		{-2, 0, 2, 0, 2, -3075.0, 0.0, -2.0, 1313.0, 0.0, -1.0},
		{3, 0, 2, 0, 2, -2904.0, 0.0, 15.0, 1233.0, 0.0, 7.0},
		{0, -1, 0, 2, 0, 4348.0, 0.0, -10.0, -81.0, 0.0, 2.0},
		{1, -1, 2, 0, 2, -2878.0, 0.0, 8.0, 1232.0, 0.0, 4.0},
		{0, 0, 0, 1, 0, -4230.0, 0.0, 5.0, -20.0, 0.0, -2.0},
		{-1, -1, 2, 2, 2, -2819.0, 0.0, 7.0, 1207.0, 0.0, 3.0},
		{-1, 0, 2, 0, 0, -4056.0, 0.0, 5.0, 40.0, 0.0, -2.0},
		{0, -1, 2, 2, 2, -2647.0, 0.0, 11.0, 1129.0, 0.0, 5.0},

		/* 61-70 */
		{-2, 0, 0, 0, 1, -2294.0, 0.0, -10.0, 1266.0, 0.0, -4.0},
		{1, 1, 2, 0, 2, 2481.0, 0.0, -7.0, -1062.0, 0.0, -3.0},
		{2, 0, 0, 0, 1, 2179.0, 0.0, -2.0, -1129.0, 0.0, -2.0},
		{-1, 1, 0, 1, 0, 3276.0, 0.0, 1.0, -9.0, 0.0, 0.0},
		{1, 1, 0, 0, 0, -3389.0, 0.0, 5.0, 35.0, 0.0, -2.0},
		{1, 0, 2, 0, 0, 3339.0, 0.0, -13.0, -107.0, 0.0, 1.0},
		{-1, 0, 2, -2, 1, -1987.0, 0.0, -6.0, 1073.0, 0.0, -2.0},
		{1, 0, 0, 0, 2, -1981.0, 0.0, 0.0, 854.0, 0.0, 0.0},
		{-1, 0, 0, 1, 0, 4026.0, 0.0, -353.0, -553.0, 0.0, -139.0},
		{0, 0, 2, 1, 2, 1660.0, 0.0, -5.0, -710.0, 0.0, -2.0},

		/* 71-77 */
		{-1, 0, 2, 4, 2, -1521.0, 0.0, 9.0, 647.0, 0.0, 4.0},
		{-1, 1, 0, 1, 1, 1314.0, 0.0, 0.0, -700.0, 0.0, 0.0},
		{0, -2, 2, -2, 1, -1283.0, 0.0, 0.0, 672.0, 0.0, 0.0},
		{1, 0, 2, 2, 1, -1331.0, 0.0, 8.0, 663.0, 0.0, 4.0},
		{-2, 0, 2, 2, 2, 1383.0, 0.0, -2.0, -594.0, 0.0, -2.0},
		{-1, 0, 0, 0, 2, 1405.0, 0.0, 4.0, -610.0, 0.0, 2.0},
		{1, 1, 2, -2, 2, 1290.0, 0.0, 0.0, -556.0, 0.0, 0.0},
	}

	/* Number of terms in the series */
	NLS := len(x)

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental epoch J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* --------------------*/
	/* LUNI-SOLAR NUTATION */
	/* --------------------*/

	/* Fundamental (Delaunay) arguments from Simon et al. (1994) */

	/* Mean anomaly of the Moon. */
	el = fmod(485868.249036+(1717915923.2178)*t, TURNAS) * DAS2R

	/* Mean anomaly of the Sun. */
	elp = fmod(1287104.79305+(129596581.0481)*t, TURNAS) * DAS2R

	/* Mean argument of the latitude of the Moon. */
	f = fmod(335779.526232+(1739527262.8478)*t, TURNAS) * DAS2R

	/* Mean elongation of the Moon from the Sun. */
	d = fmod(1072260.70369+(1602961601.2090)*t, TURNAS) * DAS2R

	/* Mean longitude of the ascending node of the Moon. */
	om = fmod(450160.398036+(-6962890.5431)*t, TURNAS) * DAS2R

	/* Initialize the nutation values. */
	dp = 0.0
	de = 0.0

	/* Summation of luni-solar nutation series (smallest terms first). */
	for i = NLS - 1; i >= 0; i-- {

		/* Argument and functions. */
		arg = fmod(float64(x[i].nl)*el+
			float64(x[i].nlp)*elp+
			float64(x[i].nf)*f+
			float64(x[i].nd)*d+
			float64(x[i].nom)*om, D2PI)
		sarg = sin(arg)
		carg = cos(arg)

		/* Term. */
		dp += (x[i].ps+x[i].pst*t)*sarg + x[i].pc*carg
		de += (x[i].ec+x[i].ect*t)*carg + x[i].es*sarg
	}

	/* Convert from 0.1 microarcsec units to radians. */
	dpsils = dp * U2R
	depsls = de * U2R

	/* ------------------------------*/
	/* IN LIEU OF PLANETARY NUTATION */
	/* ------------------------------*/

	/* Fixed offset to correct for missing terms in truncated series. */
	dpsipl = DPPLAN
	depspl = DEPLAN

	/* --------*/
	/* RESULTS */
	/* --------*/

	/* Add luni-solar and planetary components. */
	*dpsi = dpsils + dpsipl
	*deps = depsls + depspl
}

/*
Nut06a Nutation, IAU 2006/2000A
IAU 2000A nutation with adjustments to match the IAU 2006
precession.

Given:
    date1,date2   float64   TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps     float64   nutation, luni-solar + planetary (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components in longitude and obliquity are in radians
    and with respect to the mean equinox and ecliptic of date,
    2006 precession model (Hilton et al. 2006, Capitaine et al.
    2005).

 3) The function first computes the IAU 2000A nutation, then applies
    adjustments for (i) the consequences of the change in obliquity
    from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
    secular variation in the Earth's dynamical form factor J2.

 4) The present function provides classical nutation, complementing
    the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
    pole which is at current epochs accurate to a few tens of
    microarcseconds, apart from the free core nutation.

Called:
    Nut00a    nutation, IAU 2000A

References:

    Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
    Astron.Astrophys. 387, 700

    Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    Astron.Astrophys. 58, 1-16

    Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
    107, B4.  The MHB_2000 code itself was obtained on 9th September
    2002 from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

    Wallace, P.T., "Software for Implementing the IAU 2000
    Resolutions", in IERS Workshop 5.1 (2002)
*/
func Nut06a(date1, date2 float64, dpsi, deps *float64) {
	var t, fj2, dp, de float64

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Factor correcting for secular variation of J2. */
	fj2 = -2.7774e-6 * t

	/* Obtain IAU 2000A nutation. */
	Nut00a(date1, date2, &dp, &de)

	/* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
	*dpsi = dp + dp*(0.4697e-6+fj2)
	*deps = de + de*fj2

}

/*
Nut80 Nutation, IAU 1980 Nutation, IAU 1980 model.

Given:
    date1,date2   float64    TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi          float64    nutation in longitude (radians)
    deps          float64    nutation in obliquity (radians)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components are with respect to the ecliptic of
    date.

Called:
    Anpm      normalize angle into range +/- pi

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.222 (p111).
*/
func Nut80(date1, date2 float64, dpsi, deps *float64) {
	var t, el, elp, f, d, om, dp, de, arg, s, c float64
	var j int

	/* Units of 0.1 milliarcsecond to radians */
	const U2R = DAS2R / 1e4

	/* ------------------------------------------------ */
	/* Table of multiples of arguments and coefficients */
	/* ------------------------------------------------ */

	/* The units for the sine and cosine coefficients are 0.1 mas and */
	/* the same per Julian century */

	x := []struct {
		nl, nlp, nf, nd, nom int     /* coefficients of l,l',F,D,Om */
		sp, spt              float64 /* longitude sine, 1 and t coefficients */
		ce, cet              float64 /* obliquity cosine, 1 and t coefficients */
	}{

		/* 1-10 */
		{0, 0, 0, 0, 1, -171996.0, -174.2, 92025.0, 8.9},
		{0, 0, 0, 0, 2, 2062.0, 0.2, -895.0, 0.5},
		{-2, 0, 2, 0, 1, 46.0, 0.0, -24.0, 0.0},
		{2, 0, -2, 0, 0, 11.0, 0.0, 0.0, 0.0},
		{-2, 0, 2, 0, 2, -3.0, 0.0, 1.0, 0.0},
		{1, -1, 0, -1, 0, -3.0, 0.0, 0.0, 0.0},
		{0, -2, 2, -2, 1, -2.0, 0.0, 1.0, 0.0},
		{2, 0, -2, 0, 1, 1.0, 0.0, 0.0, 0.0},
		{0, 0, 2, -2, 2, -13187.0, -1.6, 5736.0, -3.1},
		{0, 1, 0, 0, 0, 1426.0, -3.4, 54.0, -0.1},

		/* 11-20 */
		{0, 1, 2, -2, 2, -517.0, 1.2, 224.0, -0.6},
		{0, -1, 2, -2, 2, 217.0, -0.5, -95.0, 0.3},
		{0, 0, 2, -2, 1, 129.0, 0.1, -70.0, 0.0},
		{2, 0, 0, -2, 0, 48.0, 0.0, 1.0, 0.0},
		{0, 0, 2, -2, 0, -22.0, 0.0, 0.0, 0.0},
		{0, 2, 0, 0, 0, 17.0, -0.1, 0.0, 0.0},
		{0, 1, 0, 0, 1, -15.0, 0.0, 9.0, 0.0},
		{0, 2, 2, -2, 2, -16.0, 0.1, 7.0, 0.0},
		{0, -1, 0, 0, 1, -12.0, 0.0, 6.0, 0.0},
		{-2, 0, 0, 2, 1, -6.0, 0.0, 3.0, 0.0},

		/* 21-30 */
		{0, -1, 2, -2, 1, -5.0, 0.0, 3.0, 0.0},
		{2, 0, 0, -2, 1, 4.0, 0.0, -2.0, 0.0},
		{0, 1, 2, -2, 1, 4.0, 0.0, -2.0, 0.0},
		{1, 0, 0, -1, 0, -4.0, 0.0, 0.0, 0.0},
		{2, 1, 0, -2, 0, 1.0, 0.0, 0.0, 0.0},
		{0, 0, -2, 2, 1, 1.0, 0.0, 0.0, 0.0},
		{0, 1, -2, 2, 0, -1.0, 0.0, 0.0, 0.0},
		{0, 1, 0, 0, 2, 1.0, 0.0, 0.0, 0.0},
		{-1, 0, 0, 1, 1, 1.0, 0.0, 0.0, 0.0},
		{0, 1, 2, -2, 0, -1.0, 0.0, 0.0, 0.0},

		/* 31-40 */
		{0, 0, 2, 0, 2, -2274.0, -0.2, 977.0, -0.5},
		{1, 0, 0, 0, 0, 712.0, 0.1, -7.0, 0.0},
		{0, 0, 2, 0, 1, -386.0, -0.4, 200.0, 0.0},
		{1, 0, 2, 0, 2, -301.0, 0.0, 129.0, -0.1},
		{1, 0, 0, -2, 0, -158.0, 0.0, -1.0, 0.0},
		{-1, 0, 2, 0, 2, 123.0, 0.0, -53.0, 0.0},
		{0, 0, 0, 2, 0, 63.0, 0.0, -2.0, 0.0},
		{1, 0, 0, 0, 1, 63.0, 0.1, -33.0, 0.0},
		{-1, 0, 0, 0, 1, -58.0, -0.1, 32.0, 0.0},
		{-1, 0, 2, 2, 2, -59.0, 0.0, 26.0, 0.0},

		/* 41-50 */
		{1, 0, 2, 0, 1, -51.0, 0.0, 27.0, 0.0},
		{0, 0, 2, 2, 2, -38.0, 0.0, 16.0, 0.0},
		{2, 0, 0, 0, 0, 29.0, 0.0, -1.0, 0.0},
		{1, 0, 2, -2, 2, 29.0, 0.0, -12.0, 0.0},
		{2, 0, 2, 0, 2, -31.0, 0.0, 13.0, 0.0},
		{0, 0, 2, 0, 0, 26.0, 0.0, -1.0, 0.0},
		{-1, 0, 2, 0, 1, 21.0, 0.0, -10.0, 0.0},
		{-1, 0, 0, 2, 1, 16.0, 0.0, -8.0, 0.0},
		{1, 0, 0, -2, 1, -13.0, 0.0, 7.0, 0.0},
		{-1, 0, 2, 2, 1, -10.0, 0.0, 5.0, 0.0},

		/* 51-60 */
		{1, 1, 0, -2, 0, -7.0, 0.0, 0.0, 0.0},
		{0, 1, 2, 0, 2, 7.0, 0.0, -3.0, 0.0},
		{0, -1, 2, 0, 2, -7.0, 0.0, 3.0, 0.0},
		{1, 0, 2, 2, 2, -8.0, 0.0, 3.0, 0.0},
		{1, 0, 0, 2, 0, 6.0, 0.0, 0.0, 0.0},
		{2, 0, 2, -2, 2, 6.0, 0.0, -3.0, 0.0},
		{0, 0, 0, 2, 1, -6.0, 0.0, 3.0, 0.0},
		{0, 0, 2, 2, 1, -7.0, 0.0, 3.0, 0.0},
		{1, 0, 2, -2, 1, 6.0, 0.0, -3.0, 0.0},
		{0, 0, 0, -2, 1, -5.0, 0.0, 3.0, 0.0},

		/* 61-70 */
		{1, -1, 0, 0, 0, 5.0, 0.0, 0.0, 0.0},
		{2, 0, 2, 0, 1, -5.0, 0.0, 3.0, 0.0},
		{0, 1, 0, -2, 0, -4.0, 0.0, 0.0, 0.0},
		{1, 0, -2, 0, 0, 4.0, 0.0, 0.0, 0.0},
		{0, 0, 0, 1, 0, -4.0, 0.0, 0.0, 0.0},
		{1, 1, 0, 0, 0, -3.0, 0.0, 0.0, 0.0},
		{1, 0, 2, 0, 0, 3.0, 0.0, 0.0, 0.0},
		{1, -1, 2, 0, 2, -3.0, 0.0, 1.0, 0.0},
		{-1, -1, 2, 2, 2, -3.0, 0.0, 1.0, 0.0},
		{-2, 0, 0, 0, 1, -2.0, 0.0, 1.0, 0.0},

		/* 71-80 */
		{3, 0, 2, 0, 2, -3.0, 0.0, 1.0, 0.0},
		{0, -1, 2, 2, 2, -3.0, 0.0, 1.0, 0.0},
		{1, 1, 2, 0, 2, 2.0, 0.0, -1.0, 0.0},
		{-1, 0, 2, -2, 1, -2.0, 0.0, 1.0, 0.0},
		{2, 0, 0, 0, 1, 2.0, 0.0, -1.0, 0.0},
		{1, 0, 0, 0, 2, -2.0, 0.0, 1.0, 0.0},
		{3, 0, 0, 0, 0, 2.0, 0.0, 0.0, 0.0},
		{0, 0, 2, 1, 2, 2.0, 0.0, -1.0, 0.0},
		{-1, 0, 0, 0, 2, 1.0, 0.0, -1.0, 0.0},
		{1, 0, 0, -4, 0, -1.0, 0.0, 0.0, 0.0},

		/* 81-90 */
		{-2, 0, 2, 2, 2, 1.0, 0.0, -1.0, 0.0},
		{-1, 0, 2, 4, 2, -2.0, 0.0, 1.0, 0.0},
		{2, 0, 0, -4, 0, -1.0, 0.0, 0.0, 0.0},
		{1, 1, 2, -2, 2, 1.0, 0.0, -1.0, 0.0},
		{1, 0, 2, 2, 1, -1.0, 0.0, 1.0, 0.0},
		{-2, 0, 2, 4, 2, -1.0, 0.0, 1.0, 0.0},
		{-1, 0, 4, 0, 2, 1.0, 0.0, 0.0, 0.0},
		{1, -1, 0, -2, 0, 1.0, 0.0, 0.0, 0.0},
		{2, 0, 2, -2, 1, 1.0, 0.0, -1.0, 0.0},
		{2, 0, 2, 2, 2, -1.0, 0.0, 0.0, 0.0},

		/* 91-100 */
		{1, 0, 0, 2, 1, -1.0, 0.0, 0.0, 0.0},
		{0, 0, 4, -2, 2, 1.0, 0.0, 0.0, 0.0},
		{3, 0, 2, -2, 2, 1.0, 0.0, 0.0, 0.0},
		{1, 0, 2, -2, 0, -1.0, 0.0, 0.0, 0.0},
		{0, 1, 2, 0, 1, 1.0, 0.0, 0.0, 0.0},
		{-1, -1, 0, 2, 1, 1.0, 0.0, 0.0, 0.0},
		{0, 0, -2, 0, 1, -1.0, 0.0, 0.0, 0.0},
		{0, 0, 2, -1, 2, -1.0, 0.0, 0.0, 0.0},
		{0, 1, 0, 2, 0, -1.0, 0.0, 0.0, 0.0},
		{1, 0, -2, -2, 0, -1.0, 0.0, 0.0, 0.0},

		/* 101-106 */
		{0, -1, 2, 0, 1, -1.0, 0.0, 0.0, 0.0},
		{1, 1, 0, -2, 1, -1.0, 0.0, 0.0, 0.0},
		{1, 0, -2, 2, 0, -1.0, 0.0, 0.0, 0.0},
		{2, 0, 0, 2, 0, 1.0, 0.0, 0.0, 0.0},
		{0, 0, 2, 4, 2, -1.0, 0.0, 0.0, 0.0},
		{0, 1, 0, 1, 0, 1.0, 0.0, 0.0, 0.0},
	}

	/* Number of terms in the series */
	NT := len(x)

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental epoch J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* --------------------- */
	/* Fundamental arguments */
	/* --------------------- */

	/* Mean longitude of Moon minus mean longitude of Moon's perigee. */
	el = Anpm((485866.733+(715922.633+(31.310+0.064*t)*t)*t)*DAS2R + fmod(1325.0*t, 1.0)*D2PI)

	/* Mean longitude of Sun minus mean longitude of Sun's perigee. */
	elp = Anpm((1287099.804+(1292581.224+(-0.577-0.012*t)*t)*t)*DAS2R + fmod(99.0*t, 1.0)*D2PI)

	/* Mean longitude of Moon minus mean longitude of Moon's node. */
	f = Anpm((335778.877+(295263.137+(-13.257+0.011*t)*t)*t)*DAS2R + fmod(1342.0*t, 1.0)*D2PI)

	/* Mean elongation of Moon from Sun. */
	d = Anpm((1072261.307+(1105601.328+(-6.891+0.019*t)*t)*t)*DAS2R + fmod(1236.0*t, 1.0)*D2PI)

	/* Longitude of the mean ascending node of the lunar orbit on the */
	/* ecliptic, measured from the mean equinox of date. */
	om = Anpm((450160.280+(-482890.539+(7.455+0.008*t)*t)*t)*DAS2R + fmod(-5.0*t, 1.0)*D2PI)

	/* --------------- */
	/* Nutation series */
	/* --------------- */

	/* Initialize nutation components. */
	dp = 0.0
	de = 0.0

	/* Sum the nutation terms, ending with the biggest. */
	for j = NT - 1; j >= 0; j-- {

		/* Form argument for current term. */
		arg = float64(x[j].nl)*el +
			float64(x[j].nlp)*elp +
			float64(x[j].nf)*f +
			float64(x[j].nd)*d +
			float64(x[j].nom)*om

		/* Accumulate current nutation term. */
		s = x[j].sp + x[j].spt*t
		c = x[j].ce + x[j].cet*t
		if s != 0.0 {
			dp += s * sin(arg)
		}
		if c != 0.0 {
			de += c * cos(arg)
		}
	}

	/* Convert results from 0.1 mas units to radians. */
	*dpsi = dp * U2R
	*deps = de * U2R
}

/*
Nutm80 Nutation matrix, IAU 1980
Form the matrix of nutation for a given date, IAU 1980 model.

Given:
    date1,date2    float64          TDB date (Note 1)

Returned:
    rmatn          [3][3]float64    nutation matrix

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(true) = rmatn * V(mean),
    where the p-vector V(true) is with respect to the true
    equatorial triad of date and the p-vector V(mean) is with
    respect to the mean equatorial triad of date.

Called:
    Nut80     nutation, IAU 1980
    Obl80     mean obliquity, IAU 1980
    Numat     form nutation matrix
*/
func Nutm80(date1, date2 float64, rmatn *[3][3]float64) {
	var dpsi, deps, epsa float64

	/* Nutation components and mean obliquity. */
	Nut80(date1, date2, &dpsi, &deps)
	epsa = Obl80(date1, date2)

	/* Build the rotation matrix. */
	Numat(epsa, dpsi, deps, rmatn)
}

/*
Obl06 Mean obliquity, IAU 2006
Mean obliquity of the ecliptic, IAU 2006 precession model.

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64   obliquity of the ecliptic (radians, Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The result is the angle between the ecliptic and mean equator of
    date date1+date2.

Reference:

    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
*/
func Obl06(date1, date2 float64) float64 {
	var t, eps0 float64

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Mean obliquity. */
	eps0 = (84381.406 + (-46.836769+(-0.0001831+(0.00200340+(-0.000000576+(-0.0000000434)*t)*t)*t)*t)*t) * DAS2R

	return eps0
}

/*
Obl80 Mean obliquity, IAU 1980 Mean obliquity of the ecliptic, IAU 1980 model.

Given:
    date1,date2   float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    obliquity of the ecliptic (radians, Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The result is the angle between the ecliptic and mean equator of
    date date1+date2.

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Expression 3.222-1 (p114).
*/
func Obl80(date1, date2 float64) float64 {
	var t, eps0 float64

	/* Interval between fundamental epoch J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Mean obliquity of date. */
	eps0 = DAS2R * (84381.448 +
		(-46.8150+
			(-0.00059+
				(0.001813)*t)*t)*t)

	return eps0
}

/*
Pb06 Zeta, z, theta precession angles, IAU 2006, including bias
This function forms three Euler angles which implement general
precession from epoch J2000.0, using the IAU 2006 model.  Frame
bias (the offset between ICRS and mean J2000.0) is included.

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned:
    bzeta        float64   1st rotation: radians cw around z
    bz           float64   3rd rotation: radians cw around z
    btheta       float64   2nd rotation: radians ccw around y

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The traditional accumulated precession angles zeta_A, z_A,
    theta_A cannot be obtained in the usual way, namely through
    polynomial expressions, because of the frame bias.  The latter
    means that two of the angles undergo rapid changes near this
    date.  They are instead the results of decomposing the
    precession-bias matrix obtained by using the Fukushima-Williams
    method, which does not suffer from the problem.  The
    decomposition returns values which can be used in the
    conventional formulation and which include frame bias.

 3) The three angles are returned in the conventional order, which
    is not the same as the order of the corresponding Euler
    rotations.  The precession-bias matrix is
    R_3(-z) x R_2(+theta) x R_3(-zeta).

 4) Should zeta_A, z_A, theta_A angles be required that do not
    contain frame bias, they are available by calling the SOFA
    function iauP06e.

Called:
    Pmat06    PB matrix, IAU 2006
    Rz        rotate around Z-axis
*/
func Pb06(date1, date2 float64, bzeta, bz, btheta *float64) {
	var r [3][3]float64
	var y, x float64

	/* Precession matrix via Fukushima-Williams angles. */
	Pmat06(date1, date2, &r)

	/* Solve for z, choosing the +/- pi alternative. */
	y = r[1][2]
	x = -r[0][2]
	if x < 0.0 {
		y = -y
		x = -x
	}

	// *bz = ( x != 0.0 || y != 0.0 ) ? - atan2(y,x) : 0.0;
	if x != 0.0 || y != 0.0 {
		*bz = -atan2(y, x)
	} else {
		*bz = 0.0
	}

	/* Derotate it out of the matrix. */
	Rz(*bz, &r)

	/* Solve for the remaining two angles. */
	y = r[0][2]
	x = r[2][2]
	// *btheta = ( x != 0.0 || y != 0.0 ) ? - atan2(y,x) : 0.0;
	if x != 0.0 || y != 0.0 {
		*btheta = -atan2(y, x)
	} else {
		*btheta = 0.0
	}

	y = -r[1][0]
	x = r[1][1]
	// *bzeta = ( x != 0.0 || y != 0.0 ) ? - atan2(y,x) : 0.0;
	if x != 0.0 || y != 0.0 {
		*bzeta = -atan2(y, x)
	} else {
		*bzeta = 0.0
	}
}

/*
Pfw06 Bias-precession Fukushima-Williams angles, IAU 2006
Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned:
    gamb         float64   F-W angle gamma_bar (radians)
    phib         float64   F-W angle phi_bar (radians)
    psib         float64   F-W angle psi_bar (radians)
    epsa         float64   F-W angle epsilon_A (radians)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole,
         E = mean ecliptic pole of date,
    and   P = mean pole of date,

    the four Fukushima-Williams angles are as follows:

      gamb = gamma_bar = epE
      phib = phi_bar = pE
      psib = psi_bar = pEP
      epsa = epsilon_A = EP

 3) The matrix representing the combined effects of frame bias and
    precession is:

      PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)

 4) The matrix representing the combined effects of frame bias,
    precession and nutation is simply:

      NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)

    where dP and dE are the nutation components with respect to the
    ecliptic of date.

Reference:

    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351

Called:
    Obl06     mean obliquity, IAU 2006
*/
func Pfw06(date1, date2 float64, gamb, phib, psib, epsa *float64) {
	var t float64

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* P03 bias+precession angles. */
	*gamb = (-0.052928 + (10.556378+(0.4932044+(-0.00031238+(-0.000002788+(0.0000000260)*t)*t)*t)*t)*t) * DAS2R
	*phib = (84381.412819 + (-46.811016+(0.0511268+(0.00053289+(-0.000000440+(-0.0000000176)*t)*t)*t)*t)*t) * DAS2R
	*psib = (-0.041775 + (5038.481484+(1.5584175+(-0.00018522+(-0.000026452+(-0.0000000148)*t)*t)*t)*t)*t) * DAS2R
	*epsa = Obl06(date1, date2)
}

/*
Pmat00 Precession matrix (including frame bias), IAU 2000
Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2000 model.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    rbp          [3][3]float64    bias-precession matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
    the p-vector V(GCRS) is with respect to the Geocentric Celestial
    Reference System (IAU, 2000) and the p-vector V(date) is with
    respect to the mean equatorial triad of the given date.

Called:
    Bp00      frame bias and precession matrices, IAU 2000

Reference:

    : Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    (2000)
*/
func Pmat00(date1, date2 float64, rbp *[3][3]float64) {
	var rb, rp [3][3]float64

	/* Obtain the required matrix (discarding others). */
	Bp00(date1, date2, &rb, &rp, rbp)
}

/*
Pmat06 Precession matrix (including frame bias), IAU 2006
Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2006 model.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    rbp          [3][3]float64    bias-precession matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
    the p-vector V(GCRS) is with respect to the Geocentric Celestial
    Reference System (IAU, 2000) and the p-vector V(date) is with
    respect to the mean equatorial triad of the given date.

Called:
    Pfw06     bias-precession F-W angles, IAU 2006
    Fw2m      F-W angles to r-matrix

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    : Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    (2000)

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Pmat06(date1, date2 float64, rbp *[3][3]float64) {
	var gamb, phib, psib, epsa float64

	/* Bias-precession Fukushima-Williams angles. */
	Pfw06(date1, date2, &gamb, &phib, &psib, &epsa)

	/* Form the matrix. */
	Fw2m(gamb, phib, psib, epsa, rbp)
}

/*
Pmat76 Precession matrix, IAU 1976
Precession matrix from J2000.0 to a specified date, IAU 1976 model.

Given:
    date1,date2 float64       ending date, TT (Note 1)

Returned:
    rmatp       [3][3]float64 precession matrix, J2000.0 -> date1+date2

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = RMATP * V(J2000),
    where the p-vector V(J2000) is with respect to the mean
    equatorial triad of epoch J2000.0 and the p-vector V(date)
    is with respect to the mean equatorial triad of the given
    date.

 3) Though the matrix method itself is rigorous, the precession
    angles are expressed through canonical polynomials which are
    valid only for a limited time span.  In addition, the IAU 1976
    precession rate is known to be imperfect.  The absolute accuracy
    of the present formulation is better than 0.1 arcsec from
    1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
    and remains below 3 arcsec for the whole of the period
    500BC to 3000AD.  The errors exceed 10 arcsec outside the
    range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
    5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.

Called:
    Prec76    accumulated precession angles, IAU 1976
    Ir        initialize r-matrix to identity
    Rz        rotate around Z-axis
    Ry        rotate around Y-axis
    Cr        copy r-matrix

References:

    Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
    equations (6) & (7), p283.

    Kaplan,G.H., 1981. USNO circular no. 163, pA2.
*/
func Pmat76(date1, date2 float64, rmatp *[3][3]float64) {
	var zeta, z, theta float64
	var wmat [3][3]float64

	/* Precession Euler angles, J2000.0 to specified date. */
	Prec76(DJ00, 0.0, date1, date2, &zeta, &z, &theta)

	/* Form the rotation matrix. */
	Ir(&wmat)
	Rz(-zeta, &wmat)
	Ry(theta, &wmat)
	Rz(-z, &wmat)
	Cr(wmat, rmatp)
}

/*
Pn00 B,P,N matrices, IAU 2000, given nutation

Precession-nutation, IAU 2000 model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)
    dpsi,deps    float64          nutation (Note 2)

Returned:
    epsa         float64          mean obliquity (Note 3)
    rb           [3][3]float64    frame bias matrix (Note 4)
    rp           [3][3]float64    precession matrix (Note 5)
    rbp          [3][3]float64    bias-precession matrix (Note 6)
    rn           [3][3]float64    nutation matrix (Note 7)
    rbpn         [3][3]float64    GCRS-to-true matrix (Note 8)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The caller is responsible for providing the nutation components;
    they are in longitude and obliquity, in radians and are with
    respect to the equinox and ecliptic of date.  For high-accuracy
    applications, free core nutation should be included as well as
    any other relevant corrections to the position of the CIP.

 3) The returned mean obliquity is consistent with the IAU 2000
    precession-nutation models.

 4) The matrix rb transforms vectors from GCRS to J2000.0 mean
    equator and equinox by applying frame bias.

 5) The matrix rp transforms vectors from J2000.0 mean equator and
    equinox to mean equator and equinox of date by applying
    precession.

 6) The matrix rbp transforms vectors from GCRS to mean equator and
    equinox of date by applying frame bias then precession.  It is
    the product rp x rb.

 7) The matrix rn transforms vectors from mean equator and equinox of
    date to true equator and equinox of date by applying the nutation
    (luni-solar + planetary).

 8) The matrix rbpn transforms vectors from GCRS to true equator and
    equinox of date.  It is the product rn x rbp, applying frame
    bias, precession and nutation in that order.

 9) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the order given.

Called:
    Pr00      IAU 2000 precession adjustments
    Obl80     mean obliquity, IAU 1980
    Bp00      frame bias and precession matrices, IAU 2000
    Cr        copy r-matrix
    Numat     form nutation matrix
    Rxr       product of two r-matrices

Reference:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.
*/
func Pn00(date1, date2 float64, dpsi, deps float64, epsa *float64, rb, rp, rbp, rn, rbpn *[3][3]float64) {
	var dpsipr, depspr float64
	var rbpw, rnw [3][3]float64

	/* IAU 2000 precession-rate adjustments. */
	Pr00(date1, date2, &dpsipr, &depspr)

	/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
	*epsa = Obl80(date1, date2) + depspr

	/* Frame bias and precession matrices and their product. */
	Bp00(date1, date2, rb, rp, &rbpw)
	Cr(rbpw, rbp)

	/* Nutation matrix. */
	Numat(*epsa, dpsi, deps, &rnw)
	Cr(rnw, rn)

	/* Bias-precession-nutation matrix (classical). */
	Rxr(rnw, rbpw, rbpn)
}

/*
Pn00a B,P,N matrices, IAU 2000A

Precession-nutation, IAU 2000A model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps    float64          nutation (Note 2)
    epsa         float64          mean obliquity (Note 3)
    rb           [3][3]float64    frame bias matrix (Note 4)
    rp           [3][3]float64    precession matrix (Note 5)
    rbp          [3][3]float64    bias-precession matrix (Note 6)
    rn           [3][3]float64    nutation matrix (Note 7)
    rbpn         [3][3]float64    GCRS-to-true matrix (Notes 8,9)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components (luni-solar + planetary, IAU 2000A) in
    longitude and obliquity are in radians and with respect to the
    equinox and ecliptic of date.  Free core nutation is omitted;
    for the utmost accuracy, use the iauPn00 function, where the
    nutation components are caller-specified.  For faster but
    slightly less accurate results, use the iauPn00b function.

 3) The mean obliquity is consistent with the IAU 2000 precession.

 4) The matrix rb transforms vectors from GCRS to J2000.0 mean
    equator and equinox by applying frame bias.

 5) The matrix rp transforms vectors from J2000.0 mean equator and
    equinox to mean equator and equinox of date by applying
    precession.

 6) The matrix rbp transforms vectors from GCRS to mean equator and
    equinox of date by applying frame bias then precession.  It is
    the product rp x rb.

 7) The matrix rn transforms vectors from mean equator and equinox
    of date to true equator and equinox of date by applying the
    nutation (luni-solar + planetary).

 8) The matrix rbpn transforms vectors from GCRS to true equator and
    equinox of date.  It is the product rn x rbp, applying frame
    bias, precession and nutation in that order.

 9) The X,Y,Z coordinates of the IAU 2000A Celestial Intermediate
    Pole are elements (3,1-3) of the GCRS-to-true matrix,
    i.e. rbpn[2][0-2].

 10)It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
    Nut00a    nutation, IAU 2000A
    Pn00      bias/precession/nutation results, IAU 2000

Reference:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.
*/
func Pn00a(date1, date2 float64, dpsi, deps, epsa *float64, rb, rp, rbp, rn, rbpn *[3][3]float64) {
	/* Nutation. */
	Nut00a(date1, date2, dpsi, deps)

	/* Remaining results. */
	Pn00(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn)
}

/*
Pn00b B,P,N matrices, IAU 2000B

Precession-nutation, IAU 2000B model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps    float64          nutation (Note 2)
    epsa         float64          mean obliquity (Note 3)
    rb           [3][3]float64    frame bias matrix (Note 4)
    rp           [3][3]float64    precession matrix (Note 5)
    rbp          [3][3]float64    bias-precession matrix (Note 6)
    rn           [3][3]float64    nutation matrix (Note 7)
    rbpn         [3][3]float64    GCRS-to-true matrix (Notes 8,9)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

       date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components (luni-solar + planetary, IAU 2000B) in
    longitude and obliquity are in radians and with respect to the
    equinox and ecliptic of date.  For more accurate results, but
    at the cost of increased computation, use the iauPn00a function.
    For the utmost accuracy, use the iauPn00 function, where the
    nutation components are caller-specified.

 3) The mean obliquity is consistent with the IAU 2000 precession.

 4) The matrix rb transforms vectors from GCRS to J2000.0 mean
    equator and equinox by applying frame bias.

 5) The matrix rp transforms vectors from J2000.0 mean equator and
    equinox to mean equator and equinox of date by applying
    precession.

 6) The matrix rbp transforms vectors from GCRS to mean equator and
    equinox of date by applying frame bias then precession.  It is
    the product rp x rb.

 7) The matrix rn transforms vectors from mean equator and equinox
    of date to true equator and equinox of date by applying the
    nutation (luni-solar + planetary).

 8) The matrix rbpn transforms vectors from GCRS to true equator and
    equinox of date.  It is the product rn x rbp, applying frame
    bias, precession and nutation in that order.

 9) The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
    Pole are elements (3,1-3) of the GCRS-to-true matrix,
    i.e. rbpn[2][0-2].

 10)It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
    Nut00b    nutation, IAU 2000B
    Pn00      bias/precession/nutation results, IAU 2000

Reference:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003).

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.
*/
func Pn00b(date1, date2 float64, dpsi, deps, epsa *float64, rb, rp, rbp, rn, rbpn *[3][3]float64) {
	/* Nutation. */
	Nut00b(date1, date2, dpsi, deps)

	/* Remaining results. */
	Pn00(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn)
}

/*
Pn06 Bias precession nutation results, IAU 2006

Precession-nutation, IAU 2006 model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)
    dpsi,deps    float64          nutation (Note 2)

Returned:
    epsa         float64          mean obliquity (Note 3)
    rb           [3][3]float64    frame bias matrix (Note 4)
    rp           [3][3]float64    precession matrix (Note 5)
    rbp          [3][3]float64    bias-precession matrix (Note 6)
    rn           [3][3]float64    nutation matrix (Note 7)
    rbpn         [3][3]float64    GCRS-to-true matrix (Notes 8,9)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

           date1          date2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The caller is responsible for providing the nutation components;
    they are in longitude and obliquity, in radians and are with
    respect to the equinox and ecliptic of date.  For high-accuracy
    applications, free core nutation should be included as well as
    any other relevant corrections to the position of the CIP.

 3) The returned mean obliquity is consistent with the IAU 2006
    precession.

 4) The matrix rb transforms vectors from GCRS to J2000.0 mean
    equator and equinox by applying frame bias.

 5) The matrix rp transforms vectors from J2000.0 mean equator and
    equinox to mean equator and equinox of date by applying
    precession.

 6) The matrix rbp transforms vectors from GCRS to mean equator and
    equinox of date by applying frame bias then precession.  It is
    the product rp x rb.

 7) The matrix rn transforms vectors from mean equator and equinox
    of date to true equator and equinox of date by applying the
    nutation (luni-solar + planetary).

 8) The matrix rbpn transforms vectors from GCRS to true equator and
    equinox of date.  It is the product rn x rbp, applying frame
    bias, precession and nutation in that order.

 9) The X,Y,Z coordinates of the Celestial Intermediate Pole are
    elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].

 10)It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
    Pfw06     bias-precession F-W angles, IAU 2006
    Fw2m      F-W angles to r-matrix
    Cr        copy r-matrix
    Tr        transpose r-matrix
    Rxr       product of two r-matrices

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Pn06(date1, date2 float64, dpsi, deps float64, epsa *float64, rb, rp, rbp, rn, rbpn *[3][3]float64) {
	var gamb, phib, psib, eps float64
	var r1, r2, rt [3][3]float64

	/* Bias-precession Fukushima-Williams angles of J2000.0 = frame bias. */
	Pfw06(DJM0, DJM00, &gamb, &phib, &psib, &eps)

	/* B matrix. */
	Fw2m(gamb, phib, psib, eps, &r1)
	Cr(r1, rb)

	/* Bias-precession Fukushima-Williams angles of date. */
	Pfw06(date1, date2, &gamb, &phib, &psib, &eps)

	/* Bias-precession matrix. */
	Fw2m(gamb, phib, psib, eps, &r2)
	Cr(r2, rbp)

	/* Solve for precession matrix. */
	Tr(r1, &rt)
	Rxr(r2, rt, rp)

	/* Equinox-based bias-precession-nutation matrix. */
	Fw2m(gamb, phib, psib+dpsi, eps+deps, &r1)
	Cr(r1, rbpn)

	/* Solve for nutation matrix. */
	Tr(r2, &rt)
	Rxr(r1, rt, rn)

	/* Obliquity, mean of date. */
	*epsa = eps

}

/*
Pn06a Bias precession nutation results, IAU 2006/2000A

Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

Given:
    date1,date2  float64          TT as a 2-part Julian Date (Note 1)

Returned:
    dpsi,deps    float64          nutation (Note 2)
    epsa         float64          mean obliquity (Note 3)
    rb           [3][3]float64    frame bias matrix (Note 4)
    rp           [3][3]float64    precession matrix (Note 5)
    rbp          [3][3]float64    bias-precession matrix (Note 6)
    rn           [3][3]float64    nutation matrix (Note 7)
    rbpn         [3][3]float64    GCRS-to-true matrix (Notes 8,9)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

           date1          date2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The nutation components (luni-solar + planetary, IAU 2000A) in
    longitude and obliquity are in radians and with respect to the
    equinox and ecliptic of date.  Free core nutation is omitted;
    for the utmost accuracy, use the iauPn06 function, where the
    nutation components are caller-specified.

 3) The mean obliquity is consistent with the IAU 2006 precession.

 4) The matrix rb transforms vectors from GCRS to mean J2000.0 by
    applying frame bias.

 5) The matrix rp transforms vectors from mean J2000.0 to mean of
    date by applying precession.

 6) The matrix rbp transforms vectors from GCRS to mean of date by
    applying frame bias then precession.  It is the product rp x rb.

 7) The matrix rn transforms vectors from mean of date to true of
    date by applying the nutation (luni-solar + planetary).

 8) The matrix rbpn transforms vectors from GCRS to true of date
    (CIP/equinox).  It is the product rn x rbp, applying frame bias,
    precession and nutation in that order.

 9) The X,Y,Z coordinates of the IAU 2006/2000A Celestial
    Intermediate Pole are elements (3,1-3) of the GCRS-to-true
    matrix, i.e. rbpn[2][0-2].

 10)It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
    Nut06a    nutation, IAU 2006/2000A
    Pn06      bias/precession/nutation results, IAU 2006

Reference:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
*/
func Pn06a(date1, date2 float64, dpsi, deps, epsa *float64, rb, rp, rbp, rn, rbpn *[3][3]float64) {
	/* Nutation. */
	Nut06a(date1, date2, dpsi, deps)

	/* Remaining results. */
	Pn06(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn)
}

/*
Pnm00a  Classical NPB matrix, IAU 2000A

Form the matrix of precession-nutation for a given date (including
frame bias), equinox based, IAU 2000A model.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rbpn        [3][3]float64 bias-precession-nutation matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways, among
    others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
    the p-vector V(date) is with respect to the true equatorial triad
    of date date1+date2 and the p-vector V(GCRS) is with respect to
    the Geocentric Celestial Reference System (IAU, 2000).

 3) A faster, but slightly less accurate, result (about 1 mas) can be
    obtained by using instead the iauPnm00b function.

Called:
    Pn00a     bias/precession/nutation, IAU 2000A

Reference:

    : Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    (2000)
*/
func Pnm00a(date1, date2 float64, rbpn *[3][3]float64) {
	var dpsi, deps, epsa float64
	var rb, rp, rbp, rn [3][3]float64

	/* Obtain the required matrix (discarding other results). */
	Pn00a(date1, date2, &dpsi, &deps, &epsa, &rb, &rp, &rbp, &rn, rbpn)
}

/*
Pnm00b Classical NPB matrix, IAU 2000B

Form the matrix of precession-nutation for a given date (including
frame bias), equinox-based, IAU 2000B model.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rbpn        [3][3]float64 bias-precession-nutation matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways, among
    others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
    the p-vector V(date) is with respect to the true equatorial triad
    of date date1+date2 and the p-vector V(GCRS) is with respect to
    the Geocentric Celestial Reference System (IAU, 2000).

 3) The present function is faster, but slightly less accurate (about
    1 mas), than the iauPnm00a function.

Called:
    Pn00b     bias/precession/nutation, IAU 2000B

Reference:

    : Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    (2000)
*/
func Pnm00b(date1, date2 float64, rbpn *[3][3]float64) {
	var dpsi, deps, epsa float64
	var rb, rp, rbp, rn [3][3]float64

	/* Obtain the required matrix (discarding other results). */
	Pn00b(date1, date2, &dpsi, &deps, &epsa, &rb, &rp, &rbp, &rn, rbpn)
}

/*
Pnm06a Classical NPB matrix, IAU 2006/2000A

Form the matrix of precession-nutation for a given date (including
frame bias), equinox based, IAU 2006 precession and IAU 2000A
nutation models.

Given:
    date1,date2 float64       TT as a 2-part Julian Date (Note 1)

Returned:
    rbpn        [3][3]float64 bias-precession-nutation matrix (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways, among
    others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
    the p-vector V(date) is with respect to the true equatorial triad
    of date date1+date2 and the p-vector V(GCRS) is with respect to
    the Geocentric Celestial Reference System (IAU, 2000).

Called:
    Pfw06     bias-precession F-W angles, IAU 2006
    Nut06a    nutation, IAU 2006/2000A
    Fw2m      F-W angles to r-matrix

Reference:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
*/
func Pnm06a(date1, date2 float64, rbpn *[3][3]float64) {
	var gamb, phib, psib, epsa, dp, de float64

	/* Fukushima-Williams angles for frame bias and precession. */
	Pfw06(date1, date2, &gamb, &phib, &psib, &epsa)

	/* Nutation components. */
	Nut06a(date1, date2, &dp, &de)

	/* Equinox based nutation x precession x bias matrix. */
	Fw2m(gamb, phib, psib+dp, epsa+de, rbpn)
}

/*
Pnm80 Precession/nutation matrix, IAU 1976/1980

Form the matrix of precession/nutation for a given date, IAU 1976
precession model, IAU 1980 nutation model.

Given:
    date1,date2     float64    TT as a 2-part Julian Date (Note 1)

Returned:
    rmatpn    [3][3]float64    combined precession/nutation matrix

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The matrix operates in the sense V(date) = rmatpn * V(J2000),
    where the p-vector V(date) is with respect to the true equatorial
    triad of date date1+date2 and the p-vector V(J2000) is with
    respect to the mean equatorial triad of epoch J2000.0.

Called:
    Pmat76    precession matrix, IAU 1976
    Nutm80    nutation matrix, IAU 1980
    Rxr       product of two r-matrices

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 3.3 (p145).
*/
func Pnm80(date1, date2 float64, rmatpn *[3][3]float64) {
	var rmatp, rmatn [3][3]float64

	/* Precession matrix, J2000.0 to date. */
	Pmat76(date1, date2, &rmatp)

	/* Nutation matrix. */
	Nutm80(date1, date2, &rmatn)

	/* Combine the matrices:  PN = N x P. */
	Rxr(rmatn, rmatp, rmatpn)
}

/*
P06e Precession angles, IAU 2006, equinox based

Given:
    date1,date2   float64   TT as a 2-part Julian Date (Note 1)

Returned (see Note 2):
    eps0          float64   epsilon_0
    psia          float64   psi_A
    oma           float64   omega_A
    bpa           float64   P_A
    bqa           float64   Q_A
    pia           float64   pi_A
    bpia          float64   Pi_A
    epsa          float64   obliquity epsilon_A
    chia          float64   chi_A
    za            float64   z_A
    zetaa         float64   zeta_A
    thetaa        float64   theta_A
    pa            float64   p_A
    gam           float64   F-W angle gamma_J2000
    phi           float64   F-W angle phi_J2000
    psi           float64   F-W angle psi_J2000

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) This function returns the set of equinox based angles for the
    Capitaine et al. "P03" precession theory, adopted by the IAU in
    2006.  The angles are set out in Table 1 of Hilton et al. (2006):

        eps0   epsilon_0   obliquity at J2000.0
        psia   psi_A       luni-solar precession
        oma    omega_A     inclination of equator wrt J2000.0 ecliptic
        bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad
        bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad
        pia    pi_A        angle between moving and J2000.0 ecliptics
        bpia   Pi_A        longitude of ascending node of the ecliptic
        epsa   epsilon_A   obliquity of the ecliptic
        chia   chi_A       planetary precession
        za     z_A         equatorial precession: -3rd 323 Euler angle
        zetaa  zeta_A      equatorial precession: -1st 323 Euler angle
        thetaa theta_A     equatorial precession: 2nd 323 Euler angle
        pa     p_A         general precession (n.b. see below)
        gam    gamma_J2000 J2000.0 RA difference of ecliptic poles
        phi    phi_J2000   J2000.0 codeclination of ecliptic pole
        psi    psi_J2000   longitude difference of equator poles, J2000.0

    The returned values are all radians.

    Note that the t^5 coefficient in the series for p_A from
    Capitaine et al. (2003) is incorrectly signed in Hilton et al.
    (2006).

 3) Hilton et al. (2006) Table 1 also contains angles that depend on
    models distinct from the P03 precession theory itself, namely the
    2000A frame bias and nutation.  The quoted polynomials are
    used in other SOFA functions:

    . iauXy06  contains the polynomial parts of the X and Y series.

    . iauS06  contains the polynomial part of the s+XY/2 series.

    . iauPfw06  implements the series for the Fukushima-Williams
      angles that are with respect to the GCRS pole (i.e. the variants
      that include frame bias).

 4) The IAU resolution stipulated that the choice of parameterization
    was left to the user, and so an IAU compliant precession
    implementation can be constructed using various combinations of
    the angles returned by the present function.

 5) The parameterization used by SOFA is the version of the Fukushima-
    Williams angles that refers directly to the GCRS pole.  These
    angles may be calculated by calling the function iauPfw06.  SOFA
    also supports the direct computation of the CIP GCRS X,Y by
    series, available by calling iauXy06.

 6) The agreement between the different parameterizations is at the
    1 microarcsecond level in the present era.

 7) When constructing a precession formulation that refers to the GCRS
    pole rather than the dynamical pole, it may (depending on the
    choice of angles) be necessary to introduce the frame bias
    explicitly.

 8) It is permissible to re-use the same variable in the returned
    arguments.  The quantities are stored in the stated order.

References:

    Capitaine, N., Wallace, P.T. & Chapront, J., 2003,
    Astron.Astrophys., 412, 567

    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351

Called:
    Obl06     mean obliquity, IAU 2006
*/
func P06e(date1, date2 float64, eps0, psia, oma, bpa *float64, bqa, pia, bpia *float64, epsa, chia, za, zetaa *float64,
	thetaa, pa *float64, gam, phi, psi *float64) {
	var t float64

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Obliquity at J2000.0. */

	*eps0 = 84381.406 * DAS2R

	/* Luni-solar precession. */

	*psia = (5038.481507 +
		(-1.0790069+
			(-0.00114045+
				(0.000132851+
					(-0.0000000951)*
						t)*t)*t)*t) * t * DAS2R

	/* Inclination of mean equator with respect to the J2000.0 ecliptic. */

	*oma = *eps0 + (-0.025754+
		(0.0512623+
			(-0.00772503+
				(-0.000000467+
					(0.0000003337)*
						t)*t)*t)*t)*t*DAS2R

	/* Ecliptic pole x, J2000.0 ecliptic triad. */

	*bpa = (4.199094 +
		(0.1939873+
			(-0.00022466+
				(-0.000000912+
					(0.0000000120)*
						t)*t)*t)*t) * t * DAS2R

	/* Ecliptic pole -y, J2000.0 ecliptic triad. */

	*bqa = (-46.811015 +
		(0.0510283+
			(0.00052413+
				(-0.000000646+
					(-0.0000000172)*
						t)*t)*t)*t) * t * DAS2R

	/* Angle between moving and J2000.0 ecliptics. */

	*pia = (46.998973 +
		(-0.0334926+
			(-0.00012559+
				(0.000000113+
					(-0.0000000022)*
						t)*t)*t)*t) * t * DAS2R

	/* Longitude of ascending node of the moving ecliptic. */

	*bpia = (629546.7936 +
		(-867.95758+
			(0.157992+
				(-0.0005371+
					(-0.00004797+
						(0.000000072)*
							t)*t)*t)*t)*t) * DAS2R

	/* Mean obliquity of the ecliptic. */

	*epsa = Obl06(date1, date2)

	/* Planetary precession. */

	*chia = (10.556403 +
		(-2.3814292+
			(-0.00121197+
				(0.000170663+
					(-0.0000000560)*
						t)*t)*t)*t) * t * DAS2R

	/* Equatorial precession: minus the third of the 323 Euler angles. */

	*za = (-2.650545 +
		(2306.077181+
			(1.0927348+
				(0.01826837+
					(-0.000028596+
						(-0.0000002904)*
							t)*t)*t)*t)*t) * DAS2R

	/* Equatorial precession: minus the first of the 323 Euler angles. */

	*zetaa = (2.650545 +
		(2306.083227+
			(0.2988499+
				(0.01801828+
					(-0.000005971+
						(-0.0000003173)*
							t)*t)*t)*t)*t) * DAS2R

	/* Equatorial precession: second of the 323 Euler angles. */

	*thetaa = (2004.191903 +
		(-0.4294934+
			(-0.04182264+
				(-0.000007089+
					(-0.0000001274)*
						t)*t)*t)*t) * t * DAS2R

	/* General precession. */

	*pa = (5028.796195 +
		(1.1054348+
			(0.00007964+
				(-0.000023857+
					(-0.0000000383)*
						t)*t)*t)*t) * t * DAS2R

	/* Fukushima-Williams angles for precession. */

	*gam = (10.556403 +
		(0.4932044+
			(-0.00031238+
				(-0.000002788+
					(0.0000000260)*
						t)*t)*t)*t) * t * DAS2R

	*phi = *eps0 + (-46.811015+
		(0.0511269+
			(0.00053289+
				(-0.000000440+
					(-0.0000000176)*
						t)*t)*t)*t)*t*DAS2R

	*psi = (5038.481507 +
		(1.5584176+
			(-0.00018522+
				(-0.000026452+
					(-0.0000000148)*
						t)*t)*t)*t) * t * DAS2R
}

/*
Pom00 Polar-motion matrix, IAU 2000

Form the matrix of polar motion for a given date, IAU 2000.

Given:
    xp,yp    float64    coordinates of the pole (radians, Note 1)
    sp       float64    the TIO locator s' (radians, Note 2)

Returned:
    rpom     [3][3]float64   polar-motion matrix (Note 3)

Notes:

 1) The arguments xp and yp are the coordinates (in radians) of the
    Celestial Intermediate Pole with respect to the International
    Terrestrial Reference System (see IERS Conventions 2003),
    measured along the meridians 0 and 90 deg west respectively.

 2) The argument sp is the TIO locator s', in radians, which
    positions the Terrestrial Intermediate Origin on the equator.  It
    is obtained from polar motion observations by numerical
    integration, and so is in essence unpredictable.  However, it is
    dominated by a secular drift of about 47 microarcseconds per
    century, and so can be taken into account by using s' = -47*t,
    where t is centuries since J2000.0.  The function iauSp00
    implements this approximation.

 3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
    that it is the final rotation when computing the pointing
    direction to a celestial source.

Called:
    Ir        initialize r-matrix to identity
    Rz        rotate around Z-axis
    Ry        rotate around Y-axis
    Rx        rotate around X-axis

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Pom00(xp, yp, sp float64, rpom *[3][3]float64) {
	/* Construct the matrix. */
	Ir(rpom)
	Rz(sp, rpom)
	Ry(-xp, rpom)
	Rx(-yp, rpom)
}

/*
Pr00 Adjustments to IAU 1976 precession, IAU 2000

Precession-rate part of the IAU 2000 precession-nutation models
(part of MHB2000).

Given:
    date1,date2    float64  TT as a 2-part Julian Date (Note 1)

Returned:
    dpsipr,depspr  float64  precession corrections (Notes 2,3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The precession adjustments are expressed as "nutation
    components", corrections in longitude and obliquity with respect
    to the J2000.0 equinox and ecliptic.

 3) Although the precession adjustments are stated to be with respect
    to Lieske et al. (1977), the MHB2000 model does not specify which
    set of Euler angles are to be used and how the adjustments are to
    be applied.  The most literal and straightforward procedure is to
    adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
    to add dpsipr to psi_A and depspr to both omega_A and eps_A.

 4) This is an implementation of one aspect of the IAU 2000A nutation
    model, formally adopted by the IAU General Assembly in 2000,
    namely MHB2000 (Mathews et al. 2002).

References:

    Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
    for the precession quantities based upon the IAU (1976) System of
    Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)

    Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
    and precession   New nutation series for nonrigid Earth and
    insights into the Earth's interior", J.Geophys.Res., 107, B4,
    2002.  The MHB2000 code itself was obtained on 9th September 2002
    from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

    Wallace, P.T., "Software for Implementing the IAU 2000
    Resolutions", in IERS Workshop 5.1 (2002).
*/
func Pr00(date1, date2 float64, dpsipr, depspr *float64) {
	var t float64

	/* Precession and obliquity corrections (radians per century) */
	const PRECOR = -0.29965 * DAS2R
	const OBLCOR = -0.02524 * DAS2R

	/* Interval between fundamental epoch J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Precession rate contributions with respect to IAU 1976/80. */
	*dpsipr = PRECOR * t
	*depspr = OBLCOR * t
}

/*
Prec76 Precession, IAU 1976

IAU 1976 precession model.

This function forms the three Euler angles which implement general
precession between two dates, using the IAU 1976 model (as for the
FK5 catalog).

Given:
    date01,date02   float64    TDB starting date (Note 1)
    date11,date12   float64    TDB ending date (Note 1)

Returned:
    zeta            float64    1st rotation: radians cw around z
    z               float64    3rd rotation: radians cw around z
    theta           float64    2nd rotation: radians ccw around y

Notes:

 1) The dates date01+date02 and date11+date12 are Julian Dates,
    apportioned in any convenient way between the arguments daten1
    and daten2.  For example, JD(TDB)=2450123.7 could be expressed in
    any of these ways, among others:

         daten1        daten2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.
    The two dates may be expressed using different methods, but at
    the risk of losing some resolution.

 2) The accumulated precession angles zeta, z, theta are expressed
    through canonical polynomials which are valid only for a limited
    time span.  In addition, the IAU 1976 precession rate is known to
    be imperfect.  The absolute accuracy of the present formulation
    is better than 0.1 arcsec from 1960AD to 2040AD, better than
    1 arcsec from 1640AD to 2360AD, and remains below 3 arcsec for
    the whole of the period 500BC to 3000AD.  The errors exceed
    10 arcsec outside the range 1200BC to 3900AD, exceed 100 arcsec
    outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC to
    8200AD.

 3) The three angles are returned in the conventional order, which
    is not the same as the order of the corresponding Euler
    rotations.  The precession matrix is
        R_3(-z) x R_2(+theta) x R_3(-zeta).

Reference:

    Lieske, J.H., 1979, Astron.Astrophys. 73, 282, equations
    (6) & (7), p283.
*/
func Prec76(date01, date02 float64, date11, date12 float64, zeta, z, theta *float64) {
	var t0, t, tas2r, w float64

	/* Interval between fundamental epoch J2000.0 and start date (JC). */
	t0 = ((date01 - DJ00) + date02) / DJC

	/* Interval over which precession required (JC). */
	t = ((date11 - date01) + (date12 - date02)) / DJC

	/* Euler angles. */
	tas2r = t * DAS2R
	w = 2306.2181 + (1.39656-0.000139*t0)*t0

	*zeta = (w + ((0.30188-0.000344*t0)+0.017998*t)*t) * tas2r

	*z = (w + ((1.09468+0.000066*t0)+0.018203*t)*t) * tas2r

	*theta = ((2004.3109 + (-0.85330-0.000217*t0)*t0) + ((-0.42665-0.000217*t0)-0.041833*t)*t) * tas2r
}

/*
S00 The CIO locator s, given X,Y, IAU 2000

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2000A precession-nutation.

Given:
    date1,date2   float64    TT as a 2-part Julian Date (Note 1)
    x,y           float64    CIP coordinates (Note 3)

Returned (function value):
    float64    the CIO locator s in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The CIO locator s is the difference between the right ascensions
    of the same point in two systems:  the two systems are the GCRS
    and the CIP,CIO, and the point is the ascending node of the
    CIP equator.  The quantity s remains below 0.1 arcsecond
    throughout 1900-2100.

 3) The series used to compute s is in fact for s+XY/2, where X and Y
    are the x and y components of the CIP unit vector;  this series
    is more compact than a direct series for s would be.  This
    function requires X,Y to be supplied by the caller, who is
    responsible for providing values that are consistent with the
    supplied date.

 4) The model is consistent with the IAU 2000A precession-nutation.

Called:
    Fal03     mean anomaly of the Moon
    Falp03    mean anomaly of the Sun
    Faf03     mean argument of the latitude of the Moon
    Fad03     mean elongation of the Moon from the Sun
    Faom03    mean longitude of the Moon's ascending node
    Fave03    mean longitude of Venus
    Fae03     mean longitude of Earth
    Fapa03    general accumulated precession in longitude

References:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func S00(date1, date2 float64, x, y float64) float64 {
	/* Time since J2000.0, in Julian centuries */
	var t float64

	/* Miscellaneous */
	var i, j int
	var a, w0, w1, w2, w3, w4, w5 float64

	/* Fundamental arguments */
	var fa [8]float64

	/* Returned value */
	var s float64

	/* --------------------- */
	/* The series for s+XY/2 */
	/* --------------------- */

	type TERM struct {
		nfa  [8]int  /* coefficients of l,l',F,D,Om,LVe,LE,pA */
		s, c float64 /* sine and cosine coefficients */
	}

	/* Polynomial coefficients */
	sp := []float64{

		/* 1-6 */
		94.00e-6,
		3808.35e-6,
		-119.94e-6,
		-72574.09e-6,
		27.70e-6,
		15.61e-6,
	}

	/* Terms of order t^0 */
	s0 := []TERM{

		/* 1-10 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, -2640.73e-6, 0.39e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -63.53e-6, 0.02e-6},
		{[8]int{0, 0, 2, -2, 3, 0, 0, 0}, -11.75e-6, -0.01e-6},
		{[8]int{0, 0, 2, -2, 1, 0, 0, 0}, -11.21e-6, -0.01e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, 4.57e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 3, 0, 0, 0}, -2.02e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 1, 0, 0, 0}, -1.98e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 3, 0, 0, 0}, 1.72e-6, 0.00e-6},
		{[8]int{0, 1, 0, 0, 1, 0, 0, 0}, 1.41e-6, 0.01e-6},
		{[8]int{0, 1, 0, 0, -1, 0, 0, 0}, 1.26e-6, 0.01e-6},

		/* 11-20 */
		{[8]int{1, 0, 0, 0, -1, 0, 0, 0}, 0.63e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, 1, 0, 0, 0}, 0.63e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 3, 0, 0, 0}, -0.46e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 1, 0, 0, 0}, -0.45e-6, 0.00e-6},
		{[8]int{0, 0, 4, -4, 4, 0, 0, 0}, -0.36e-6, 0.00e-6},
		{[8]int{0, 0, 1, -1, 1, -8, 12, 0}, 0.24e-6, 0.12e-6},
		{[8]int{0, 0, 2, 0, 0, 0, 0, 0}, -0.32e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, -0.28e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 3, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 1, 0, 0, 0}, -0.26e-6, 0.00e-6},

		/* 21-30 */
		{[8]int{0, 0, 2, -2, 0, 0, 0, 0}, 0.21e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -3, 0, 0, 0}, -0.19e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -1, 0, 0, 0}, -0.18e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 0, 8, -13, -1}, 0.10e-6, -0.05e-6},
		{[8]int{0, 0, 0, 2, 0, 0, 0, 0}, -0.15e-6, 0.00e-6},
		{[8]int{2, 0, -2, 0, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 2, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, 1, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{0, 0, 4, -2, 4, 0, 0, 0}, -0.13e-6, 0.00e-6},

		/* 31-33 */
		{[8]int{0, 0, 2, -2, 4, 0, 0, 0}, 0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -3, 0, 0, 0}, -0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -1, 0, 0, 0}, -0.11e-6, 0.00e-6},
	}

	/* Terms of order t^1 */
	s1 := []TERM{

		/* 1-3 */
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -0.07e-6, 3.57e-6},
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 1.71e-6, -0.03e-6},
		{[8]int{0, 0, 2, -2, 3, 0, 0, 0}, 0.00e-6, 0.48e-6},
	}

	/* Terms of order t^2 */
	s2 := []TERM{

		/* 1-10 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 743.53e-6, -0.17e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, 56.91e-6, 0.06e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, 9.84e-6, -0.01e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -8.85e-6, 0.01e-6},
		{[8]int{0, 1, 0, 0, 0, 0, 0, 0}, -6.38e-6, -0.05e-6},
		{[8]int{1, 0, 0, 0, 0, 0, 0, 0}, -3.07e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 2, 0, 0, 0}, 2.23e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 1, 0, 0, 0}, 1.67e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 2, 0, 0, 0}, 1.30e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -2, 0, 0, 0}, 0.93e-6, 0.00e-6},

		/* 11-20 */
		{[8]int{1, 0, 0, -2, 0, 0, 0, 0}, 0.68e-6, 0.00e-6},
		{[8]int{0, 0, 2, -2, 1, 0, 0, 0}, -0.55e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -2, 0, 0, 0}, 0.53e-6, 0.00e-6},
		{[8]int{0, 0, 0, 2, 0, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, 1, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, -2, -2, -2, 0, 0, 0}, -0.26e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, -1, 0, 0, 0}, -0.25e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 1, 0, 0, 0}, 0.22e-6, 0.00e-6},
		{[8]int{2, 0, 0, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
		{[8]int{2, 0, -2, 0, -1, 0, 0, 0}, 0.20e-6, 0.00e-6},

		/* 21-25 */
		{[8]int{0, 0, 2, 2, 2, 0, 0, 0}, 0.17e-6, 0.00e-6},
		{[8]int{2, 0, 2, 0, 2, 0, 0, 0}, 0.13e-6, 0.00e-6},
		{[8]int{2, 0, 0, 0, 0, 0, 0, 0}, -0.13e-6, 0.00e-6},
		{[8]int{1, 0, 2, -2, 2, 0, 0, 0}, -0.12e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 0, 0, 0, 0}, -0.11e-6, 0.00e-6},
	}

	/* Terms of order t^3 */
	s3 := []TERM{

		/* 1-4 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 0.30e-6, -23.51e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, -0.03e-6, -1.39e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, -0.01e-6, -0.24e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, 0.00e-6, 0.22e-6},
	}

	/* Terms of order t^4 */
	s4 := []TERM{

		/* 1-1 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, -0.26e-6, -0.01e-6},
	}

	/* Number of terms in the series */
	NS0 := len(s0)
	NS1 := len(s1)
	NS2 := len(s2)
	NS3 := len(s3)
	NS4 := len(s4)

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental epoch J2000.0 and current date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Fundamental Arguments (from IERS Conventions 2003) */

	/* Mean anomaly of the Moon. */
	fa[0] = Fal03(t)

	/* Mean anomaly of the Sun. */
	fa[1] = Falp03(t)

	/* Mean longitude of the Moon minus that of the ascending node. */
	fa[2] = Faf03(t)

	/* Mean elongation of the Moon from the Sun. */
	fa[3] = Fad03(t)

	/* Mean longitude of the ascending node of the Moon. */
	fa[4] = Faom03(t)

	/* Mean longitude of Venus. */
	fa[5] = Fave03(t)

	/* Mean longitude of Earth. */
	fa[6] = Fae03(t)

	/* General precession in longitude. */
	fa[7] = Fapa03(t)

	/* Evaluate s. */
	w0 = sp[0]
	w1 = sp[1]
	w2 = sp[2]
	w3 = sp[3]
	w4 = sp[4]
	w5 = sp[5]

	for i = NS0 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s0[i].nfa[j]) * fa[j]
		}
		w0 += s0[i].s*sin(a) + s0[i].c*cos(a)
	}

	for i = NS1 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s1[i].nfa[j]) * fa[j]
		}
		w1 += s1[i].s*sin(a) + s1[i].c*cos(a)
	}

	for i = NS2 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s2[i].nfa[j]) * fa[j]
		}
		w2 += s2[i].s*sin(a) + s2[i].c*cos(a)
	}

	for i = NS3 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s3[i].nfa[j]) * fa[j]
		}
		w3 += s3[i].s*sin(a) + s3[i].c*cos(a)
	}

	for i = NS4 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s4[i].nfa[j]) * fa[j]
		}
		w4 += s4[i].s*sin(a) + s4[i].c*cos(a)
	}

	s = (w0+
		(w1+
			(w2+
				(w3+
					(w4+
						w5*t)*t)*t)*t)*t)*DAS2R - x*y/2.0

	return s
}

/*
S00a The CIO locator s, IAU 2000A

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000A
precession-nutation model.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    the CIO locator s in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The CIO locator s is the difference between the right ascensions
    of the same point in two systems.  The two systems are the GCRS
    and the CIP,CIO, and the point is the ascending node of the
    CIP equator.  The CIO locator s remains a small fraction of
    1 arcsecond throughout 1900-2100.

 3) The series used to compute s is in fact for s+XY/2, where X and Y
    are the x and y components of the CIP unit vector;  this series
    is more compact than a direct series for s would be.  The present
    function uses the full IAU 2000A nutation model when predicting
    the CIP position.  Faster results, with no significant loss of
    accuracy, can be obtained via the function iauS00b, which uses
    instead the IAU 2000B truncated model.

Called:
    Pnm00a    classical NPB matrix, IAU 2000A
    Bnp2xy    extract CIP X,Y from the BPN matrix
    S00       the CIO locator s, given X,Y, IAU 2000A

References:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func S00a(date1, date2 float64) float64 {
	var rbpn [3][3]float64
	var x, y, s float64

	/* Bias-precession-nutation-matrix, IAU 2000A. */
	Pnm00a(date1, date2, &rbpn)

	/* Extract the CIP coordinates. */
	Bpn2xy(rbpn, &x, &y)

	/* Compute the CIO locator s, given the CIP coordinates. */
	s = S00(date1, date2, x, y)

	return s
}

/*
S00b The CIO locator s, IAU 2000B

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000B
precession-nutation model.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    the CIO locator s in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The CIO locator s is the difference between the right ascensions
    of the same point in two systems.  The two systems are the GCRS
    and the CIP,CIO, and the point is the ascending node of the
    CIP equator.  The CIO locator s remains a small fraction of
    1 arcsecond throughout 1900-2100.

 3) The series used to compute s is in fact for s+XY/2, where X and Y
    are the x and y components of the CIP unit vector;  this series
    is more compact than a direct series for s would be.  The present
    function uses the IAU 2000B truncated nutation model when
    predicting the CIP position.  The function iauS00a uses instead
    the full IAU 2000A model, but with no significant increase in
    accuracy and at some cost in speed.

Called:
    Pnm00b    classical NPB matrix, IAU 2000B
    Bnp2xy    extract CIP X,Y from the BPN matrix
    S00       the CIO locator s, given X,Y, IAU 2000A

References:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func S00b(date1, date2 float64) float64 {
	var rbpn [3][3]float64
	var x, y, s float64

	/* Bias-precession-nutation-matrix, IAU 2000B. */
	Pnm00b(date1, date2, &rbpn)

	/* Extract the CIP coordinates. */
	Bpn2xy(rbpn, &x, &y)

	/* Compute the CIO locator s, given the CIP coordinates. */
	s = S00(date1, date2, x, y)

	return s
}

/*
S06 The CIO locator s, given X,Y, IAU 2006

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2006/2000A precession-nutation.

Given:
    date1,date2   float64    TT as a 2-part Julian Date (Note 1)
    x,y           float64    CIP coordinates (Note 3)

Returned (function value):
    float64    the CIO locator s in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The CIO locator s is the difference between the right ascensions
    of the same point in two systems:  the two systems are the GCRS
    and the CIP,CIO, and the point is the ascending node of the
    CIP equator.  The quantity s remains below 0.1 arcsecond
    throughout 1900-2100.

 3) The series used to compute s is in fact for s+XY/2, where X and Y
    are the x and y components of the CIP unit vector;  this series
    is more compact than a direct series for s would be.  This
    function requires X,Y to be supplied by the caller, who is
    responsible for providing values that are consistent with the
    supplied date.

 4) The model is consistent with the "P03" precession (Capitaine et
    al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
    2000A nutation (with P03 adjustments).

Called:
    Fal03     mean anomaly of the Moon
    Falp03    mean anomaly of the Sun
    Faf03     mean argument of the latitude of the Moon
    Fad03     mean elongation of the Moon from the Sun
    Faom03    mean longitude of the Moon's ascending node
    Fave03    mean longitude of Venus
    Fae03     mean longitude of Earth
    Fapa03    general accumulated precession in longitude

References:

    Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
    Astrophys. 432, 355

    McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG
*/
func S06(date1, date2 float64, x, y float64) float64 {
	/* Time since J2000.0, in Julian centuries */
	var t float64

	/* Miscellaneous */
	var i, j int
	var a, w0, w1, w2, w3, w4, w5 float64

	/* Fundamental arguments */
	var fa [8]float64

	/* Returned value */
	var s float64

	/* --------------------- */
	/* The series for s+XY/2 */
	/* --------------------- */

	type TERM struct {
		nfa  [8]int  /* coefficients of l,l',F,D,Om,LVe,LE,pA */
		s, c float64 /* sine and cosine coefficients */
	}

	/* Polynomial coefficients */
	sp := []float64{

		/* 1-6 */
		94.00e-6,
		3808.65e-6,
		-122.68e-6,
		-72574.11e-6,
		27.98e-6,
		15.62e-6,
	}

	/* Terms of order t^0 */
	s0 := []TERM{

		/* 1-10 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, -2640.73e-6, 0.39e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -63.53e-6, 0.02e-6},
		{[8]int{0, 0, 2, -2, 3, 0, 0, 0}, -11.75e-6, -0.01e-6},
		{[8]int{0, 0, 2, -2, 1, 0, 0, 0}, -11.21e-6, -0.01e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, 4.57e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 3, 0, 0, 0}, -2.02e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 1, 0, 0, 0}, -1.98e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 3, 0, 0, 0}, 1.72e-6, 0.00e-6},
		{[8]int{0, 1, 0, 0, 1, 0, 0, 0}, 1.41e-6, 0.01e-6},
		{[8]int{0, 1, 0, 0, -1, 0, 0, 0}, 1.26e-6, 0.01e-6},

		/* 11-20 */
		{[8]int{1, 0, 0, 0, -1, 0, 0, 0}, 0.63e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, 1, 0, 0, 0}, 0.63e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 3, 0, 0, 0}, -0.46e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 1, 0, 0, 0}, -0.45e-6, 0.00e-6},
		{[8]int{0, 0, 4, -4, 4, 0, 0, 0}, -0.36e-6, 0.00e-6},
		{[8]int{0, 0, 1, -1, 1, -8, 12, 0}, 0.24e-6, 0.12e-6},
		{[8]int{0, 0, 2, 0, 0, 0, 0, 0}, -0.32e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, -0.28e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 3, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 1, 0, 0, 0}, -0.26e-6, 0.00e-6},

		/* 21-30 */
		{[8]int{0, 0, 2, -2, 0, 0, 0, 0}, 0.21e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -3, 0, 0, 0}, -0.19e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -1, 0, 0, 0}, -0.18e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 0, 8, -13, -1}, 0.10e-6, -0.05e-6},
		{[8]int{0, 0, 0, 2, 0, 0, 0, 0}, -0.15e-6, 0.00e-6},
		{[8]int{2, 0, -2, 0, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 2, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, 1, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{0, 0, 4, -2, 4, 0, 0, 0}, -0.13e-6, 0.00e-6},

		/* 31-33 */
		{[8]int{0, 0, 2, -2, 4, 0, 0, 0}, 0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -3, 0, 0, 0}, -0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -1, 0, 0, 0}, -0.11e-6, 0.00e-6},
	}

	/* Terms of order t^1 */
	s1 := []TERM{

		/* 1 - 3 */
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -0.07e-6, 3.57e-6},
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 1.73e-6, -0.03e-6},
		{[8]int{0, 0, 2, -2, 3, 0, 0, 0}, 0.00e-6, 0.48e-6},
	}

	/* Terms of order t^2 */
	s2 := []TERM{

		/* 1-10 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 743.52e-6, -0.17e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, 56.91e-6, 0.06e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, 9.84e-6, -0.01e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, -8.85e-6, 0.01e-6},
		{[8]int{0, 1, 0, 0, 0, 0, 0, 0}, -6.38e-6, -0.05e-6},
		{[8]int{1, 0, 0, 0, 0, 0, 0, 0}, -3.07e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 2, 0, 0, 0}, 2.23e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 1, 0, 0, 0}, 1.67e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 2, 0, 0, 0}, 1.30e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -2, 0, 0, 0}, 0.93e-6, 0.00e-6},

		/* 11-20 */
		{[8]int{1, 0, 0, -2, 0, 0, 0, 0}, 0.68e-6, 0.00e-6},
		{[8]int{0, 0, 2, -2, 1, 0, 0, 0}, -0.55e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -2, 0, 0, 0}, 0.53e-6, 0.00e-6},
		{[8]int{0, 0, 0, 2, 0, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, 1, 0, 0, 0}, -0.27e-6, 0.00e-6},
		{[8]int{1, 0, -2, -2, -2, 0, 0, 0}, -0.26e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, -1, 0, 0, 0}, -0.25e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 1, 0, 0, 0}, 0.22e-6, 0.00e-6},
		{[8]int{2, 0, 0, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
		{[8]int{2, 0, -2, 0, -1, 0, 0, 0}, 0.20e-6, 0.00e-6},

		/* 21-25 */
		{[8]int{0, 0, 2, 2, 2, 0, 0, 0}, 0.17e-6, 0.00e-6},
		{[8]int{2, 0, 2, 0, 2, 0, 0, 0}, 0.13e-6, 0.00e-6},
		{[8]int{2, 0, 0, 0, 0, 0, 0, 0}, -0.13e-6, 0.00e-6},
		{[8]int{1, 0, 2, -2, 2, 0, 0, 0}, -0.12e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 0, 0, 0, 0}, -0.11e-6, 0.00e-6},
	}

	/* Terms of order t^3 */
	s3 := []TERM{

		/* 1-4 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 0.30e-6, -23.42e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, -0.03e-6, -1.46e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, -0.01e-6, -0.25e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, 0.00e-6, 0.23e-6},
	}

	/* Terms of order t^4 */
	s4 := []TERM{

		/* 1-1 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, -0.26e-6, -0.01e-6},
	}

	/* Number of terms in the series */
	NS0 := len(s0)
	NS1 := len(s1)
	NS2 := len(s2)
	NS3 := len(s3)
	NS4 := len(s4)

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental epoch J2000.0 and current date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Fundamental Arguments (from IERS Conventions 2003) */

	/* Mean anomaly of the Moon. */
	fa[0] = Fal03(t)

	/* Mean anomaly of the Sun. */
	fa[1] = Falp03(t)

	/* Mean longitude of the Moon minus that of the ascending node. */
	fa[2] = Faf03(t)

	/* Mean elongation of the Moon from the Sun. */
	fa[3] = Fad03(t)

	/* Mean longitude of the ascending node of the Moon. */
	fa[4] = Faom03(t)

	/* Mean longitude of Venus. */
	fa[5] = Fave03(t)

	/* Mean longitude of Earth. */
	fa[6] = Fae03(t)

	/* General precession in longitude. */
	fa[7] = Fapa03(t)

	/* Evaluate s. */
	w0 = sp[0]
	w1 = sp[1]
	w2 = sp[2]
	w3 = sp[3]
	w4 = sp[4]
	w5 = sp[5]

	for i = NS0 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s0[i].nfa[j]) * fa[j]
		}
		w0 += s0[i].s*sin(a) + s0[i].c*cos(a)
	}

	for i = NS1 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s1[i].nfa[j]) * fa[j]
		}
		w1 += s1[i].s*sin(a) + s1[i].c*cos(a)
	}

	for i = NS2 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s2[i].nfa[j]) * fa[j]
		}
		w2 += s2[i].s*sin(a) + s2[i].c*cos(a)
	}

	for i = NS3 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s3[i].nfa[j]) * fa[j]
		}
		w3 += s3[i].s*sin(a) + s3[i].c*cos(a)
	}

	for i = NS4 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(s4[i].nfa[j]) * fa[j]
		}
		w4 += s4[i].s*sin(a) + s4[i].c*cos(a)
	}

	s = (w0+(w1+(w2+(w3+(w4+w5*t)*t)*t)*t)*t)*DAS2R - x*y/2.0

	return s
}

/*
S06a The CIO locator s, IAU 2006/2000A

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2006
precession and IAU 2000A nutation models.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    the CIO locator s in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The CIO locator s is the difference between the right ascensions
    of the same point in two systems.  The two systems are the GCRS
    and the CIP,CIO, and the point is the ascending node of the
    CIP equator.  The CIO locator s remains a small fraction of
    1 arcsecond throughout 1900-2100.

 3) The series used to compute s is in fact for s+XY/2, where X and Y
    are the x and y components of the CIP unit vector;  this series is
    more compact than a direct series for s would be.  The present
    function uses the full IAU 2000A nutation model when predicting
    the CIP position.

Called:
    Pnm06a    classical NPB matrix, IAU 2006/2000A
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S06       the CIO locator s, given X,Y, IAU 2006

References:

    Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    "Expressions for the Celestial Intermediate Pole and Celestial
    Ephemeris Origin consistent with the IAU 2000A precession-
    nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

    n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func S06a(date1, date2 float64) float64 {
	var rnpb [3][3]float64
	var x, y, s float64

	/* Bias-precession-nutation-matrix, IAU 20006/2000A. */
	Pnm06a(date1, date2, &rnpb)

	/* Extract the CIP coordinates. */
	Bpn2xy(rnpb, &x, &y)

	/* Compute the CIO locator s, given the CIP coordinates. */
	s = S06(date1, date2, x, y)

	return s
}

/*
Sp00 The TIO locator s', IERS 2003

The TIO locator s', positioning the Terrestrial Intermediate Origin
on the equator of the Celestial Intermediate Pole.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    the TIO locator s' in radians (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The TIO locator s' is obtained from polar motion observations by
    numerical integration, and so is in essence unpredictable.
    However, it is dominated by a secular drift of about
    47 microarcseconds per century, which is the approximation
    evaluated by the present function.

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Sp00(date1, date2 float64) float64 {
	var t, sp float64

	/* Interval between fundamental epoch J2000.0 and current date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Approximate s'. */
	sp = -47e-6 * t * DAS2R

	return sp
}

/*
Xy06 CIP, IAU 2006/2000A, from series

X,Y coordinates of celestial intermediate pole from series based
on IAU 2006 precession and IAU 2000A nutation.

Given:
    date1,date2  float64     TT as a 2-part Julian Date (Note 1)

Returned:
    x,y          float64     CIP X,Y coordinates (Note 2)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The X,Y coordinates are those of the unit vector towards the
    celestial intermediate pole.  They represent the combined effects
    of frame bias, precession and nutation.

 3) The fundamental arguments used are as adopted in IERS Conventions
    (2003) and are from Simon et al. (1994) and Souchay et al.
    (1999).

 4) This is an alternative to the angles-based method, via the SOFA
    function iauFw2xy and as used in iauXys06a for example.  The two
    methods agree at the 1 microarcsecond level (at present), a
    negligible amount compared with the intrinsic accuracy of the
    models.  However, it would be unwise to mix the two methods
    (angles-based and series-based) in a single application.

Called:
    Fal03     mean anomaly of the Moon
    Falp03    mean anomaly of the Sun
    Faf03     mean argument of the latitude of the Moon
    Fad03     mean elongation of the Moon from the Sun
    Faom03    mean longitude of the Moon's ascending node
    Fame03    mean longitude of Mercury
    Fave03    mean longitude of Venus
    Fae03     mean longitude of Earth
    Fama03    mean longitude of Mars
    Faju03    mean longitude of Jupiter
    Fasa03    mean longitude of Saturn
    Faur03    mean longitude of Uranus
    Fane03    mean longitude of Neptune
    Fapa03    general accumulated precession in longitude

References:

    Capitaine, N., Wallace, P.T. & Chapront, J., 2003,
    Astron.Astrophys., 412, 567

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG

    Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Xy06(date1, date2 float64, x, y *float64) {
	/* Maximum power of T in the polynomials for X and Y */
	const MAXPT = 5

	/* Polynomial coefficients (arcsec, X then Y). */
	xyp := [2][MAXPT + 1]float64{
		{-0.016617,
			2004.191898,
			-0.4297829,
			-0.19861834,
			0.000007578,
			0.0000059285,
		},
		{-0.006951,
			-0.025896,
			-22.4072747,
			0.00190059,
			0.001112526,
			0.0000001358,
		},
	}

	/* Fundamental-argument multipliers:  luni-solar terms */
	mfals := [][5]int{

		/* 1-10 */
		{0, 0, 0, 0, 1},
		{0, 0, 2, -2, 2},
		{0, 0, 2, 0, 2},
		{0, 0, 0, 0, 2},
		{0, 1, 0, 0, 0},
		{0, 1, 2, -2, 2},
		{1, 0, 0, 0, 0},
		{0, 0, 2, 0, 1},
		{1, 0, 2, 0, 2},
		{0, 1, -2, 2, -2},

		/* 11-20 */
		{0, 0, 2, -2, 1},
		{1, 0, -2, 0, -2},
		{1, 0, 0, -2, 0},
		{1, 0, 0, 0, 1},
		{1, 0, 0, 0, -1},
		{1, 0, -2, -2, -2},
		{1, 0, 2, 0, 1},
		{2, 0, -2, 0, -1},
		{0, 0, 0, 2, 0},
		{0, 0, 2, 2, 2},

		/* 21-30 */
		{2, 0, 0, -2, 0},
		{0, 2, -2, 2, -2},
		{2, 0, 2, 0, 2},
		{1, 0, 2, -2, 2},
		{1, 0, -2, 0, -1},
		{2, 0, 0, 0, 0},
		{0, 0, 2, 0, 0},
		{0, 1, 0, 0, 1},
		{1, 0, 0, -2, -1},
		{0, 2, 2, -2, 2},

		/* 31-40 */
		{0, 0, 2, -2, 0},
		{1, 0, 0, -2, 1},
		{0, 1, 0, 0, -1},
		{0, 2, 0, 0, 0},
		{1, 0, -2, -2, -1},
		{1, 0, 2, 2, 2},
		{0, 1, 2, 0, 2},
		{2, 0, -2, 0, 0},
		{0, 0, 2, 2, 1},
		{0, 1, -2, 0, -2},

		/* 41-50 */
		{0, 0, 0, 2, 1},
		{1, 0, 2, -2, 1},
		{2, 0, 0, -2, -1},
		{2, 0, 2, -2, 2},
		{2, 0, 2, 0, 1},
		{0, 0, 0, 2, -1},
		{0, 1, -2, 2, -1},
		{1, 1, 0, -2, 0},
		{2, 0, 0, -2, 1},
		{1, 0, 0, 2, 0},

		/* 51-60 */
		{0, 1, 2, -2, 1},
		{1, -1, 0, 0, 0},
		{0, 1, -1, 1, -1},
		{2, 0, -2, 0, -2},
		{0, 1, 0, -2, 0},
		{1, 0, 0, -1, 0},
		{3, 0, 2, 0, 2},
		{0, 0, 0, 1, 0},
		{1, -1, 2, 0, 2},
		{1, 1, -2, -2, -2},

		/* 61-70 */
		{1, 0, -2, 0, 0},
		{2, 0, 0, 0, -1},
		{0, 1, -2, -2, -2},
		{1, 1, 2, 0, 2},
		{2, 0, 0, 0, 1},
		{1, 1, 0, 0, 0},
		{1, 0, -2, 2, -1},
		{1, 0, 2, 0, 0},
		{1, -1, 0, -1, 0},
		{1, 0, 0, 0, 2},

		/* 71-80 */
		{1, 0, -1, 0, -1},
		{0, 0, 2, 1, 2},
		{1, 0, -2, -4, -2},
		{1, -1, 0, -1, -1},
		{1, 0, 2, 2, 1},
		{0, 2, -2, 2, -1},
		{1, 0, 0, 0, -2},
		{2, 0, -2, -2, -2},
		{1, 1, 2, -2, 2},
		{2, 0, -2, -4, -2},

		/* 81-90 */
		{1, 0, -4, 0, -2},
		{2, 0, 2, -2, 1},
		{1, 0, 0, -1, -1},
		{2, 0, 2, 2, 2},
		{3, 0, 0, 0, 0},
		{1, 0, 0, 2, 1},
		{0, 0, 2, -2, -1},
		{3, 0, 2, -2, 2},
		{0, 0, 4, -2, 2},
		{1, 0, 0, -4, 0},

		/* 91-100 */
		{0, 1, 2, 0, 1},
		{2, 0, 0, -4, 0},
		{1, 1, 0, -2, -1},
		{2, 0, -2, 0, 1},
		{0, 0, 2, 0, -1},
		{0, 1, -2, 0, -1},
		{0, 1, 0, 0, 2},
		{0, 0, 2, -1, 2},
		{0, 0, 2, 4, 2},
		{2, 1, 0, -2, 0},

		/* 101-110 */
		{1, 1, 0, -2, 1},
		{1, -1, 0, -2, 0},
		{1, -1, 0, -1, -2},
		{1, -1, 0, 0, 1},
		{0, 1, -2, 2, 0},
		{0, 1, 0, 0, -2},
		{1, -1, 2, 2, 2},
		{1, 0, 0, 2, -1},
		{1, -1, -2, -2, -2},
		{3, 0, 2, 0, 1},

		/* 111-120 */
		{0, 1, 2, 2, 2},
		{1, 0, 2, -2, 0},
		{1, 1, -2, -2, -1},
		{1, 0, 2, -4, 1},
		{0, 1, -2, -2, -1},
		{2, -1, 2, 0, 2},
		{0, 0, 0, 2, 2},
		{1, -1, 2, 0, 1},
		{1, -1, -2, 0, -2},
		{0, 1, 0, 2, 0},

		/* 121-130 */
		{0, 1, 2, -2, 0},
		{0, 0, 0, 1, 1},
		{1, 0, -2, -2, 0},
		{0, 3, 2, -2, 2},
		{2, 1, 2, 0, 2},
		{1, 1, 0, 0, 1},
		{2, 0, 0, 2, 0},
		{1, 1, 2, 0, 1},
		{1, 0, 0, -2, -2},
		{1, 0, -2, 2, 0},

		/* 131-140 */
		{1, 0, -1, 0, -2},
		{0, 1, 0, -2, 1},
		{0, 1, 0, 1, 0},
		{0, 0, 0, 1, -1},
		{1, 0, -2, 2, -2},
		{1, -1, 0, 0, -1},
		{0, 0, 0, 4, 0},
		{1, -1, 0, 2, 0},
		{1, 0, 2, 1, 2},
		{1, 0, 2, -1, 2},

		/* 141-150 */
		{0, 0, 2, 1, 1},
		{1, 0, 0, -2, 2},
		{1, 0, -2, 0, 1},
		{1, 0, -2, -4, -1},
		{0, 0, 2, 2, 0},
		{1, 1, 2, -2, 1},
		{1, 0, -2, 1, -1},
		{0, 0, 1, 0, 1},
		{2, 0, -2, -2, -1},
		{4, 0, 2, 0, 2},

		/* 151-160 */
		{2, -1, 0, 0, 0},
		{2, 1, 2, -2, 2},
		{0, 1, 2, 1, 2},
		{1, 0, 4, -2, 2},
		{1, 1, 0, 0, -1},
		{2, 0, 2, 0, 0},
		{2, 0, -2, -4, -1},
		{1, 0, -1, 0, 0},
		{1, 0, 0, 1, 0},
		{0, 1, 0, 2, 1},

		/* 161-170 */
		{1, 0, -4, 0, -1},
		{1, 0, 0, -4, -1},
		{2, 0, 2, 2, 1},
		{2, 1, 0, 0, 0},
		{0, 0, 2, -3, 2},
		{1, 2, 0, -2, 0},
		{0, 3, 0, 0, 0},
		{0, 0, 4, 0, 2},
		{0, 0, 2, -4, 1},
		{2, 0, 0, -2, -2},

		/* 171-180 */
		{1, 1, -2, -4, -2},
		{0, 1, 0, -2, -1},
		{0, 0, 0, 4, 1},
		{3, 0, 2, -2, 1},
		{1, 0, 2, 4, 2},
		{1, 1, -2, 0, -2},
		{0, 0, 4, -2, 1},
		{2, -2, 0, -2, 0},
		{2, 1, 0, -2, -1},
		{0, 2, 0, -2, 0},

		/* 181-190 */
		{1, 0, 0, -1, 1},
		{1, 1, 2, 2, 2},
		{3, 0, 0, 0, -1},
		{2, 0, 0, -4, -1},
		{3, 0, 2, 2, 2},
		{0, 0, 2, 4, 1},
		{0, 2, -2, -2, -2},
		{1, -1, 0, -2, -1},
		{0, 0, 2, -1, 1},
		{2, 0, 0, 2, 1},

		/* 191-200 */
		{1, -1, -2, 2, -1},
		{0, 0, 0, 2, -2},
		{2, 0, 0, -4, 1},
		{1, 0, 0, -4, 1},
		{2, 0, 2, -4, 1},
		{4, 0, 2, -2, 2},
		{2, 1, -2, 0, -1},
		{2, 1, -2, -4, -2},
		{3, 0, 0, -4, 0},
		{1, -1, 2, 2, 1},

		/* 201-210 */
		{1, -1, -2, 0, -1},
		{0, 2, 0, 0, 1},
		{1, 2, -2, -2, -2},
		{1, 1, 0, -4, 0},
		{2, 0, 0, -2, 2},
		{0, 2, 2, -2, 1},
		{1, 0, 2, 0, -1},
		{2, 1, 0, -2, 1},
		{2, -1, -2, 0, -1},
		{1, -1, -2, -2, -1},

		/* 211-220 */
		{0, 1, -2, 1, -2},
		{1, 0, -4, 2, -2},
		{0, 1, 2, 2, 1},
		{3, 0, 0, 0, 1},
		{2, -1, 2, 2, 2},
		{0, 1, -2, -4, -2},
		{1, 0, -2, -3, -2},
		{2, 0, 0, 0, 2},
		{1, -1, 0, -2, -2},
		{2, 0, -2, 2, -1},

		/* 221-230 */
		{0, 2, -2, 0, -2},
		{3, 0, -2, 0, -1},
		{2, -1, 2, 0, 1},
		{1, 0, -2, -1, -2},
		{0, 0, 2, 0, 3},
		{2, 0, -4, 0, -2},
		{2, 1, 0, -4, 0},
		{1, 1, -2, 1, -1},
		{0, 2, 2, 0, 2},
		{1, -1, 2, -2, 2},

		/* 231-240 */
		{1, -1, 0, -2, 1},
		{2, 1, 2, 0, 1},
		{1, 0, 2, -4, 2},
		{1, 1, -2, 0, -1},
		{1, 1, 0, 2, 0},
		{1, 0, 0, -3, 0},
		{2, 0, 2, -1, 2},
		{0, 2, 0, 0, -1},
		{2, -1, 0, -2, 0},
		{4, 0, 0, 0, 0},

		/* 241-250 */
		{2, 1, -2, -2, -2},
		{0, 2, -2, 2, 0},
		{1, 0, 2, 1, 1},
		{1, 0, -1, 0, -3},
		{3, -1, 2, 0, 2},
		{2, 0, 2, -2, 0},
		{1, -2, 0, 0, 0},
		{2, 0, 0, 0, -2},
		{1, 0, 0, 4, 0},
		{0, 1, 0, 1, 1},

		/* 251-260 */
		{1, 0, 2, 2, 0},
		{0, 1, 0, 2, -1},
		{0, 1, 0, 1, -1},
		{0, 0, 2, -2, 3},
		{3, 1, 2, 0, 2},
		{1, 1, 2, 1, 2},
		{1, 1, -2, 2, -1},
		{2, -1, 2, -2, 2},
		{1, -2, 2, 0, 2},
		{1, 0, 2, -4, 0},

		/* 261-270 */
		{0, 0, 1, 0, 0},
		{1, 0, 2, -3, 1},
		{1, -2, 0, -2, 0},
		{2, 0, 0, 2, -1},
		{1, 1, 2, -4, 1},
		{4, 0, 2, 0, 1},
		{0, 1, 2, 1, 1},
		{1, 2, 2, -2, 2},
		{2, 0, 2, 1, 2},
		{2, 1, 2, -2, 1},

		/* 271-280 */
		{1, 0, 2, -1, 1},
		{1, 0, 4, -2, 1},
		{1, -1, 2, -2, 1},
		{0, 1, 0, -4, 0},
		{3, 0, -2, -2, -2},
		{0, 0, 4, -4, 2},
		{2, 0, -4, -2, -2},
		{2, -2, 0, -2, -1},
		{1, 0, 2, -2, -1},
		{2, 0, -2, -6, -2},

		/* 281-290 */
		{1, 0, -2, 1, -2},
		{1, 0, -2, 2, 1},
		{1, -1, 0, 2, -1},
		{1, 0, -2, 1, 0},
		{2, -1, 0, -2, 1},
		{1, -1, 0, 2, 1},
		{2, 0, -2, -2, 0},
		{1, 0, 2, -3, 2},
		{0, 0, 0, 4, -1},
		{2, -1, 0, 0, 1},

		/* 291-300 */
		{2, 0, 4, -2, 2},
		{0, 0, 2, 3, 2},
		{0, 1, 4, -2, 2},
		{0, 1, -2, 2, 1},
		{1, 1, 0, 2, 1},
		{1, 0, 0, 4, 1},
		{0, 0, 4, 0, 1},
		{2, 0, 0, -3, 0},
		{1, 0, 0, -1, -2},
		{1, -2, -2, -2, -2},

		/* 301-310 */
		{3, 0, 0, 2, 0},
		{2, 0, 2, -4, 2},
		{1, 1, -2, -4, -1},
		{1, 0, -2, -6, -2},
		{2, -1, 0, 0, -1},
		{2, -1, 0, 2, 0},
		{0, 1, 2, -2, -1},
		{1, 1, 0, 1, 0},
		{1, 2, 0, -2, -1},
		{1, 0, 0, 1, -1},

		/* 311-320 */
		{0, 0, 1, 0, 2},
		{3, 1, 2, -2, 2},
		{1, 0, -4, -2, -2},
		{1, 0, 2, 4, 1},
		{1, -2, 2, 2, 2},
		{1, -1, -2, -4, -2},
		{0, 0, 2, -4, 2},
		{0, 0, 2, -3, 1},
		{2, 1, -2, 0, 0},
		{3, 0, -2, -2, -1},

		/* 321-330 */
		{2, 0, 2, 4, 2},
		{0, 0, 0, 0, 3},
		{2, -1, -2, -2, -2},
		{2, 0, 0, -1, 0},
		{3, 0, 2, -4, 2},
		{2, 1, 2, 2, 2},
		{0, 0, 3, 0, 3},
		{1, 1, 2, 2, 1},
		{2, 1, 0, 0, -1},
		{1, 2, 0, -2, 1},

		/* 331-340 */
		{3, 0, 2, 2, 1},
		{1, -1, -2, 2, -2},
		{1, 1, 0, -1, 0},
		{1, 2, 0, 0, 0},
		{1, 0, 4, 0, 2},
		{1, -1, 2, 4, 2},
		{2, 1, 0, 0, 1},
		{1, 0, 0, 2, 2},
		{1, -1, -2, 2, 0},
		{0, 2, -2, -2, -1},

		/* 341-350 */
		{2, 0, -2, 0, 2},
		{5, 0, 2, 0, 2},
		{3, 0, -2, -6, -2},
		{1, -1, 2, -1, 2},
		{3, 0, 0, -4, -1},
		{1, 0, 0, 1, 1},
		{1, 0, -4, 2, -1},
		{0, 1, 2, -4, 1},
		{1, 2, 2, 0, 2},
		{0, 1, 0, -2, -2},

		/* 351-360 */
		{0, 0, 2, -1, 0},
		{1, 0, 1, 0, 1},
		{0, 2, 0, -2, 1},
		{3, 0, 2, 0, 0},
		{1, 1, -2, 1, 0},
		{2, 1, -2, -4, -1},
		{3, -1, 0, 0, 0},
		{2, -1, -2, 0, 0},
		{4, 0, 2, -2, 1},
		{2, 0, -2, 2, 0},

		/* 361-370 */
		{1, 1, 2, -2, 0},
		{1, 0, -2, 4, -1},
		{1, 0, -2, -2, 1},
		{2, 0, 2, -4, 0},
		{1, 1, 0, -2, -2},
		{1, 1, -2, -2, 0},
		{1, 0, 1, -2, 1},
		{2, -1, -2, -4, -2},
		{3, 0, -2, 0, -2},
		{0, 1, -2, -2, 0},

		/* 371-380 */
		{3, 0, 0, -2, -1},
		{1, 0, -2, -3, -1},
		{0, 1, 0, -4, -1},
		{1, -2, 2, -2, 1},
		{0, 1, -2, 1, -1},
		{1, -1, 0, 0, 2},
		{2, 0, 0, 1, 0},
		{1, -2, 0, 2, 0},
		{1, 2, -2, -2, -1},
		{0, 0, 4, -4, 1},

		/* 381-390 */
		{0, 1, 2, 4, 2},
		{0, 1, -4, 2, -2},
		{3, 0, -2, 0, 0},
		{2, -1, 2, 2, 1},
		{0, 1, -2, -4, -1},
		{4, 0, 2, 2, 2},
		{2, 0, -2, -3, -2},
		{2, 0, 0, -6, 0},
		{1, 0, 2, 0, 3},
		{3, 1, 0, 0, 0},

		/* 391-400 */
		{3, 0, 0, -4, 1},
		{1, -1, 2, 0, 0},
		{1, -1, 0, -4, 0},
		{2, 0, -2, 2, -2},
		{1, 1, 0, -2, 2},
		{4, 0, 0, -2, 0},
		{2, 2, 0, -2, 0},
		{0, 1, 2, 0, 0},
		{1, 1, 0, -4, 1},
		{1, 0, 0, -4, -2},

		/* 401-410 */
		{0, 0, 0, 1, 2},
		{3, 0, 0, 2, 1},
		{1, 1, 0, -4, -1},
		{0, 0, 2, 2, -1},
		{1, 1, 2, 0, 0},
		{1, -1, 2, -4, 1},
		{1, 1, 0, 0, 2},
		{0, 0, 2, 6, 2},
		{4, 0, -2, -2, -1},
		{2, 1, 0, -4, -1},

		/* 411-420 */
		{0, 0, 0, 3, 1},
		{1, -1, -2, 0, 0},
		{0, 0, 2, 1, 0},
		{1, 0, 0, 2, -2},
		{3, -1, 2, 2, 2},
		{3, -1, 2, -2, 2},
		{1, 0, 0, -1, 2},
		{1, -2, 2, -2, 2},
		{0, 1, 0, 2, 2},
		{0, 1, -2, -1, -2},

		/* 421-430 */
		{1, 1, -2, 0, 0},
		{0, 2, 2, -2, 0},
		{3, -1, -2, -1, -2},
		{1, 0, 0, -6, 0},
		{1, 0, -2, -4, 0},
		{2, 1, 0, -4, 1},
		{2, 0, 2, 0, -1},
		{2, 0, -4, 0, -1},
		{0, 0, 3, 0, 2},
		{2, 1, -2, -2, -1},

		/* 431-440 */
		{1, -2, 0, 0, 1},
		{2, -1, 0, -4, 0},
		{0, 0, 0, 3, 0},
		{5, 0, 2, -2, 2},
		{1, 2, -2, -4, -2},
		{1, 0, 4, -4, 2},
		{0, 0, 4, -1, 2},
		{3, 1, 0, -4, 0},
		{3, 0, 0, -6, 0},
		{2, 0, 0, 2, 2},

		/* 441-450 */
		{2, -2, 2, 0, 2},
		{1, 0, 0, -3, 1},
		{1, -2, -2, 0, -2},
		{1, -1, -2, -3, -2},
		{0, 0, 2, -2, -2},
		{2, 0, -2, -4, 0},
		{1, 0, -4, 0, 0},
		{0, 1, 0, -1, 0},
		{4, 0, 0, 0, -1},
		{3, 0, 2, -1, 2},

		/* 451-460 */
		{3, -1, 2, 0, 1},
		{2, 0, 2, -1, 1},
		{1, 2, 2, -2, 1},
		{1, 1, 0, 2, -1},
		{0, 2, 2, 0, 1},
		{3, 1, 2, 0, 1},
		{1, 1, 2, 1, 1},
		{1, 1, 0, -1, 1},
		{1, -2, 0, -2, -1},
		{4, 0, 0, -4, 0},

		/* 461-470 */
		{2, 1, 0, 2, 0},
		{1, -1, 0, 4, 0},
		{0, 1, 0, -2, 2},
		{0, 0, 2, 0, -2},
		{1, 0, -1, 0, 1},
		{3, 0, 2, -2, 0},
		{2, 0, 2, 2, 0},
		{1, 2, 0, -4, 0},
		{1, -1, 0, -3, 0},
		{0, 1, 0, 4, 0},

		/* 471 - 480 */
		{0, 1, -2, 0, 0},
		{2, 2, 2, -2, 2},
		{0, 0, 0, 1, -2},
		{0, 2, -2, 0, -1},
		{4, 0, 2, -4, 2},
		{2, 0, -4, 2, -2},
		{2, -1, -2, 0, -2},
		{1, 1, 4, -2, 2},
		{1, 1, 2, -4, 2},
		{1, 0, 2, 3, 2},

		/* 481-490 */
		{1, 0, 0, 4, -1},
		{0, 0, 0, 4, 2},
		{2, 0, 0, 4, 0},
		{1, 1, -2, 2, 0},
		{2, 1, 2, 1, 2},
		{2, 1, 2, -4, 1},
		{2, 0, 2, 1, 1},
		{2, 0, -4, -2, -1},
		{2, 0, -2, -6, -1},
		{2, -1, 2, -1, 2},

		/* 491-500 */
		{1, -2, 2, 0, 1},
		{1, -2, 0, -2, 1},
		{1, -1, 0, -4, -1},
		{0, 2, 2, 2, 2},
		{0, 2, -2, -4, -2},
		{0, 1, 2, 3, 2},
		{0, 1, 0, -4, 1},
		{3, 0, 0, -2, 1},
		{2, 1, -2, 0, 1},
		{2, 0, 4, -2, 1},

		/* 501-510 */
		{2, 0, 0, -3, -1},
		{2, -2, 0, -2, 1},
		{2, -1, 2, -2, 1},
		{1, 0, 0, -6, -1},
		{1, -2, 0, 0, -1},
		{1, -2, -2, -2, -1},
		{0, 1, 4, -2, 1},
		{0, 0, 2, 3, 1},
		{2, -1, 0, -1, 0},
		{1, 3, 0, -2, 0},

		/* 511-520 */
		{0, 3, 0, -2, 0},
		{2, -2, 2, -2, 2},
		{0, 0, 4, -2, 0},
		{4, -1, 2, 0, 2},
		{2, 2, -2, -4, -2},
		{4, 1, 2, 0, 2},
		{4, -1, -2, -2, -2},
		{2, 1, 0, -2, -2},
		{2, 1, -2, -6, -2},
		{2, 0, 0, -1, 1},

		/* 521-530 */
		{2, -1, -2, 2, -1},
		{1, 1, -2, 2, -2},
		{1, 1, -2, -3, -2},
		{1, 0, 3, 0, 3},
		{1, 0, -2, 1, 1},
		{1, 0, -2, 0, 2},
		{1, -1, 2, 1, 2},
		{1, -1, 0, 0, -2},
		{1, -1, -4, 2, -2},
		{0, 3, -2, -2, -2},

		/* 531-540 */
		{0, 1, 0, 4, 1},
		{0, 0, 4, 2, 2},
		{3, 0, -2, -2, 0},
		{2, -2, 0, 0, 0},
		{1, 1, 2, -4, 0},
		{1, 1, 0, -3, 0},
		{1, 0, 2, -3, 0},
		{1, -1, 2, -2, 0},
		{0, 2, 0, 2, 0},
		{0, 0, 2, 4, 0},

		/* 541-550 */
		{1, 0, 1, 0, 0},
		{3, 1, 2, -2, 1},
		{3, 0, 4, -2, 2},
		{3, 0, 2, 1, 2},
		{3, 0, 0, 2, -1},
		{3, 0, 0, 0, 2},
		{3, 0, -2, 2, -1},
		{2, 0, 4, -4, 2},
		{2, 0, 2, -3, 2},
		{2, 0, 0, 4, 1},

		/* 551-560 */
		{2, 0, 0, -3, 1},
		{2, 0, -4, 2, -1},
		{2, 0, -2, -2, 1},
		{2, -2, 2, 2, 2},
		{2, -2, 0, -2, -2},
		{2, -1, 0, 2, 1},
		{2, -1, 0, 2, -1},
		{1, 1, 2, 4, 2},
		{1, 1, 0, 1, 1},
		{1, 1, 0, 1, -1},

		/* 561-570 */
		{1, 1, -2, -6, -2},
		{1, 0, 0, -3, -1},
		{1, 0, -4, -2, -1},
		{1, 0, -2, -6, -1},
		{1, -2, 2, 2, 1},
		{1, -2, -2, 2, -1},
		{1, -1, -2, -4, -1},
		{0, 2, 0, 0, 2},
		{0, 1, 2, -4, 2},
		{0, 1, -2, 4, -1},

		/* 571-580 */
		{5, 0, 0, 0, 0},
		{3, 0, 0, -3, 0},
		{2, 2, 0, -4, 0},
		{1, -1, 2, 2, 0},
		{0, 1, 0, 3, 0},
		{4, 0, -2, 0, -1},
		{3, 0, -2, -6, -1},
		{3, 0, -2, -1, -1},
		{2, 1, 2, 2, 1},
		{2, 1, 0, 2, 1},

		/* 581-590 */
		{2, 0, 2, 4, 1},
		{2, 0, 2, -6, 1},
		{2, 0, 2, -2, -1},
		{2, 0, 0, -6, -1},
		{2, -1, -2, -2, -1},
		{1, 2, 2, 0, 1},
		{1, 2, 0, 0, 1},
		{1, 0, 4, 0, 1},
		{1, 0, 2, -6, 1},
		{1, 0, 2, -4, -1},

		/* 591-600 */
		{1, 0, -1, -2, -1},
		{1, -1, 2, 4, 1},
		{1, -1, 2, -3, 1},
		{1, -1, 0, 4, 1},
		{1, -1, -2, 1, -1},
		{0, 1, 2, -2, 3},
		{3, 0, 0, -2, 0},
		{1, 0, 1, -2, 0},
		{0, 2, 0, -4, 0},
		{0, 0, 2, -4, 0},

		/* 601-610 */
		{0, 0, 1, -1, 0},
		{0, 0, 0, 6, 0},
		{0, 2, 0, 0, -2},
		{0, 1, -2, 2, -3},
		{4, 0, 0, 2, 0},
		{3, 0, 0, -1, 0},
		{3, -1, 0, 2, 0},
		{2, 1, 0, 1, 0},
		{2, 1, 0, -6, 0},
		{2, -1, 2, 0, 0},

		/* 611-620 */
		{1, 0, 2, -1, 0},
		{1, -1, 0, 1, 0},
		{1, -1, -2, -2, 0},
		{0, 1, 2, 2, 0},
		{0, 0, 2, -3, 0},
		{2, 2, 0, -2, -1},
		{2, -1, -2, 0, 1},
		{1, 2, 2, -4, 1},
		{0, 1, 4, -4, 2},
		{0, 0, 0, 3, 2},

		/* 621-630 */
		{5, 0, 2, 0, 1},
		{4, 1, 2, -2, 2},
		{4, 0, -2, -2, 0},
		{3, 1, 2, 2, 2},
		{3, 1, 0, -2, 0},
		{3, 1, -2, -6, -2},
		{3, 0, 0, 0, -2},
		{3, 0, -2, -4, -2},
		{3, -1, 0, -3, 0},
		{3, -1, 0, -2, 0},

		/* 631-640 */
		{2, 1, 2, 0, 0},
		{2, 1, 2, -4, 2},
		{2, 1, 2, -2, 0},
		{2, 1, 0, -3, 0},
		{2, 1, -2, 0, -2},
		{2, 0, 0, -4, 2},
		{2, 0, 0, -4, -2},
		{2, 0, -2, -5, -2},
		{2, -1, 2, 4, 2},
		{2, -1, 0, -2, 2},

		/* 641-650 */
		{1, 3, -2, -2, -2},
		{1, 1, 0, 0, -2},
		{1, 1, 0, -6, 0},
		{1, 1, -2, 1, -2},
		{1, 1, -2, -1, -2},
		{1, 0, 2, 1, 0},
		{1, 0, 0, 3, 0},
		{1, 0, 0, -4, 2},
		{1, 0, -2, 4, -2},
		{1, -2, 0, -1, 0},

		/* 651-NFLS */
		{0, 1, -4, 2, -1},
		{1, 0, -2, 0, -3},
		{0, 0, 4, -4, 4},
	}

	/* Number of frequencies:  luni-solar */
	NFLS := len(mfals)

	/* Fundamental-argument multipliers:  planetary terms */
	mfapl := [][14]int{

		/* 1-10 */
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -1, 2, 0, 0, 0, 0, 0},

		/* 11-20 */
		{0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 2, -5, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, -8, 3, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 6, -8, 3, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 0},

		/* 21-30 */
		{0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2},
		{0, 0, 1, -1, 1, 0, 0, 0, -2, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, -1, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 2, 0, 0, 0, 0},

		/* 31-40 */
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 1},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0},

		/* 41-50 */
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, -2, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 2},
		{1, 0, 0, 0, 0, 0, -18, 16, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2},

		/* 51-60 */
		{0, 0, 1, -1, 1, 0, -5, 7, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, -10, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, 0, -5, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 2},
		{1, 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1},
		{1, 0, -2, 0, -2, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 2, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0},

		/* 61-70 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, 0, 3, -8, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -11, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0},

		/* 71-80 */
		{0, 0, 0, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 2},
		{0, 0, 1, -1, 1, 0, 0, -5, 8, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 0},

		/* 81-90 */
		{2, 0, 0, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, -1},
		{2, 0, 0, -2, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 8, -13, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 5, 0, 0, 0},
		{1, 0, 0, -1, 0, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2},
		{1, 0, 0, 0, -1, 0, -18, 16, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 2, -5, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0},

		/* 91-100 */
		{1, 0, 0, -2, 0, 0, 19, -21, 3, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -8, 13, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -9, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2},
		{1, 0, 0, 0, 1, 0, -18, 16, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 6, -16, 4, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 3, -7, 0, 0, 0, 0, 0, -2},

		/* 101-110 */
		{0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
		{2, 0, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, -1, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2},

		/* 111-120 */
		{0, 0, 0, 0, 1, 0, 0, 1, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2},
		{0, 0, 2, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, -6, 8, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0},

		/* 121-130 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -10, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, -2},
		{1, 0, 0, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0},

		/* 131-140 */
		{0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, -3, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0, 2, -3, 0, 0, 0, 0, 0, 0},

		/* 141-150 */
		{1, 0, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 9, -11, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -4, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, -1, 0, 0, 0, 2},

		/* 151-160 */
		{1, 0, 0, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -4, 10, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, -2},
		{0, 0, 2, -2, 1, 0, -4, 4, 0, 0, 0, 0, 0, 0},

		/* 161-170 */
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 2},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 2, 0},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -9, 13, 0, 0, 0, 0, 0},
		{2, 0, 2, 0, 2, 0, 0, 2, 0, -3, 0, 0, 0, 0},

		/* 171-180 */
		{0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, 0, 2, 0, 0, 0},
		{1, 0, 0, -1, -1, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1},
		{1, 0, 2, 0, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{1, 0, -2, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -2, 4, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0},

		/* 181-190 */
		{0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1},
		{0, 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, -8, 3, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 6, -10, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 7, -8, 3, 0, 0, 0, 2},
		{0, 0, 0, 0, 1, 0, -3, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, -5, 7, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 1},

		/* 191-200 */
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -10, 0, 0, 0, 0, 0, -2},
		{1, 0, 0, -2, 0, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 2, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, 1, 0, 0, -9, 15, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0},

		/* 201-210 */
		{0, 0, 0, 0, 0, 0, 0, 1, -4, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, 2},
		{2, 0, 0, -2, 1, 0, -6, 8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, 1, 0, 3, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 8, -14, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},

		/* 211-220 */
		{0, 0, 0, 0, 1, 0, 0, 8, -15, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2},
		{2, 0, -1, -1, 0, 0, 0, 3, -7, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0},

		/* 221-230 */
		{2, 0, 0, -2, 0, 0, 0, -6, 8, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -5, 6, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -9, 4, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, -2},

		/* 231-240 */
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -4, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1},
		{0, 0, 0, 0, 0, 0, 7, -11, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 3, -5, 4, 0, 0, 0, 0, 2},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, -1, 1, 0, 0, 0},
		{2, 0, 0, 0, 0, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, -1},

		/* 241-250 */
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -8, 3, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 2, -4, 0, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 3, -5, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 2},
		{0, 0, 2, -2, 2, 0, -8, 11, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 0, 0, 0},

		/* 251-260 */
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 7, -9, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, -2, -2, -2, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 1},

		/* 261-270 */
		{0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 2, -5, 0, 0, 2},
		{2, 0, 0, -2, -1, 0, 0, -2, 0, 0, 5, 0, 0, 0},
		{2, 0, 0, -2, -1, 0, -6, 8, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 2, -5, 0, 0, 2},
		{0, 0, 0, 0, 1, 0, 3, -7, 4, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0},

		/* 271-280 */
		{0, 0, 1, -1, 0, 0, 0, -1, 0, -2, 5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 6, -15, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 2},
		{1, 0, 0, -1, 0, 0, 0, -3, 4, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -3, 7, -4, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, 0, -2, 0, 0, 0, 2},

		/* 281-290 */
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 1},
		{0, 0, 2, -2, 2, 0, -5, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 6, -11, 0, 0, 0, 0, -2},

		/* 291-300 */
		{0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 3, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 9, -12, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 1},
		{0, 0, 1, -1, 0, 0, -8, 12, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, -1},

		/* 301-310 */
		{0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, -1, 0, 0, 0, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, -5, 0, 0, 0, 0, -2},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 3, -1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2},

		/* 311-320 */
		{0, 0, 0, 0, 0, 0, 9, -9, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 3, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 2, -4, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2},
		{0, 0, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0},

		/* 321-330 */
		{0, 0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, 0, -3, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
		{2, 0, -1, -1, -1, 0, 0, -1, 0, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 5, -10, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 1},
		{0, 0, 2, -2, 1, -1, 0, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 2, 0, 0},

		/* 331-340 */
		{0, 0, 0, 0, 1, 0, 3, -5, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 9, -9, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -8, 11, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -2, 0, 0, 2, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 2, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 2, -6, 0, 0, 0, 0, 0, -2},

		/* 341-350 */
		{0, 0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 5, -2, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 7, -13, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2},
		{0, 0, 2, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -8, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 8, -10, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 1},

		/* 351-360 */
		{0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -4, 0, 0, 0, 0},
		{2, 0, 0, -2, -1, 0, 0, -5, 6, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, -2},
		{2, 0, -1, -1, -1, 0, 0, 3, -7, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0},

		/* 361-370 */
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 4, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 6, -11, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 1, 0, 0, -6, 8, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 1, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 6, -5, 0, 0, 0, 0, 2},
		{1, 0, -2, -2, -2, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{0, 0, 0, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1},

		/* 371-380 */
		{0, 0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, 0, -2, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -2, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 1, -6, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2},
		{0, 0, 0, 0, 0, 0, 3, -5, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 7, -13, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 2},

		/* 381-390 */
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 2, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -8, 15, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, -2, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{2, 0, -1, -1, -1, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{1, 0, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{1, 0, -1, 1, -1, 0, -18, 17, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 2, -2, -1, 0, -5, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0},

		/* 391-400 */
		{0, 0, 0, 0, 1, 0, 2, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -16, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2},
		{0, 0, 0, 0, 2, 0, 0, -1, 2, 0, 0, 0, 0, 0},
		{2, 0, -1, -1, -2, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -10, 0, 0, 0, 0, 0, -1},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 4, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 2},
		{2, 0, 0, -2, -1, 0, 0, -2, 0, 4, -5, 0, 0, 0},

		/* 401-410 */
		{2, 0, 0, -2, -1, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{2, 0, -1, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -1, -1, 0, 0, -2, 2, 0, 0, 0, 0, 0},
		{1, 0, -1, -1, -1, 0, 20, -20, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 1, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -2, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 5, -8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0},

		/* 411-420 */
		{0, 0, 0, 0, 0, 0, 9, -11, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1},
		{0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -2},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, -2, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, -2, 5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0},

		/* 421-430 */
		{0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, -6, 0, 0, 0, 0, -2},
		{1, 0, 0, -2, 0, 0, 20, -21, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -12, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 8, -12, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 9, -17, 0, 0, 0, 0, 0},

		/* 431-440 */
		{0, 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 1, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, -7, 0, 0, 0, 0, -2},
		{1, 0, 0, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0},
		{1, 0, -2, 0, -2, 0, -10, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -9, 17, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -4, 0, 0, 0, 0, 0, -2},
		{1, 0, -2, -2, -2, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{1, 0, -1, 1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0},

		/* 441-450 */
		{0, 0, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, 0, -1, 0, 0, 1, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, -5, 7, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -10, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, -5, 0, 0, 0, -2},

		/* 451-460 */
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -5, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 1},
		{1, 0, 0, -2, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -7, 4, 0, 0, 0, 0, 0},
		{2, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 1, 0, -2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 6, -10, 0, 0, 0, 0, -2},

		/* 461-470 */
		{1, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -3, 0, 3, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, -5, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 1, -3, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -4, 6, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, -1, 0, 0},
		{0, 0, 1, -1, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0},

		/* 471-480 */
		{0, 0, 0, 0, 1, 0, 3, -4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -10, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 3, -8, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 7, -9, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 2},

		/* 481-490 */
		{0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 3, -8, 3, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1},
		{2, 0, 0, -2, -1, 0, 0, -6, 8, 0, 0, 0, 0, 0},
		{2, 0, -1, -1, 1, 0, 0, 3, -7, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -7, 9, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, -1},

		/* 491-500 */
		{0, 0, 1, -1, 2, 0, -8, 12, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 1, 0, 0, -5, 6, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, -1, 0, 0, -2, 0, 3, -1, 0, 0, 0},
		{1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{1, 0, 0, -2, -1, 0, 0, -2, 0, 2, 0, 0, 0, 0},

		/* 501-510 */
		{1, 0, 0, -1, -1, 0, 0, -3, 4, 0, 0, 0, 0, 0},
		{1, 0, -1, 0, -1, 0, -3, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -4, 4, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, -8, 11, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, 0, 0, -9, 13, 0, 0, 0, 0, 0},
		{0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, 1, -4, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, 0, -1, 0, 1, -3, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 7, -13, 0, 0, 0, 0, 0},

		/* 511-520 */
		{0, 0, 0, 0, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 7, -11, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 1},

		/* 521-530 */
		{0, 0, 0, 0, 0, 0, 1, -4, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 9, -17, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 4, -8, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -4, 8, -3, 0, 0, 0, 0},

		/* 531-540 */
		{2, 0, 0, -2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 17, -16, 0, -2, 0, 0, 0, 0},
		{1, 0, 0, -1, 0, 0, 0, -2, 2, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -4, 0, 0, 0, 0},

		/* 541-550 */
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2},
		{2, 0, 0, -2, 0, 0, 0, -4, 4, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 2, 2, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 0, -4, 8, -3, 0, 0, 0, 0},

		/* 551-560 */
		{1, 0, 0, -2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, 0, -4, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -2, 2, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, -4, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 2, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0},

		/* 561-570 */
		{0, 0, 0, 0, 0, 0, 8, -9, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0},
		{2, 0, -2, -2, -2, 0, 0, -2, 0, 2, 0, 0, 0, 0},
		{1, 0, 0, 0, 1, 0, -10, 3, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, -1, 0, -10, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, 2, -3, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 0, 0, 0},

		/* 571-580 */
		{0, 0, 2, 0, 2, 0, -2, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, 0, 2, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0},
		{2, 0, 2, -2, 2, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{2, 0, 1, -3, 1, 0, -6, 7, 0, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 2, -5, 0, 0, 0, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 5, -5, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 1, 5, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 0, 5, 0, 0, 0},

		/* 581-590 */
		{2, 0, 0, -2, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0},
		{2, 0, 0, -2, 0, 0, -4, 4, 0, 0, 0, 0, 0, 0},
		{2, 0, -2, 0, -2, 0, 0, 5, -9, 0, 0, 0, 0, 0},
		{2, 0, -1, -1, 0, 0, 0, -1, 0, 3, 0, 0, 0, 0},
		{1, 0, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{1, 0, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0},
		{1, 0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0},
		{1, 0, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0},
		{1, 0, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0},

		/* 591-600 */
		{1, 0, 0, 0, 0, 0, 0, -2, 0, 3, 0, 0, 0, 0},
		{1, 0, 0, -2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0},
		{1, 0, -2, -2, -2, 0, 0, 1, 0, -1, 0, 0, 0, 0},
		{1, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{1, 0, -1, -1, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0},
		{0, 0, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -2, 0, 1, 0, 0, 0, 0},
		{0, 0, 2, -2, 1, 0, 0, -10, 15, 0, 0, 0, 0, 0},
		{0, 0, 2, -2, 0, -1, 0, 2, 0, 0, 0, 0, 0, 0},

		/* 601-610 */
		{0, 0, 1, -1, 2, 0, 0, -1, 0, 0, -1, 0, 0, 0},
		{0, 0, 1, -1, 2, 0, -3, 4, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -4, 6, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, 0, -1, 0, 0, -2, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, -1, -1, 0, -5, 7, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0},

		/* 611-620 */
		{0, 0, 0, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 2, 0, -3, 5, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 9, -13, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 8, -14, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 8, -11, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, -2},

		/* 621-630 */
		{0, 0, 0, 0, 0, 0, 5, -6, -4, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 4, -8, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 3, -3, 0, 2, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 7, -12, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, -2},

		/* 631-640 */
		{0, 0, 0, 0, 0, 0, 0, 6, -8, 1, 5, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 6, -10, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, 0, -4, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, -1},
		{0, 0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 5, -16, 4, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 5, -13, 0, 0, 0, 0, -2},

		/* 641-650 */
		{0, 0, 0, 0, 0, 0, 0, 3, 0, -5, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, -9, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 3, -7, 0, 0, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, -8, 1, 5, 0, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 1, -5, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -3, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 5, 0, 0, 0},

		/* 651-NFPL */
		{0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -6, 3, 0, -2},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
	}

	/* Number of frequencies:  planetary */
	NFPL := len(mfapl)

	/* Pointers into amplitudes array, one pointer per frequency */
	nc := []int{

		/* 1-100 */
		1, 21, 37, 51, 65, 79, 91, 103, 115, 127,
		139, 151, 163, 172, 184, 196, 207, 219, 231, 240,
		252, 261, 273, 285, 297, 309, 318, 327, 339, 351,
		363, 372, 384, 396, 405, 415, 423, 435, 444, 452,
		460, 467, 474, 482, 490, 498, 506, 513, 521, 528,
		536, 543, 551, 559, 566, 574, 582, 590, 597, 605,
		613, 620, 628, 636, 644, 651, 658, 666, 674, 680,
		687, 695, 702, 710, 717, 725, 732, 739, 746, 753,
		760, 767, 774, 782, 790, 798, 805, 812, 819, 826,
		833, 840, 846, 853, 860, 867, 874, 881, 888, 895,

		/* 101-200 */
		901, 908, 914, 921, 928, 934, 941, 948, 955, 962,
		969, 976, 982, 989, 996, 1003, 1010, 1017, 1024, 1031,
		1037, 1043, 1050, 1057, 1064, 1071, 1078, 1084, 1091, 1098,
		1104, 1112, 1118, 1124, 1131, 1138, 1145, 1151, 1157, 1164,
		1171, 1178, 1185, 1192, 1199, 1205, 1212, 1218, 1226, 1232,
		1239, 1245, 1252, 1259, 1266, 1272, 1278, 1284, 1292, 1298,
		1304, 1310, 1316, 1323, 1329, 1335, 1341, 1347, 1353, 1359,
		1365, 1371, 1377, 1383, 1389, 1396, 1402, 1408, 1414, 1420,
		1426, 1434, 1440, 1446, 1452, 1459, 1465, 1471, 1477, 1482,
		1488, 1493, 1499, 1504, 1509, 1514, 1520, 1527, 1532, 1538,

		/* 201-300 */
		1543, 1548, 1553, 1558, 1564, 1569, 1574, 1579, 1584, 1589,
		1594, 1596, 1598, 1600, 1602, 1605, 1608, 1610, 1612, 1617,
		1619, 1623, 1625, 1627, 1629, 1632, 1634, 1640, 1642, 1644,
		1646, 1648, 1650, 1652, 1654, 1658, 1660, 1662, 1664, 1668,
		1670, 1672, 1673, 1675, 1679, 1681, 1683, 1684, 1686, 1688,
		1690, 1693, 1695, 1697, 1701, 1703, 1705, 1707, 1709, 1711,
		1712, 1715, 1717, 1721, 1723, 1725, 1727, 1729, 1731, 1733,
		1735, 1737, 1739, 1741, 1743, 1745, 1747, 1749, 1751, 1753,
		1755, 1757, 1759, 1761, 1762, 1764, 1766, 1768, 1769, 1771,
		1773, 1775, 1777, 1779, 1781, 1783, 1785, 1787, 1788, 1790,

		/* 301-400 */
		1792, 1794, 1796, 1798, 1800, 1802, 1804, 1806, 1807, 1809,
		1811, 1815, 1817, 1819, 1821, 1823, 1825, 1827, 1829, 1831,
		1833, 1835, 1837, 1839, 1840, 1842, 1844, 1848, 1850, 1852,
		1854, 1856, 1858, 1859, 1860, 1862, 1864, 1866, 1868, 1869,
		1871, 1873, 1875, 1877, 1879, 1881, 1883, 1885, 1887, 1889,
		1891, 1892, 1896, 1898, 1900, 1901, 1903, 1905, 1907, 1909,
		1910, 1911, 1913, 1915, 1919, 1921, 1923, 1927, 1929, 1931,
		1933, 1935, 1937, 1939, 1943, 1945, 1947, 1948, 1949, 1951,
		1953, 1955, 1957, 1958, 1960, 1962, 1964, 1966, 1968, 1970,
		1971, 1973, 1974, 1975, 1977, 1979, 1980, 1981, 1982, 1984,

		/* 401-500 */
		1986, 1988, 1990, 1992, 1994, 1995, 1997, 1999, 2001, 2003,
		2005, 2007, 2008, 2009, 2011, 2013, 2015, 2017, 2019, 2021,
		2023, 2024, 2025, 2027, 2029, 2031, 2033, 2035, 2037, 2041,
		2043, 2045, 2046, 2047, 2049, 2051, 2053, 2055, 2056, 2057,
		2059, 2061, 2063, 2065, 2067, 2069, 2070, 2071, 2072, 2074,
		2076, 2078, 2080, 2082, 2084, 2086, 2088, 2090, 2092, 2094,
		2095, 2096, 2097, 2099, 2101, 2105, 2106, 2107, 2108, 2109,
		2110, 2111, 2113, 2115, 2119, 2121, 2123, 2125, 2127, 2129,
		2131, 2133, 2135, 2136, 2137, 2139, 2141, 2143, 2145, 2147,
		2149, 2151, 2153, 2155, 2157, 2159, 2161, 2163, 2165, 2167,

		/* 501-600 */
		2169, 2171, 2173, 2175, 2177, 2179, 2181, 2183, 2185, 2186,
		2187, 2188, 2192, 2193, 2195, 2197, 2199, 2201, 2203, 2205,
		2207, 2209, 2211, 2213, 2217, 2219, 2221, 2223, 2225, 2227,
		2229, 2231, 2233, 2234, 2235, 2236, 2237, 2238, 2239, 2240,
		2241, 2244, 2246, 2248, 2250, 2252, 2254, 2256, 2258, 2260,
		2262, 2264, 2266, 2268, 2270, 2272, 2274, 2276, 2278, 2280,
		2282, 2284, 2286, 2288, 2290, 2292, 2294, 2296, 2298, 2300,
		2302, 2303, 2304, 2305, 2306, 2307, 2309, 2311, 2313, 2315,
		2317, 2319, 2321, 2323, 2325, 2327, 2329, 2331, 2333, 2335,
		2337, 2341, 2343, 2345, 2347, 2349, 2351, 2352, 2355, 2356,

		/* 601-700 */
		2357, 2358, 2359, 2361, 2363, 2364, 2365, 2366, 2367, 2368,
		2369, 2370, 2371, 2372, 2373, 2374, 2376, 2378, 2380, 2382,
		2384, 2385, 2386, 2387, 2388, 2389, 2390, 2391, 2392, 2393,
		2394, 2395, 2396, 2397, 2398, 2399, 2400, 2401, 2402, 2403,
		2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412, 2413,
		2414, 2415, 2417, 2418, 2430, 2438, 2445, 2453, 2460, 2468,
		2474, 2480, 2488, 2496, 2504, 2512, 2520, 2527, 2535, 2543,
		2550, 2558, 2566, 2574, 2580, 2588, 2596, 2604, 2612, 2619,
		2627, 2634, 2642, 2648, 2656, 2664, 2671, 2679, 2685, 2693,
		2701, 2709, 2717, 2725, 2733, 2739, 2747, 2753, 2761, 2769,

		/* 701-800 */
		2777, 2785, 2793, 2801, 2809, 2817, 2825, 2833, 2841, 2848,
		2856, 2864, 2872, 2878, 2884, 2892, 2898, 2906, 2914, 2922,
		2930, 2938, 2944, 2952, 2958, 2966, 2974, 2982, 2988, 2996,
		3001, 3009, 3017, 3025, 3032, 3039, 3045, 3052, 3059, 3067,
		3069, 3076, 3083, 3090, 3098, 3105, 3109, 3111, 3113, 3120,
		3124, 3128, 3132, 3136, 3140, 3144, 3146, 3150, 3158, 3161,
		3165, 3166, 3168, 3172, 3176, 3180, 3182, 3185, 3189, 3193,
		3194, 3197, 3200, 3204, 3208, 3212, 3216, 3219, 3221, 3222,
		3226, 3230, 3234, 3238, 3242, 3243, 3247, 3251, 3254, 3258,
		3262, 3266, 3270, 3274, 3275, 3279, 3283, 3287, 3289, 3293,

		/* 801-900 */
		3296, 3300, 3303, 3307, 3311, 3315, 3319, 3321, 3324, 3327,
		3330, 3334, 3338, 3340, 3342, 3346, 3350, 3354, 3358, 3361,
		3365, 3369, 3373, 3377, 3381, 3385, 3389, 3393, 3394, 3398,
		3402, 3406, 3410, 3413, 3417, 3421, 3425, 3429, 3433, 3435,
		3439, 3443, 3446, 3450, 3453, 3457, 3458, 3461, 3464, 3468,
		3472, 3476, 3478, 3481, 3485, 3489, 3493, 3497, 3501, 3505,
		3507, 3511, 3514, 3517, 3521, 3524, 3525, 3527, 3529, 3533,
		3536, 3540, 3541, 3545, 3548, 3551, 3555, 3559, 3563, 3567,
		3569, 3570, 3574, 3576, 3578, 3582, 3586, 3590, 3593, 3596,
		3600, 3604, 3608, 3612, 3616, 3620, 3623, 3626, 3630, 3632,

		/* 901-1000 */
		3636, 3640, 3643, 3646, 3648, 3652, 3656, 3660, 3664, 3667,
		3669, 3671, 3675, 3679, 3683, 3687, 3689, 3693, 3694, 3695,
		3699, 3703, 3705, 3707, 3710, 3713, 3717, 3721, 3725, 3729,
		3733, 3736, 3740, 3744, 3748, 3752, 3754, 3757, 3759, 3763,
		3767, 3770, 3773, 3777, 3779, 3783, 3786, 3790, 3794, 3798,
		3801, 3805, 3809, 3813, 3817, 3821, 3825, 3827, 3831, 3835,
		3836, 3837, 3840, 3844, 3848, 3852, 3856, 3859, 3863, 3867,
		3869, 3871, 3875, 3879, 3883, 3887, 3890, 3894, 3898, 3901,
		3905, 3909, 3913, 3917, 3921, 3922, 3923, 3924, 3926, 3930,
		3932, 3936, 3938, 3940, 3944, 3948, 3952, 3956, 3959, 3963,

		/* 1001-1100 */
		3965, 3969, 3973, 3977, 3979, 3981, 3982, 3986, 3989, 3993,
		3997, 4001, 4004, 4006, 4009, 4012, 4016, 4020, 4024, 4026,
		4028, 4032, 4036, 4040, 4044, 4046, 4050, 4054, 4058, 4060,
		4062, 4063, 4064, 4068, 4071, 4075, 4077, 4081, 4083, 4087,
		4089, 4091, 4095, 4099, 4101, 4103, 4105, 4107, 4111, 4115,
		4119, 4123, 4127, 4129, 4131, 4135, 4139, 4141, 4143, 4145,
		4149, 4153, 4157, 4161, 4165, 4169, 4173, 4177, 4180, 4183,
		4187, 4191, 4195, 4198, 4201, 4205, 4209, 4212, 4213, 4216,
		4217, 4221, 4223, 4226, 4230, 4234, 4236, 4240, 4244, 4248,
		4252, 4256, 4258, 4262, 4264, 4266, 4268, 4270, 4272, 4276,

		/* 1101-1200 */
		4279, 4283, 4285, 4287, 4289, 4293, 4295, 4299, 4300, 4301,
		4305, 4309, 4313, 4317, 4319, 4323, 4325, 4329, 4331, 4333,
		4335, 4337, 4341, 4345, 4349, 4351, 4353, 4357, 4361, 4365,
		4367, 4369, 4373, 4377, 4381, 4383, 4387, 4389, 4391, 4395,
		4399, 4403, 4407, 4411, 4413, 4414, 4415, 4418, 4419, 4421,
		4423, 4427, 4429, 4431, 4433, 4435, 4437, 4439, 4443, 4446,
		4450, 4452, 4456, 4458, 4460, 4462, 4466, 4469, 4473, 4477,
		4481, 4483, 4487, 4489, 4491, 4493, 4497, 4499, 4501, 4504,
		4506, 4510, 4513, 4514, 4515, 4518, 4521, 4522, 4525, 4526,
		4527, 4530, 4533, 4534, 4537, 4541, 4542, 4543, 4544, 4545,

		/* 1201-1300 */
		4546, 4547, 4550, 4553, 4554, 4555, 4558, 4561, 4564, 4567,
		4568, 4571, 4574, 4575, 4578, 4581, 4582, 4585, 4586, 4588,
		4590, 4592, 4596, 4598, 4602, 4604, 4608, 4612, 4613, 4616,
		4619, 4622, 4623, 4624, 4625, 4626, 4629, 4632, 4633, 4636,
		4639, 4640, 4641, 4642, 4643, 4644, 4645, 4648, 4649, 4650,
		4651, 4652, 4653, 4656, 4657, 4660, 4661, 4664, 4667, 4670,
		4671, 4674, 4675, 4676, 4677, 4678, 4681, 4682, 4683, 4684,
		4687, 4688, 4689, 4692, 4693, 4696, 4697, 4700, 4701, 4702,
		4703, 4704, 4707, 4708, 4711, 4712, 4715, 4716, 4717, 4718,
		4719, 4720, 4721, 4722, 4723, 4726, 4729, 4730, 4733, 4736,

		/* 1301-(NFLS+NFPL) */
		4737, 4740, 4741, 4742, 4745, 4746, 4749, 4752, 4753,
	}

	/* Amplitude coefficients (microarcsec);  indexed using the nc array. */
	a := []float64{

		/* 1-105 */
		-6844318.44, 9205236.26, 1328.67, 1538.18, 205833.11,
		153041.79, -3309.73, 853.32, 2037.98, -2301.27,
		81.46, 120.56, -20.39, -15.22, 1.73, -1.61, -0.10, 0.11,
		-0.02, -0.02, -523908.04, 573033.42, -544.75, -458.66,
		12814.01, 11714.49, 198.97, -290.91, 155.74, -143.27,
		-2.75, -1.03, -1.27, -1.16, 0.00, -0.01, -90552.22,
		97846.69, 111.23, 137.41, 2187.91, 2024.68, 41.44, -51.26,
		26.92, -24.46, -0.46, -0.28, -0.22, -0.20, 82168.76,
		-89618.24, -27.64, -29.05, -2004.36, -1837.32,
		-36.07, 48.00, -24.43, 22.41, 0.47, 0.24, 0.20, 0.18,
		58707.02, 7387.02, 470.05, -192.40, 164.33, -1312.21,
		-179.73, -28.93, -17.36, -1.83, -0.50, 3.57, 0.00, 0.13,
		-20557.78, 22438.42, -20.84, -17.40, 501.82, 459.68,
		59.20, -67.30, 6.08, -5.61, -1.36, -1.19, 28288.28,
		-674.99, -34.69, 35.80, -15.07, -632.54, -11.19, 0.78, -8.41,
		0.17, 0.01, 0.07, -15406.85, 20069.50, 15.12,

		/* 106-219 */
		31.80, 448.76, 344.50, -5.77, 1.41, 4.59, -5.02, 0.17,
		0.24, -11991.74, 12902.66, 32.46, 36.70, 288.49,
		268.14, 5.70, -7.06, 3.57, -3.23, -0.06, -0.04,
		-8584.95, -9592.72, 4.42, -13.20, -214.50, 192.06,
		23.87, 29.83, 2.54, 2.40, 0.60, -0.48, 5095.50,
		-6918.22, 7.19, 3.92, -154.91, -113.94, 2.86, -1.04,
		-1.52, 1.73, -0.07, -0.10, -4910.93, -5331.13,
		0.76, 0.40, -119.21, 109.81, 2.16, 3.20, 1.46, 1.33,
		0.04, -0.02, -6245.02, -123.48, -6.68, -8.20, -2.76,
		139.64, 2.71, 0.15, 1.86, 2511.85, -3323.89, 1.07,
		-0.90, -74.33, -56.17, 1.16, -0.01, -0.75, 0.83, -0.02,
		-0.04, 2307.58, 3143.98, -7.52, 7.50, 70.31, -51.60, 1.46,
		0.16, -0.69, -0.79, 0.02, -0.05, 2372.58, 2554.51, 5.93,
		-6.60, 57.12, -53.05, -0.96, -1.24, -0.71, -0.64, -0.01,
		-2053.16, 2636.13, 5.13, 7.80, 58.94, 45.91, -0.42,
		-0.12, 0.61, -0.66, 0.02, 0.03, -1825.49,

		/* 220-339 */
		-2423.59, 1.23, -2.00, -54.19, 40.82, -1.07, -1.02,
		0.54, 0.61, -0.04, 0.04, 2521.07, -122.28, -5.97, 2.90,
		-2.73, -56.37, -0.82, 0.13, -0.75, -1534.09, 1645.01,
		6.29, 6.80, 36.78, 34.30, 0.92, -1.25, 0.46, -0.41,
		-0.02, -0.01, 1898.27, 47.70, -0.72, 2.50, 1.07, -42.45,
		-0.94, 0.02, -0.56, -1292.02, -1387.00, 0.00,
		0.00, -31.01, 28.89, 0.68, 0.00, 0.38, 0.35, -0.01,
		-0.01, -1234.96, 1323.81, 5.21, 5.90, 29.60, 27.61,
		0.74, -1.22, 0.37, -0.33, -0.02, -0.01, 1137.48,
		-1233.89, -0.04, -0.30, -27.59, -25.43, -0.61, 1.00,
		-0.34, 0.31, 0.01, 0.01, -813.13, -1075.60, 0.40,
		0.30, -24.05, 18.18, -0.40, -0.01, 0.24, 0.27, -0.01,
		0.01, 1163.22, -60.90, -2.94, 1.30, -1.36, -26.01, -0.58,
		0.07, -0.35, 1029.70, -55.55, -2.63, 1.10, -1.25, -23.02,
		-0.52, 0.06, -0.31, -556.26, 852.85, 3.16, -4.48, 19.06,
		12.44, -0.81, -0.27, 0.17, -0.21, 0.00, 0.02, -603.52,

		/* 340-467 */
		-800.34, 0.44, 0.10, -17.90, 13.49, -0.08, -0.01, 0.18,
		0.20, -0.01, 0.01, -628.24, 684.99, -0.64, -0.50, 15.32,
		14.05, 3.18, -4.19, 0.19, -0.17, -0.09, -0.07, -866.48,
		-16.26, 0.52, -1.30, -0.36, 19.37, 0.43, -0.01, 0.26,
		-512.37, 695.54, -1.47, -1.40, 15.55, 11.46, -0.16, 0.03,
		0.15, -0.17, 0.01, 0.01, 506.65, 643.75, 2.54, -2.62,
		14.40, -11.33, -0.77, -0.06, -0.15, -0.16, 0.00, 0.01,
		664.57, 16.81, -0.40, 1.00, 0.38, -14.86, -3.71, -0.09,
		-0.20, 405.91, 522.11, 0.99, -1.50, 11.67, -9.08, -0.25,
		-0.02, -0.12, -0.13, -305.78, 326.60, 1.75, 1.90, 7.30,
		6.84, 0.20, -0.04, 300.99, -325.03, -0.44, -0.50, -7.27,
		-6.73, -1.01, 0.01, 0.00, 0.08, 0.00, 0.02, 438.51,
		10.47, -0.56, -0.20, 0.24, -9.81, -0.24, 0.01, -0.13,
		-264.02, 335.24, 0.99, 1.40, 7.49, 5.90, -0.27, -0.02,
		284.09, 307.03, 0.32, -0.40, 6.87, -6.35, -0.99, -0.01,
		-250.54, 327.11, 0.08, 0.40, 7.31, 5.60, -0.30, 230.72,

		/* 468-595 */
		-304.46, 0.08, -0.10, -6.81, -5.16, 0.27, 229.78, 304.17,
		-0.60, 0.50, 6.80, -5.14, 0.33, 0.01, 256.30, -276.81,
		-0.28, -0.40, -6.19, -5.73, -0.14, 0.01, -212.82, 269.45,
		0.84, 1.20, 6.02, 4.76, 0.14, -0.02, 196.64, 272.05,
		-0.84, 0.90, 6.08, -4.40, 0.35, 0.02, 188.95, 272.22,
		-0.12, 0.30, 6.09, -4.22, 0.34, -292.37, -5.10, -0.32,
		-0.40, -0.11, 6.54, 0.14, 0.01, 161.79, -220.67, 0.24,
		0.10, -4.93, -3.62, -0.08, 261.54, -19.94, -0.95, 0.20,
		-0.45, -5.85, -0.13, 0.02, 142.16, -190.79, 0.20, 0.10,
		-4.27, -3.18, -0.07, 187.95, -4.11, -0.24, 0.30, -0.09,
		-4.20, -0.09, 0.01, 0.00, 0.00, -79.08, 167.90, 0.04,
		0.00, 3.75, 1.77, 121.98, 131.04, -0.08, 0.10, 2.93,
		-2.73, -0.06, -172.95, -8.11, -0.40, -0.20, -0.18, 3.87,
		0.09, 0.01, -160.15, -55.30, -14.04, 13.90, -1.23, 3.58,
		0.40, 0.31, -115.40, 123.20, 0.60, 0.70, 2.75, 2.58,
		0.08, -0.01, -168.26, -2.00, 0.20, -0.20, -0.04, 3.76,

		/* 596-723 */
		0.08, -114.49, 123.20, 0.32, 0.40, 2.75, 2.56, 0.07,
		-0.01, 112.14, 120.70, 0.28, -0.30, 2.70, -2.51, -0.07,
		-0.01, 161.34, 4.03, 0.20, 0.20, 0.09, -3.61, -0.08,
		91.31, 126.64, -0.40, 0.40, 2.83, -2.04, -0.04, 0.01,
		105.29, 112.90, 0.44, -0.50, 2.52, -2.35, -0.07, -0.01,
		98.69, -106.20, -0.28, -0.30, -2.37, -2.21, -0.06, 0.01,
		86.74, -112.94, -0.08, -0.20, -2.53, -1.94, -0.05, -134.81,
		3.51, 0.20, -0.20, 0.08, 3.01, 0.07, 79.03, 107.31,
		-0.24, 0.20, 2.40, -1.77, -0.04, 0.01, 132.81, -10.77,
		-0.52, 0.10, -0.24, -2.97, -0.07, 0.01, -130.31, -0.90,
		0.04, 0.00, 0.00, 2.91, -78.56, 85.32, 0.00, 0.00,
		1.91, 1.76, 0.04, 0.00, 0.00, -41.53, 89.10, 0.02,
		0.00, 1.99, 0.93, 66.03, -71.00, -0.20, -0.20, -1.59,
		-1.48, -0.04, 60.50, 64.70, 0.36, -0.40, 1.45, -1.35,
		-0.04, -0.01, -52.27, -70.01, 0.00, 0.00, -1.57, 1.17,
		0.03, -52.95, 66.29, 0.32, 0.40, 1.48, 1.18, 0.04,

		/* 724-851 */
		-0.01, 51.02, 67.25, 0.00, 0.00, 1.50, -1.14, -0.03,
		-55.66, -60.92, 0.16, -0.20, -1.36, 1.24, 0.03, -54.81,
		-59.20, -0.08, 0.20, -1.32, 1.23, 0.03, 51.32, -55.60,
		0.00, 0.00, -1.24, -1.15, -0.03, 48.29, 51.80, 0.20,
		-0.20, 1.16, -1.08, -0.03, -45.59, -49.00, -0.12, 0.10,
		-1.10, 1.02, 0.03, 40.54, -52.69, -0.04, -0.10, -1.18,
		-0.91, -0.02, -40.58, -49.51, -1.00, 1.00, -1.11, 0.91,
		0.04, 0.02, -43.76, 46.50, 0.36, 0.40, 1.04, 0.98,
		0.03, -0.01, 62.65, -5.00, -0.24, 0.00, -0.11, -1.40,
		-0.03, 0.01, -38.57, 49.59, 0.08, 0.10, 1.11, 0.86,
		0.02, -33.22, -44.04, 0.08, -0.10, -0.98, 0.74, 0.02,
		37.15, -39.90, -0.12, -0.10, -0.89, -0.83, -0.02, 36.68,
		-39.50, -0.04, -0.10, -0.88, -0.82, -0.02, -53.22, -3.91,
		-0.20, 0.00, -0.09, 1.19, 0.03, 32.43, -42.19, -0.04,
		-0.10, -0.94, -0.73, -0.02, -51.00, -2.30, -0.12, -0.10,
		0.00, 1.14, -29.53, -39.11, 0.04, 0.00, -0.87, 0.66,

		/* 852-979 */
		0.02, 28.50, -38.92, -0.08, -0.10, -0.87, -0.64, -0.02,
		26.54, 36.95, -0.12, 0.10, 0.83, -0.59, -0.01, 26.54,
		34.59, 0.04, -0.10, 0.77, -0.59, -0.02, 28.35, -32.55,
		-0.16, 0.20, -0.73, -0.63, -0.01, -28.00, 30.40, 0.00,
		0.00, 0.68, 0.63, 0.01, -27.61, 29.40, 0.20, 0.20,
		0.66, 0.62, 0.02, 40.33, 0.40, -0.04, 0.10, 0.00,
		-0.90, -23.28, 31.61, -0.08, -0.10, 0.71, 0.52, 0.01,
		37.75, 0.80, 0.04, 0.10, 0.00, -0.84, 23.66, 25.80,
		0.00, 0.00, 0.58, -0.53, -0.01, 21.01, -27.91, 0.00,
		0.00, -0.62, -0.47, -0.01, -34.81, 2.89, 0.04, 0.00,
		0.00, 0.78, -23.49, -25.31, 0.00, 0.00, -0.57, 0.53,
		0.01, -23.47, 25.20, 0.16, 0.20, 0.56, 0.52, 0.02,
		19.58, 27.50, -0.12, 0.10, 0.62, -0.44, -0.01, -22.67,
		-24.40, -0.08, 0.10, -0.55, 0.51, 0.01, -19.97, 25.00,
		0.12, 0.20, 0.56, 0.45, 0.01, 21.28, -22.80, -0.08,
		-0.10, -0.51, -0.48, -0.01, -30.47, 0.91, 0.04, 0.00,

		/* 980-1107 */
		0.00, 0.68, 18.58, 24.00, 0.04, -0.10, 0.54, -0.42,
		-0.01, -18.02, 24.40, -0.04, -0.10, 0.55, 0.40, 0.01,
		17.74, 22.50, 0.08, -0.10, 0.50, -0.40, -0.01, -19.41,
		20.70, 0.08, 0.10, 0.46, 0.43, 0.01, -18.64, 20.11,
		0.00, 0.00, 0.45, 0.42, 0.01, -16.75, 21.60, 0.04,
		0.10, 0.48, 0.37, 0.01, -18.42, -20.00, 0.00, 0.00,
		-0.45, 0.41, 0.01, -26.77, 1.41, 0.08, 0.00, 0.00,
		0.60, -26.17, -0.19, 0.00, 0.00, 0.00, 0.59, -15.52,
		20.51, 0.00, 0.00, 0.46, 0.35, 0.01, -25.42, -1.91,
		-0.08, 0.00, -0.04, 0.57, 0.45, -17.42, 18.10, 0.00,
		0.00, 0.40, 0.39, 0.01, 16.39, -17.60, -0.08, -0.10,
		-0.39, -0.37, -0.01, -14.37, 18.91, 0.00, 0.00, 0.42,
		0.32, 0.01, 23.39, -2.40, -0.12, 0.00, 0.00, -0.52,
		14.32, -18.50, -0.04, -0.10, -0.41, -0.32, -0.01, 15.69,
		17.08, 0.00, 0.00, 0.38, -0.35, -0.01, -22.99, 0.50,
		0.04, 0.00, 0.00, 0.51, 0.00, 0.00, 14.47, -17.60,

		/* 1108-1235 */
		-0.01, 0.00, -0.39, -0.32, -13.33, 18.40, -0.04, -0.10,
		0.41, 0.30, 22.47, -0.60, -0.04, 0.00, 0.00, -0.50,
		-12.78, -17.41, 0.04, 0.00, -0.39, 0.29, 0.01, -14.10,
		-15.31, 0.04, 0.00, -0.34, 0.32, 0.01, 11.98, 16.21,
		-0.04, 0.00, 0.36, -0.27, -0.01, 19.65, -1.90, -0.08,
		0.00, 0.00, -0.44, 19.61, -1.50, -0.08, 0.00, 0.00,
		-0.44, 13.41, -14.30, -0.04, -0.10, -0.32, -0.30, -0.01,
		-13.29, 14.40, 0.00, 0.00, 0.32, 0.30, 0.01, 11.14,
		-14.40, -0.04, 0.00, -0.32, -0.25, -0.01, 12.24, -13.38,
		0.04, 0.00, -0.30, -0.27, -0.01, 10.07, -13.81, 0.04,
		0.00, -0.31, -0.23, -0.01, 10.46, 13.10, 0.08, -0.10,
		0.29, -0.23, -0.01, 16.55, -1.71, -0.08, 0.00, 0.00,
		-0.37, 9.75, -12.80, 0.00, 0.00, -0.29, -0.22, -0.01,
		9.11, 12.80, 0.00, 0.00, 0.29, -0.20, 0.00, 0.00,
		-6.44, -13.80, 0.00, 0.00, -0.31, 0.14, -9.19, -12.00,
		0.00, 0.00, -0.27, 0.21, -10.30, 10.90, 0.08, 0.10,

		/* 1236-1363 */
		0.24, 0.23, 0.01, 14.92, -0.80, -0.04, 0.00, 0.00,
		-0.33, 10.02, -10.80, 0.00, 0.00, -0.24, -0.22, -0.01,
		-9.75, 10.40, 0.04, 0.00, 0.23, 0.22, 0.01, 9.67,
		-10.40, -0.04, 0.00, -0.23, -0.22, -0.01, -8.28, -11.20,
		0.04, 0.00, -0.25, 0.19, 13.32, -1.41, -0.08, 0.00,
		0.00, -0.30, 8.27, 10.50, 0.04, 0.00, 0.23, -0.19,
		0.00, 0.00, 13.13, 0.00, 0.00, 0.00, 0.00, -0.29,
		-12.93, 0.70, 0.04, 0.00, 0.00, 0.29, 7.91, -10.20,
		0.00, 0.00, -0.23, -0.18, -7.84, -10.00, -0.04, 0.00,
		-0.22, 0.18, 7.44, 9.60, 0.00, 0.00, 0.21, -0.17,
		-7.64, 9.40, 0.08, 0.10, 0.21, 0.17, 0.01, -11.38,
		0.60, 0.04, 0.00, 0.00, 0.25, -7.48, 8.30, 0.00,
		0.00, 0.19, 0.17, -10.98, -0.20, 0.00, 0.00, 0.00,
		0.25, 10.98, 0.20, 0.00, 0.00, 0.00, -0.25, 7.40,
		-7.90, -0.04, 0.00, -0.18, -0.17, -6.09, 8.40, -0.04,
		0.00, 0.19, 0.14, -6.94, -7.49, 0.00, 0.00, -0.17,

		/* 1364-1491 */
		0.16, 6.92, 7.50, 0.04, 0.00, 0.17, -0.15, 6.20,
		8.09, 0.00, 0.00, 0.18, -0.14, -6.12, 7.80, 0.04,
		0.00, 0.17, 0.14, 5.85, -7.50, 0.00, 0.00, -0.17,
		-0.13, -6.48, 6.90, 0.08, 0.10, 0.15, 0.14, 0.01,
		6.32, 6.90, 0.00, 0.00, 0.15, -0.14, 5.61, -7.20,
		0.00, 0.00, -0.16, -0.13, 9.07, 0.00, 0.00, 0.00,
		0.00, -0.20, 5.25, 6.90, 0.00, 0.00, 0.15, -0.12,
		-8.47, -0.40, 0.00, 0.00, 0.00, 0.19, 6.32, -5.39,
		-1.11, 1.10, -0.12, -0.14, 0.02, 0.02, 5.73, -6.10,
		-0.04, 0.00, -0.14, -0.13, 4.70, 6.60, -0.04, 0.00,
		0.15, -0.11, -4.90, -6.40, 0.00, 0.00, -0.14, 0.11,
		-5.33, 5.60, 0.04, 0.10, 0.13, 0.12, 0.01, -4.81,
		6.00, 0.04, 0.00, 0.13, 0.11, 5.13, 5.50, 0.04,
		0.00, 0.12, -0.11, 4.50, 5.90, 0.00, 0.00, 0.13,
		-0.10, -4.22, 6.10, 0.00, 0.00, 0.14, -4.53, 5.70,
		0.00, 0.00, 0.13, 0.10, 4.18, 5.70, 0.00, 0.00,

		/* 1492-1619 */
		0.13, -4.75, -5.19, 0.00, 0.00, -0.12, 0.11, -4.06,
		5.60, 0.00, 0.00, 0.13, -3.98, 5.60, -0.04, 0.00,
		0.13, 4.02, -5.40, 0.00, 0.00, -0.12, 4.49, -4.90,
		-0.04, 0.00, -0.11, -0.10, -3.62, -5.40, -0.16, 0.20,
		-0.12, 0.00, 0.01, 4.38, 4.80, 0.00, 0.00, 0.11,
		-6.40, -0.10, 0.00, 0.00, 0.00, 0.14, -3.98, 5.00,
		0.04, 0.00, 0.11, -3.82, -5.00, 0.00, 0.00, -0.11,
		-3.71, 5.07, 0.00, 0.00, 0.11, 4.14, 4.40, 0.00,
		0.00, 0.10, -6.01, -0.50, -0.04, 0.00, 0.00, 0.13,
		-4.04, 4.39, 0.00, 0.00, 0.10, 3.45, -4.72, 0.00,
		0.00, -0.11, 3.31, 4.71, 0.00, 0.00, 0.11, 3.26,
		-4.50, 0.00, 0.00, -0.10, -3.26, -4.50, 0.00, 0.00,
		-0.10, -3.34, -4.40, 0.00, 0.00, -0.10, -3.74, -4.00,
		3.70, 4.00, 3.34, -4.30, 3.30, -4.30, -3.66, 3.90,
		0.04, 3.66, 3.90, 0.04, -3.62, -3.90, -3.61, 3.90,
		-0.20, 5.30, 0.00, 0.00, 0.12, 3.06, 4.30, 3.30,

		/* 1620-1747 */
		4.00, 0.40, 0.20, 3.10, 4.10, -3.06, 3.90, -3.30,
		-3.60, -3.30, 3.36, 0.01, 3.14, 3.40, -4.57, -0.20,
		0.00, 0.00, 0.00, 0.10, -2.70, -3.60, 2.94, -3.20,
		-2.90, 3.20, 2.47, -3.40, 2.55, -3.30, 2.80, -3.08,
		2.51, 3.30, -4.10, 0.30, -0.12, -0.10, 4.10, 0.20,
		-2.74, 3.00, 2.46, 3.23, -3.66, 1.20, -0.20, 0.20,
		3.74, -0.40, -2.51, -2.80, -3.74, 2.27, -2.90, 0.00,
		0.00, -2.50, 2.70, -2.51, 2.60, -3.50, 0.20, 3.38,
		-2.22, -2.50, 3.26, -0.40, 1.95, -2.60, 3.22, -0.40,
		-0.04, -1.79, -2.60, 1.91, 2.50, 0.74, 3.05, -0.04,
		0.08, 2.11, -2.30, -2.11, 2.20, -1.87, -2.40, 2.03,
		-2.20, -2.03, 2.20, 2.98, 0.00, 0.00, 2.98, -1.71,
		2.40, 2.94, -0.10, -0.12, 0.10, 1.67, 2.40, -1.79,
		2.30, -1.79, 2.20, -1.67, 2.20, 1.79, -2.00, 1.87,
		-1.90, 1.63, -2.10, -1.59, 2.10, 1.55, -2.10, -1.55,
		2.10, -2.59, -0.20, -1.75, -1.90, -1.75, 1.90, -1.83,

		/* 1748-1875 */
		-1.80, 1.51, 2.00, -1.51, -2.00, 1.71, 1.80, 1.31,
		2.10, -1.43, 2.00, 1.43, 2.00, -2.43, -1.51, 1.90,
		-1.47, 1.90, 2.39, 0.20, -2.39, 1.39, 1.90, 1.39,
		-1.80, 1.47, -1.60, 1.47, -1.60, 1.43, -1.50, -1.31,
		1.60, 1.27, -1.60, -1.27, 1.60, 1.27, -1.60, 2.03,
		1.35, 1.50, -1.39, -1.40, 1.95, -0.20, -1.27, 1.49,
		1.19, 1.50, 1.27, 1.40, 1.15, 1.50, 1.87, -0.10,
		-1.12, -1.50, 1.87, -1.11, -1.50, -1.11, -1.50, 0.00,
		0.00, 1.19, 1.40, 1.27, -1.30, -1.27, -1.30, -1.15,
		1.40, -1.23, 1.30, -1.23, -1.30, 1.22, -1.29, 1.07,
		-1.40, 1.75, -0.20, -1.03, -1.40, -1.07, 1.20, -1.03,
		1.15, 1.07, 1.10, 1.51, -1.03, 1.10, 1.03, -1.10,
		0.00, 0.00, -1.03, -1.10, 0.91, -1.20, -0.88, -1.20,
		-0.88, 1.20, -0.95, 1.10, -0.95, -1.10, 1.43, -1.39,
		0.95, -1.00, -0.95, 1.00, -0.80, 1.10, 0.91, -1.00,
		-1.35, 0.88, 1.00, -0.83, 1.00, -0.91, 0.90, 0.91,

		/* 1876-2003 */
		0.90, 0.88, -0.90, -0.76, -1.00, -0.76, 1.00, 0.76,
		1.00, -0.72, 1.00, 0.84, -0.90, 0.84, 0.90, 1.23,
		0.00, 0.00, -0.52, -1.10, -0.68, 1.00, 1.19, -0.20,
		1.19, 0.76, 0.90, 1.15, -0.10, 1.15, -0.10, 0.72,
		-0.90, -1.15, -1.15, 0.68, 0.90, -0.68, 0.90, -1.11,
		0.00, 0.00, 0.20, 0.79, 0.80, -1.11, -0.10, 0.00,
		0.00, -0.48, -1.00, -0.76, -0.80, -0.72, -0.80, -1.07,
		-0.10, 0.64, 0.80, -0.64, -0.80, 0.64, 0.80, 0.40,
		0.60, 0.52, -0.50, -0.60, -0.80, -0.71, 0.70, -0.99,
		0.99, 0.56, 0.80, -0.56, 0.80, 0.68, -0.70, 0.68,
		0.70, -0.95, -0.64, 0.70, 0.64, 0.70, -0.60, 0.70,
		-0.60, -0.70, -0.91, -0.10, -0.51, 0.76, -0.91, -0.56,
		0.70, 0.88, 0.88, -0.63, -0.60, 0.55, -0.60, -0.80,
		0.80, -0.80, -0.52, 0.60, 0.52, 0.60, 0.52, -0.60,
		-0.48, 0.60, 0.48, 0.60, 0.48, 0.60, -0.76, 0.44,
		-0.60, 0.52, -0.50, -0.52, 0.50, 0.40, 0.60, -0.40,

		/* 2004-2131 */
		-0.60, 0.40, -0.60, 0.72, -0.72, -0.51, -0.50, -0.48,
		0.50, 0.48, -0.50, -0.48, 0.50, -0.48, 0.50, 0.48,
		-0.50, -0.48, -0.50, -0.68, -0.68, 0.44, 0.50, -0.64,
		-0.10, -0.64, -0.10, -0.40, 0.50, 0.40, 0.50, 0.40,
		0.50, 0.00, 0.00, -0.40, -0.50, -0.36, -0.50, 0.36,
		-0.50, 0.60, -0.60, 0.40, -0.40, 0.40, 0.40, -0.40,
		0.40, -0.40, 0.40, -0.56, -0.56, 0.36, -0.40, -0.36,
		0.40, 0.36, -0.40, -0.36, -0.40, 0.36, 0.40, 0.36,
		0.40, -0.52, 0.52, 0.52, 0.32, 0.40, -0.32, 0.40,
		-0.32, 0.40, -0.32, 0.40, 0.32, -0.40, -0.32, -0.40,
		0.32, -0.40, 0.28, -0.40, -0.28, 0.40, 0.28, -0.40,
		0.28, 0.40, 0.48, -0.48, 0.48, 0.36, -0.30, -0.36,
		-0.30, 0.00, 0.00, 0.20, 0.40, -0.44, 0.44, -0.44,
		-0.44, -0.44, -0.44, 0.32, -0.30, 0.32, 0.30, 0.24,
		0.30, -0.12, -0.10, -0.28, 0.30, 0.28, 0.30, 0.28,
		0.30, 0.28, -0.30, 0.28, -0.30, 0.28, -0.30, 0.28,

		/* 2132-2259 */
		0.30, -0.28, 0.30, 0.40, 0.40, -0.24, 0.30, 0.24,
		-0.30, 0.24, -0.30, -0.24, -0.30, 0.24, 0.30, 0.24,
		-0.30, -0.24, 0.30, 0.24, -0.30, -0.24, -0.30, 0.24,
		-0.30, 0.24, 0.30, -0.24, 0.30, -0.24, 0.30, 0.20,
		-0.30, 0.20, -0.30, 0.20, -0.30, 0.20, 0.30, 0.20,
		-0.30, 0.20, -0.30, 0.20, 0.30, 0.20, 0.30, -0.20,
		-0.30, 0.20, -0.30, 0.20, -0.30, -0.36, -0.36, -0.36,
		-0.04, 0.30, 0.12, -0.10, -0.32, -0.24, 0.20, 0.24,
		0.20, 0.20, -0.20, -0.20, -0.20, -0.20, -0.20, 0.20,
		0.20, 0.20, -0.20, 0.20, 0.20, 0.20, 0.20, -0.20,
		-0.20, 0.00, 0.00, -0.20, -0.20, -0.20, 0.20, -0.20,
		0.20, 0.20, -0.20, -0.20, -0.20, 0.20, 0.20, 0.20,
		0.20, 0.20, -0.20, 0.20, -0.20, 0.28, 0.28, 0.28,
		0.28, 0.28, 0.28, -0.28, 0.28, 0.12, 0.00, 0.24,
		0.16, -0.20, 0.16, -0.20, 0.16, -0.20, 0.16, 0.20,
		-0.16, 0.20, 0.16, 0.20, -0.16, 0.20, -0.16, 0.20,

		/* 2260-2387 */
		-0.16, 0.20, 0.16, -0.20, 0.16, 0.20, 0.16, -0.20,
		-0.16, 0.20, -0.16, -0.20, -0.16, 0.20, 0.16, 0.20,
		0.16, -0.20, 0.16, -0.20, 0.16, 0.20, 0.16, 0.20,
		0.16, 0.20, -0.16, -0.20, 0.16, 0.20, -0.16, 0.20,
		0.16, 0.20, -0.16, -0.20, 0.16, -0.20, 0.16, -0.20,
		-0.16, -0.20, 0.24, -0.24, -0.24, 0.24, 0.24, 0.12,
		0.20, 0.12, 0.20, -0.12, -0.20, 0.12, -0.20, 0.12,
		-0.20, -0.12, 0.20, -0.12, 0.20, -0.12, -0.20, 0.12,
		0.20, 0.12, 0.20, 0.12, -0.20, -0.12, 0.20, 0.12,
		-0.20, -0.12, 0.20, 0.12, 0.20, 0.00, 0.00, -0.12,
		0.20, -0.12, 0.20, 0.12, -0.20, -0.12, 0.20, 0.12,
		0.20, 0.00, -0.21, -0.20, 0.00, 0.00, 0.20, -0.20,
		-0.20, -0.20, 0.20, -0.16, -0.10, 0.00, 0.17, 0.16,
		0.16, 0.16, 0.16, -0.16, 0.16, 0.16, -0.16, 0.16,
		-0.16, 0.16, 0.12, 0.10, 0.12, -0.10, -0.12, 0.10,
		-0.12, 0.10, 0.12, -0.10, -0.12, 0.12, -0.12, 0.12,

		/* 2388-2515 */
		-0.12, 0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12,
		-0.12, 0.12, 0.12, 0.12, 0.12, -0.12, -0.12, 0.12,
		0.12, 0.12, -0.12, 0.12, -0.12, -0.12, -0.12, 0.12,
		-0.12, -0.12, 0.12, 0.00, 0.11, 0.11, -122.67, 164.70,
		203.78, 273.50, 3.58, 2.74, 6.18, -4.56, 0.00, -0.04,
		0.00, -0.07, 57.44, -77.10, 95.82, 128.60, -1.77, -1.28,
		2.85, -2.14, 82.14, 89.50, 0.00, 0.00, 2.00, -1.84,
		-0.04, 47.73, -64.10, 23.79, 31.90, -1.45, -1.07, 0.69,
		-0.53, -46.38, 50.50, 0.00, 0.00, 1.13, 1.04, 0.02,
		-18.38, 0.00, 63.80, 0.00, 0.00, 0.41, 0.00, -1.43,
		59.07, 0.00, 0.00, 0.00, 0.00, -1.32, 57.28, 0.00,
		0.00, 0.00, 0.00, -1.28, -48.65, 0.00, -1.15, 0.00,
		0.00, 1.09, 0.00, 0.03, -18.30, 24.60, -17.30, -23.20,
		0.56, 0.41, -0.51, 0.39, -16.91, 26.90, 8.43, 13.30,
		0.60, 0.38, 0.31, -0.19, 1.23, -1.70, -19.13, -25.70,
		-0.03, -0.03, -0.58, 0.43, -0.72, 0.90, -17.34, -23.30,

		/* 2516-2643 */
		0.03, 0.02, -0.52, 0.39, -19.49, -21.30, 0.00, 0.00,
		-0.48, 0.44, 0.01, 20.57, -20.10, 0.64, 0.70, -0.45,
		-0.46, 0.00, -0.01, 4.89, 5.90, -16.55, 19.90, 0.14,
		-0.11, 0.44, 0.37, 18.22, 19.80, 0.00, 0.00, 0.44,
		-0.41, -0.01, 4.89, -5.30, -16.51, -18.00, -0.11, -0.11,
		-0.41, 0.37, -17.86, 0.00, 17.10, 0.00, 0.00, 0.40,
		0.00, -0.38, 0.32, 0.00, 24.42, 0.00, 0.00, -0.01,
		0.00, -0.55, -23.79, 0.00, 0.00, 0.00, 0.00, 0.53,
		14.72, -16.00, -0.32, 0.00, -0.36, -0.33, -0.01, 0.01,
		3.34, -4.50, 11.86, 15.90, -0.11, -0.07, 0.35, -0.27,
		-3.26, 4.40, 11.62, 15.60, 0.09, 0.07, 0.35, -0.26,
		-19.53, 0.00, 5.09, 0.00, 0.00, 0.44, 0.00, -0.11,
		-13.48, 14.70, 0.00, 0.00, 0.33, 0.30, 0.01, 10.86,
		-14.60, 3.18, 4.30, -0.33, -0.24, 0.09, -0.07, -11.30,
		-15.10, 0.00, 0.00, -0.34, 0.25, 0.01, 2.03, -2.70,
		10.82, 14.50, -0.07, -0.05, 0.32, -0.24, 17.46, 0.00,

		/* 2644-2771 */
		0.00, 0.00, 0.00, -0.39, 16.43, 0.00, 0.52, 0.00,
		0.00, -0.37, 0.00, -0.01, 9.35, 0.00, 13.29, 0.00,
		0.00, -0.21, 0.00, -0.30, -10.42, 11.40, 0.00, 0.00,
		0.25, 0.23, 0.01, 0.44, 0.50, -10.38, 11.30, 0.02,
		-0.01, 0.25, 0.23, -14.64, 0.00, 0.00, 0.00, 0.00,
		0.33, 0.56, 0.80, -8.67, 11.70, 0.02, -0.01, 0.26,
		0.19, 13.88, 0.00, -2.47, 0.00, 0.00, -0.31, 0.00,
		0.06, -1.99, 2.70, 7.72, 10.30, 0.06, 0.04, 0.23,
		-0.17, -0.20, 0.00, 13.05, 0.00, 0.00, 0.00, 0.00,
		-0.29, 6.92, -9.30, 3.34, 4.50, -0.21, -0.15, 0.10,
		-0.07, -6.60, 0.00, 10.70, 0.00, 0.00, 0.15, 0.00,
		-0.24, -8.04, -8.70, 0.00, 0.00, -0.19, 0.18, -10.58,
		0.00, -3.10, 0.00, 0.00, 0.24, 0.00, 0.07, -7.32,
		8.00, -0.12, -0.10, 0.18, 0.16, 1.63, 1.70, 6.96,
		-7.60, 0.03, -0.04, -0.17, -0.16, -3.62, 0.00, 9.86,
		0.00, 0.00, 0.08, 0.00, -0.22, 0.20, -0.20, -6.88,

		/* 2772-2899 */
		-7.50, 0.00, 0.00, -0.17, 0.15, -8.99, 0.00, 4.02,
		0.00, 0.00, 0.20, 0.00, -0.09, -1.07, 1.40, -5.69,
		-7.70, 0.03, 0.02, -0.17, 0.13, 6.48, -7.20, -0.48,
		-0.50, -0.16, -0.14, -0.01, 0.01, 5.57, -7.50, 1.07,
		1.40, -0.17, -0.12, 0.03, -0.02, 8.71, 0.00, 3.54,
		0.00, 0.00, -0.19, 0.00, -0.08, 0.40, 0.00, 9.27,
		0.00, 0.00, -0.01, 0.00, -0.21, -6.13, 6.70, -1.19,
		-1.30, 0.15, 0.14, -0.03, 0.03, 5.21, -5.70, -2.51,
		-2.60, -0.13, -0.12, -0.06, 0.06, 5.69, -6.20, -0.12,
		-0.10, -0.14, -0.13, -0.01, 2.03, -2.70, 4.53, 6.10,
		-0.06, -0.05, 0.14, -0.10, 5.01, 5.50, -2.51, 2.70,
		0.12, -0.11, 0.06, 0.06, -1.91, 2.60, -4.38, -5.90,
		0.06, 0.04, -0.13, 0.10, 4.65, -6.30, 0.00, 0.00,
		-0.14, -0.10, -5.29, 5.70, 0.00, 0.00, 0.13, 0.12,
		-2.23, -4.00, -4.65, 4.20, -0.09, 0.05, 0.10, 0.10,
		-4.53, 6.10, 0.00, 0.00, 0.14, 0.10, 2.47, 2.70,

		/* 2900-3027 */
		-4.46, 4.90, 0.06, -0.06, 0.11, 0.10, -5.05, 5.50,
		0.84, 0.90, 0.12, 0.11, 0.02, -0.02, 4.97, -5.40,
		-1.71, 0.00, -0.12, -0.11, 0.00, 0.04, -0.99, -1.30,
		4.22, -5.70, -0.03, 0.02, -0.13, -0.09, 0.99, 1.40,
		4.22, -5.60, 0.03, -0.02, -0.13, -0.09, -4.69, -5.20,
		0.00, 0.00, -0.12, 0.10, -3.42, 0.00, 6.09, 0.00,
		0.00, 0.08, 0.00, -0.14, -4.65, -5.10, 0.00, 0.00,
		-0.11, 0.10, 0.00, 0.00, -4.53, -5.00, 0.00, 0.00,
		-0.11, 0.10, -2.43, -2.70, -3.82, 4.20, -0.06, 0.05,
		0.10, 0.09, 0.00, 0.00, -4.53, 4.90, 0.00, 0.00,
		0.11, 0.10, -4.49, -4.90, 0.00, 0.00, -0.11, 0.10,
		2.67, -2.90, -3.62, -3.90, -0.06, -0.06, -0.09, 0.08,
		3.94, -5.30, 0.00, 0.00, -0.12, -3.38, 3.70, -2.78,
		-3.10, 0.08, 0.08, -0.07, 0.06, 3.18, -3.50, -2.82,
		-3.10, -0.08, -0.07, -0.07, 0.06, -5.77, 0.00, 1.87,
		0.00, 0.00, 0.13, 0.00, -0.04, 3.54, -4.80, -0.64,

		/* 3028-3155 */
		-0.90, -0.11, 0.00, -0.02, -3.50, -4.70, 0.68, -0.90,
		-0.11, 0.00, -0.02, 5.49, 0.00, 0.00, 0.00, 0.00,
		-0.12, 1.83, -2.50, 2.63, 3.50, -0.06, 0.00, 0.08,
		3.02, -4.10, 0.68, 0.90, -0.09, 0.00, 0.02, 0.00,
		0.00, 5.21, 0.00, 0.00, 0.00, 0.00, -0.12, -3.54,
		3.80, 2.70, 3.60, -1.35, 1.80, 0.08, 0.00, 0.04,
		-2.90, 3.90, 0.68, 0.90, 0.09, 0.00, 0.02, 0.80,
		-1.10, -2.78, -3.70, -0.02, 0.00, -0.08, 4.10, 0.00,
		-2.39, 0.00, 0.00, -0.09, 0.00, 0.05, -1.59, 2.10,
		2.27, 3.00, 0.05, 0.00, 0.07, -2.63, 3.50, -0.48,
		-0.60, -2.94, -3.20, -2.94, 3.20, 2.27, -3.00, -1.11,
		-1.50, -0.07, 0.00, -0.03, -0.56, -0.80, -2.35, 3.10,
		0.00, -0.60, -3.42, 1.90, -0.12, -0.10, 2.63, -2.90,
		2.51, 2.80, -0.64, 0.70, -0.48, -0.60, 2.19, -2.90,
		0.24, -0.30, 2.15, 2.90, 2.15, -2.90, 0.52, 0.70,
		2.07, -2.80, -3.10, 0.00, 1.79, 0.00, 0.00, 0.07,

		/* 3156-3283 */
		0.00, -0.04, 0.88, 0.00, -3.46, 2.11, 2.80, -0.36,
		0.50, 3.54, -0.20, -3.50, -1.39, 1.50, -1.91, -2.10,
		-1.47, 2.00, 1.39, 1.90, 2.07, -2.30, 0.91, 1.00,
		1.99, -2.70, 3.30, 0.00, 0.60, -0.44, -0.70, -1.95,
		2.60, 2.15, -2.40, -0.60, -0.70, 3.30, 0.84, 0.00,
		-3.10, -3.10, 0.00, -0.72, -0.32, 0.40, -1.87, -2.50,
		1.87, -2.50, 0.32, 0.40, -0.24, 0.30, -1.87, -2.50,
		-0.24, -0.30, 1.87, -2.50, -2.70, 0.00, 1.55, 2.03,
		2.20, -2.98, -1.99, -2.20, 0.12, -0.10, -0.40, 0.50,
		1.59, 2.10, 0.00, 0.00, -1.79, 2.00, -1.03, 1.40,
		-1.15, -1.60, 0.32, 0.50, 1.39, -1.90, 2.35, -1.27,
		1.70, 0.60, 0.80, -0.32, -0.40, 1.35, -1.80, 0.44,
		0.00, 2.23, -0.84, 0.90, -1.27, -1.40, -1.47, 1.60,
		-0.28, -0.30, -0.28, 0.40, -1.27, -1.70, 0.28, -0.40,
		-1.43, -1.50, 0.00, 0.00, -1.27, -1.70, 2.11, -0.32,
		-0.40, -1.23, 1.60, 1.19, -1.30, -0.72, -0.80, 0.72,

		/* 3284-3411 */
		-0.80, -1.15, -1.30, -1.35, -1.50, -1.19, -1.60, -0.12,
		0.20, 1.79, 0.00, -0.88, -0.28, 0.40, 1.11, 1.50,
		-1.83, 0.00, 0.56, -0.12, 0.10, -1.27, -1.40, 0.00,
		0.00, 1.15, 1.50, -0.12, 0.20, 1.11, 1.50, 0.36,
		-0.50, -1.07, -1.40, -1.11, 1.50, 1.67, 0.00, 0.80,
		-1.11, 0.00, 1.43, 1.23, -1.30, -0.24, -1.19, -1.30,
		-0.24, 0.20, -0.44, -0.90, -0.95, 1.10, 1.07, -1.40,
		1.15, -1.30, 1.03, -1.10, -0.56, -0.60, -0.68, 0.90,
		-0.76, -1.00, -0.24, -0.30, 0.95, -1.30, 0.56, 0.70,
		0.84, -1.10, -0.56, 0.00, -1.55, 0.91, -1.30, 0.28,
		0.30, 0.16, -0.20, 0.95, 1.30, 0.40, -0.50, -0.88,
		-1.20, 0.95, -1.10, -0.48, -0.50, 0.00, 0.00, -1.07,
		1.20, 0.44, -0.50, 0.95, 1.10, 0.00, 0.00, 0.92,
		-1.30, 0.95, 1.00, -0.52, 0.60, 1.59, 0.24, -0.40,
		0.91, 1.20, 0.84, -1.10, -0.44, -0.60, 0.84, 1.10,
		-0.44, 0.60, -0.44, 0.60, -0.84, -1.10, -0.80, 0.00,

		/* 3412-3539 */
		1.35, 0.76, 0.20, -0.91, -1.00, 0.20, -0.30, -0.91,
		-1.20, -0.95, 1.00, -0.48, -0.50, 0.88, 1.00, 0.48,
		-0.50, -0.95, -1.10, 0.20, -0.20, -0.99, 1.10, -0.84,
		1.10, -0.24, -0.30, 0.20, -0.30, 0.84, 1.10, -1.39,
		0.00, -0.28, -0.16, 0.20, 0.84, 1.10, 0.00, 0.00,
		1.39, 0.00, 0.00, -0.95, 1.00, 1.35, -0.99, 0.00,
		0.88, -0.52, 0.00, -1.19, 0.20, 0.20, 0.76, -1.00,
		0.00, 0.00, 0.76, 1.00, 0.00, 0.00, 0.76, 1.00,
		-0.76, 1.00, 0.00, 0.00, 1.23, 0.76, 0.80, -0.32,
		0.40, -0.72, 0.80, -0.40, -0.40, 0.00, 0.00, -0.80,
		-0.90, -0.68, 0.90, -0.16, -0.20, -0.16, -0.20, 0.68,
		-0.90, -0.36, 0.50, -0.56, -0.80, 0.72, -0.90, 0.44,
		-0.60, -0.48, -0.70, -0.16, 0.00, -1.11, 0.32, 0.00,
		-1.07, 0.60, -0.80, -0.28, -0.40, -0.64, 0.00, 0.91,
		1.11, 0.64, -0.90, 0.76, -0.80, 0.00, 0.00, -0.76,
		-0.80, 1.03, 0.00, -0.36, -0.64, -0.70, 0.36, -0.40,

		/* 3540-3667 */
		1.07, 0.36, -0.50, -0.52, -0.70, 0.60, 0.00, 0.88,
		0.95, 0.00, 0.48, 0.16, -0.20, 0.60, 0.80, 0.16,
		-0.20, -0.60, -0.80, 0.00, -1.00, 0.12, 0.20, 0.16,
		-0.20, 0.68, 0.70, 0.59, -0.80, -0.99, -0.56, -0.60,
		0.36, -0.40, -0.68, -0.70, -0.68, -0.70, -0.36, -0.50,
		-0.44, 0.60, 0.64, 0.70, -0.12, 0.10, -0.52, 0.60,
		0.36, 0.40, 0.00, 0.00, 0.95, -0.84, 0.00, 0.44,
		0.56, 0.60, 0.32, -0.30, 0.00, 0.00, 0.60, 0.70,
		0.00, 0.00, 0.60, 0.70, -0.12, -0.20, 0.52, -0.70,
		0.00, 0.00, 0.56, 0.70, -0.12, 0.10, -0.52, -0.70,
		0.00, 0.00, 0.88, -0.76, 0.00, -0.44, 0.00, 0.00,
		-0.52, -0.70, 0.52, -0.70, 0.36, -0.40, -0.44, -0.50,
		0.00, 0.00, 0.60, 0.60, 0.84, 0.00, 0.12, -0.24,
		0.00, 0.80, -0.56, 0.60, -0.32, -0.30, 0.48, -0.50,
		0.28, -0.30, -0.48, -0.50, 0.12, 0.20, 0.48, -0.60,
		0.48, 0.60, -0.12, 0.20, 0.24, 0.00, 0.76, -0.52,

		/* 3668-3795 */
		-0.60, -0.52, 0.60, 0.48, -0.50, -0.24, -0.30, 0.12,
		-0.10, 0.48, 0.60, 0.52, -0.20, 0.36, 0.40, -0.44,
		0.50, -0.24, -0.30, -0.48, -0.60, -0.44, -0.60, -0.12,
		0.10, 0.76, 0.76, 0.20, -0.20, 0.48, 0.50, 0.40,
		-0.50, -0.24, -0.30, 0.44, -0.60, 0.44, -0.60, 0.36,
		0.00, -0.64, 0.72, 0.00, -0.12, 0.00, -0.10, -0.40,
		-0.60, -0.20, -0.20, -0.44, 0.50, -0.44, 0.50, 0.20,
		0.20, -0.44, -0.50, 0.20, -0.20, -0.20, 0.20, -0.44,
		-0.50, 0.64, 0.00, 0.32, -0.36, 0.50, -0.20, -0.30,
		0.12, -0.10, 0.48, 0.50, -0.12, 0.30, -0.36, -0.50,
		0.00, 0.00, 0.48, 0.50, -0.48, 0.50, 0.68, 0.00,
		-0.12, 0.56, -0.40, 0.44, -0.50, -0.12, -0.10, 0.24,
		0.30, -0.40, 0.40, 0.64, 0.00, -0.24, 0.64, 0.00,
		-0.20, 0.00, 0.00, 0.44, -0.50, 0.44, 0.50, -0.12,
		0.20, -0.36, -0.50, 0.12, 0.00, 0.64, -0.40, 0.50,
		0.00, 0.10, 0.00, 0.00, -0.40, 0.50, 0.00, 0.00,

		/* 3796-3923 */
		-0.40, -0.50, 0.56, 0.00, 0.28, 0.00, 0.10, 0.36,
		0.50, 0.00, -0.10, 0.36, -0.50, 0.36, 0.50, 0.00,
		-0.10, 0.24, -0.20, -0.36, -0.40, 0.16, 0.20, 0.40,
		-0.40, 0.00, 0.00, -0.36, -0.50, -0.36, -0.50, -0.32,
		-0.50, -0.12, 0.10, 0.20, 0.20, -0.36, 0.40, -0.60,
		0.60, 0.28, 0.00, 0.52, 0.12, -0.10, 0.40, 0.40,
		0.00, -0.50, 0.20, -0.20, -0.32, 0.40, 0.16, 0.20,
		-0.16, 0.20, 0.32, 0.40, 0.56, 0.00, -0.12, 0.32,
		-0.40, -0.16, -0.20, 0.00, 0.00, 0.40, 0.40, -0.40,
		-0.40, -0.40, 0.40, -0.36, 0.40, 0.12, 0.10, 0.00,
		0.10, 0.36, 0.40, 0.00, -0.10, 0.36, 0.40, -0.36,
		0.40, 0.00, 0.10, 0.32, 0.00, 0.44, 0.12, 0.20,
		0.28, -0.40, 0.00, 0.00, 0.36, 0.40, 0.32, -0.40,
		-0.16, 0.12, 0.10, 0.32, -0.40, 0.20, 0.30, -0.24,
		0.30, 0.00, 0.10, 0.32, 0.40, 0.00, -0.10, -0.32,
		-0.40, -0.32, 0.40, 0.00, 0.10, -0.52, -0.52, 0.52,

		/* 3924-4051 */
		0.32, -0.40, 0.00, 0.00, 0.32, 0.40, 0.32, -0.40,
		0.00, 0.00, -0.32, -0.40, -0.32, 0.40, 0.32, 0.40,
		0.00, 0.00, 0.32, 0.40, 0.00, 0.00, -0.32, -0.40,
		0.00, 0.00, 0.32, 0.40, 0.16, 0.20, 0.32, -0.30,
		-0.16, 0.00, -0.48, -0.20, 0.20, -0.28, -0.30, 0.28,
		-0.40, 0.00, 0.00, 0.28, -0.40, 0.00, 0.00, 0.28,
		-0.40, 0.00, 0.00, -0.28, -0.40, 0.28, 0.40, -0.28,
		-0.40, -0.48, -0.20, 0.20, 0.24, 0.30, 0.44, 0.00,
		0.16, 0.24, 0.30, 0.16, -0.20, 0.24, 0.30, -0.12,
		0.20, 0.20, 0.30, -0.16, 0.20, 0.00, 0.00, 0.44,
		-0.32, 0.30, 0.24, 0.00, -0.36, 0.36, 0.00, 0.24,
		0.12, -0.20, 0.20, 0.30, -0.12, 0.00, -0.28, 0.30,
		-0.24, 0.30, 0.12, 0.10, -0.28, -0.30, -0.28, 0.30,
		0.00, 0.00, -0.28, -0.30, 0.00, 0.00, -0.28, -0.30,
		0.00, 0.00, 0.28, 0.30, 0.00, 0.00, -0.28, -0.30,
		-0.28, 0.30, 0.00, 0.00, -0.28, -0.30, 0.00, 0.00,

		/* 4052-4179 */
		0.28, 0.30, 0.00, 0.00, -0.28, 0.30, 0.28, -0.30,
		-0.28, 0.30, 0.40, 0.40, -0.24, 0.30, 0.00, -0.10,
		0.16, 0.00, 0.36, -0.20, 0.30, -0.12, -0.10, -0.24,
		-0.30, 0.00, 0.00, -0.24, 0.30, -0.24, 0.30, 0.00,
		0.00, -0.24, 0.30, -0.24, 0.30, 0.24, -0.30, 0.00,
		0.00, 0.24, -0.30, 0.00, 0.00, 0.24, 0.30, 0.24,
		-0.30, 0.24, 0.30, -0.24, 0.30, -0.24, 0.30, -0.20,
		0.20, -0.16, -0.20, 0.00, 0.00, -0.32, 0.20, 0.00,
		0.10, 0.20, -0.30, 0.20, -0.20, 0.12, 0.20, -0.16,
		0.20, 0.16, 0.20, 0.20, 0.30, 0.20, 0.30, 0.00,
		0.00, -0.20, 0.30, 0.00, 0.00, 0.20, 0.30, -0.20,
		-0.30, -0.20, -0.30, 0.20, -0.30, 0.00, 0.00, 0.20,
		0.30, 0.00, 0.00, 0.20, 0.30, 0.00, 0.00, 0.20,
		0.30, 0.00, 0.00, 0.20, 0.30, 0.00, 0.00, 0.20,
		-0.30, 0.00, 0.00, -0.20, -0.30, 0.00, 0.00, -0.20,
		0.30, 0.00, 0.00, -0.20, 0.30, 0.00, 0.00, 0.36,

		/* 4180-4307 */
		0.00, 0.00, 0.36, 0.12, 0.10, -0.24, 0.20, 0.12,
		-0.20, -0.16, -0.20, -0.13, 0.10, 0.22, 0.21, 0.20,
		0.00, -0.28, 0.32, 0.00, -0.12, -0.20, -0.20, 0.12,
		-0.10, 0.12, 0.10, -0.20, 0.20, 0.00, 0.00, -0.32,
		0.32, 0.00, 0.00, 0.32, 0.32, 0.00, 0.00, -0.24,
		-0.20, 0.24, 0.20, 0.20, 0.00, -0.24, 0.00, 0.00,
		-0.24, -0.20, 0.00, 0.00, 0.24, 0.20, -0.24, -0.20,
		0.00, 0.00, -0.24, 0.20, 0.16, -0.20, 0.12, 0.10,
		0.20, 0.20, 0.00, -0.10, -0.12, 0.10, -0.16, -0.20,
		-0.12, -0.10, -0.16, 0.20, 0.20, 0.20, 0.00, 0.00,
		-0.20, 0.20, -0.20, 0.20, -0.20, 0.20, -0.20, 0.20,
		0.20, -0.20, -0.20, -0.20, 0.00, 0.00, -0.20, 0.20,
		0.20, 0.00, -0.20, 0.00, 0.00, -0.20, 0.20, -0.20,
		0.20, -0.20, -0.20, -0.20, -0.20, 0.00, 0.00, 0.20,
		0.20, 0.20, 0.20, 0.12, -0.20, -0.12, -0.10, 0.28,
		-0.28, 0.16, -0.20, 0.00, -0.10, 0.00, 0.10, -0.16,

		/* 4308-4435 */
		0.20, 0.00, -0.10, -0.16, -0.20, 0.00, -0.10, 0.16,
		-0.20, 0.16, -0.20, 0.00, 0.00, 0.16, 0.20, -0.16,
		0.20, 0.00, 0.00, 0.16, 0.20, 0.16, -0.20, 0.16,
		-0.20, -0.16, 0.20, 0.16, -0.20, 0.00, 0.00, 0.16,
		0.20, 0.00, 0.00, 0.16, 0.20, 0.00, 0.00, -0.16,
		-0.20, 0.16, -0.20, -0.16, -0.20, 0.00, 0.00, -0.16,
		-0.20, 0.00, 0.00, -0.16, 0.20, 0.00, 0.00, 0.16,
		-0.20, 0.16, 0.20, 0.16, 0.20, 0.00, 0.00, -0.16,
		-0.20, 0.00, 0.00, -0.16, -0.20, 0.00, 0.00, 0.16,
		0.20, 0.16, 0.20, 0.00, 0.00, 0.16, 0.20, 0.16,
		-0.20, 0.16, 0.20, 0.00, 0.00, -0.16, 0.20, 0.00,
		0.10, 0.12, -0.20, 0.12, -0.20, 0.00, -0.10, 0.00,
		-0.10, 0.12, 0.20, 0.00, -0.10, -0.12, 0.20, -0.15,
		0.20, -0.24, 0.24, 0.00, 0.00, 0.24, 0.24, 0.12,
		-0.20, -0.12, -0.20, 0.00, 0.00, 0.12, 0.20, 0.12,
		-0.20, 0.12, 0.20, 0.12, 0.20, 0.12, 0.20, 0.12,

		/* 4436-4563 */
		-0.20, -0.12, 0.20, 0.00, 0.00, 0.12, 0.20, 0.12,
		0.00, -0.20, 0.00, 0.00, -0.12, -0.20, 0.12, -0.20,
		0.00, 0.00, 0.12, 0.20, -0.12, 0.20, -0.12, 0.20,
		0.12, -0.20, 0.00, 0.00, 0.12, 0.20, 0.20, 0.00,
		0.12, 0.00, 0.00, -0.12, 0.20, 0.00, 0.00, -0.12,
		-0.20, 0.00, 0.00, -0.12, -0.20, -0.12, -0.20, 0.00,
		0.00, 0.12, -0.20, 0.12, -0.20, 0.12, 0.20, -0.12,
		-0.20, 0.00, 0.00, 0.12, -0.20, 0.12, -0.20, 0.12,
		0.20, 0.12, 0.00, 0.20, -0.12, -0.20, 0.00, 0.00,
		0.12, 0.20, -0.16, 0.00, 0.16, -0.20, 0.20, 0.00,
		0.00, -0.20, 0.00, 0.00, -0.20, 0.20, 0.00, 0.00,
		0.20, 0.20, -0.20, 0.00, 0.00, -0.20, 0.12, 0.00,
		-0.16, 0.20, 0.00, 0.00, 0.20, 0.12, -0.10, 0.00,
		0.10, 0.16, -0.16, -0.16, -0.16, -0.16, -0.16, 0.00,
		0.00, -0.16, 0.00, 0.00, -0.16, -0.16, -0.16, 0.00,
		0.00, -0.16, 0.00, 0.00, 0.16, 0.00, 0.00, 0.16,

		/* 4564-4691 */
		0.00, 0.00, 0.16, 0.16, 0.00, 0.00, -0.16, 0.00,
		0.00, -0.16, -0.16, 0.00, 0.00, 0.16, 0.00, 0.00,
		-0.16, -0.16, 0.00, 0.00, -0.16, -0.16, 0.12, 0.10,
		0.12, -0.10, 0.12, 0.10, 0.00, 0.00, 0.12, 0.10,
		-0.12, 0.10, 0.00, 0.00, 0.12, 0.10, 0.12, -0.10,
		0.00, 0.00, -0.12, -0.10, 0.00, 0.00, 0.12, 0.10,
		0.12, 0.00, 0.00, 0.12, 0.00, 0.00, -0.12, 0.00,
		0.00, 0.12, 0.12, 0.12, 0.12, 0.12, 0.00, 0.00,
		0.12, 0.00, 0.00, 0.12, 0.12, 0.00, 0.00, 0.12,
		0.00, 0.00, 0.12, -0.12, -0.12, 0.12, 0.12, -0.12,
		-0.12, 0.00, 0.00, 0.12, -0.12, 0.12, 0.12, -0.12,
		-0.12, 0.00, 0.00, -0.12, -0.12, 0.00, 0.00, -0.12,
		0.12, 0.00, 0.00, 0.12, 0.00, 0.00, 0.12, 0.00,
		0.00, 0.12, -0.12, 0.00, 0.00, -0.12, 0.12, -0.12,
		-0.12, 0.12, 0.00, 0.00, 0.12, 0.12, 0.12, -0.12,
		0.00, 0.00, -0.12, -0.12, -0.12, 0.00, 0.00, -0.12,

		/* 4692-NA */
		-0.12, 0.00, 0.00, 0.12, 0.12, 0.00, 0.00, -0.12,
		-0.12, -0.12, -0.12, 0.12, 0.00, 0.00, 0.12, -0.12,
		0.00, 0.00, -0.12, -0.12, 0.00, 0.00, 0.12, -0.12,
		-0.12, -0.12, -0.12, 0.12, 0.12, -0.12, -0.12, 0.00,
		0.00, -0.12, 0.00, 0.00, -0.12, 0.12, 0.00, 0.00,
		0.12, 0.00, 0.00, -0.12, -0.12, 0.00, 0.00, -0.12,
		-0.12, 0.12, 0.00, 0.00, 0.12, 0.12, 0.00, 0.00,
		0.12, 0.00, 0.00, 0.12, 0.12, 0.08, 0.00, 0.04,
	}

	/* Number of amplitude coefficients */
	NA := len(a)

	/* Amplitude usage: X or Y, sin or cos, power of T. */
	jaxy := []int{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}
	jasc := []int{0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0}
	japt := []int{0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4}

	/* Miscellaneous */
	var t, w float64
	var pt [MAXPT + 1]float64
	var fa [14]float64
	var xypr, xypl, xyls, sc [2]float64
	var arg float64

	var jpt, i, j, jxy, ialast, ifreq, m, ia, jsc int

	/* ------------------------------------------------------------------ */

	/* Interval between fundamental date J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Powers of T. */
	w = 1.0
	for jpt = 0; jpt <= MAXPT; jpt++ {
		pt[jpt] = w
		w *= t
	}

	/* Initialize totals in X and Y:  polynomial, luni-solar, planetary. */
	for jxy = 0; jxy < 2; jxy++ {
		xypr[jxy] = 0.0
		xyls[jxy] = 0.0
		xypl[jxy] = 0.0
	}

	/* --------------------------------- */
	/* Fundamental arguments (IERS 2003) */
	/* --------------------------------- */

	/* Mean anomaly of the Moon. */
	fa[0] = Fal03(t)

	/* Mean anomaly of the Sun. */
	fa[1] = Falp03(t)

	/* Mean argument of the latitude of the Moon. */
	fa[2] = Faf03(t)

	/* Mean elongation of the Moon from the Sun. */
	fa[3] = Fad03(t)

	/* Mean longitude of the ascending node of the Moon. */
	fa[4] = Faom03(t)

	/* Planetary longitudes, Mercury through Neptune. */
	fa[5] = Fame03(t)
	fa[6] = Fave03(t)
	fa[7] = Fae03(t)
	fa[8] = Fama03(t)
	fa[9] = Faju03(t)
	fa[10] = Fasa03(t)
	fa[11] = Faur03(t)
	fa[12] = Fane03(t)

	/* General accumulated precession in longitude. */
	fa[13] = Fapa03(t)

	/* -------------------------------------- */
	/* Polynomial part of precession-nutation */
	/* -------------------------------------- */

	for jxy = 0; jxy < 2; jxy++ {
		for j = MAXPT; j >= 0; j-- {
			xypr[jxy] += xyp[jxy][j] * pt[j]
		}
	}

	/* ---------------------------------- */
	/* Nutation periodic terms, planetary */
	/* ---------------------------------- */

	/* Work backwards through the coefficients per frequency list. */
	ialast = NA
	for ifreq = NFPL - 1; ifreq >= 0; ifreq-- {

		/* Obtain the argument functions. */
		arg = 0.0
		for i = 0; i < 14; i++ {
			m = mfapl[ifreq][i]
			if m != 0 {
				arg += float64(m) * fa[i]
			}
		}
		sc[0] = sin(arg)
		sc[1] = cos(arg)

		/* Work backwards through the amplitudes at this frequency. */
		ia = nc[ifreq+NFLS]
		for i = ialast; i >= ia; i-- {

			/* Coefficient number (0 = 1st). */
			j = i - ia

			/* X or Y. */
			jxy = jaxy[j]

			/* Sin or cos. */
			jsc = jasc[j]

			/* Power of T. */
			jpt = japt[j]

			/* Accumulate the component. */
			xypl[jxy] += a[i-1] * sc[jsc] * pt[jpt]
		}
		ialast = ia - 1
	}

	/* ----------------------------------- */
	/* Nutation periodic terms, luni-solar */
	/* ----------------------------------- */

	/* Continue working backwards through the number of coefficients list. */
	for ifreq = NFLS - 1; ifreq >= 0; ifreq-- {

		/* Obtain the argument functions. */
		arg = 0.0
		for i = 0; i < 5; i++ {
			m = mfals[ifreq][i]
			if m != 0 {
				arg += float64(m) * fa[i]
			}
		}
		sc[0] = sin(arg)
		sc[1] = cos(arg)

		/* Work backwards through the amplitudes at this frequency. */
		ia = nc[ifreq]
		for i = ialast; i >= ia; i-- {

			/* Coefficient number (0 = 1st). */
			j = i - ia

			/* X or Y. */
			jxy = jaxy[j]

			/* Sin or cos. */
			jsc = jasc[j]

			/* Power of T. */
			jpt = japt[j]

			/* Accumulate the component. */
			xyls[jxy] += a[i-1] * sc[jsc] * pt[jpt]
		}
		ialast = ia - 1
	}

	/* ------------------------------------ */
	/* Results:  CIP unit vector components */
	/* ------------------------------------ */

	*x = DAS2R * (xypr[0] + (xyls[0]+xypl[0])/1e6)
	*y = DAS2R * (xypr[1] + (xyls[1]+xypl[1])/1e6)
}

/*
Xys00a CIP and s, IAU 2000A

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000A
precession-nutation model.

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned:
    x,y          float64   Celestial Intermediate Pole (Note 2)
    s            float64   the CIO locator s (Note 3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The Celestial Intermediate Pole coordinates are the x,y
    components of the unit vector in the Geocentric Celestial
    Reference System.

 3) The CIO locator s (in radians) positions the Celestial
    Intermediate Origin on the equator of the CIP.

 4) A faster, but slightly less accurate result (about 1 mas for
    X,Y), can be obtained by using instead the iauXys00b function.

Called:
    Pnm00a    classical NPB matrix, IAU 2000A
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S00       the CIO locator s, given X,Y, IAU 2000A

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Xys00a(date1, date2 float64, x, y, s *float64) {
	var rbpn [3][3]float64

	/* Form the bias-precession-nutation matrix, IAU 2000A. */
	Pnm00a(date1, date2, &rbpn)

	/* Extract X,Y. */
	Bpn2xy(rbpn, x, y)

	/* Obtain s. */
	*s = S00(date1, date2, *x, *y)
}

/*
Xys00b CIP and s, IAU 2000B

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000B
precession-nutation model.

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned:
    x,y          float64   Celestial Intermediate Pole (Note 2)
    s            float64   the CIO locator s (Note 3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The Celestial Intermediate Pole coordinates are the x,y
    components of the unit vector in the Geocentric Celestial
    Reference System.

 3) The CIO locator s (in radians) positions the Celestial
    Intermediate Origin on the equator of the CIP.

 4) The present function is faster, but slightly less accurate (about
    1 mas in X,Y), than the iauXys00a function.

Called:
    Pnm00b    classical NPB matrix, IAU 2000B
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S00       the CIO locator s, given X,Y, IAU 2000A

Reference:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Xys00b(date1, date2 float64, x, y, s *float64) {
	var rbpn [3][3]float64

	/* Form the bias-precession-nutation matrix, IAU 2000A. */
	Pnm00b(date1, date2, &rbpn)

	/* Extract X,Y. */
	Bpn2xy(rbpn, x, y)

	/* Obtain s. */
	*s = S00(date1, date2, *x, *y)
}

/*
Xys06a CIP and s, IAU 2006/2000A

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2006
precession and IAU 2000A nutation models.

Given:
    date1,date2  float64  TT as a 2-part Julian Date (Note 1)

Returned:
    x,y          float64  Celestial Intermediate Pole (Note 2)
    s            float64  the CIO locator s (Note 3)

Notes:

 1) The TT date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TT)=2450123.7 could be expressed in any of these ways,
    among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 method is best matched to the way
    the argument is handled internally and will deliver the
    optimum resolution.  The MJD method and the date & time methods
    are both good compromises between resolution and convenience.

 2) The Celestial Intermediate Pole coordinates are the x,y components
    of the unit vector in the Geocentric Celestial Reference System.

 3) The CIO locator s (in radians) positions the Celestial
    Intermediate Origin on the equator of the CIP.

 4) Series-based solutions for generating X and Y are also available:
    see Capitaine & Wallace (2006) and iauXy06.

Called:
    Pnm06a    classical NPB matrix, IAU 2006/2000A
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S06       the CIO locator s, given X,Y, IAU 2006

References:

    Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Xys06a(date1, date2 float64, x, y, s *float64) {
	var rbpn [3][3]float64

	/* Form the bias-precession-nutation matrix, IAU 2006/2000A. */
	Pnm06a(date1, date2, &rbpn)

	/* Extract X,Y. */
	Bpn2xy(rbpn, x, y)

	/* Obtain s. */
	*s = S06(date1, date2, *x, *y)
}
