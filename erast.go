// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Earth Rotation and Sidereal Time

/*
Ee00 Equation of the equinoxes, IAU 2000

The equation of the equinoxes, compatible with IAU 2000 resolutions,
given the nutation in longitude and the mean obliquity.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)
    epsa         float64    mean obliquity (Note 2)
    dpsi         float64    nutation in longitude (Note 3)

Returned (function value):
                 float64    equation of the equinoxes (Note 4)

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

 2) The obliquity, in radians, is mean of date.

 3) The result, which is in radians, operates in the following sense:

       Greenwich apparent ST = GMST + equation of the equinoxes

 4) The result is compatible with the IAU 2000 resolutions.  For
    further details, see IERS Conventions 2003 and Capitaine et al.
    (2002).

Called:
    Eect00    equation of the equinoxes complementary terms

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Ee00(date1, date2 float64, epsa, dpsi float64) float64 {
	var ee float64

	/* Equation of the equinoxes. */
	ee = dpsi*cos(epsa) + Eect00(date1, date2)

	return ee
}

/*
Ee00a Equation of the equinoxes, IAU 2000A

Equation of the equinoxes, compatible with IAU 2000 resolutions.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    equation of the equinoxes (Note 2)

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

 2) The result, which is in radians, operates in the following sense:

       Greenwich apparent ST = GMST + equation of the equinoxes

 3) The result is compatible with the IAU 2000 resolutions.  For
    further details, see IERS Conventions 2003 and Capitaine et al.
    (2002).

Called:
    Pr00      IAU 2000 precession adjustments
    Obl80     mean obliquity, IAU 1980
    Nut00a    nutation, IAU 2000A
    Ee00      equation of the equinoxes, IAU 2000

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003).

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004).
*/
func Ee00a(date1, date2 float64) float64 {
	var dpsipr, depspr, epsa, dpsi, deps, ee float64

	/* IAU 2000 precession-rate adjustments. */
	Pr00(date1, date2, &dpsipr, &depspr)

	/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
	epsa = Obl80(date1, date2) + depspr

	/* Nutation in longitude. */
	Nut00a(date1, date2, &dpsi, &deps)

	/* Equation of the equinoxes. */
	ee = Ee00(date1, date2, epsa, dpsi)

	return ee
}

/*
Ee00b	Equation of the equinoxes, IAU 2000B

Equation of the equinoxes, compatible with IAU 2000 resolutions but
using the truncated nutation model IAU 2000B.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                 float64    equation of the equinoxes (Note 2)

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

 2) The result, which is in radians, operates in the following sense:

       Greenwich apparent ST = GMST + equation of the equinoxes

 3) The result is compatible with the IAU 2000 resolutions except
    that accuracy has been compromised (1 mas) for the sake of speed.
    For further details, see McCarthy & Luzum (2003), IERS
    Conventions 2003 and Capitaine et al. (2003).

Called:
    Pr00      IAU 2000 precession adjustments
    Obl80     mean obliquity, IAU 1980
    Nut00b    nutation, IAU 2000B
    Ee00      equation of the equinoxes, IAU 2000

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003)

    McCarthy, D.D. & Luzum, B.J., "An abridged model of the
    precession-nutation of the celestial pole", Celestial Mechanics &
    Dynamical Astronomy, 85, 37-49 (2003)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Ee00b(date1, date2 float64) float64 {
	var dpsipr, depspr, epsa, dpsi, deps, ee float64

	/* IAU 2000 precession-rate adjustments. */
	Pr00(date1, date2, &dpsipr, &depspr)

	/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
	epsa = Obl80(date1, date2) + depspr

	/* Nutation in longitude. */
	Nut00b(date1, date2, &dpsi, &deps)

	/* Equation of the equinoxes. */
	ee = Ee00(date1, date2, epsa, dpsi)

	return ee
}

/*
Ee06a	Equation of the equinoxes, IAU 2006/2000A

Equation of the equinoxes, compatible with IAU 2000 resolutions and
IAU 2006/2000A precession-nutation.

Given:
    date1,date2  float64    TT as a 2-part Julian Date (Note 1)

Returned (function value):
    float64    equation of the equinoxes (Note 2)

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

 2) The result, which is in radians, operates in the following sense:

       Greenwich apparent ST = GMST + equation of the equinoxes

Called:
    Anpm      normalize angle into range +/- pi
    Gst06a    Greenwich apparent sidereal time, IAU 2006/2000A
    Gmst06    Greenwich mean sidereal time, IAU 2006

Reference:

    McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
    IERS Technical Note No. 32, BKG
*/
func Ee06a(date1, date2 float64) float64 {
	var gst06a, gmst06, ee float64

	/* Apparent and mean sidereal times. */
	gst06a = Gst06a(0.0, 0.0, date1, date2)
	gmst06 = Gmst06(0.0, 0.0, date1, date2)

	/* Equation of the equinoxes. */
	ee = Anpm(gst06a - gmst06)

	return ee

}

/*
Eect00	Equation of the equinoxes complementary terms, consistent with
IAU 2000 resolutions.

Given:
    date1,date2  float64   TT as a 2-part Julian Date (Note 1)

Returned (function value):
                 float64   complementary terms (Note 2)

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

 2) The "complementary terms" are part of the equation of the
    equinoxes (EE), classically the difference between apparent and
    mean Sidereal Time:

       GAST = GMST + EE

    with:

       EE = dpsi * cos(eps)

    where dpsi is the nutation in longitude and eps is the obliquity
    of date.  However, if the rotation of the Earth were constant in
    an inertial frame the classical formulation would lead to
    apparent irregularities in the UT1 timescale traceable to side-
    effects of precession-nutation.  In order to eliminate these
    effects from UT1, "complementary terms" were introduced in 1994
    (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
    1993):

       GAST = GMST + CT + EE

    By convention, the complementary terms are included as part of
    the equation of the equinoxes rather than as part of the mean
    Sidereal Time.  This slightly compromises the "geometrical"
    interpretation of mean sidereal time but is otherwise
    inconsequential.

    The present function computes CT in the above expression,
    compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
    IERS Conventions 2003).

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

    Capitaine, N. & Gontier, A.-M., Astron.Astrophys., 275,
    645-650 (1993)

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astron.Astrophys., 406,
    1135-1149 (2003)

    IAU Resolution C7, Recommendation 3 (1994)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Eect00(date1, date2 float64) float64 {
	/* Time since J2000.0, in Julian centuries */
	var t float64

	/* Miscellaneous */
	var i, j int
	var a, s0, s1 float64

	/* Fundamental arguments */
	var fa [14]float64

	/* Returned value. */
	var eect float64

	/* ----------------------------------------- */
	/* The series for the EE complementary terms */
	/* ----------------------------------------- */

	type TERM struct {
		nfa  [8]int  /* coefficients of l,l',F,D,Om,LVe,LE,pA */
		s, c float64 /* sine and cosine coefficients */
	}

	/* Terms of order t^0 */
	e0 := []TERM{

		/* 1-10 */
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, 2640.96e-6, -0.39e-6},
		{[8]int{0, 0, 0, 0, 2, 0, 0, 0}, 63.52e-6, -0.02e-6},
		{[8]int{0, 0, 2, -2, 3, 0, 0, 0}, 11.75e-6, 0.01e-6},
		{[8]int{0, 0, 2, -2, 1, 0, 0, 0}, 11.21e-6, 0.01e-6},
		{[8]int{0, 0, 2, -2, 2, 0, 0, 0}, -4.55e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 3, 0, 0, 0}, 2.02e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 1, 0, 0, 0}, 1.98e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 3, 0, 0, 0}, -1.72e-6, 0.00e-6},
		{[8]int{0, 1, 0, 0, 1, 0, 0, 0}, -1.41e-6, -0.01e-6},
		{[8]int{0, 1, 0, 0, -1, 0, 0, 0}, -1.26e-6, -0.01e-6},

		/* 11-20 */
		{[8]int{1, 0, 0, 0, -1, 0, 0, 0}, -0.63e-6, 0.00e-6},
		{[8]int{1, 0, 0, 0, 1, 0, 0, 0}, -0.63e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 3, 0, 0, 0}, 0.46e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 1, 0, 0, 0}, 0.45e-6, 0.00e-6},
		{[8]int{0, 0, 4, -4, 4, 0, 0, 0}, 0.36e-6, 0.00e-6},
		{[8]int{0, 0, 1, -1, 1, -8, 12, 0}, -0.24e-6, -0.12e-6},
		{[8]int{0, 0, 2, 0, 0, 0, 0, 0}, 0.32e-6, 0.00e-6},
		{[8]int{0, 0, 2, 0, 2, 0, 0, 0}, 0.28e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 3, 0, 0, 0}, 0.27e-6, 0.00e-6},
		{[8]int{1, 0, 2, 0, 1, 0, 0, 0}, 0.26e-6, 0.00e-6},

		/* 21-30 */
		{[8]int{0, 0, 2, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -3, 0, 0, 0}, 0.19e-6, 0.00e-6},
		{[8]int{0, 1, -2, 2, -1, 0, 0, 0}, 0.18e-6, 0.00e-6},
		{[8]int{0, 0, 0, 0, 0, 8, -13, -1}, -0.10e-6, 0.05e-6},
		{[8]int{0, 0, 0, 2, 0, 0, 0, 0}, 0.15e-6, 0.00e-6},
		{[8]int{2, 0, -2, 0, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, 1, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{0, 1, 2, -2, 2, 0, 0, 0}, -0.14e-6, 0.00e-6},
		{[8]int{1, 0, 0, -2, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
		{[8]int{0, 0, 4, -2, 4, 0, 0, 0}, 0.13e-6, 0.00e-6},

		/* 31-33 */
		{[8]int{0, 0, 2, -2, 4, 0, 0, 0}, -0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -3, 0, 0, 0}, 0.11e-6, 0.00e-6},
		{[8]int{1, 0, -2, 0, -1, 0, 0, 0}, 0.11e-6, 0.00e-6},
	}

	/* Terms of order t^1 */
	e1 := []TERM{
		{[8]int{0, 0, 0, 0, 1, 0, 0, 0}, -0.87e-6, 0.00e-6},
	}

	/* Number of terms in the series */
	NE0 := len(e0)
	NE1 := len(e1)

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

	/* Evaluate the EE complementary terms. */
	s0 = 0.0
	s1 = 0.0

	for i = NE0 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(e0[i].nfa[j]) * fa[j]
		}
		s0 += e0[i].s*sin(a) + e0[i].c*cos(a)
	}

	for i = NE1 - 1; i >= 0; i-- {
		a = 0.0
		for j = 0; j < 8; j++ {
			a += float64(e1[i].nfa[j]) * fa[j]
		}
		s1 += e1[i].s*sin(a) + e1[i].c*cos(a)
	}

	eect = (s0 + s1*t) * DAS2R

	return eect
}

/*
Eqeq94	Equation of the equinoxes, IAU 1994

Equation of the equinoxes, IAU 1994 model.

Given:
    date1,date2   float64     TDB date (Note 1)

Returned (function value):
    float64     equation of the equinoxes (Note 2)

Notes:

 1) The date date1+date2 is a Julian Date, apportioned in any
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

 2) The result, which is in radians, operates in the following sense:

       Greenwich apparent ST = GMST + equation of the equinoxes

Called:
    Anpm      normalize angle into range +/- pi
    Nut80     nutation, IAU 1980
    Obl80     mean obliquity, IAU 1980

References:

    IAU Resolution C7, Recommendation 3 (1994).

    Capitaine, N. & Gontier, A.-M., 1993, Astron.Astrophys., 275,
    645-650.
*/
func Eqeq94(date1, date2 float64) float64 {
	var t, om, dpsi, deps, eps0, ee float64

	/* Interval between fundamental epoch J2000.0 and given date (JC). */
	t = ((date1 - DJ00) + date2) / DJC

	/* Longitude of the mean ascending node of the lunar orbit on the */
	/* ecliptic, measured from the mean equinox of date. */
	om = Anpm((450160.280+(-482890.539+
		(7.455+0.008*t)*t)*t)*DAS2R +
		fmod(-5.0*t, 1.0)*D2PI)

	/* Nutation components and mean obliquity. */
	Nut80(date1, date2, &dpsi, &deps)
	eps0 = Obl80(date1, date2)

	/* Equation of the equinoxes. */
	ee = dpsi*cos(eps0) + DAS2R*(0.00264*sin(om)+0.000063*sin(om+om))

	return ee
}

/*
Era00 Earth Rotation Angle, IAU 2000

Given:
    dj1,dj2   float64    UT1 as a 2-part Julian Date (see note)

Returned (function value):
    float64    Earth rotation angle (radians), range 0-2pi

Notes:

 1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
    convenient way between the arguments dj1 and dj2.  For example,
    JD(UT1)=2450123.7 could be expressed in any of these ways,
    among others:

            dj1            dj2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  The date & time method is
    best matched to the algorithm used:  maximum precision is
    delivered when the dj1 argument is for 0hrs UT1 on the day in
    question and the dj2 argument lies in the range 0 to 1, or vice
    versa.

 2) The algorithm is adapted from Expression 22 of Capitaine et al.
    2000.  The time argument has been expressed in days directly,
    and, to retain precision, integer contributions have been
    eliminated.  The same formulation is given in IERS Conventions
    (2003), Chap. 5, Eq. 14.

Called:
    Anp       normalize angle into range 0 to 2pi

References:

    Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
    Astrophys., 355, 398-405.

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Era00(dj1, dj2 float64) float64 {
	var d1, d2, t, f, theta float64

	/* Days since fundamental epoch. */
	if dj1 < dj2 {
		d1 = dj1
		d2 = dj2
	} else {
		d1 = dj2
		d2 = dj1
	}
	t = d1 + (d2 - DJ00)

	/* Fractional part of T (days). */
	f = fmod(d1, 1.0) + fmod(d2, 1.0)

	/* Earth rotation angle at this UT1. */
	theta = Anp(D2PI * (f + 0.7790572732640 + 0.00273781191135448*t))

	return theta

}

/*
Gmst00	Greenwich Mean Sidereal Time, IAU 2000

Greenwich mean sidereal time (model consistent with IAU 2000
resolutions).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)
    tta,ttb    float64    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich mean sidereal time (radians)

Notes:

 1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    Julian Dates, apportioned in any convenient way between the
    argument pairs.  For example, JD(UT1)=2450123.7 could be
    expressed in any of these ways, among others:

           Part A         Part B

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable (in the case of UT;  the TT is not at all critical
    in this respect).  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    Rotation Angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    and TT to predict the effects of precession.  If UT1 is used for
    both purposes, errors of order 100 microarcseconds result.

 3) This GMST is compatible with the IAU 2000 resolutions and must be
    used only in conjunction with other IAU 2000 compatible
    components such as precession-nutation and equation of the
    equinoxes.

 4) The result is returned in the range 0 to 2pi.

 5) The algorithm is from Capitaine et al. (2003) and IERS
    Conventions 2003.

Called:
    Era00     Earth rotation angle, IAU 2000
    Anp       normalize angle into range 0 to 2pi

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Gmst00(uta, utb float64, tta, ttb float64) float64 {
	var t, gmst float64

	/* TT Julian centuries since J2000.0. */
	t = ((tta - DJ00) + ttb) / DJC

	/* Greenwich Mean Sidereal Time, IAU 2000. */
	gmst = Anp(Era00(uta, utb) +
		(0.014506+
			(4612.15739966+
				(1.39667721+
					(-0.00009344+
						(0.00001882)*t)*t)*t)*t)*DAS2R)

	return gmst
}

/*
Gmst06	Greenwich Mean Sidereal Time, IAU 2006

Greenwich mean sidereal time (consistent with IAU 2006 precession).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)
    tta,ttb    float64    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich mean sidereal time (radians)

Notes:

 1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    Julian Dates, apportioned in any convenient way between the
    argument pairs.  For example, JD=2450123.7 could be expressed in
    any of these ways, among others:

           Part A        Part B

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable (in the case of UT;  the TT is not at all critical
    in this respect).  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    rotation angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    and TT to predict the effects of precession.  If UT1 is used for
    both purposes, errors of order 100 microarcseconds result.

 3) This GMST is compatible with the IAU 2006 precession and must not
    be used with other precession models.

 4) The result is returned in the range 0 to 2pi.

Called:
    Era00     Earth rotation angle, IAU 2000
    Anp       normalize angle into range 0 to 2pi

Reference:

    Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
    Astron.Astrophys. 432, 355
*/
func Gmst06(uta, utb float64, tta, ttb float64) float64 {
	var t, gmst float64

	/* TT Julian centuries since J2000.0. */
	t = ((tta - DJ00) + ttb) / DJC

	/* Greenwich mean sidereal time, IAU 2006. */
	gmst = Anp(Era00(uta, utb) +
		(0.014506+
			(4612.156534+
				(1.3915817+
					(-0.00000044+
						(-0.000029956+
							(-0.0000000368)*t)*t)*t)*t)*t)*DAS2R)

	return gmst
}

/*
Gmst82	Greenwich Mean Sidereal Time, IAU 1982

Universal Time to Greenwich mean sidereal time (IAU 1982 model).

Given:
    dj1,dj2    float64    UT1 Julian Date (see note)

Returned (function value):
    float64    Greenwich mean sidereal time (radians)

Notes:

 1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
    convenient way between the arguments dj1 and dj2.  For example,
    JD(UT1)=2450123.7 could be expressed in any of these ways,
    among others:

            dj1            dj2

        2450123.7          0          (JD method)
         2451545        -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5         0.2         (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  The date & time method is
    best matched to the algorithm used:  maximum accuracy (or, at
    least, minimum noise) is delivered when the dj1 argument is for
    0hrs UT1 on the day in question and the dj2 argument lies in the
    range 0 to 1, or vice versa.

 2) The algorithm is based on the IAU 1982 expression.  This is
    always described as giving the GMST at 0 hours UT1.  In fact, it
    gives the difference between the GMST and the UT, the steady
    4-minutes-per-day drawing-ahead of ST with respect to UT.  When
    whole days are ignored, the expression happens to equal the GMST
    at 0 hours UT1 each day.

 3) In this function, the entire UT1 (the sum of the two arguments
    dj1 and dj2) is used directly as the argument for the standard
    formula, the constant term of which is adjusted by 12 hours to
    take account of the noon phasing of Julian Date.  The UT1 is then
    added, but omitting whole days to conserve accuracy.

Called:
    Anp       normalize angle into range 0 to 2pi

References:

    Transactions of the International Astronomical Union,
    XVIII B, 67 (1983).

    Aoki et al., Astron.Astrophys., 105, 359-361 (1982).
*/
func Gmst82(dj1, dj2 float64) float64 {
	/* Coefficients of IAU 1982 GMST-UT1 model */
	A := 24110.54841 - DAYSEC/2.0
	B := 8640184.812866
	C := 0.093104
	D := -6.2e-6

	/* The first constant, A, has to be adjusted by 12 hours because the */
	/* UT1 is supplied as a Julian date, which begins at noon.           */

	var d1, d2, t, f, gmst float64

	/* Julian centuries since fundamental epoch. */
	if dj1 < dj2 {
		d1 = dj1
		d2 = dj2
	} else {
		d1 = dj2
		d2 = dj1
	}
	t = (d1 + (d2 - DJ00)) / DJC

	/* Fractional part of JD(UT1), in seconds. */
	f = DAYSEC * (fmod(d1, 1.0) + fmod(d2, 1.0))

	/* GMST at this UT1. */
	gmst = Anp(DS2R * ((A + (B+(C+D*t)*t)*t) + f))

	return gmst
}

/*
Gst00a	Greenwich Apparent Sidereal Time, IAU 2000A

Greenwich apparent sidereal time (consistent with IAU 2000
resolutions).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)
    tta,ttb    float64    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich apparent sidereal time (radians)

Notes:

 1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    Julian Dates, apportioned in any convenient way between the
    argument pairs.  For example, JD(UT1)=2450123.7 could be
    expressed in any of these ways, among others:

            uta            utb

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable (in the case of UT;  the TT is not at all critical
    in this respect).  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    Rotation Angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    and TT to predict the effects of precession-nutation.  If UT1 is
    used for both purposes, errors of order 100 microarcseconds
    result.

 3) This GAST is compatible with the IAU 2000 resolutions and must be
    used only in conjunction with other IAU 2000 compatible
    components such as precession-nutation.

 4) The result is returned in the range 0 to 2pi.

 5) The algorithm is from Capitaine et al. (2003) and IERS
    Conventions 2003.

Called:
    Gmst00    Greenwich mean sidereal time, IAU 2000
    Ee00a     equation of the equinoxes, IAU 2000A
    Anp       normalize angle into range 0 to 2pi

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Gst00a(uta, utb float64, tta, ttb float64) float64 {
	var gmst00, ee00a, gst float64

	gmst00 = Gmst00(uta, utb, tta, ttb)
	ee00a = Ee00a(tta, ttb)
	gst = Anp(gmst00 + ee00a)

	return gst
}

/*
Gst00b	Greenwich Apparent Sidereal Time, IAU 2000B

Greenwich apparent sidereal time (consistent with IAU 2000
resolutions but using the truncated nutation model IAU 2000B).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich apparent sidereal time (radians)

Notes:

 1) The UT1 date uta+utb is a Julian Date, apportioned in any
    convenient way between the argument pair.  For example,
    JD(UT1)=2450123.7 could be expressed in any of these ways,
    among others:

            uta            utb

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    Rotation Angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) The result is compatible with the IAU 2000 resolutions, except
    that accuracy has been compromised for the sake of speed and
    convenience in two respects:

    . UT is used instead of TDB (or TT) to compute the precession
      component of GMST and the equation of the equinoxes.  This
      results in errors of order 0.1 mas at present.

    . The IAU 2000B abridged nutation model (McCarthy & Luzum, 2003)
      is used, introducing errors of up to 1 mas.

 3) This GAST is compatible with the IAU 2000 resolutions and must be
    used only in conjunction with other IAU 2000 compatible
    components such as precession-nutation.

 4) The result is returned in the range 0 to 2pi.

 5) The algorithm is from Capitaine et al. (2003) and IERS
    Conventions 2003.

Called:
    Gmst00    Greenwich mean sidereal time, IAU 2000
    Ee00b     equation of the equinoxes, IAU 2000B
    Anp       normalize angle into range 0 to 2pi

References:

    Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    implement the IAU 2000 definition of UT1", Astronomy &
    Astrophysics, 406, 1135-1149 (2003)

    McCarthy, D.D. & Luzum, B.J., "An abridged model of the
    precession-nutation of the celestial pole", Celestial Mechanics &
    Dynamical Astronomy, 85, 37-49 (2003)

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Gst00b(uta, utb float64) float64 {
	var gmst00, ee00b, gst float64

	gmst00 = Gmst00(uta, utb, uta, utb)
	ee00b = Ee00b(uta, utb)
	gst = Anp(gmst00 + ee00b)

	return gst
}

/*
Gst06	Greenwich Apparent Sidereal Time, IAU 2006 given NPB matrix

Given:
    uta,utb  float64        UT1 as a 2-part Julian Date (Notes 1,2)
    tta,ttb  float64        TT as a 2-part Julian Date (Notes 1,2)
    rnpb     [3][3]float64  nutation x precession x bias matrix

Returned (function value):
    float64        Greenwich apparent sidereal time (radians)

Notes:

 1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    Julian Dates, apportioned in any convenient way between the
    argument pairs.  For example, JD(UT1)=2450123.7 could be
    expressed in any of these ways, among others:

            uta            utb

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable (in the case of UT;  the TT is not at all critical
    in this respect).  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    rotation angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    and TT to predict the effects of precession-nutation.  If UT1 is
    used for both purposes, errors of order 100 microarcseconds
    result.

 3) Although the function uses the IAU 2006 series for s+XY/2, it is
    otherwise independent of the precession-nutation model and can in
    practice be used with any equinox-based NPB matrix.

 4) The result is returned in the range 0 to 2pi.

Called:
    Bpn2xy    extract CIP X,Y coordinates from NPB matrix
    S06       the CIO locator s, given X,Y, IAU 2006
    Anp       normalize angle into range 0 to 2pi
    Era00     Earth rotation angle, IAU 2000
    Eors      equation of the origins, given NPB matrix and s

Reference:

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Gst06(uta, utb float64, tta, ttb float64, rnpb [3][3]float64) float64 {
	var x, y, s, era, eors, gst float64

	/* Extract CIP coordinates. */
	Bpn2xy(rnpb, &x, &y)

	/* The CIO locator, s. */
	s = S06(tta, ttb, x, y)

	/* Greenwich apparent sidereal time. */
	era = Era00(uta, utb)
	eors = Eors(rnpb, s)
	gst = Anp(era - eors)

	return gst
}

/*
Gst06a	Greenwich Apparent Sidereal Time, IAU 2006/2000A

Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
resolutions).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)
    tta,ttb    float64    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich apparent sidereal time (radians)

Notes:

 1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    Julian Dates, apportioned in any convenient way between the
    argument pairs.  For example, JD(UT1)=2450123.7 could be
    expressed in any of these ways, among others:

            uta            utb

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in
    cases where the loss of several decimal digits of resolution
    is acceptable (in the case of UT;  the TT is not at all critical
    in this respect).  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    rotation angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    and TT to predict the effects of precession-nutation.  If UT1 is
    used for both purposes, errors of order 100 microarcseconds
    result.

 3) This GAST is compatible with the IAU 2000/2006 resolutions and
    must be used only in conjunction with IAU 2006 precession and
    IAU 2000A nutation.

 4) The result is returned in the range 0 to 2pi.

Called:
    Pnm06a    classical NPB matrix, IAU 2006/2000A
    Gst06     Greenwich apparent ST, IAU 2006, given NPB matrix

Reference:

    Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*/
func Gst06a(uta, utb float64, tta, ttb float64) float64 {
	var rnpb [3][3]float64
	var gst float64

	/* Classical nutation x precession x bias matrix, IAU 2000A. */
	Pnm06a(tta, ttb, &rnpb)

	/* Greenwich apparent sidereal time. */
	gst = Gst06(uta, utb, tta, ttb, rnpb)

	return gst
}

/*
Gst94	Greenwich Apparent Sidereal Time, IAU 1994

Greenwich apparent sidereal time (consistent with IAU 1982/94
resolutions).

Given:
    uta,utb    float64    UT1 as a 2-part Julian Date (Notes 1,2)

Returned (function value):
    float64    Greenwich apparent sidereal time (radians)

Notes:

 1) The UT1 date uta+utb is a Julian Date, apportioned in any
    convenient way between the argument pair.  For example,
    JD(UT1)=2450123.7 could be expressed in any of these ways, among
    others:

            uta            utb

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  For UT, the date & time
    method is best matched to the algorithm that is used by the Earth
    Rotation Angle function, called internally:  maximum precision is
    delivered when the uta argument is for 0hrs UT1 on the day in
    question and the utb argument lies in the range 0 to 1, or vice
    versa.

 2) The result is compatible with the IAU 1982 and 1994 resolutions,
    except that accuracy has been compromised for the sake of
    convenience in that UT is used instead of TDB (or TT) to compute
    the equation of the equinoxes.

 3) This GAST must be used only in conjunction with contemporaneous
    IAU standards such as 1976 precession, 1980 obliquity and 1982
    nutation.  It is not compatible with the IAU 2000 resolutions.

 4) The result is returned in the range 0 to 2pi.

Called:
    Gmst82    Greenwich mean sidereal time, IAU 1982
    Eqeq94    equation of the equinoxes, IAU 1994
    Anp       normalize angle into range 0 to 2pi

References:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

    IAU Resolution C7, Recommendation 3 (1994)
*/
func Gst94(uta, utb float64) float64 {
	var gmst82, eqeq94, gst float64

	gmst82 = Gmst82(uta, utb)
	eqeq94 = Eqeq94(uta, utb)
	gst = Anp(gmst82 + eqeq94)

	return gst
}
