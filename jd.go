// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Calendar

/*
Cal2jd Gregorian Calendar to Julian Date.

Given:
   iy,im,id  int     year, month, day in Gregorian calendar (Note 1)

Returned:
   djm0      float64  MJD zero-point: always 2400000.5
   djm       float64  Modified Julian Date for 0 hrs

Returned (function value):
   int   status:
		 0 = OK
		-1 = bad year   (Note 3: JD not computed)
		-2 = bad month  (JD not computed)
		-3 = bad day    (JD computed)

Notes:

 1. The algorithm used is valid from -4800 March 1, but this
    implementation rejects dates before -4799 January 1.

 2. The Julian Date is returned in two pieces, in the usual SOFA
    manner, which is designed to preserve time resolution.  The
    Julian Date is available as a single number by adding djm0 and
    djm.

 3. In early eras the conversion is from the "Proleptic Gregorian
    Calendar";  no account is taken of the date(s) of adoption of
    the Gregorian Calendar, nor is the AD/BC numbering convention
    observed.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).
*/
func Cal2jd(iy, im, id int, djm0, djm *float64) int {
	var j, ly, my int
	var iypmy int

	/* Earliest year allowed (4800BC) */
	const IYMIN = -4799

	/* Month lengths in days */
	mtab := []int{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}

	/* Preset status. */
	j = 0

	/* Validate year and month. */
	if iy < IYMIN {
		return -1
	}
	if im < 1 || im > 12 {
		return -2
	}

	/* If February in a leap year, 1, otherwise 0. */
	// ly = ((im == 2) && !(iy%4) && (iy%100 || !(iy%400)));
	if (im == 2) && (iy%4 == 0) && ((iy%100 != 0) || (iy%400 == 0)) {
		ly = 1
	} else {
		ly = 0
	}

	/* Validate day, taking into account leap years. */
	if (id < 1) || (id > (mtab[im-1] + ly)) {
		j = -3
	}

	/* Return result. */
	my = (im - 14) / 12
	iypmy = (iy + my)
	*djm0 = DJM0
	*djm = float64((1461*(iypmy+4800))/4 + (367*(im-2-12*my))/12 - (3*((iypmy+4900)/100))/4 + id - 2432076)

	/* Return status. */
	return j

}

/*
Jd2cal Julian Date to Gregorian year, month, day, and fraction of a day.

Given:
    dj1,dj2   float64   Julian Date (Notes 1, 2)

Returned (arguments):
    iy        int       year
    im        int       month
    id        int       day
    fd        float64   fraction of day

Returned (function value):
    int     status:
            0 = OK
           -1 = unacceptable date (Note 1)

Notes:

 1) The earliest valid date is -68569.5 (-4900 March 1).  The
    largest value accepted is 1e9.

 2) The Julian Date is apportioned in any convenient way between
    the arguments dj1 and dj2.  For example, JD=2450123.7 could
    be expressed in any of these ways, among others:

        dj1             dj2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

    Separating integer and fraction uses the "compensated summation"
    algorithm of Kahan-Neumaier to preserve as much precision as
    possible irrespective of the jd1+jd2 apportionment.

 3) In early eras the conversion is from the "proleptic Gregorian
    calendar";  no account is taken of the date(s) of adoption of
    the Gregorian calendar, nor is the AD/BC numbering convention
    observed.

References:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 12.92 (p604).

    Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
    Computing, 76, 279-293 (2006), Section 3.
*/
func Jd2cal(dj1, dj2 float64, iy, im, id *int, fd *float64) int {
	/* Minimum and maximum allowed JD */
	const DJMIN = -68569.5
	const DJMAX = 1e9

	var jd, i, l, n, k int64
	var dj, f1, f2, d, s, cs, x, t, f float64
	var v [2]float64

	/* Verify date is acceptable. */
	dj = dj1 + dj2
	if dj < DJMIN || dj > DJMAX {
		return -1
	}

	/* Separate day and fraction (where -0.5 <= fraction < 0.5). */
	d = dnint(dj1)
	f1 = dj1 - d
	jd = int64(d)
	d = dnint(dj2)
	f2 = dj2 - d
	jd += int64(d)

	/* Compute f1+f2+0.5 using compensated summation (Klein 2006). */
	s = 0.5
	cs = 0.0
	v[0] = f1
	v[1] = f2
	for i := 0; i < 2; i++ {
		x = v[i]
		t = s + x
		if fabs(s) >= fabs(x) {
			cs += (s - t) + x
		} else {
			cs += (x - t) + s
		}

		s = t
		if s >= 1.0 {
			jd++
			s -= 1.0
		}
	}
	f = s + cs
	cs = f - s

	/* Deal with negative f. */
	if f < 0.0 {

		/* Compensated summation: assume that |s| <= 1.0. */
		f = s + 1.0
		cs += (1.0 - f) + s
		s = f
		f = s + cs
		cs = f - s
		jd--
	}

	/* Deal with f that is 1.0 or more (when rounded to double). */
	if (f - 1.0) >= -DBL_EPSILON/4.0 {

		/* Compensated summation: assume that |s| <= 1.0. */
		t = s - 1.0
		cs += (s - t) - 1.0
		s = t
		f = s + cs
		if -DBL_EPSILON/2.0 < f {
			jd++
			f = gmax(f, 0.0)
		}
	}

	/* Express day in Gregorian calendar. */
	l = jd + 68569
	n = (4 * l) / 146097
	l -= (146097*n + 3) / 4
	i = (4000 * (l + 1)) / 1461001
	l -= (1461*i)/4 - 31
	k = (80 * l) / 2447
	*id = int(l - (2447*k)/80)
	l = k / 11
	*im = int(k + 2 - 12*l)
	*iy = int(100*(n-49) + i + l)
	*fd = f

	/* Success. */
	return 0
}

/*
Jdcalf  Julian Date to Gregorian Calendar, expressed in a form convenient
for formatting messages:  rounded to a specified precision.

Given:
    ndp       int       number of decimal places of days in fraction
    dj1,dj2   float64   dj1+dj2 = Julian Date (Note 1)

Returned:
    iymdf     [4]int   year, month, day, fraction in Gregorian calendar

Returned (function value):
    int     status:
            -1 = date out of range
             0 = OK
            +1 = NDP not 0-9 (interpreted as 0)

Notes:

 1) The Julian Date is apportioned in any convenient way between
    the arguments dj1 and dj2.  For example, JD=2450123.7 could
    be expressed in any of these ways, among others:

            dj1            dj2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

 2) In early eras the conversion is from the "Proleptic Gregorian
    Calendar";  no account is taken of the date(s) of adoption of
    the Gregorian Calendar, nor is the AD/BC numbering convention
    observed.

 3) See also the function Jd2cal.

 4) The number of decimal places ndp should be 4 or less if internal
    overflows are to be avoided on platforms which use 16-bit
    integers.

Called:
    Jd2cal    JD to Gregorian calendar

Reference:

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 12.92 (p604).

*/
func Jdcalf(ndp int, dj1, dj2 float64, iymdf *[4]int) int {
	var j, js int
	var denom, d1, d2, f1, f2, d, djd, f, rf float64

	/* Denominator of fraction (e.g. 100 for 2 decimal places). */
	if (ndp >= 0) && (ndp <= 9) {
		j = 0
		denom = pow10(ndp)
	} else {
		j = 1
		denom = 1.0
	}

	/* Copy the date, big then small. */
	if fabs(dj1) >= fabs(dj2) {
		d1 = dj1
		d2 = dj2
	} else {
		d1 = dj2
		d2 = dj1
	}

	/* Realign to midnight (without rounding error). */
	d1 -= 0.5

	/* Separate day and fraction (as precisely as possible). */
	d = dnint(d1)
	f1 = d1 - d
	djd = d
	d = dnint(d2)
	f2 = d2 - d
	djd += d
	d = dnint(f1 + f2)
	f = (f1 - d) + f2
	if f < 0.0 {
		f += 1.0
		d -= 1.0
	}
	djd += d

	/* Round the total fraction to the specified number of places. */
	rf = dnint(f*denom) / denom

	/* Re-align to noon. */
	djd += 0.5

	/* Convert to Gregorian calendar. */
	js = Jd2cal(djd, rf, &iymdf[0], &iymdf[1], &iymdf[2], &f)
	if js == 0 {
		iymdf[3] = int(dnint(f * denom))
	} else {
		j = js
	}

	/* Return the status. */
	return j
}

/*
Epb  Julian Date to Besselian Epoch.

Given:
	dj1,dj2    float64     Julian Date (see note)

Returned (function value):
	float64    Besselian Epoch.

Note:

	The Julian Date is supplied in two pieces, in the usual SOFA
	manner, which is designed to preserve time resolution.  The
	Julian Date is available as a single number by adding dj1 and
	dj2.  The maximum resolution is achieved if dj1 is 2451545.0
	(J2000.0).

Reference:

	Lieske, J.H., 1979. Astron.Astrophys., 73, 282.
*/
func Epb(dj1, dj2 float64) float64 {
	/* J2000.0-B1900.0 (2415019.81352) in days */
	const D1900 = 36524.68648

	return 1900.0 + ((dj1-DJ00)+(dj2+D1900))/DTY
}

/*
Epb2jd Besselian Epoch to Julian Date.

Given:
	epb      float64    Besselian Epoch (e.g. 1957.3)

Returned:
	djm0     float64    MJD zero-point: always 2400000.5
	djm      float64    Modified Julian Date

Note:

	The Julian Date is returned in two pieces, in the usual SOFA
	manner, which is designed to preserve time resolution.  The
	Julian Date is available as a single number by adding djm0 and
	djm.

Reference:
	Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
*/
func Epb2jd(epb float64, djm0, djm *float64) {
	*djm0 = DJM0
	*djm = 15019.81352 + (epb-1900.0)*DTY
}

/*
Epj Julian Date to Julian Epoch.

Given:
	dj1,dj2    float64     Julian Date (see note)

Returned (function value):
	float64    Julian Epoch

Note:

	The Julian Date is supplied in two pieces, in the usual SOFA
	manner, which is designed to preserve time resolution.  The
	Julian Date is available as a single number by adding dj1 and
	dj2.  The maximum resolution is achieved if dj1 is 2451545.0
	(J2000.0).

Reference:
	Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
*/
func Epj(dj1, dj2 float64) float64 {
	epj := 2000.0 + ((dj1-DJ00)+dj2)/DJY

	return epj
}

/*
Epj2jd Julian Epoch to Julian Date.

Given:
	epj      float64    Julian Epoch (e.g. 1996.8)

Returned:
	djm0     float64    MJD zero-point: always 2400000.5
	djm      float64    Modified Julian Date

Note:

	The Julian Date is returned in two pieces, in the usual SOFA
	manner, which is designed to preserve time resolution.  The
	Julian Date is available as a single number by adding djm0 and
	djm.

Reference:
	Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
*/
func Epj2jd(epj float64, djm0, djm *float64) {
	*djm0 = DJM0
	*djm = DJM00 + (epj-2000.0)*365.25
}
