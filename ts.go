// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

import (
	"math"
	"strings"
)

// SOFA Time Scale and Calendar Tools

// Time Scales (20)

/*
D2dtf Format 2-part JD for output

Format for output a 2-part Julian Date (or in the case of UTC a
quasi-JD form that includes special provision for leap seconds).

Given:

	scale     string  time scale ID (Note 1)
	ndp       int     resolution (Note 2)
	d1,d2     float64  time as a 2-part Julian Date (Notes 3,4)

Returned:

	iy,im,id  int     year, month, day in Gregorian calendar (Note 5)
	ihmsf     [4]int  hours, minutes, seconds, fraction (Note 1)

Returned (function value):

	int     status: +1 = dubious year (Note 5)
	                 0 = OK
	                -1 = unacceptable date (Note 6)

Notes:

 1. scale identifies the time scale.  Only the value "UTC" (in upper
    case) is significant, and enables handling of leap seconds (see
    Note 4).

 2. ndp is the number of decimal places in the seconds field, and can
    have negative as well as positive values, such as:

    ndp         resolution
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001

    The limits are platform dependent, but a safe range is -5 to +9.

 3. d1+d2 is Julian Date, apportioned in any convenient way between
    the two arguments, for example where d1 is the Julian Day Number
    and d2 is the fraction of a day.  In the case of UTC, where the
    use of JD is problematical, special conventions apply:  see the
    next note.

 4. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The SOFA internal convention is that
    the quasi-JD day represents UTC days whether the length is 86399,
    86400 or 86401 SI seconds.  In the 1960-1972 era there were
    smaller jumps (in either direction) each time the linear UTC(TAI)
    expression was changed, and these "mini-leaps" are also included
    in the SOFA convention.

 5. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

 6. For calendar conventions and limitations, see Cal2jd.

Called:

	Jd2cal    JD to Gregorian calendar
	D2tf      decompose days to hms
	Dat       delta(AT) = TAI-UTC
*/
func D2dtf(scale string, ndp int, d1, d2 float64, iy, im, id *int, ihmsf *[4]int) int {
	var leap bool
	var s byte
	var iy1, im1, id1, js, iy2, im2, id2, i int
	var ihmsf1 [4]int
	var a1, b1, fd, dat0, dat12, w, dat24, dleap float64

	/* The two-part JD. */
	a1 = d1
	b1 = d2

	/* Provisional calendar date. */
	js = Jd2cal(a1, b1, &iy1, &im1, &id1, &fd)
	if js != 0 {
		return -1
	}

	/* Is this a leap second day? */
	leap = false
	//if ( ! strcmp(scale,"UTC") ) {
	if strings.Compare(scale, "UTC") == 0 {
		/* TAI-UTC at 0h today. */
		js = Dat(iy1, im1, id1, 0.0, &dat0)
		if js < 0 {
			return -1
		}

		/* TAI-UTC at 12h today (to detect drift). */
		js = Dat(iy1, im1, id1, 0.5, &dat12)
		if js < 0 {
			return -1
		}

		/* TAI-UTC at 0h tomorrow (to detect jumps). */
		js = Jd2cal(a1+1.5, b1-fd, &iy2, &im2, &id2, &w)
		if js != 0 {
			return -1
		}
		js = Dat(iy2, im2, id2, 0.0, &dat24)
		if js < 0 {
			return -1
		}

		/* Any sudden change in TAI-UTC (seconds). */
		dleap = dat24 - (2.0*dat12 - dat0)

		/* If leap second day, scale the fraction of a day into SI. */
		// leap = (math.Abs(dleap) > 0.5)

		// if math.Abs(dleap) > 0.5 {
		// 	leap = 1
		// }
		leap = (math.Abs(dleap) > 0.5)

		if leap {
			fd += fd * dleap / DAYSEC
		}
	}

	/* Provisional time of day. */
	D2tf(ndp, fd, &s, &ihmsf1)

	/* Has the (rounded) time gone past 24h? */
	if ihmsf1[0] > 23 {

		/* Yes.  We probably need tomorrow's calendar date. */
		js = Jd2cal(a1+1.5, b1-fd, &iy2, &im2, &id2, &w)
		if js != 0 {
			return -1
		}

		/* Is today a leap second day? */
		if !leap {

			/* No.  Use 0h tomorrow. */
			iy1 = iy2
			im1 = im2
			id1 = id2
			ihmsf1[0] = 0
			ihmsf1[1] = 0
			ihmsf1[2] = 0

		} else {

			/* Yes.  Are we past the leap second itself? */
			if ihmsf1[2] > 0 {

				/* Yes.  Use tomorrow but allow for the leap second. */
				iy1 = iy2
				im1 = im2
				id1 = id2
				ihmsf1[0] = 0
				ihmsf1[1] = 0
				ihmsf1[2] = 0

			} else {
				/* No.  Use 23 59 60... today. */
				ihmsf1[0] = 23
				ihmsf1[1] = 59
				ihmsf1[2] = 60
			}

			/* If rounding to 10s or coarser always go up to new day. */
			if ndp < 0 && ihmsf1[2] == 60 {
				iy1 = iy2
				im1 = im2
				id1 = id2
				ihmsf1[0] = 0
				ihmsf1[1] = 0
				ihmsf1[2] = 0
			}
		}
	}

	/* Results. */
	*iy = iy1
	*im = im1
	*id = id1
	for i = 0; i < 4; i++ {
		ihmsf[i] = ihmsf1[i]
	}

	/* Status. */
	return js
}

/*
Dtf2d	Encode time and date fields into 2-part JD

Encode date and time fields into 2-part Julian Date (or in the case
of UTC a quasi-JD form that includes special provision for leap
seconds).

Given:

	scale     string   time scale ID (Note 1)
	iy,im,id  int      year, month, day in Gregorian calendar (Note 2)
	ihr,imn   int      hour, minute
	sec       float64  seconds

Returned:

	d1,d2     float64  2-part Julian Date (Notes 3,4)

Returned (function value):

	int     status: +3 = both of next two
	                +2 = time is after end of day (Note 5)
	                +1 = dubious year (Note 6)
	                 0 = OK
	                -1 = bad year
	                -2 = bad month
	                -3 = bad day
	                -4 = bad hour
	                -5 = bad minute
	                -6 = bad second (<0)

Notes:

 1. scale identifies the time scale.  Only the value "UTC" (in upper
    case) is significant, and enables handling of leap seconds (see
    Note 4).

 2. For calendar conventions and limitations, see Cal2jd.

 3. The sum of the results, d1+d2, is Julian Date, where normally d1
    is the Julian Day Number and d2 is the fraction of a day.  In the
    case of UTC, where the use of JD is problematical, special
    conventions apply:  see the next note.

 4. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The SOFA internal convention is that
    the quasi-JD day represents UTC days whether the length is 86399,
    86400 or 86401 SI seconds.  In the 1960-1972 era there were
    smaller jumps (in either direction) each time the linear UTC(TAI)
    expression was changed, and these "mini-leaps" are also included
    in the SOFA convention.

 5. The warning status "time is after end of day" usually means that
    the sec argument is greater than 60.0.  However, in a day ending
    in a leap second the limit changes to 61.0 (or 59.0 in the case
    of a negative leap second).

 6. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

 7. Only in the case of continuous and regular time scales (TAI, TT,
    TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
    speaking.  In the other cases (UT1 and UTC) the result must be
    used with circumspection;  in particular the difference between
    two such results cannot be interpreted as a precise time
    interval.

Called:

	Cal2jd    Gregorian calendar to JD
	Dat       delta(AT) = TAI-UTC
	Jd2cal    JD to Gregorian calendar
*/
func Dtf2d(scale string, iy, im, id, ihr, imn int, sec float64, d1, d2 *float64) int {
	var js, iy2, im2, id2 int
	var dj, w, day, seclim, dat0, dat12, dat24, dleap, time float64

	/* Today's Julian Day Number. */
	js = Cal2jd(iy, im, id, &dj, &w)
	if js != 0 {
		return js
	}
	dj += w

	/* Day length and final minute length in seconds (provisional). */
	day = DAYSEC
	seclim = 60.0

	/* Deal with the UTC leap second case. */
	if strings.Compare(scale, "UTC") == 0 {

		/* TAI-UTC at 0h today. */
		js = Dat(iy, im, id, 0.0, &dat0)
		if js < 0 {
			return js
		}

		/* TAI-UTC at 12h today (to detect drift). */
		js = Dat(iy, im, id, 0.5, &dat12)
		if js < 0 {
			return js
		}

		/* TAI-UTC at 0h tomorrow (to detect jumps). */
		js = Jd2cal(dj, 1.5, &iy2, &im2, &id2, &w)
		if js != 0 {
			return js
		}
		js = Dat(iy2, im2, id2, 0.0, &dat24)
		if js < 0 {
			return js
		}

		/* Any sudden change in TAI-UTC between today and tomorrow. */
		dleap = dat24 - (2.0*dat12 - dat0)

		/* If leap second day, correct the day and final minute lengths. */
		day += dleap
		if ihr == 23 && imn == 59 {
			seclim += dleap
		}

		/* End of UTC-specific actions. */
	}

	/* Validate the time. */
	if ihr >= 0 && ihr <= 23 {
		if imn >= 0 && imn <= 59 {
			if sec >= 0 {
				if sec >= seclim {
					js += 2
				}
			} else {
				js = -6
			}
		} else {
			js = -5
		}
	} else {
		js = -4
	}
	if js < 0 {
		return js
	}

	/* The time in days. */
	time = (60.0*(float64(60*ihr+imn)) + sec) / float64(day)

	/* Return the date and time. */
	*d1 = dj
	*d2 = time

	/* Status. */
	return js
}

/*
Dat		Delta(AT) (=TAI-UTC) for a given UTC date

For a given UTC date, calculate Delta(AT) = TAI-UTC.

	:------------------------------------------:
	:                                          :
	:                 IMPORTANT                :
	:                                          :
	:  A new version of this function must be  :
	:  produced whenever a new leap second is  :
	:  announced.  There are four items to     :
	:  change on each such occasion:           :
	:                                          :
	:  1) A new line must be added to the set  :
	:     of statements that initialize the    :
	:     array "changes".                     :
	:                                          :
	:  2) The constant IYV must be set to the  :
	:     current year.                        :
	:                                          :
	:  3) The "Latest leap second" comment     :
	:     below must be set to the new leap    :
	:     second date.                         :
	:                                          :
	:  4) The "This revision" comment, later,  :
	:     must be set to the current date.     :
	:                                          :
	:  Change (2) must also be carried out     :
	:  whenever the function is re-issued,     :
	:  even if no leap seconds have been       :
	:  added.                                  :
	:                                          :
	:  Latest leap second:  2016 December 31   :
	:                                          :
	:__________________________________________:

Given:

	iy     int       UTC:  year (Notes 1 and 2)
	im     int             month (Note 2)
	id     int             day (Notes 2 and 3)
	fd     float64         fraction of day (Note 4)

Returned:

	deltat float64   TAI minus UTC, seconds

Returned (function value):

	int      status (Note 5):
	          1 = dubious year (Note 1)
	          0 = OK
	         -1 = bad year
	         -2 = bad month
	         -3 = bad day (Note 3)
	         -4 = bad fraction (Note 4)
	         -5 = internal error (Note 5)

Notes:

 1. UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
    to call the function with an earlier date.  If this is attempted,
    zero is returned together with a warning status.

    Because leap seconds cannot, in principle, be predicted in
    advance, a reliable check for dates beyond the valid range is
    impossible.  To guard against gross errors, a year five or more
    after the release year of the present function (see the constant
    IYV) is considered dubious.  In this case a warning status is
    returned but the result is computed in the normal way.

    For both too-early and too-late years, the warning status is +1.
    This is distinct from the error status -1, which signifies a year
    so early that JD could not be computed.

 2. If the specified date is for a day which ends with a leap second,
    the TAI-UTC value returned is for the period leading up to the
    leap second.  If the date is for a day which begins as a leap
    second ends, the TAI-UTC returned is for the period following the
    leap second.

 3. The day number must be in the normal calendar range, for example
    1 through 30 for April.  The "almanac" convention of allowing
    such dates as January 0 and December 32 is not supported in this
    function, in order to avoid confusion near leap seconds.

 4. The fraction of day is used only for dates before the
    introduction of leap seconds, the first of which occurred at the
    end of 1971.  It is tested for validity (0 to 1 is the valid
    range) even if not used;  if invalid, zero is used and status -4
    is returned.  For many applications, setting fd to zero is
    acceptable;  the resulting error is always less than 3 ms (and
    occurs only pre-1972).

 5. The status value returned in the case where there are multiple
    errors refers to the first error detected.  For example, if the
    month and day are 13 and 32 respectively, status -2 (bad month)
    will be returned.  The "internal error" status refers to a
    case that is impossible but causes some compilers to issue a
    warning.

 6. In cases where a valid result is not available, zero is returned.

References:

 1. For dates from 1961 January 1 onwards, the expressions from the
    file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.

 2. The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
    the 1992 Explanatory Supplement.

Called:

	Cal2jd    Gregorian calendar to JD
*/
func Dat(iy, im, id int, fd float64, deltat *float64) int {
	/* Release year for this version of Dat */
	const IYV = 2023

	/* Reference dates (MJD) and drift rates (s/day), pre leap seconds */
	drift := [...][2]float64{
		{37300.0, 0.0012960},
		{37300.0, 0.0012960},
		{37300.0, 0.0012960},
		{37665.0, 0.0011232},
		{37665.0, 0.0011232},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{38761.0, 0.0012960},
		{39126.0, 0.0025920},
		{39126.0, 0.0025920},
	}

	/* Number of Delta(AT) expressions before leap seconds were introduced */
	// enum { NERA1 = (int) (sizeof drift / sizeof (double) / 2) };
	NERA1 := len(drift)

	/* Dates and Delta(AT)s */
	changes := [...]struct {
		iyear int
		month int
		delat float64
	}{
		{1960, 1, 1.4178180},
		{1961, 1, 1.4228180},
		{1961, 8, 1.3728180},
		{1962, 1, 1.8458580},
		{1963, 11, 1.9458580},
		{1964, 1, 3.2401300},
		{1964, 4, 3.3401300},
		{1964, 9, 3.4401300},
		{1965, 1, 3.5401300},
		{1965, 3, 3.6401300},
		{1965, 7, 3.7401300},
		{1965, 9, 3.8401300},
		{1966, 1, 4.3131700},
		{1968, 2, 4.2131700},
		{1972, 1, 10.0},
		{1972, 7, 11.0},
		{1973, 1, 12.0},
		{1974, 1, 13.0},
		{1975, 1, 14.0},
		{1976, 1, 15.0},
		{1977, 1, 16.0},
		{1978, 1, 17.0},
		{1979, 1, 18.0},
		{1980, 1, 19.0},
		{1981, 7, 20.0},
		{1982, 7, 21.0},
		{1983, 7, 22.0},
		{1985, 7, 23.0},
		{1988, 1, 24.0},
		{1990, 1, 25.0},
		{1991, 1, 26.0},
		{1992, 7, 27.0},
		{1993, 7, 28.0},
		{1994, 7, 29.0},
		{1996, 1, 30.0},
		{1997, 7, 31.0},
		{1999, 1, 32.0},
		{2006, 1, 33.0},
		{2009, 1, 34.0},
		{2012, 7, 35.0},
		{2015, 7, 36.0},
		{2017, 1, 37.0},
	}

	/* Number of Delta(AT) changes */
	// enum { NDAT = (int) (sizeof changes / sizeof changes[0]) };
	NDAT := len(changes)

	/* Miscellaneous local variables */
	var j, i, m int
	var da, djm0, djm float64

	/* Initialize the result to zero. */
	// *deltat = da = 0.0;
	*deltat = 0.0
	// da = 0.0

	/* If invalid fraction of a day, set error status and give up. */
	if fd < 0.0 || fd > 1.0 {
		return -4
	}

	/* Convert the date into an MJD. */
	j = Cal2jd(iy, im, id, &djm0, &djm)

	/* If invalid year, month, or day, give up. */
	if j < 0 {
		return j
	}

	/* If pre-UTC year, set warning status and give up. */
	if iy < changes[0].iyear {
		return 1
	}

	/* If suspiciously late year, set warning status but proceed. */
	if iy > IYV+5 {
		j = 1
	}

	/* Combine year and month to form a date-ordered integer... */
	m = 12*iy + im

	/* ...and use it to find the preceding table entry. */
	for i = NDAT - 1; i >= 0; i-- {
		if m >= (12*changes[i].iyear + changes[i].month) {
			break
		}
	}

	/* Prevent underflow warnings. */
	if i < 0 {
		return -5
	}

	/* Get the Delta(AT). */
	da = changes[i].delat

	/* If pre-1972, adjust for drift. */
	if i < NERA1 {
		da += (djm + fd - drift[i][0]) * drift[i][1]
	}

	/* Return the Delta(AT) value. */
	*deltat = da

	/* Return the status. */
	return j
}

/*
Dtdb		TDB-TT

An approximation to TDB-TT, the difference between barycentric
dynamical time and terrestrial time, for an observer on the Earth.

The different time scales - proper, coordinate and realized - are
related to each other:

	        TAI             <-  physically realized
	         :
	      offset            <-  observed (nominally +32.184s)
	         :
	        TT              <-  terrestrial time
	         :
	rate adjustment (L_G)   <-  definition of TT
	         :
	        TCG             <-  time scale for GCRS
	         :
	  "periodic" terms      <-  Dtdb  is an implementation
	         :
	rate adjustment (L_C)   <-  function of solar-system ephemeris
	         :
	        TCB             <-  time scale for BCRS
	         :
	rate adjustment (-L_B)  <-  definition of TDB
	         :
	        TDB             <-  TCB scaled to track TT
	         :
	  "periodic" terms      <-  -Dtdb is an approximation
	         :
	        TT              <-  terrestrial time

Adopted values for the various constants can be found in the IERS
Conventions (McCarthy & Petit 2003).

Given:

	date1,date2   float64  date, TDB (Notes 1-3)
	ut            float64  universal time (UT1, fraction of one day)
	elong         float64  longitude (east positive, radians)
	u             float64  distance from Earth spin axis (km)
	v             float64  distance north of equatorial plane (km)

Returned (function value):

	float64  TDB-TT (seconds)

Notes:

 1. The date date1+date2 is a Julian Date, apportioned in any
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

    Although the date is, formally, barycentric dynamical time (TDB),
    the terrestrial dynamical time (TT) can be used with no practical
    effect on the accuracy of the prediction.

 2. TT can be regarded as a coordinate time that is realized as an
    offset of 32.184s from International Atomic Time, TAI.  TT is a
    specific linear transformation of geocentric coordinate time TCG,
    which is the time scale for the Geocentric Celestial Reference
    System, GCRS.

 3. TDB is a coordinate time, and is a specific linear transformation
    of barycentric coordinate time TCB, which is the time scale for
    the Barycentric Celestial Reference System, BCRS.

 4. The difference TCG-TCB depends on the masses and positions of the
    bodies of the solar system and the velocity of the Earth.  It is
    dominated by a rate difference, the residual being of a periodic
    character.  The latter, which is modeled by the present function,
    comprises a main (annual) sinusoidal term of amplitude
    approximately 0.00166 seconds, plus planetary terms up to about
    20 microseconds, and lunar and diurnal terms up to 2 microseconds.
    These effects come from the changing transverse Doppler effect
    and gravitational red-shift as the observer (on the Earth's
    surface) experiences variations in speed (with respect to the
    BCRS) and gravitational potential.

 5. TDB can be regarded as the same as TCB but with a rate adjustment
    to keep it close to TT, which is convenient for many applications.
    The history of successive attempts to define TDB is set out in
    Resolution 3 adopted by the IAU General Assembly in 2006, which
    defines a fixed TDB(TCB) transformation that is consistent with
    contemporary solar-system ephemerides.  Future ephemerides will
    imply slightly changed transformations between TCG and TCB, which
    could introduce a linear drift between TDB and TT;  however, any
    such drift is unlikely to exceed 1 nanosecond per century.

 6. The geocentric TDB-TT model used in the present function is that of
    Fairhead & Bretagnon (1990), in its full form.  It was originally
    supplied by Fairhead (private communications with P.T.Wallace,

 1990. as a Fortran subroutine.  The present C function contains an
    adaptation of the Fairhead code.  The numerical results are
    essentially unaffected by the changes, the differences with
    respect to the Fairhead & Bretagnon original being at the 1e-20 s
    level.

    The topocentric part of the model is from Moyer (1981) and
    Murray (1983), with fundamental arguments adapted from
    Simon et al. 1994.  It is an approximation to the expression
    ( v / c ) . ( r / c ), where v is the barycentric velocity of
    the Earth, r is the geocentric position of the observer and
    c is the speed of light.

    By supplying zeroes for u and v, the topocentric part of the
    model can be nullified, and the function will return the Fairhead
    & Bretagnon result alone.

 7. During the interval 1950-2050, the absolute accuracy is better
    than +/- 3 nanoseconds relative to time ephemerides obtained by
    direct numerical integrations based on the JPL DE405 solar system
    ephemeris.

 8. It must be stressed that the present function is merely a model,
    and that numerical integration of solar-system ephemerides is the
    definitive method for predicting the relationship between TCG and
    TCB and hence between TT and TDB.

References:

	Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
	(1990).

	IAU 2006 Resolution 3.

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Moyer, T.D., Cel.Mech., 23, 33 (1981).

	Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).

	Seidelmann, P.K. et al., Explanatory Supplement to the
	Astronomical Almanac, Chapter 2, University Science Books (1992).

	Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
	Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
*/
func Dtdb(date1, date2 float64, ut, elong, u, v float64) float64 {
	var t, tsol, w, elsun, emsun, d, elj, els, wt, w0, w1, w2, w3, w4,
		wf, wj float64
	var j int

	/*
	 =====================
	 Fairhead et al. model
	 =====================

	 787 sets of three coefficients.

	 Each set is
	    amplitude (microseconds)
	      frequency (radians per Julian millennium since J2000.0)
	      phase (radians)

	 Sets   1-474 are the T0 terms
	  "   475-679  "   "  T1
	  "   680-764  "   "  T2
	  "   765-784  "   "  T3
	  "   785-787  "   "  T4
	*/

	fairhd := [787][3]float64{
		/* 1, 10 */
		{1656.674564e-6, 6283.075849991, 6.240054195},
		{22.417471e-6, 5753.384884897, 4.296977442},
		{13.839792e-6, 12566.151699983, 6.196904410},
		{4.770086e-6, 529.690965095, 0.444401603},
		{4.676740e-6, 6069.776754553, 4.021195093},
		{2.256707e-6, 213.299095438, 5.543113262},
		{1.694205e-6, -3.523118349, 5.025132748},
		{1.554905e-6, 77713.771467920, 5.198467090},
		{1.276839e-6, 7860.419392439, 5.988822341},
		{1.193379e-6, 5223.693919802, 3.649823730},

		/* 11, 20 */
		{1.115322e-6, 3930.209696220, 1.422745069},
		{0.794185e-6, 11506.769769794, 2.322313077},
		{0.447061e-6, 26.298319800, 3.615796498},
		{0.435206e-6, -398.149003408, 4.349338347},
		{0.600309e-6, 1577.343542448, 2.678271909},
		{0.496817e-6, 6208.294251424, 5.696701824},
		{0.486306e-6, 5884.926846583, 0.520007179},
		{0.432392e-6, 74.781598567, 2.435898309},
		{0.468597e-6, 6244.942814354, 5.866398759},
		{0.375510e-6, 5507.553238667, 4.103476804},
		/* 21, 30 */
		{0.243085e-6, -775.522611324, 3.651837925},
		{0.173435e-6, 18849.227549974, 6.153743485},
		{0.230685e-6, 5856.477659115, 4.773852582},
		{0.203747e-6, 12036.460734888, 4.333987818},
		{0.143935e-6, -796.298006816, 5.957517795},
		{0.159080e-6, 10977.078804699, 1.890075226},
		{0.119979e-6, 38.133035638, 4.551585768},
		{0.118971e-6, 5486.777843175, 1.914547226},
		{0.116120e-6, 1059.381930189, 0.873504123},
		{0.137927e-6, 11790.629088659, 1.135934669},
		/* 31, 40 */
		{0.098358e-6, 2544.314419883, 0.092793886},
		{0.101868e-6, -5573.142801634, 5.984503847},
		{0.080164e-6, 206.185548437, 2.095377709},
		{0.079645e-6, 4694.002954708, 2.949233637},
		{0.062617e-6, 20.775395492, 2.654394814},
		{0.075019e-6, 2942.463423292, 4.980931759},
		{0.064397e-6, 5746.271337896, 1.280308748},
		{0.063814e-6, 5760.498431898, 4.167901731},
		{0.048042e-6, 2146.165416475, 1.495846011},
		{0.048373e-6, 155.420399434, 2.251573730},
		/* 41, 50 */
		{0.058844e-6, 426.598190876, 4.839650148},
		{0.046551e-6, -0.980321068, 0.921573539},
		{0.054139e-6, 17260.154654690, 3.411091093},
		{0.042411e-6, 6275.962302991, 2.869567043},
		{0.040184e-6, -7.113547001, 3.565975565},
		{0.036564e-6, 5088.628839767, 3.324679049},
		{0.040759e-6, 12352.852604545, 3.981496998},
		{0.036507e-6, 801.820931124, 6.248866009},
		{0.036955e-6, 3154.687084896, 5.071801441},
		{0.042732e-6, 632.783739313, 5.720622217},
		/* 51, 60 */
		{0.042560e-6, 161000.685737473, 1.270837679},
		{0.040480e-6, 15720.838784878, 2.546610123},
		{0.028244e-6, -6286.598968340, 5.069663519},
		{0.033477e-6, 6062.663207553, 4.144987272},
		{0.034867e-6, 522.577418094, 5.210064075},
		{0.032438e-6, 6076.890301554, 0.749317412},
		{0.030215e-6, 7084.896781115, 3.389610345},
		{0.029247e-6, -71430.695617928, 4.183178762},
		{0.033529e-6, 9437.762934887, 2.404714239},
		{0.032423e-6, 8827.390269875, 5.541473556},
		/* 61, 70 */
		{0.027567e-6, 6279.552731642, 5.040846034},
		{0.029862e-6, 12139.553509107, 1.770181024},
		{0.022509e-6, 10447.387839604, 1.460726241},
		{0.020937e-6, 8429.241266467, 0.652303414},
		{0.020322e-6, 419.484643875, 3.735430632},
		{0.024816e-6, -1194.447010225, 1.087136918},
		{0.025196e-6, 1748.016413067, 2.901883301},
		{0.021691e-6, 14143.495242431, 5.952658009},
		{0.017673e-6, 6812.766815086, 3.186129845},
		{0.022567e-6, 6133.512652857, 3.307984806},
		/* 71, 80 */
		{0.016155e-6, 10213.285546211, 1.331103168},
		{0.014751e-6, 1349.867409659, 4.308933301},
		{0.015949e-6, -220.412642439, 4.005298270},
		{0.015974e-6, -2352.866153772, 6.145309371},
		{0.014223e-6, 17789.845619785, 2.104551349},
		{0.017806e-6, 73.297125859, 3.475975097},
		{0.013671e-6, -536.804512095, 5.971672571},
		{0.011942e-6, 8031.092263058, 2.053414715},
		{0.014318e-6, 16730.463689596, 3.016058075},
		{0.012462e-6, 103.092774219, 1.737438797},
		/* 81, 90 */
		{0.010962e-6, 3.590428652, 2.196567739},
		{0.015078e-6, 19651.048481098, 3.969480770},
		{0.010396e-6, 951.718406251, 5.717799605},
		{0.011707e-6, -4705.732307544, 2.654125618},
		{0.010453e-6, 5863.591206116, 1.913704550},
		{0.012420e-6, 4690.479836359, 4.734090399},
		{0.011847e-6, 5643.178563677, 5.489005403},
		{0.008610e-6, 3340.612426700, 3.661698944},
		{0.011622e-6, 5120.601145584, 4.863931876},
		{0.010825e-6, 553.569402842, 0.842715011},
		/* 91, 100 */
		{0.008666e-6, -135.065080035, 3.293406547},
		{0.009963e-6, 149.563197135, 4.870690598},
		{0.009858e-6, 6309.374169791, 1.061816410},
		{0.007959e-6, 316.391869657, 2.465042647},
		{0.010099e-6, 283.859318865, 1.942176992},
		{0.007147e-6, -242.728603974, 3.661486981},
		{0.007505e-6, 5230.807466803, 4.920937029},
		{0.008323e-6, 11769.853693166, 1.229392026},
		{0.007490e-6, -6256.777530192, 3.658444681},
		{0.009370e-6, 149854.400134205, 0.673880395},
		/* 101, 110 */
		{0.007117e-6, 38.027672636, 5.294249518},
		{0.007857e-6, 12168.002696575, 0.525733528},
		{0.007019e-6, 6206.809778716, 0.837688810},
		{0.006056e-6, 955.599741609, 4.194535082},
		{0.008107e-6, 13367.972631107, 3.793235253},
		{0.006731e-6, 5650.292110678, 5.639906583},
		{0.007332e-6, 36.648562930, 0.114858677},
		{0.006366e-6, 4164.311989613, 2.262081818},
		{0.006858e-6, 5216.580372801, 0.642063318},
		{0.006919e-6, 6681.224853400, 6.018501522},
		/* 111, 120 */
		{0.006826e-6, 7632.943259650, 3.458654112},
		{0.005308e-6, -1592.596013633, 2.500382359},
		{0.005096e-6, 11371.704689758, 2.547107806},
		{0.004841e-6, 5333.900241022, 0.437078094},
		{0.005582e-6, 5966.683980335, 2.246174308},
		{0.006304e-6, 11926.254413669, 2.512929171},
		{0.006603e-6, 23581.258177318, 5.393136889},
		{0.005123e-6, -1.484472708, 2.999641028},
		{0.004648e-6, 1589.072895284, 1.275847090},
		{0.005119e-6, 6438.496249426, 1.486539246},
		/* 121, 130 */
		{0.004521e-6, 4292.330832950, 6.140635794},
		{0.005680e-6, 23013.539539587, 4.557814849},
		{0.005488e-6, -3.455808046, 0.090675389},
		{0.004193e-6, 7234.794256242, 4.869091389},
		{0.003742e-6, 7238.675591600, 4.691976180},
		{0.004148e-6, -110.206321219, 3.016173439},
		{0.004553e-6, 11499.656222793, 5.554998314},
		{0.004892e-6, 5436.993015240, 1.475415597},
		{0.004044e-6, 4732.030627343, 1.398784824},
		{0.004164e-6, 12491.370101415, 5.650931916},
		/* 131, 140 */
		{0.004349e-6, 11513.883316794, 2.181745369},
		{0.003919e-6, 12528.018664345, 5.823319737},
		{0.003129e-6, 6836.645252834, 0.003844094},
		{0.004080e-6, -7058.598461315, 3.690360123},
		{0.003270e-6, 76.266071276, 1.517189902},
		{0.002954e-6, 6283.143160294, 4.447203799},
		{0.002872e-6, 28.449187468, 1.158692983},
		{0.002881e-6, 735.876513532, 0.349250250},
		{0.003279e-6, 5849.364112115, 4.893384368},
		{0.003625e-6, 6209.778724132, 1.473760578},
		/* 141, 150 */
		{0.003074e-6, 949.175608970, 5.185878737},
		{0.002775e-6, 9917.696874510, 1.030026325},
		{0.002646e-6, 10973.555686350, 3.918259169},
		{0.002575e-6, 25132.303399966, 6.109659023},
		{0.003500e-6, 263.083923373, 1.892100742},
		{0.002740e-6, 18319.536584880, 4.320519510},
		{0.002464e-6, 202.253395174, 4.698203059},
		{0.002409e-6, 2.542797281, 5.325009315},
		{0.003354e-6, -90955.551694697, 1.942656623},
		{0.002296e-6, 6496.374945429, 5.061810696},
		/* 151, 160 */
		{0.003002e-6, 6172.869528772, 2.797822767},
		{0.003202e-6, 27511.467873537, 0.531673101},
		{0.002954e-6, -6283.008539689, 4.533471191},
		{0.002353e-6, 639.897286314, 3.734548088},
		{0.002401e-6, 16200.772724501, 2.605547070},
		{0.003053e-6, 233141.314403759, 3.029030662},
		{0.003024e-6, 83286.914269554, 2.355556099},
		{0.002863e-6, 17298.182327326, 5.240963796},
		{0.002103e-6, -7079.373856808, 5.756641637},
		{0.002303e-6, 83996.847317911, 2.013686814},
		/* 161, 170 */
		{0.002303e-6, 18073.704938650, 1.089100410},
		{0.002381e-6, 63.735898303, 0.759188178},
		{0.002493e-6, 6386.168624210, 0.645026535},
		{0.002366e-6, 3.932153263, 6.215885448},
		{0.002169e-6, 11015.106477335, 4.845297676},
		{0.002397e-6, 6243.458341645, 3.809290043},
		{0.002183e-6, 1162.474704408, 6.179611691},
		{0.002353e-6, 6246.427287062, 4.781719760},
		{0.002199e-6, -245.831646229, 5.956152284},
		{0.001729e-6, 3894.181829542, 1.264976635},
		/* 171, 180 */
		{0.001896e-6, -3128.388765096, 4.914231596},
		{0.002085e-6, 35.164090221, 1.405158503},
		{0.002024e-6, 14712.317116458, 2.752035928},
		{0.001737e-6, 6290.189396992, 5.280820144},
		{0.002229e-6, 491.557929457, 1.571007057},
		{0.001602e-6, 14314.168113050, 4.203664806},
		{0.002186e-6, 454.909366527, 1.402101526},
		{0.001897e-6, 22483.848574493, 4.167932508},
		{0.001825e-6, -3738.761430108, 0.545828785},
		{0.001894e-6, 1052.268383188, 5.817167450},
		/* 181, 190 */
		{0.001421e-6, 20.355319399, 2.419886601},
		{0.001408e-6, 10984.192351700, 2.732084787},
		{0.001847e-6, 10873.986030480, 2.903477885},
		{0.001391e-6, -8635.942003763, 0.593891500},
		{0.001388e-6, -7.046236698, 1.166145902},
		{0.001810e-6, -88860.057071188, 0.487355242},
		{0.001288e-6, -1990.745017041, 3.913022880},
		{0.001297e-6, 23543.230504682, 3.063805171},
		{0.001335e-6, -266.607041722, 3.995764039},
		{0.001376e-6, 10969.965257698, 5.152914309},
		/* 191, 200 */
		{0.001745e-6, 244287.600007027, 3.626395673},
		{0.001649e-6, 31441.677569757, 1.952049260},
		{0.001416e-6, 9225.539273283, 4.996408389},
		{0.001238e-6, 4804.209275927, 5.503379738},
		{0.001472e-6, 4590.910180489, 4.164913291},
		{0.001169e-6, 6040.347246017, 5.841719038},
		{0.001039e-6, 5540.085789459, 2.769753519},
		{0.001004e-6, -170.672870619, 0.755008103},
		{0.001284e-6, 10575.406682942, 5.306538209},
		{0.001278e-6, 71.812653151, 4.713486491},
		/* 201, 210 */
		{0.001321e-6, 18209.330263660, 2.624866359},
		{0.001297e-6, 21228.392023546, 0.382603541},
		{0.000954e-6, 6282.095528923, 0.882213514},
		{0.001145e-6, 6058.731054289, 1.169483931},
		{0.000979e-6, 5547.199336460, 5.448375984},
		{0.000987e-6, -6262.300454499, 2.656486959},
		{0.001070e-6, -154717.609887482, 1.827624012},
		{0.000991e-6, 4701.116501708, 4.387001801},
		{0.001155e-6, -14.227094002, 3.042700750},
		{0.001176e-6, 277.034993741, 3.335519004},
		/* 211, 220 */
		{0.000890e-6, 13916.019109642, 5.601498297},
		{0.000884e-6, -1551.045222648, 1.088831705},
		{0.000876e-6, 5017.508371365, 3.969902609},
		{0.000806e-6, 15110.466119866, 5.142876744},
		{0.000773e-6, -4136.910433516, 0.022067765},
		{0.001077e-6, 175.166059800, 1.844913056},
		{0.000954e-6, -6284.056171060, 0.968480906},
		{0.000737e-6, 5326.786694021, 4.923831588},
		{0.000845e-6, -433.711737877, 4.749245231},
		{0.000819e-6, 8662.240323563, 5.991247817},
		/* 221, 230 */
		{0.000852e-6, 199.072001436, 2.189604979},
		{0.000723e-6, 17256.631536341, 6.068719637},
		{0.000940e-6, 6037.244203762, 6.197428148},
		{0.000885e-6, 11712.955318231, 3.280414875},
		{0.000706e-6, 12559.038152982, 2.824848947},
		{0.000732e-6, 2379.164473572, 2.501813417},
		{0.000764e-6, -6127.655450557, 2.236346329},
		{0.000908e-6, 131.541961686, 2.521257490},
		{0.000907e-6, 35371.887265976, 3.370195967},
		{0.000673e-6, 1066.495477190, 3.876512374},
		/* 231, 240 */
		{0.000814e-6, 17654.780539750, 4.627122566},
		{0.000630e-6, 36.027866677, 0.156368499},
		{0.000798e-6, 515.463871093, 5.151962502},
		{0.000798e-6, 148.078724426, 5.909225055},
		{0.000806e-6, 309.278322656, 6.054064447},
		{0.000607e-6, -39.617508346, 2.839021623},
		{0.000601e-6, 412.371096874, 3.984225404},
		{0.000646e-6, 11403.676995575, 3.852959484},
		{0.000704e-6, 13521.751441591, 2.300991267},
		{0.000603e-6, -65147.619767937, 4.140083146},
		/* 241, 250 */
		{0.000609e-6, 10177.257679534, 0.437122327},
		{0.000631e-6, 5767.611978898, 4.026532329},
		{0.000576e-6, 11087.285125918, 4.760293101},
		{0.000674e-6, 14945.316173554, 6.270510511},
		{0.000726e-6, 5429.879468239, 6.039606892},
		{0.000710e-6, 28766.924424484, 5.672617711},
		{0.000647e-6, 11856.218651625, 3.397132627},
		{0.000678e-6, -5481.254918868, 6.249666675},
		{0.000618e-6, 22003.914634870, 2.466427018},
		{0.000738e-6, 6134.997125565, 2.242668890},
		/* 251, 260 */
		{0.000660e-6, 625.670192312, 5.864091907},
		{0.000694e-6, 3496.032826134, 2.668309141},
		{0.000531e-6, 6489.261398429, 1.681888780},
		{0.000611e-6, -143571.324284214, 2.424978312},
		{0.000575e-6, 12043.574281889, 4.216492400},
		{0.000553e-6, 12416.588502848, 4.772158039},
		{0.000689e-6, 4686.889407707, 6.224271088},
		{0.000495e-6, 7342.457780181, 3.817285811},
		{0.000567e-6, 3634.621024518, 1.649264690},
		{0.000515e-6, 18635.928454536, 3.945345892},
		/* 261, 270 */
		{0.000486e-6, -323.505416657, 4.061673868},
		{0.000662e-6, 25158.601719765, 1.794058369},
		{0.000509e-6, 846.082834751, 3.053874588},
		{0.000472e-6, -12569.674818332, 5.112133338},
		{0.000461e-6, 6179.983075773, 0.513669325},
		{0.000641e-6, 83467.156352816, 3.210727723},
		{0.000520e-6, 10344.295065386, 2.445597761},
		{0.000493e-6, 18422.629359098, 1.676939306},
		{0.000478e-6, 1265.567478626, 5.487314569},
		{0.000472e-6, -18.159247265, 1.999707589},
		/* 271, 280 */
		{0.000559e-6, 11190.377900137, 5.783236356},
		{0.000494e-6, 9623.688276691, 3.022645053},
		{0.000463e-6, 5739.157790895, 1.411223013},
		{0.000432e-6, 16858.482532933, 1.179256434},
		{0.000574e-6, 72140.628666286, 1.758191830},
		{0.000484e-6, 17267.268201691, 3.290589143},
		{0.000550e-6, 4907.302050146, 0.864024298},
		{0.000399e-6, 14.977853527, 2.094441910},
		{0.000491e-6, 224.344795702, 0.878372791},
		{0.000432e-6, 20426.571092422, 6.003829241},
		/* 281, 290 */
		{0.000481e-6, 5749.452731634, 4.309591964},
		{0.000480e-6, 5757.317038160, 1.142348571},
		{0.000485e-6, 6702.560493867, 0.210580917},
		{0.000426e-6, 6055.549660552, 4.274476529},
		{0.000480e-6, 5959.570433334, 5.031351030},
		{0.000466e-6, 12562.628581634, 4.959581597},
		{0.000520e-6, 39302.096962196, 4.788002889},
		{0.000458e-6, 12132.439962106, 1.880103788},
		{0.000470e-6, 12029.347187887, 1.405611197},
		{0.000416e-6, -7477.522860216, 1.082356330},
		/* 291, 300 */
		{0.000449e-6, 11609.862544012, 4.179989585},
		{0.000465e-6, 17253.041107690, 0.353496295},
		{0.000362e-6, -4535.059436924, 1.583849576},
		{0.000383e-6, 21954.157609398, 3.747376371},
		{0.000389e-6, 17.252277143, 1.395753179},
		{0.000331e-6, 18052.929543158, 0.566790582},
		{0.000430e-6, 13517.870106233, 0.685827538},
		{0.000368e-6, -5756.908003246, 0.731374317},
		{0.000330e-6, 10557.594160824, 3.710043680},
		{0.000332e-6, 20199.094959633, 1.652901407},
		/* 301, 310 */
		{0.000384e-6, 11933.367960670, 5.827781531},
		{0.000387e-6, 10454.501386605, 2.541182564},
		{0.000325e-6, 15671.081759407, 2.178850542},
		{0.000318e-6, 138.517496871, 2.253253037},
		{0.000305e-6, 9388.005909415, 0.578340206},
		{0.000352e-6, 5749.861766548, 3.000297967},
		{0.000311e-6, 6915.859589305, 1.693574249},
		{0.000297e-6, 24072.921469776, 1.997249392},
		{0.000363e-6, -640.877607382, 5.071820966},
		{0.000323e-6, 12592.450019783, 1.072262823},
		/* 311, 320 */
		{0.000341e-6, 12146.667056108, 4.700657997},
		{0.000290e-6, 9779.108676125, 1.812320441},
		{0.000342e-6, 6132.028180148, 4.322238614},
		{0.000329e-6, 6268.848755990, 3.033827743},
		{0.000374e-6, 17996.031168222, 3.388716544},
		{0.000285e-6, -533.214083444, 4.687313233},
		{0.000338e-6, 6065.844601290, 0.877776108},
		{0.000276e-6, 24.298513841, 0.770299429},
		{0.000336e-6, -2388.894020449, 5.353796034},
		{0.000290e-6, 3097.883822726, 4.075291557},
		/* 321, 330 */
		{0.000318e-6, 709.933048357, 5.941207518},
		{0.000271e-6, 13095.842665077, 3.208912203},
		{0.000331e-6, 6073.708907816, 4.007881169},
		{0.000292e-6, 742.990060533, 2.714333592},
		{0.000362e-6, 29088.811415985, 3.215977013},
		{0.000280e-6, 12359.966151546, 0.710872502},
		{0.000267e-6, 10440.274292604, 4.730108488},
		{0.000262e-6, 838.969287750, 1.327720272},
		{0.000250e-6, 16496.361396202, 0.898769761},
		{0.000325e-6, 20597.243963041, 0.180044365},
		/* 331, 340 */
		{0.000268e-6, 6148.010769956, 5.152666276},
		{0.000284e-6, 5636.065016677, 5.655385808},
		{0.000301e-6, 6080.822454817, 2.135396205},
		{0.000294e-6, -377.373607916, 3.708784168},
		{0.000236e-6, 2118.763860378, 1.733578756},
		{0.000234e-6, 5867.523359379, 5.575209112},
		{0.000268e-6, -226858.238553767, 0.069432392},
		{0.000265e-6, 167283.761587465, 4.369302826},
		{0.000280e-6, 28237.233459389, 5.304829118},
		{0.000292e-6, 12345.739057544, 4.096094132},
		/* 341, 350 */
		{0.000223e-6, 19800.945956225, 3.069327406},
		{0.000301e-6, 43232.306658416, 6.205311188},
		{0.000264e-6, 18875.525869774, 1.417263408},
		{0.000304e-6, -1823.175188677, 3.409035232},
		{0.000301e-6, 109.945688789, 0.510922054},
		{0.000260e-6, 813.550283960, 2.389438934},
		{0.000299e-6, 316428.228673312, 5.384595078},
		{0.000211e-6, 5756.566278634, 3.789392838},
		{0.000209e-6, 5750.203491159, 1.661943545},
		{0.000240e-6, 12489.885628707, 5.684549045},
		/* 351, 360 */
		{0.000216e-6, 6303.851245484, 3.862942261},
		{0.000203e-6, 1581.959348283, 5.549853589},
		{0.000200e-6, 5642.198242609, 1.016115785},
		{0.000197e-6, -70.849445304, 4.690702525},
		{0.000227e-6, 6287.008003254, 2.911891613},
		{0.000197e-6, 533.623118358, 1.048982898},
		{0.000205e-6, -6279.485421340, 1.829362730},
		{0.000209e-6, -10988.808157535, 2.636140084},
		{0.000208e-6, -227.526189440, 4.127883842},
		{0.000191e-6, 415.552490612, 4.401165650},
		/* 361, 370 */
		{0.000190e-6, 29296.615389579, 4.175658539},
		{0.000264e-6, 66567.485864652, 4.601102551},
		{0.000256e-6, -3646.350377354, 0.506364778},
		{0.000188e-6, 13119.721102825, 2.032195842},
		{0.000185e-6, -209.366942175, 4.694756586},
		{0.000198e-6, 25934.124331089, 3.832703118},
		{0.000195e-6, 4061.219215394, 3.308463427},
		{0.000234e-6, 5113.487598583, 1.716090661},
		{0.000188e-6, 1478.866574064, 5.686865780},
		{0.000222e-6, 11823.161639450, 1.942386641},
		/* 371, 380 */
		{0.000181e-6, 10770.893256262, 1.999482059},
		{0.000171e-6, 6546.159773364, 1.182807992},
		{0.000206e-6, 70.328180442, 5.934076062},
		{0.000169e-6, 20995.392966449, 2.169080622},
		{0.000191e-6, 10660.686935042, 5.405515999},
		{0.000228e-6, 33019.021112205, 4.656985514},
		{0.000184e-6, -4933.208440333, 3.327476868},
		{0.000220e-6, -135.625325010, 1.765430262},
		{0.000166e-6, 23141.558382925, 3.454132746},
		{0.000191e-6, 6144.558353121, 5.020393445},
		/* 381, 390 */
		{0.000180e-6, 6084.003848555, 0.602182191},
		{0.000163e-6, 17782.732072784, 4.960593133},
		{0.000225e-6, 16460.333529525, 2.596451817},
		{0.000222e-6, 5905.702242076, 3.731990323},
		{0.000204e-6, 227.476132789, 5.636192701},
		{0.000159e-6, 16737.577236597, 3.600691544},
		{0.000200e-6, 6805.653268085, 0.868220961},
		{0.000187e-6, 11919.140866668, 2.629456641},
		{0.000161e-6, 127.471796607, 2.862574720},
		{0.000205e-6, 6286.666278643, 1.742882331},
		/* 391, 400 */
		{0.000189e-6, 153.778810485, 4.812372643},
		{0.000168e-6, 16723.350142595, 0.027860588},
		{0.000149e-6, 11720.068865232, 0.659721876},
		{0.000189e-6, 5237.921013804, 5.245313000},
		{0.000143e-6, 6709.674040867, 4.317625647},
		{0.000146e-6, 4487.817406270, 4.815297007},
		{0.000144e-6, -664.756045130, 5.381366880},
		{0.000175e-6, 5127.714692584, 4.728443327},
		{0.000162e-6, 6254.626662524, 1.435132069},
		{0.000187e-6, 47162.516354635, 1.354371923},
		/* 401, 410 */
		{0.000146e-6, 11080.171578918, 3.369695406},
		{0.000180e-6, -348.924420448, 2.490902145},
		{0.000148e-6, 151.047669843, 3.799109588},
		{0.000157e-6, 6197.248551160, 1.284375887},
		{0.000167e-6, 146.594251718, 0.759969109},
		{0.000133e-6, -5331.357443741, 5.409701889},
		{0.000154e-6, 95.979227218, 3.366890614},
		{0.000148e-6, -6418.140930027, 3.384104996},
		{0.000128e-6, -6525.804453965, 3.803419985},
		{0.000130e-6, 11293.470674356, 0.939039445},
		/* 411, 420 */
		{0.000152e-6, -5729.506447149, 0.734117523},
		{0.000138e-6, 210.117701700, 2.564216078},
		{0.000123e-6, 6066.595360816, 4.517099537},
		{0.000140e-6, 18451.078546566, 0.642049130},
		{0.000126e-6, 11300.584221356, 3.485280663},
		{0.000119e-6, 10027.903195729, 3.217431161},
		{0.000151e-6, 4274.518310832, 4.404359108},
		{0.000117e-6, 6072.958148291, 0.366324650},
		{0.000165e-6, -7668.637425143, 4.298212528},
		{0.000117e-6, -6245.048177356, 5.379518958},
		/* 421, 430 */
		{0.000130e-6, -5888.449964932, 4.527681115},
		{0.000121e-6, -543.918059096, 6.109429504},
		{0.000162e-6, 9683.594581116, 5.720092446},
		{0.000141e-6, 6219.339951688, 0.679068671},
		{0.000118e-6, 22743.409379516, 4.881123092},
		{0.000129e-6, 1692.165669502, 0.351407289},
		{0.000126e-6, 5657.405657679, 5.146592349},
		{0.000114e-6, 728.762966531, 0.520791814},
		{0.000120e-6, 52.596639600, 0.948516300},
		{0.000115e-6, 65.220371012, 3.504914846},
		/* 431, 440 */
		{0.000126e-6, 5881.403728234, 5.577502482},
		{0.000158e-6, 163096.180360983, 2.957128968},
		{0.000134e-6, 12341.806904281, 2.598576764},
		{0.000151e-6, 16627.370915377, 3.985702050},
		{0.000109e-6, 1368.660252845, 0.014730471},
		{0.000131e-6, 6211.263196841, 0.085077024},
		{0.000146e-6, 5792.741760812, 0.708426604},
		{0.000146e-6, -77.750543984, 3.121576600},
		{0.000107e-6, 5341.013788022, 0.288231904},
		{0.000138e-6, 6281.591377283, 2.797450317},
		/* 441, 450 */
		{0.000113e-6, -6277.552925684, 2.788904128},
		{0.000115e-6, -525.758811831, 5.895222200},
		{0.000138e-6, 6016.468808270, 6.096188999},
		{0.000139e-6, 23539.707386333, 2.028195445},
		{0.000146e-6, -4176.041342449, 4.660008502},
		{0.000107e-6, 16062.184526117, 4.066520001},
		{0.000142e-6, 83783.548222473, 2.936315115},
		{0.000128e-6, 9380.959672717, 3.223844306},
		{0.000135e-6, 6205.325306007, 1.638054048},
		{0.000101e-6, 2699.734819318, 5.481603249},
		/* 451, 460 */
		{0.000104e-6, -568.821874027, 2.205734493},
		{0.000103e-6, 6321.103522627, 2.440421099},
		{0.000119e-6, 6321.208885629, 2.547496264},
		{0.000138e-6, 1975.492545856, 2.314608466},
		{0.000121e-6, 137.033024162, 4.539108237},
		{0.000123e-6, 19402.796952817, 4.538074405},
		{0.000119e-6, 22805.735565994, 2.869040566},
		{0.000133e-6, 64471.991241142, 6.056405489},
		{0.000129e-6, -85.827298831, 2.540635083},
		{0.000131e-6, 13613.804277336, 4.005732868},
		/* 461, 470 */
		{0.000104e-6, 9814.604100291, 1.959967212},
		{0.000112e-6, 16097.679950283, 3.589026260},
		{0.000123e-6, 2107.034507542, 1.728627253},
		{0.000121e-6, 36949.230808424, 6.072332087},
		{0.000108e-6, -12539.853380183, 3.716133846},
		{0.000113e-6, -7875.671863624, 2.725771122},
		{0.000109e-6, 4171.425536614, 4.033338079},
		{0.000101e-6, 6247.911759770, 3.441347021},
		{0.000113e-6, 7330.728427345, 0.656372122},
		{0.000113e-6, 51092.726050855, 2.791483066},
		/* 471, 480 */
		{0.000106e-6, 5621.842923210, 1.815323326},
		{0.000101e-6, 111.430161497, 5.711033677},
		{0.000103e-6, 909.818733055, 2.812745443},
		{0.000101e-6, 1790.642637886, 1.965746028},

		/* T */
		{102.156724e-6, 6283.075849991, 4.249032005},
		{1.706807e-6, 12566.151699983, 4.205904248},
		{0.269668e-6, 213.299095438, 3.400290479},
		{0.265919e-6, 529.690965095, 5.836047367},
		{0.210568e-6, -3.523118349, 6.262738348},
		{0.077996e-6, 5223.693919802, 4.670344204},
		/* 481, 490 */
		{0.054764e-6, 1577.343542448, 4.534800170},
		{0.059146e-6, 26.298319800, 1.083044735},
		{0.034420e-6, -398.149003408, 5.980077351},
		{0.032088e-6, 18849.227549974, 4.162913471},
		{0.033595e-6, 5507.553238667, 5.980162321},
		{0.029198e-6, 5856.477659115, 0.623811863},
		{0.027764e-6, 155.420399434, 3.745318113},
		{0.025190e-6, 5746.271337896, 2.980330535},
		{0.022997e-6, -796.298006816, 1.174411803},
		{0.024976e-6, 5760.498431898, 2.467913690},
		/* 491, 500 */
		{0.021774e-6, 206.185548437, 3.854787540},
		{0.017925e-6, -775.522611324, 1.092065955},
		{0.013794e-6, 426.598190876, 2.699831988},
		{0.013276e-6, 6062.663207553, 5.845801920},
		{0.011774e-6, 12036.460734888, 2.292832062},
		{0.012869e-6, 6076.890301554, 5.333425680},
		{0.012152e-6, 1059.381930189, 6.222874454},
		{0.011081e-6, -7.113547001, 5.154724984},
		{0.010143e-6, 4694.002954708, 4.044013795},
		{0.009357e-6, 5486.777843175, 3.416081409},
		/* 501, 510 */
		{0.010084e-6, 522.577418094, 0.749320262},
		{0.008587e-6, 10977.078804699, 2.777152598},
		{0.008628e-6, 6275.962302991, 4.562060226},
		{0.008158e-6, -220.412642439, 5.806891533},
		{0.007746e-6, 2544.314419883, 1.603197066},
		{0.007670e-6, 2146.165416475, 3.000200440},
		{0.007098e-6, 74.781598567, 0.443725817},
		{0.006180e-6, -536.804512095, 1.302642751},
		{0.005818e-6, 5088.628839767, 4.827723531},
		{0.004945e-6, -6286.598968340, 0.268305170},
		/* 511, 520 */
		{0.004774e-6, 1349.867409659, 5.808636673},
		{0.004687e-6, -242.728603974, 5.154890570},
		{0.006089e-6, 1748.016413067, 4.403765209},
		{0.005975e-6, -1194.447010225, 2.583472591},
		{0.004229e-6, 951.718406251, 0.931172179},
		{0.005264e-6, 553.569402842, 2.336107252},
		{0.003049e-6, 5643.178563677, 1.362634430},
		{0.002974e-6, 6812.766815086, 1.583012668},
		{0.003403e-6, -2352.866153772, 2.552189886},
		{0.003030e-6, 419.484643875, 5.286473844},
		/* 521, 530 */
		{0.003210e-6, -7.046236698, 1.863796539},
		{0.003058e-6, 9437.762934887, 4.226420633},
		{0.002589e-6, 12352.852604545, 1.991935820},
		{0.002927e-6, 5216.580372801, 2.319951253},
		{0.002425e-6, 5230.807466803, 3.084752833},
		{0.002656e-6, 3154.687084896, 2.487447866},
		{0.002445e-6, 10447.387839604, 2.347139160},
		{0.002990e-6, 4690.479836359, 6.235872050},
		{0.002890e-6, 5863.591206116, 0.095197563},
		{0.002498e-6, 6438.496249426, 2.994779800},
		/* 531, 540 */
		{0.001889e-6, 8031.092263058, 3.569003717},
		{0.002567e-6, 801.820931124, 3.425611498},
		{0.001803e-6, -71430.695617928, 2.192295512},
		{0.001782e-6, 3.932153263, 5.180433689},
		{0.001694e-6, -4705.732307544, 4.641779174},
		{0.001704e-6, -1592.596013633, 3.997097652},
		{0.001735e-6, 5849.364112115, 0.417558428},
		{0.001643e-6, 8429.241266467, 2.180619584},
		{0.001680e-6, 38.133035638, 4.164529426},
		{0.002045e-6, 7084.896781115, 0.526323854},
		/* 541, 550 */
		{0.001458e-6, 4292.330832950, 1.356098141},
		{0.001437e-6, 20.355319399, 3.895439360},
		{0.001738e-6, 6279.552731642, 0.087484036},
		{0.001367e-6, 14143.495242431, 3.987576591},
		{0.001344e-6, 7234.794256242, 0.090454338},
		{0.001438e-6, 11499.656222793, 0.974387904},
		{0.001257e-6, 6836.645252834, 1.509069366},
		{0.001358e-6, 11513.883316794, 0.495572260},
		{0.001628e-6, 7632.943259650, 4.968445721},
		{0.001169e-6, 103.092774219, 2.838496795},
		/* 551, 560 */
		{0.001162e-6, 4164.311989613, 3.408387778},
		{0.001092e-6, 6069.776754553, 3.617942651},
		{0.001008e-6, 17789.845619785, 0.286350174},
		{0.001008e-6, 639.897286314, 1.610762073},
		{0.000918e-6, 10213.285546211, 5.532798067},
		{0.001011e-6, -6256.777530192, 0.661826484},
		{0.000753e-6, 16730.463689596, 3.905030235},
		{0.000737e-6, 11926.254413669, 4.641956361},
		{0.000694e-6, 3340.612426700, 2.111120332},
		{0.000701e-6, 3894.181829542, 2.760823491},
		/* 561, 570 */
		{0.000689e-6, -135.065080035, 4.768800780},
		{0.000700e-6, 13367.972631107, 5.760439898},
		{0.000664e-6, 6040.347246017, 1.051215840},
		{0.000654e-6, 5650.292110678, 4.911332503},
		{0.000788e-6, 6681.224853400, 4.699648011},
		{0.000628e-6, 5333.900241022, 5.024608847},
		{0.000755e-6, -110.206321219, 4.370971253},
		{0.000628e-6, 6290.189396992, 3.660478857},
		{0.000635e-6, 25132.303399966, 4.121051532},
		{0.000534e-6, 5966.683980335, 1.173284524},
		/* 571, 580 */
		{0.000543e-6, -433.711737877, 0.345585464},
		{0.000517e-6, -1990.745017041, 5.414571768},
		{0.000504e-6, 5767.611978898, 2.328281115},
		{0.000485e-6, 5753.384884897, 1.685874771},
		{0.000463e-6, 7860.419392439, 5.297703006},
		{0.000604e-6, 515.463871093, 0.591998446},
		{0.000443e-6, 12168.002696575, 4.830881244},
		{0.000570e-6, 199.072001436, 3.899190272},
		{0.000465e-6, 10969.965257698, 0.476681802},
		{0.000424e-6, -7079.373856808, 1.112242763},
		/* 581, 590 */
		{0.000427e-6, 735.876513532, 1.994214480},
		{0.000478e-6, -6127.655450557, 3.778025483},
		{0.000414e-6, 10973.555686350, 5.441088327},
		{0.000512e-6, 1589.072895284, 0.107123853},
		{0.000378e-6, 10984.192351700, 0.915087231},
		{0.000402e-6, 11371.704689758, 4.107281715},
		{0.000453e-6, 9917.696874510, 1.917490952},
		{0.000395e-6, 149.563197135, 2.763124165},
		{0.000371e-6, 5739.157790895, 3.112111866},
		{0.000350e-6, 11790.629088659, 0.440639857},
		/* 591, 600 */
		{0.000356e-6, 6133.512652857, 5.444568842},
		{0.000344e-6, 412.371096874, 5.676832684},
		{0.000383e-6, 955.599741609, 5.559734846},
		{0.000333e-6, 6496.374945429, 0.261537984},
		{0.000340e-6, 6055.549660552, 5.975534987},
		{0.000334e-6, 1066.495477190, 2.335063907},
		{0.000399e-6, 11506.769769794, 5.321230910},
		{0.000314e-6, 18319.536584880, 2.313312404},
		{0.000424e-6, 1052.268383188, 1.211961766},
		{0.000307e-6, 63.735898303, 3.169551388},
		/* 601, 610 */
		{0.000329e-6, 29.821438149, 6.106912080},
		{0.000357e-6, 6309.374169791, 4.223760346},
		{0.000312e-6, -3738.761430108, 2.180556645},
		{0.000301e-6, 309.278322656, 1.499984572},
		{0.000268e-6, 12043.574281889, 2.447520648},
		{0.000257e-6, 12491.370101415, 3.662331761},
		{0.000290e-6, 625.670192312, 1.272834584},
		{0.000256e-6, 5429.879468239, 1.913426912},
		{0.000339e-6, 3496.032826134, 4.165930011},
		{0.000283e-6, 3930.209696220, 4.325565754},
		/* 611, 620 */
		{0.000241e-6, 12528.018664345, 3.832324536},
		{0.000304e-6, 4686.889407707, 1.612348468},
		{0.000259e-6, 16200.772724501, 3.470173146},
		{0.000238e-6, 12139.553509107, 1.147977842},
		{0.000236e-6, 6172.869528772, 3.776271728},
		{0.000296e-6, -7058.598461315, 0.460368852},
		{0.000306e-6, 10575.406682942, 0.554749016},
		{0.000251e-6, 17298.182327326, 0.834332510},
		{0.000290e-6, 4732.030627343, 4.759564091},
		{0.000261e-6, 5884.926846583, 0.298259862},
		/* 621, 630 */
		{0.000249e-6, 5547.199336460, 3.749366406},
		{0.000213e-6, 11712.955318231, 5.415666119},
		{0.000223e-6, 4701.116501708, 2.703203558},
		{0.000268e-6, -640.877607382, 0.283670793},
		{0.000209e-6, 5636.065016677, 1.238477199},
		{0.000193e-6, 10177.257679534, 1.943251340},
		{0.000182e-6, 6283.143160294, 2.456157599},
		{0.000184e-6, -227.526189440, 5.888038582},
		{0.000182e-6, -6283.008539689, 0.241332086},
		{0.000228e-6, -6284.056171060, 2.657323816},
		/* 631, 640 */
		{0.000166e-6, 7238.675591600, 5.930629110},
		{0.000167e-6, 3097.883822726, 5.570955333},
		{0.000159e-6, -323.505416657, 5.786670700},
		{0.000154e-6, -4136.910433516, 1.517805532},
		{0.000176e-6, 12029.347187887, 3.139266834},
		{0.000167e-6, 12132.439962106, 3.556352289},
		{0.000153e-6, 202.253395174, 1.463313961},
		{0.000157e-6, 17267.268201691, 1.586837396},
		{0.000142e-6, 83996.847317911, 0.022670115},
		{0.000152e-6, 17260.154654690, 0.708528947},
		/* 641, 650 */
		{0.000144e-6, 6084.003848555, 5.187075177},
		{0.000135e-6, 5756.566278634, 1.993229262},
		{0.000134e-6, 5750.203491159, 3.457197134},
		{0.000144e-6, 5326.786694021, 6.066193291},
		{0.000160e-6, 11015.106477335, 1.710431974},
		{0.000133e-6, 3634.621024518, 2.836451652},
		{0.000134e-6, 18073.704938650, 5.453106665},
		{0.000134e-6, 1162.474704408, 5.326898811},
		{0.000128e-6, 5642.198242609, 2.511652591},
		{0.000160e-6, 632.783739313, 5.628785365},
		/* 651, 660 */
		{0.000132e-6, 13916.019109642, 0.819294053},
		{0.000122e-6, 14314.168113050, 5.677408071},
		{0.000125e-6, 12359.966151546, 5.251984735},
		{0.000121e-6, 5749.452731634, 2.210924603},
		{0.000136e-6, -245.831646229, 1.646502367},
		{0.000120e-6, 5757.317038160, 3.240883049},
		{0.000134e-6, 12146.667056108, 3.059480037},
		{0.000137e-6, 6206.809778716, 1.867105418},
		{0.000141e-6, 17253.041107690, 2.069217456},
		{0.000129e-6, -7477.522860216, 2.781469314},
		/* 661, 670 */
		{0.000116e-6, 5540.085789459, 4.281176991},
		{0.000116e-6, 9779.108676125, 3.320925381},
		{0.000129e-6, 5237.921013804, 3.497704076},
		{0.000113e-6, 5959.570433334, 0.983210840},
		{0.000122e-6, 6282.095528923, 2.674938860},
		{0.000140e-6, -11.045700264, 4.957936982},
		{0.000108e-6, 23543.230504682, 1.390113589},
		{0.000106e-6, -12569.674818332, 0.429631317},
		{0.000110e-6, -266.607041722, 5.501340197},
		{0.000115e-6, 12559.038152982, 4.691456618},
		/* 671, 680 */
		{0.000134e-6, -2388.894020449, 0.577313584},
		{0.000109e-6, 10440.274292604, 6.218148717},
		{0.000102e-6, -543.918059096, 1.477842615},
		{0.000108e-6, 21228.392023546, 2.237753948},
		{0.000101e-6, -4535.059436924, 3.100492232},
		{0.000103e-6, 76.266071276, 5.594294322},
		{0.000104e-6, 949.175608970, 5.674287810},
		{0.000101e-6, 13517.870106233, 2.196632348},
		{0.000100e-6, 11933.367960670, 4.056084160},

		/* T^2 */
		{4.322990e-6, 6283.075849991, 2.642893748},
		/* 681, 690 */
		{0.406495e-6, 0.000000000, 4.712388980},
		{0.122605e-6, 12566.151699983, 2.438140634},
		{0.019476e-6, 213.299095438, 1.642186981},
		{0.016916e-6, 529.690965095, 4.510959344},
		{0.013374e-6, -3.523118349, 1.502210314},
		{0.008042e-6, 26.298319800, 0.478549024},
		{0.007824e-6, 155.420399434, 5.254710405},
		{0.004894e-6, 5746.271337896, 4.683210850},
		{0.004875e-6, 5760.498431898, 0.759507698},
		{0.004416e-6, 5223.693919802, 6.028853166},
		/* 691, 700 */
		{0.004088e-6, -7.113547001, 0.060926389},
		{0.004433e-6, 77713.771467920, 3.627734103},
		{0.003277e-6, 18849.227549974, 2.327912542},
		{0.002703e-6, 6062.663207553, 1.271941729},
		{0.003435e-6, -775.522611324, 0.747446224},
		{0.002618e-6, 6076.890301554, 3.633715689},
		{0.003146e-6, 206.185548437, 5.647874613},
		{0.002544e-6, 1577.343542448, 6.232904270},
		{0.002218e-6, -220.412642439, 1.309509946},
		{0.002197e-6, 5856.477659115, 2.407212349},
		/* 701, 710 */
		{0.002897e-6, 5753.384884897, 5.863842246},
		{0.001766e-6, 426.598190876, 0.754113147},
		{0.001738e-6, -796.298006816, 2.714942671},
		{0.001695e-6, 522.577418094, 2.629369842},
		{0.001584e-6, 5507.553238667, 1.341138229},
		{0.001503e-6, -242.728603974, 0.377699736},
		{0.001552e-6, -536.804512095, 2.904684667},
		{0.001370e-6, -398.149003408, 1.265599125},
		{0.001889e-6, -5573.142801634, 4.413514859},
		{0.001722e-6, 6069.776754553, 2.445966339},
		/* 711, 720 */
		{0.001124e-6, 1059.381930189, 5.041799657},
		{0.001258e-6, 553.569402842, 3.849557278},
		{0.000831e-6, 951.718406251, 2.471094709},
		{0.000767e-6, 4694.002954708, 5.363125422},
		{0.000756e-6, 1349.867409659, 1.046195744},
		{0.000775e-6, -11.045700264, 0.245548001},
		{0.000597e-6, 2146.165416475, 4.543268798},
		{0.000568e-6, 5216.580372801, 4.178853144},
		{0.000711e-6, 1748.016413067, 5.934271972},
		{0.000499e-6, 12036.460734888, 0.624434410},
		/* 721, 730 */
		{0.000671e-6, -1194.447010225, 4.136047594},
		{0.000488e-6, 5849.364112115, 2.209679987},
		{0.000621e-6, 6438.496249426, 4.518860804},
		{0.000495e-6, -6286.598968340, 1.868201275},
		{0.000456e-6, 5230.807466803, 1.271231591},
		{0.000451e-6, 5088.628839767, 0.084060889},
		{0.000435e-6, 5643.178563677, 3.324456609},
		{0.000387e-6, 10977.078804699, 4.052488477},
		{0.000547e-6, 161000.685737473, 2.841633844},
		{0.000522e-6, 3154.687084896, 2.171979966},
		/* 731, 740 */
		{0.000375e-6, 5486.777843175, 4.983027306},
		{0.000421e-6, 5863.591206116, 4.546432249},
		{0.000439e-6, 7084.896781115, 0.522967921},
		{0.000309e-6, 2544.314419883, 3.172606705},
		{0.000347e-6, 4690.479836359, 1.479586566},
		{0.000317e-6, 801.820931124, 3.553088096},
		{0.000262e-6, 419.484643875, 0.606635550},
		{0.000248e-6, 6836.645252834, 3.014082064},
		{0.000245e-6, -1592.596013633, 5.519526220},
		{0.000225e-6, 4292.330832950, 2.877956536},
		/* 741, 750 */
		{0.000214e-6, 7234.794256242, 1.605227587},
		{0.000205e-6, 5767.611978898, 0.625804796},
		{0.000180e-6, 10447.387839604, 3.499954526},
		{0.000229e-6, 199.072001436, 5.632304604},
		{0.000214e-6, 639.897286314, 5.960227667},
		{0.000175e-6, -433.711737877, 2.162417992},
		{0.000209e-6, 515.463871093, 2.322150893},
		{0.000173e-6, 6040.347246017, 2.556183691},
		{0.000184e-6, 6309.374169791, 4.732296790},
		{0.000227e-6, 149854.400134205, 5.385812217},
		/* 751, 760 */
		{0.000154e-6, 8031.092263058, 5.120720920},
		{0.000151e-6, 5739.157790895, 4.815000443},
		{0.000197e-6, 7632.943259650, 0.222827271},
		{0.000197e-6, 74.781598567, 3.910456770},
		{0.000138e-6, 6055.549660552, 1.397484253},
		{0.000149e-6, -6127.655450557, 5.333727496},
		{0.000137e-6, 3894.181829542, 4.281749907},
		{0.000135e-6, 9437.762934887, 5.979971885},
		{0.000139e-6, -2352.866153772, 4.715630782},
		{0.000142e-6, 6812.766815086, 0.513330157},
		/* 761, 770 */
		{0.000120e-6, -4705.732307544, 0.194160689},
		{0.000131e-6, -71430.695617928, 0.000379226},
		{0.000124e-6, 6279.552731642, 2.122264908},
		{0.000108e-6, -6256.777530192, 0.883445696},

		/* T^3 */
		{0.143388e-6, 6283.075849991, 1.131453581},
		{0.006671e-6, 12566.151699983, 0.775148887},
		{0.001480e-6, 155.420399434, 0.480016880},
		{0.000934e-6, 213.299095438, 6.144453084},
		{0.000795e-6, 529.690965095, 2.941595619},
		{0.000673e-6, 5746.271337896, 0.120415406},
		/* 771, 780 */
		{0.000672e-6, 5760.498431898, 5.317009738},
		{0.000389e-6, -220.412642439, 3.090323467},
		{0.000373e-6, 6062.663207553, 3.003551964},
		{0.000360e-6, 6076.890301554, 1.918913041},
		{0.000316e-6, -21.340641002, 5.545798121},
		{0.000315e-6, -242.728603974, 1.884932563},
		{0.000278e-6, 206.185548437, 1.266254859},
		{0.000238e-6, -536.804512095, 4.532664830},
		{0.000185e-6, 522.577418094, 4.578313856},
		{0.000245e-6, 18849.227549974, 0.587467082},
		/* 781, 787 */
		{0.000180e-6, 426.598190876, 5.151178553},
		{0.000200e-6, 553.569402842, 5.355983739},
		{0.000141e-6, 5223.693919802, 1.336556009},
		{0.000104e-6, 5856.477659115, 4.239842759},

		/* T^4 */
		{0.003826e-6, 6283.075849991, 5.705257275},
		{0.000303e-6, 12566.151699983, 5.407132842},
		{0.000209e-6, 155.420399434, 1.989815753},
	}

	/* Time since J2000.0 in Julian millennia. */
	t = ((date1 - DJ00) + date2) / DJM

	/* ================= */
	/* Topocentric terms */
	/* ================= */

	/* Convert UT to local solar time in radians. */
	tsol = math.Mod(ut, 1.0)*D2PI + elong

	/* FUNDAMENTAL ARGUMENTS:  Simon et al. 1994. */

	/* Combine time argument (millennia) with deg/arcsec factor. */
	w = t / 3600.0

	/* Sun Mean Longitude. */
	elsun = math.Mod(280.46645683+1296027711.03429*w, 360.0) * DD2R

	/* Sun Mean Anomaly. */
	emsun = math.Mod(357.52910918+1295965810.481*w, 360.0) * DD2R

	/* Mean Elongation of Moon from Sun. */
	d = math.Mod(297.85019547+16029616012.090*w, 360.0) * DD2R

	/* Mean Longitude of Jupiter. */
	elj = math.Mod(34.35151874+109306899.89453*w, 360.0) * DD2R

	/* Mean Longitude of Saturn. */
	els = math.Mod(50.07744430+44046398.47038*w, 360.0) * DD2R

	/* TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983. */
	wt = +0.00029e-10*u*math.Sin(tsol+elsun-els) +
		+0.00100e-10*u*math.Sin(tsol-2.0*emsun) +
		+0.00133e-10*u*math.Sin(tsol-d) +
		+0.00133e-10*u*math.Sin(tsol+elsun-elj) +
		-0.00229e-10*u*math.Sin(tsol+2.0*elsun+emsun) +
		-0.02200e-10*v*math.Cos(elsun+emsun) +
		+0.05312e-10*u*math.Sin(tsol-emsun) +
		-0.13677e-10*u*math.Sin(tsol+2.0*elsun) +
		-1.31840e-10*v*math.Cos(elsun) +
		+3.17679e-10*u*math.Sin(tsol)

	/* ===================== */
	/* Fairhead et al. model */
	/* ===================== */

	/* T0 */
	w0 = 0
	for j = 473; j >= 0; j-- {
		w0 += fairhd[j][0] * math.Sin(fairhd[j][1]*t+fairhd[j][2])
	}

	/* T1 */
	w1 = 0
	for j = 678; j >= 474; j-- {
		w1 += fairhd[j][0] * math.Sin(fairhd[j][1]*t+fairhd[j][2])
	}

	/* T2 */
	w2 = 0
	for j = 763; j >= 679; j-- {
		w2 += fairhd[j][0] * math.Sin(fairhd[j][1]*t+fairhd[j][2])
	}

	/* T3 */
	w3 = 0
	for j = 783; j >= 764; j-- {
		w3 += fairhd[j][0] * math.Sin(fairhd[j][1]*t+fairhd[j][2])
	}

	/* T4 */
	w4 = 0
	for j = 786; j >= 784; j-- {
		w4 += fairhd[j][0] * math.Sin(fairhd[j][1]*t+fairhd[j][2])
	}

	/* Multiply by powers of T and combine. */
	wf = t*(t*(t*(t*w4+w3)+w2)+w1) + w0

	/* Adjustments to use JPL planetary masses instead of IAU. */
	wj = 0.00065e-6*math.Sin(6069.776754*t+4.021194) +
		0.00033e-6*math.Sin(213.299095*t+5.543132) +
		(-0.00196e-6 * math.Sin(6208.294251*t+5.696701)) +
		(-0.00173e-6 * math.Sin(74.781599*t+2.435900)) +
		0.03638e-6*t*t

	/* ============ */
	/* Final result */
	/* ============ */

	/* TDB-TT in seconds. */
	w = wt + wf + wj

	return w
}

/*
Taitt	TAI to TT.

Time scale transformation:  International Atomic Time, TAI, to
Terrestrial Time, TT.

Given:

	tai1,tai2  float64    TAI as a 2-part Julian Date

Returned:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Note:

	tai1+tai2 is Julian Date, apportioned in any convenient way
	between the two arguments, for example where tai1 is the Julian
	Day Number and tai2 is the fraction of a day.  The returned
	tt1,tt2 follow suit.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Taitt(tai1, tai2 float64, tt1, tt2 *float64) int {
	/* TT minus TAI (days). */
	const dtat = TTMTAI / DAYSEC

	/* Result, safeguarding precision. */
	if math.Abs(tai1) > math.Abs(tai2) {
		*tt1 = tai1
		*tt2 = tai2 + dtat
	} else {
		*tt1 = tai1 + dtat
		*tt2 = tai2
	}

	/* Status (always OK). */
	return 0
}

/*
Taiut1	TAI to UT1

Time scale transformation:  International Atomic Time, TAI, to
Universal Time, UT1.

Given:

	tai1,tai2  float64    TAI as a 2-part Julian Date
	dta        float64    UT1-TAI in seconds

Returned:

	ut11,ut12  float64    UT1 as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tai1+tai2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tai1 is the Julian
    Day Number and tai2 is the fraction of a day.  The returned
    UT11,UT12 follow suit.

 2. The argument dta, i.e. UT1-TAI, is an observed quantity, and is
    available from IERS tabulations.

Reference:

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Taiut1(tai1, tai2, dta float64, ut11, ut12 *float64) int {
	var dtad float64

	/* Result, safeguarding precision. */
	dtad = dta / DAYSEC
	if math.Abs(tai1) > math.Abs(tai2) {
		*ut11 = tai1
		*ut12 = tai2 + dtad
	} else {
		*ut11 = tai1 + dtad
		*ut12 = tai2
	}

	/* Status (always OK). */
	return 0
}

/*
Taiutc	TAI to UTC

Time scale transformation:  International Atomic Time, TAI, to
Coordinated Universal Time, UTC.

Given:

	tai1,tai2  float64   TAI as a 2-part Julian Date (Note 1)

Returned:

	utc1,utc2  float64   UTC as a 2-part quasi Julian Date (Notes 1-3)

Returned (function value):

	int      status: +1 = dubious year (Note 4)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. tai1+tai2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tai1 is the Julian
    Day Number and tai2 is the fraction of a day.  The returned utc1
    and utc2 form an analogous pair, except that a special convention
    is used, to deal with the problem of leap seconds - see the next
    note.

 2. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the JD day represents UTC days whether the
    length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
    there were smaller jumps (in either direction) each time the
    linear UTC(TAI) expression was changed, and these "mini-leaps"
    are also included in the SOFA convention.

 3. The function D2dtf can be used to transform the UTC quasi-JD
    into calendar date and clock time, including UTC leap second
    handling.

 4. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

Called:

	Utctai    UTC to TAI

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Taiutc(tai1, tai2 float64, utc1, utc2 *float64) int {
	var big1 bool
	var i, j int
	var a1, a2, u1, u2, g1, g2 float64

	/* Put the two parts of the TAI into big-first order. */
	big1 = (math.Abs(tai1) >= math.Abs(tai2))
	if big1 {
		a1 = tai1
		a2 = tai2
	} else {
		a1 = tai2
		a2 = tai1
	}

	/* Initial guess for UTC. */
	u1 = a1
	u2 = a2

	/* Iterate (though in most cases just once is enough). */
	for i = 0; i < 3; i++ {

		/* Guessed UTC to TAI. */
		j = Utctai(u1, u2, &g1, &g2)
		if j < 0 {
			return j
		}

		/* Adjust guessed UTC. */
		u2 += a1 - g1
		u2 += a2 - g2
	}

	/* Return the UTC result, preserving the TAI order. */
	if big1 {
		*utc1 = u1
		*utc2 = u2
	} else {
		*utc1 = u2
		*utc2 = u1
	}

	/* Status. */
	return j
}

/*
Tcbtdb	TCB to TDB

Time scale transformation:  Barycentric Coordinate Time, TCB, to
Barycentric Dynamical Time, TDB.

Given:

	tcb1,tcb2  float64    TCB as a 2-part Julian Date

Returned:

	tdb1,tdb2  float64    TDB as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tcb1+tcb2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tcb1 is the Julian
    Day Number and tcb2 is the fraction of a day.  The returned
    tdb1,tdb2 follow suit.

 2. The 2006 IAU General Assembly introduced a conventional linear
    transformation between TDB and TCB.  This transformation
    compensates for the drift between TCB and terrestrial time TT,
    and keeps TDB approximately centered on TT.  Because the
    relationship between TT and TCB depends on the adopted solar
    system ephemeris, the degree of alignment between TDB and TT over
    long intervals will vary according to which ephemeris is used.
    Former definitions of TDB attempted to avoid this problem by
    stipulating that TDB and TT should differ only by periodic
    effects.  This is a good description of the nature of the
    relationship but eluded precise mathematical formulation.  The
    conventional linear relationship adopted in 2006 sidestepped
    these difficulties whilst delivering a TDB that in practice was
    consistent with values before that date.

 3. TDB is essentially the same as Teph, the time argument for the
    JPL solar system ephemerides.

Reference:

	IAU 2006 Resolution B3
*/
func Tcbtdb(tcb1, tcb2 float64, tdb1, tdb2 *float64) int {
	/* 1977 Jan 1 00:00:32.184 TT, as two-part JD */
	const t77td = DJM0 + DJM77
	const t77tf = TTMTAI / DAYSEC

	/* TDB (days) at TAI 1977 Jan 1.0 */
	const tdb0 = TDB0 / DAYSEC

	var d float64

	/* Result, safeguarding precision. */
	if math.Abs(tcb1) > math.Abs(tcb2) {
		d = tcb1 - t77td
		*tdb1 = tcb1
		*tdb2 = tcb2 + tdb0 - (d+(tcb2-t77tf))*ELB
	} else {
		d = tcb2 - t77td
		*tdb1 = tcb1 + tdb0 - (d+(tcb1-t77tf))*ELB
		*tdb2 = tcb2
	}

	/* Status (always OK). */
	return 0
}

/*
Tcgtt	TCG to TT

Time scale transformation:  Geocentric Coordinate Time, TCG, to
Terrestrial Time, TT.

Given:

	tcg1,tcg2  float64    TCG as a 2-part Julian Date

Returned:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Note:

	tcg1+tcg2 is Julian Date, apportioned in any convenient way
	between the two arguments, for example where tcg1 is the Julian
	Day Number and tcg22 is the fraction of a day.  The returned
	tt1,tt2 follow suit.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	IAU 2000 Resolution B1.9
*/
func Tcgtt(tcg1, tcg2 float64, tt1, tt2 *float64) int {
	/* 1977 Jan 1 00:00:32.184 TT, as MJD */
	const t77t = DJM77 + TTMTAI/DAYSEC

	/* Result, safeguarding precision. */
	if math.Abs(tcg1) > math.Abs(tcg2) {
		*tt1 = tcg1
		*tt2 = tcg2 - ((tcg1-DJM0)+(tcg2-t77t))*ELG
	} else {
		*tt1 = tcg1 - ((tcg2-DJM0)+(tcg1-t77t))*ELG
		*tt2 = tcg2
	}

	/* OK status. */
	return 0
}

/*
Tdbtcb	TDB to TCB

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Barycentric Coordinate Time, TCB.

Given:

	tdb1,tdb2  float64    TDB as a 2-part Julian Date

Returned:

	tcb1,tcb2  float64    TCB as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tdb1 is the Julian
    Day Number and tdb2 is the fraction of a day.  The returned
    tcb1,tcb2 follow suit.

 2. The 2006 IAU General Assembly introduced a conventional linear
    transformation between TDB and TCB.  This transformation
    compensates for the drift between TCB and terrestrial time TT,
    and keeps TDB approximately centered on TT.  Because the
    relationship between TT and TCB depends on the adopted solar
    system ephemeris, the degree of alignment between TDB and TT over
    long intervals will vary according to which ephemeris is used.
    Former definitions of TDB attempted to avoid this problem by
    stipulating that TDB and TT should differ only by periodic
    effects.  This is a good description of the nature of the
    relationship but eluded precise mathematical formulation.  The
    conventional linear relationship adopted in 2006 sidestepped
    these difficulties whilst delivering a TDB that in practice was
    consistent with values before that date.

 3. TDB is essentially the same as Teph, the time argument for the
    JPL solar system ephemerides.

Reference:

	IAU 2006 Resolution B3
*/
func Tdbtcb(tdb1, tdb2 float64, tcb1, tcb2 *float64) int {
	/* 1977 Jan 1 00:00:32.184 TT, as two-part JD */
	const t77td = DJM0 + DJM77
	const t77tf = TTMTAI / DAYSEC

	/* TDB (days) at TAI 1977 Jan 1.0 */
	const tdb0 = TDB0 / DAYSEC

	/* TDB to TCB rate */
	const elbb = ELB / (1.0 - ELB)

	var d, f float64

	/* Result, preserving date format but safeguarding precision. */
	if math.Abs(tdb1) > math.Abs(tdb2) {
		d = t77td - tdb1
		f = tdb2 - tdb0
		*tcb1 = tdb1
		*tcb2 = f - (d-(f-t77tf))*elbb
	} else {
		d = t77td - tdb2
		f = tdb1 - tdb0
		*tcb1 = f - (d-(f-t77tf))*elbb
		*tcb2 = tdb2
	}

	/* Status (always OK). */
	return 0
}

/*
Tdbtt	TDB to TT

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Terrestrial Time, TT.

Given:

	tdb1,tdb2  float64    TDB as a 2-part Julian Date
	dtr        float64    TDB-TT in seconds

Returned:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tdb1 is the Julian
    Day Number and tdb2 is the fraction of a day.  The returned
    tt1,tt2 follow suit.

 2. The argument dtr represents the quasi-periodic component of the
    GR transformation between TT and TCB.  It is dependent upon the
    adopted solar-system ephemeris, and can be obtained by numerical
    integration, by interrogating a precomputed time ephemeris or by
    evaluating a model such as that implemented in the SOFA function
    Dtdb.   The quantity is dominated by an annual term of 1.7 ms
    amplitude.

 3. TDB is essentially the same as Teph, the time argument for the
    JPL solar system ephemerides.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	IAU 2006 Resolution 3
*/
func Tdbtt(tdb1, tdb2, dtr float64, tt1, tt2 *float64) int {
	var dtrd float64

	/* Result, safeguarding precision. */
	dtrd = dtr / DAYSEC
	if math.Abs(tdb1) > math.Abs(tdb2) {
		*tt1 = tdb1
		*tt2 = tdb2 - dtrd
	} else {
		*tt1 = tdb1 - dtrd
		*tt2 = tdb2
	}

	/* Status (always OK). */
	return 0
}

/*
Tttai	TT to TAI

Time scale transformation:  Terrestrial Time, TT, to International
Atomic Time, TAI.

Given:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned:

	tai1,tai2  float64    TAI as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Note:

	tt1+tt2 is Julian Date, apportioned in any convenient way between
	the two arguments, for example where tt1 is the Julian Day Number
	and tt2 is the fraction of a day.  The returned tai1,tai2 follow
	suit.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Tttai(tt1, tt2 float64, tai1, tai2 *float64) int {
	/* TT minus TAI (days). */
	const dtat = TTMTAI / DAYSEC

	/* Result, safeguarding precision. */
	if math.Abs(tt1) > math.Abs(tt2) {
		*tai1 = tt1
		*tai2 = tt2 - dtat
	} else {
		*tai1 = tt1 - dtat
		*tai2 = tt2
	}

	/* Status (always OK). */
	return 0

}

/*
Tttcg	TT to TCG

Time scale transformation:  Terrestrial Time, TT, to Geocentric
Coordinate Time, TCG.

Given:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned:

	tcg1,tcg2  float64    TCG as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Note:

	tt1+tt2 is Julian Date, apportioned in any convenient way between
	the two arguments, for example where tt1 is the Julian Day Number
	and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
	suit.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	IAU 2000 Resolution B1.9
*/
func Tttcg(tt1, tt2 float64, tcg1, tcg2 *float64) int {
	/* 1977 Jan 1 00:00:32.184 TT, as MJD */
	const t77t = DJM77 + TTMTAI/DAYSEC

	/* TT to TCG rate */
	const elgg = ELG / (1.0 - ELG)

	/* Result, safeguarding precision. */
	if math.Abs(tt1) > math.Abs(tt2) {
		*tcg1 = tt1
		*tcg2 = tt2 + ((tt1-DJM0)+(tt2-t77t))*elgg
	} else {
		*tcg1 = tt1 + ((tt2-DJM0)+(tt1-t77t))*elgg
		*tcg2 = tt2
	}

	/* Status (always OK). */
	return 0
}

/*
Tttdb	TT to TDB

Time scale transformation:  Terrestrial Time, TT, to Barycentric
Dynamical Time, TDB.

Given:

	tt1,tt2    float64    TT as a 2-part Julian Date
	dtr        float64    TDB-TT in seconds

Returned:

	tdb1,tdb2  float64    TDB as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tt1+tt2 is Julian Date, apportioned in any convenient way between
    the two arguments, for example where tt1 is the Julian Day Number
    and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
    suit.

 2. The argument dtr represents the quasi-periodic component of the
    GR transformation between TT and TCB.  It is dependent upon the
    adopted solar-system ephemeris, and can be obtained by numerical
    integration, by interrogating a precomputed time ephemeris or by
    evaluating a model such as that implemented in the SOFA function
    Dtdb.   The quantity is dominated by an annual term of 1.7 ms
    amplitude.

 3. TDB is essentially the same as Teph, the time argument for the JPL
    solar system ephemerides.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	IAU 2006 Resolution 3
*/
func Tttdb(tt1, tt2, dtr float64, tdb1, tdb2 *float64) int {
	var dtrd float64

	/* Result, safeguarding precision. */
	dtrd = dtr / DAYSEC
	if math.Abs(tt1) > math.Abs(tt2) {
		*tdb1 = tt1
		*tdb2 = tt2 + dtrd
	} else {
		*tdb1 = tt1 + dtrd
		*tdb2 = tt2
	}

	/* Status (always OK). */
	return 0
}

/*
Ttut1	TT to UT1

Time scale transformation:  Terrestrial Time, TT, to Universal Time,
UT1.

Given:

	tt1,tt2    float64    TT as a 2-part Julian Date
	dt         float64    TT-UT1 in seconds

Returned:

	ut11,ut12  float64    UT1 as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. tt1+tt2 is Julian Date, apportioned in any convenient way between
    the two arguments, for example where tt1 is the Julian Day Number
    and tt2 is the fraction of a day.  The returned ut11,ut12 follow
    suit.

 2. The argument dt is classical Delta T.

Reference:

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Ttut1(tt1, tt2, dt float64, ut11, ut12 *float64) int {
	var dtd float64

	/* Result, safeguarding precision. */
	dtd = dt / DAYSEC
	if math.Abs(tt1) > math.Abs(tt2) {
		*ut11 = tt1
		*ut12 = tt2 - dtd
	} else {
		*ut11 = tt1 - dtd
		*ut12 = tt2
	}

	/* Status (always OK). */
	return 0
}

/*
Ut1tai	UT1 to TAI

Time scale transformation:  Universal Time, UT1, to International
Atomic Time, TAI.

Given:

	ut11,ut12  float64    UT1 as a 2-part Julian Date
	dta        float64    UT1-TAI in seconds

Returned:

	tai1,tai2  float64    TAI as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. ut11+ut12 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where ut11 is the Julian
    Day Number and ut12 is the fraction of a day.  The returned
    tai1,tai2 follow suit.

 2. The argument dta, i.e. UT1-TAI, is an observed quantity, and is
    available from IERS tabulations.

Reference:

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Ut1tai(ut11, ut12, dta float64, tai1, tai2 *float64) int {
	var dtad float64

	/* Result, safeguarding precision. */
	dtad = dta / DAYSEC
	if math.Abs(ut11) > math.Abs(ut12) {
		*tai1 = ut11
		*tai2 = ut12 - dtad
	} else {
		*tai1 = ut11 - dtad
		*tai2 = ut12
	}

	/* Status (always OK). */
	return 0
}

/*
Ut1tt	UT1 to TT

Time scale transformation:  Universal Time, UT1, to Terrestrial
Time, TT.

Given:

	ut11,ut12  float64    UT1 as a 2-part Julian Date
	dt         float64    TT-UT1 in seconds

Returned:

	tt1,tt2    float64    TT as a 2-part Julian Date

Returned (function value):

	int       status:  0 = OK

Notes:

 1. ut11+ut12 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where ut11 is the Julian
    Day Number and ut12 is the fraction of a day.  The returned
    tt1,tt2 follow suit.

 2. The argument dt is classical Delta T.

Reference:

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Ut1tt(ut11, ut12, dt float64, tt1, tt2 *float64) int {
	var dtd float64

	/* Result, safeguarding precision. */
	dtd = dt / DAYSEC
	if math.Abs(ut11) > math.Abs(ut12) {
		*tt1 = ut11
		*tt2 = ut12 + dtd
	} else {
		*tt1 = ut11 + dtd
		*tt2 = ut12
	}

	/* Status (always OK). */
	return 0
}

/*
Ut1utc	UT1 to UTC

Time scale transformation:  Universal Time, UT1, to Coordinated
Universal Time, UTC.

Given:

	ut11,ut12  float64   UT1 as a 2-part Julian Date (Note 1)
	dut1       float64   Delta UT1: UT1-UTC in seconds (Note 2)

Returned:

	utc1,utc2  float64   UTC as a 2-part quasi Julian Date (Notes 3,4)

Returned (function value):

	int      status: +1 = dubious year (Note 5)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. ut11+ut12 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where ut11 is the Julian
    Day Number and ut12 is the fraction of a day.  The returned utc1
    and utc2 form an analogous pair, except that a special convention
    is used, to deal with the problem of leap seconds - see Note 3.

 2. Delta UT1 can be obtained from tabulations provided by the
    International Earth Rotation and Reference Systems Service.  The
    value changes abruptly by 1s at a leap second;  however, close to
    a leap second the algorithm used here is tolerant of the "wrong"
    choice of value being made.

 3. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the returned quasi-JD UTC1+UTC2 represents UTC
    days whether the length is 86399, 86400 or 86401 SI seconds.

 4. The function D2dtf can be used to transform the UTC quasi-JD
    into calendar date and clock time, including UTC leap second
    handling.

 5. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

Called:

	Jd2cal    JD to Gregorian calendar
	Dat       delta(AT) = TAI-UTC
	Cal2jd    Gregorian calendar to JD

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Ut1utc(ut11, ut12, dut1 float64, utc1, utc2 *float64) int {
	var big1 bool
	var i, iy, im, id, js int
	var duts, u1, u2, d1, dats1, d2, fd, dats2, ddats, us1, us2, du float64

	/* UT1-UTC in seconds. */
	duts = dut1

	/* Put the two parts of the UT1 into big-first order. */
	big1 = (math.Abs(ut11) >= math.Abs(ut12))
	if big1 {
		u1 = ut11
		u2 = ut12
	} else {
		u1 = ut12
		u2 = ut11
	}

	/* See if the UT1 can possibly be in a leap-second day. */
	d1 = u1
	dats1 = 0
	for i = -1; i <= 3; i++ {
		d2 = u2 + float64(i)
		if Jd2cal(d1, d2, &iy, &im, &id, &fd) != 0 {
			return -1
		}
		js = Dat(iy, im, id, 0.0, &dats2)
		if js < 0 {
			return -1
		}
		if i == -1 {
			dats1 = dats2
		}
		ddats = dats2 - dats1
		if math.Abs(ddats) >= 0.5 {

			/* Yes, leap second nearby: ensure UT1-UTC is "before" value. */
			if ddats*duts >= 0 {
				duts -= ddats
			}

			/* UT1 for the start of the UTC day that ends in a leap. */
			if Cal2jd(iy, im, id, &d1, &d2) != 0 {
				return -1
			}
			us1 = d1
			us2 = d2 - 1.0 + duts/DAYSEC

			/* Is the UT1 after this point? */
			du = u1 - us1
			du += u2 - us2
			if du > 0 {

				/* Yes:  fraction of the current UTC day that has elapsed. */
				fd = du * DAYSEC / (DAYSEC + ddats)

				/* Ramp UT1-UTC to bring about SOFA's JD(UTC) convention. */
				//duts += ddats * ( fd <= 1.0 ? fd : 1.0 );
				if fd <= 1.0 {
					duts += ddats * fd
				} else {
					duts += ddats
				}
			}

			/* Done. */
			break
		}
		dats1 = dats2
	}

	/* Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC. */
	u2 -= duts / DAYSEC

	/* Result, safeguarding precision. */
	if big1 {
		*utc1 = u1
		*utc2 = u2
	} else {
		*utc1 = u2
		*utc2 = u1
	}

	/* Status. */
	return js
}

/*
Utctai	UTC to TAI

Time scale transformation:  Coordinated Universal Time, UTC, to
International Atomic Time, TAI.

Given:

	utc1,utc2  float64   UTC as a 2-part quasi Julian Date (Notes 1-4)

Returned:

	tai1,tai2  float64   TAI as a 2-part Julian Date (Note 5)

Returned (function value):

	int      status: +1 = dubious year (Note 3)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

 2. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the JD day represents UTC days whether the
    length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
    there were smaller jumps (in either direction) each time the
    linear UTC(TAI) expression was changed, and these "mini-leaps"
    are also included in the SOFA convention.

 3. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

 4. The function Dtf2d converts from calendar date and time of day
    into 2-part Julian Date, and in the case of UTC implements the
    leap-second-ambiguity convention described above.

 5. The returned TAI1,TAI2 are such that their sum is the TAI Julian
    Date.

Called:

	Jd2cal    JD to Gregorian calendar
	Dat       delta(AT) = TAI-UTC
	Cal2jd    Gregorian calendar to JD

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Explanatory Supplement to the Astronomical Almanac,
	P. Kenneth Seidelmann (ed), University Science Books (1992)
*/
func Utctai(utc1, utc2 float64, tai1, tai2 *float64) int {
	var big1 bool
	var iy, im, id, j, iyt, imt, idt int
	var u1, u2, fd, dat0, dat12, w, dat24, dlod, dleap, z1, z2, a2 float64

	/* Put the two parts of the UTC into big-first order. */
	big1 = (math.Abs(utc1) >= math.Abs(utc2))
	if big1 {
		u1 = utc1
		u2 = utc2
	} else {
		u1 = utc2
		u2 = utc1
	}

	/* Get TAI-UTC at 0h today. */
	j = Jd2cal(u1, u2, &iy, &im, &id, &fd)
	if j != 0 {
		return j
	}
	j = Dat(iy, im, id, 0.0, &dat0)
	if j < 0 {
		return j
	}

	/* Get TAI-UTC at 12h today (to detect drift). */
	j = Dat(iy, im, id, 0.5, &dat12)
	if j < 0 {
		return j
	}

	/* Get TAI-UTC at 0h tomorrow (to detect jumps). */
	j = Jd2cal(u1+1.5, u2-fd, &iyt, &imt, &idt, &w)
	if j != 0 {
		return j
	}
	j = Dat(iyt, imt, idt, 0.0, &dat24)
	if j < 0 {
		return j
	}

	/* Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP). */
	dlod = 2.0 * (dat12 - dat0)
	dleap = dat24 - (dat0 + dlod)

	/* Remove any scaling applied to spread leap into preceding day. */
	fd *= (DAYSEC + dleap) / DAYSEC

	/* Scale from (pre-1972) UTC seconds to SI seconds. */
	fd *= (DAYSEC + dlod) / DAYSEC

	/* Today's calendar date to 2-part JD. */
	if Cal2jd(iy, im, id, &z1, &z2) != 0 {
		return -1
	}

	/* Assemble the TAI result, preserving the UTC split and order. */
	a2 = z1 - u1
	a2 += z2
	a2 += fd + dat0/DAYSEC
	if big1 {
		*tai1 = u1
		*tai2 = a2
	} else {
		*tai1 = a2
		*tai2 = u1
	}

	/* Status. */
	return j
}

/*
Utcut1	UTC to UT1

Time scale transformation:  Coordinated Universal Time, UTC, to
Universal Time, UT1.

Given:

	utc1,utc2  float64   UTC as a 2-part quasi Julian Date (Notes 1-4)
	dut1       float64   Delta UT1 = UT1-UTC in seconds (Note 5)

Returned:

	ut11,ut12  float64   UT1 as a 2-part Julian Date (Note 6)

Returned (function value):

	int      status: +1 = dubious year (Note 3)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

 2. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the JD day represents UTC days whether the
    length is 86399, 86400 or 86401 SI seconds.

 3. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

 4. The function Dtf2d converts from calendar date and time of
    day into 2-part Julian Date, and in the case of UTC implements
    the leap-second-ambiguity convention described above.

 5. Delta UT1 can be obtained from tabulations provided by the
    International Earth Rotation and Reference Systems Service.
    It is the caller's responsibility to supply a dut1 argument
    containing the UT1-UTC value that matches the given UTC.

 6. The returned ut11,ut12 are such that their sum is the UT1 Julian
    Date.

References:

	   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	   IERS Technical Note No. 32, BKG (2004)

	   Explanatory Supplement to the Astronomical Almanac,
	   P. Kenneth Seidelmann (ed), University Science Books (1992)

	Called:
	   Jd2cal    JD to Gregorian calendar
	   Dat       delta(AT) = TAI-UTC
	   Utctai    UTC to TAI
	   Taiut1    TAI to UT1
*/
func Utcut1(utc1, utc2, dut1 float64, ut11, ut12 *float64) int {
	var iy, im, id, js, jw int
	var w, dat, dta, tai1, tai2 float64

	/* Look up TAI-UTC. */
	if Jd2cal(utc1, utc2, &iy, &im, &id, &w) != 0 {
		return -1
	}
	js = Dat(iy, im, id, 0.0, &dat)
	if js < 0 {
		return -1
	}

	/* Form UT1-TAI. */
	dta = dut1 - dat

	/* UTC to TAI to UT1. */
	jw = Utctai(utc1, utc2, &tai1, &tai2)
	if jw < 0 {
		return -1
	} else if jw > 0 {
		js = jw
	}
	if Taiut1(tai1, tai2, dta, ut11, ut12) != 0 {
		return -1
	}

	/* Status. */
	return js

}
