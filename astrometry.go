// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// SOFA Astrometry Tools

/*
Ab Apply stellar aberration

Apply aberration to transform natural direction into proper
direction.

Given:

	pnat    [3]float64   natural direction to the source (unit vector)
	v       [3]float64   observer barycentric velocity in units of c
	s       float64      distance between the Sun and the observer (au)
	bm1     float64      sqrt(1-|v|^2): reciprocal of Lorenz factor

Returned:

	ppr     [3]float64   proper direction to source (unit vector)

Notes:

 1. The algorithm is based on Expr. (7.40) in the Explanatory
    Supplement (Urban & Seidelmann 2013), but with the following
    changes:

    o  Rigorous rather than approximate normalization is applied.

    o  The gravitational potential term from Expr. (7) in
    Klioner (2003) is added, taking into account only the Sun's
    contribution.  This has a maximum effect of about
    0.4 microarcsecond.

 2. In almost all cases, the maximum accuracy will be limited by the
    supplied velocity.  For example, if the SOFA Epv00 function is
    used, errors of up to 5 microarcseconds could occur.

References:

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013).

	Klioner, Sergei A., "A practical relativistic model for micro-
	arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

Called:

	Pdp       scalar product of two p-vectors
*/
func Ab(pnat, v [3]float64, s, bm1 float64, ppr *[3]float64) {
	var i int
	var pdv, w1, w2, r2, w, r float64
	var p [3]float64

	pdv = Pdp(pnat, v)
	w1 = 1.0 + pdv/(1.0+bm1)
	w2 = SRS / s
	r2 = 0.0
	for i = 0; i < 3; i++ {
		w = pnat[i]*bm1 + w1*v[i] + w2*(v[i]-pdv*pnat[i])
		p[i] = w
		r2 = r2 + w*w
	}
	r = sqrt(r2)
	for i = 0; i < 3; i++ {
		ppr[i] = p[i] / r
	}
}

/*
Apcg Prepare for ICRS <-> GCRS, geocentric, special

For a geocentric observer, prepare star-independent astrometry
parameters for transformations between ICRS and GCRS coordinates.
The Earth ephemeris is supplied by the caller.

The parameters produced by this function are required in the
parallax, light deflection and aberration parts of the astrometric
transformation chain.

Given:

	date1  float64       TDB as a 2-part...
	date2  float64       ...Julian Date (Note 1)
	ebpv   [2][3]float64 Earth barycentric pos/vel (au, au/day)
	ehp    [3]float64    Earth heliocentric position (au)

Returned:

	astrom  ASTROM      star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       unchanged
	xpl    float64       unchanged
	ypl    float64       unchanged
	sphi   float64       unchanged
	cphi   float64       unchanged
	diurab float64       unchanged
	eral   float64       unchanged
	refa   float64       unchanged
	refb   float64       unchanged

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 4. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Apcs      astrometry parameters, ICRS-GCRS, space observer
*/
func Apcg(date1, date2 float64, ebpv [2][3]float64, ehp [3]float64, astrom *ASTROM) {
	/* Geocentric observer */
	pv := [2][3]float64{
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
	}

	/* Compute the star-independent astrometry parameters. */
	Apcs(date1, date2, pv, ebpv, ehp, astrom)
}

/*
Apcg13 Prepare for ICRS <-> GCRS, geocentric

For a geocentric observer, prepare star-independent astrometry
parameters for transformations between ICRS and GCRS coordinates.
The caller supplies the date, and SOFA models are used to predict
the Earth ephemeris.

The parameters produced by this function are required in the
parallax, light deflection and aberration parts of the astrometric
transformation chain.

Given:

	date1  float64     TDB as a 2-part...
	date2  float64     ...Julian Date (Note 1)

Returned:

	astrom  ASTROM     star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       unchanged
	xpl    float64       unchanged
	ypl    float64       unchanged
	sphi   float64       unchanged
	cphi   float64       unchanged
	diurab float64       unchanged
	eral   float64       unchanged
	refa   float64       unchanged
	refb   float64       unchanged

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. In cases where the caller wishes to supply his own Earth
    ephemeris, the function Apcg can be used instead of the present
    function.

 4. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 5. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Epv00     Earth position and velocity
	Apcg      astrometry parameters, ICRS-GCRS, geocenter
*/
func Apcg13(date1, date2 float64, astrom *ASTROM) {
	var ehpv, ebpv [2][3]float64

	/* Earth barycentric & heliocentric position/velocity (au, au/d). */
	Epv00(date1, date2, &ehpv, &ebpv)

	/* Compute the star-independent astrometry parameters. */
	Apcg(date1, date2, ebpv, ehpv[0], astrom)
}

/*
Apci Prepare for ICRS <-> CIRS, terrestrial, special

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
caller.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the astrometric transformation chain.

Given:

	date1  float64       TDB as a 2-part...
	date2  float64       ...Julian Date (Note 1)
	ebpv   [2][3]float64 Earth barycentric position/velocity (au, au/day)
	ehp    [3]float64    Earth heliocentric position (au)
	x,y    float64       CIP X,Y (components of unit vector)
	s      float64       the CIO locator s (radians)

Returned:

	astrom ASTROM     star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       unchanged
	xpl    float64       unchanged
	ypl    float64       unchanged
	sphi   float64       unchanged
	cphi   float64       unchanged
	diurab float64       unchanged
	eral   float64       unchanged
	refa   float64       unchanged
	refb   float64       unchanged

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. In cases where the caller does not wish to provide the Earth
    ephemeris and CIP/CIO, the function Apci13 can be used instead
    of the present function.  This computes the required quantities
    using other SOFA functions.

 4. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 5. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Apcg      astrometry parameters, ICRS-GCRS, geocenter
	C2ixys    celestial-to-intermediate matrix, given X,Y and s
*/
func Apci(date1, date2 float64, ebpv [2][3]float64, ehp [3]float64, x, y, s float64, astrom *ASTROM) {
	/* Star-independent astrometry parameters for geocenter. */
	Apcg(date1, date2, ebpv, ehp, astrom)

	/* CIO based BPN matrix. */
	C2ixys(x, y, s, &astrom.Bpn)
}

/*
Apci13 Prepare for ICRS <-> CIRS, terrestrial

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The caller supplies the date, and SOFA models are used
to predict the Earth ephemeris and CIP/CIO.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the astrometric transformation chain.

Given:

	date1  float64      TDB as a 2-part...
	date2  float64      ...Julian Date (Note 1)

Returned:

	astrom  ASTROM    star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       unchanged
	 xpl    float64       unchanged
	 ypl    float64       unchanged
	 sphi   float64       unchanged
	 cphi   float64       unchanged
	 diurab float64       unchanged
	 eral   float64       unchanged
	 refa   float64       unchanged
	 refb   float64       unchanged
	eo      float64       equation of the origins (ERA-GST)

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. In cases where the caller wishes to supply his own Earth
    ephemeris and CIP/CIO, the function Apci can be used instead
    of the present function.

 4. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 5. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Epv00     Earth position and velocity
	Pnm06a    classical NPB matrix, IAU 2006/2000A
	Bpn2xy    extract CIP X,Y coordinates from NPB matrix
	S06       the CIO locator s, given X,Y, IAU 2006
	Apci      astrometry parameters, ICRS-CIRS
	Eors      equation of the origins, given NPB matrix and s
*/
func Apci13(date1, date2 float64, astrom *ASTROM, eo *float64) {
	var ehpv, ebpv [2][3]float64
	var r [3][3]float64
	var x, y, s float64

	/* Earth barycentric & heliocentric position/velocity (au, au/d). */
	Epv00(date1, date2, &ehpv, &ebpv)

	/* Form the equinox based BPN matrix, IAU 2006/2000A. */
	Pnm06a(date1, date2, &r)

	/* Extract CIP X,Y. */
	Bpn2xy(r, &x, &y)

	/* Obtain CIO locator s. */
	s = S06(date1, date2, x, y)

	/* Compute the star-independent astrometry parameters. */
	Apci(date1, date2, ebpv, ehpv[0], x, y, s, astrom)

	/* Equation of the origins. */
	*eo = Eors(r, s)
}

/*
Apco Prepare for ICRS <-> observed, terrestrial, special

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed
coordinates.  The caller supplies the Earth ephemeris, the Earth
rotation information and the refraction constants as well as the
site coordinates.

Given:

	date1  float64       TDB as a 2-part...
	date2  float64       ...Julian Date (Note 1)
	ebpv   [2][3]float64 Earth barycentric PV (au, au/day, Note 2)
	ehp    [3]float64    Earth heliocentric P (au, Note 2)
	x,y    float64       CIP X,Y (components of unit vector)
	s      float64       the CIO locator s (radians)
	theta  float64       Earth rotation angle (radians)
	elong  float64       longitude (radians, east +ve, Note 3)
	phi    float64       latitude (geodetic, radians, Note 3)
	hm     float64       height above ellipsoid (m, geodetic, Note 3)
	xp,yp  float64       polar motion coordinates (radians, Note 4)
	sp     float64       the TIO locator s' (radians, Note 4)
	refa   float64       refraction constant A (radians, Note 5)
	refb   float64       refraction constant B (radians, Note 5)

Returned:

	astrom ASTROM     star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       adjusted longitude (radians)
	xpl    float64       polar motion xp wrt local meridian (radians)
	ypl    float64       polar motion yp wrt local meridian (radians)
	sphi   float64       sine of geodetic latitude
	cphi   float64       cosine of geodetic latitude
	diurab float64       magnitude of diurnal aberration vector
	eral   float64       "local" Earth rotation angle (radians)
	refa   float64       refraction constant A (radians)
	refb   float64       refraction constant B (radians)

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. The vectors eb, eh, and all the astrom vectors, are with respect
    to BCRS axes.

 3. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
    CONVENTION:  the longitude required by the present function is
    right-handed, i.e. east-positive, in accordance with geographical
    convention.

    The adjusted longitude stored in the astrom array takes into
    account the TIO locator and polar motion.

 4. xp and yp are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions), measured along the
    meridians 0 and 90 deg west respectively.  sp is the TIO locator
    s', in radians, which positions the Terrestrial Intermediate
    Origin on the equator.  For many applications, xp, yp and
    (especially) sp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto the
    local meridian.

 5. The refraction constants refa and refb are for use in a
    dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
    (i.e. refracted) zenith distance and dZ is the amount of
    refraction.

 6. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

 7. In cases where the caller does not wish to provide the Earth
    Ephemeris, the Earth rotation information and refraction
    constants, the function Apco13 can be used instead of the
    present function.  This starts from UTC and weather readings etc.
    and computes suitable values using other SOFA functions.

 8. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 9. The context structure astrom produced by this function is used by
    Atioq, Atoiq, Atciq* and Aticq*.

Called:

	Ir        initialize r-matrix to identity
	Rz        rotate around Z-axis
	Ry        rotate around Y-axis
	Rx        rotate around X-axis
	Anpm      normalize angle into range +/- pi
	C2ixys    celestial-to-intermediate matrix, given X,Y and s
	Pvtob     position/velocity of terrestrial station
	Trxpv     product of transpose of r-matrix and pv-vector
	Apcs      astrometry parameters, ICRS-GCRS, space observer
	Cr        copy r-matrix
*/
func Apco(date1, date2 float64, ebpv [2][3]float64, ehp [3]float64, x, y, s float64, theta,
	elong, phi float64, hm, xp, yp, sp float64, refa, refb float64, astrom *ASTROM) {
	var r [3][3]float64
	var a, b, eral, c float64
	var pvc, pv [2][3]float64

	/* Form the rotation matrix, CIRS to apparent [HA,Dec]. */
	Ir(&r)
	Rz(theta+sp, &r)
	Ry(-xp, &r)
	Rx(-yp, &r)
	Rz(elong, &r)

	/* Solve for local Earth rotation angle. */
	a = r[0][0]
	b = r[0][1]
	// eral = ( a != 0.0 || b != 0.0 ) ?  atan2(b, a) : 0.0;
	if a != 0.0 || b != 0.0 {
		eral = atan2(b, a)
	} else {
		eral = 0.0
	}
	astrom.Eral = eral

	/* Solve for polar motion [X,Y] with respect to local meridian. */
	a = r[0][0]
	c = r[0][2]
	astrom.Xpl = atan2(c, sqrt(a*a+b*b))
	a = r[1][2]
	b = r[2][2]
	// astrom.Ypl = ( a != 0.0 || b != 0.0 ) ? -atan2(a, b) : 0.0;
	if a != 0.0 || b != 0.0 {
		astrom.Ypl = -atan2(a, b)
	} else {
		astrom.Ypl = 0.0
	}

	/* Adjusted longitude. */
	astrom.Along = Anpm(eral - theta)

	/* Functions of latitude. */
	astrom.Sphi = sin(phi)
	astrom.Cphi = cos(phi)

	/* Refraction constants. */
	astrom.Refa = refa
	astrom.Refb = refb

	/* Disable the (redundant) diurnal aberration step. */
	astrom.Diurab = 0.0

	/* CIO based BPN matrix. */
	C2ixys(x, y, s, &r)

	/* Observer's geocentric position and velocity (m, m/s, CIRS). */
	Pvtob(elong, phi, hm, xp, yp, sp, theta, &pvc)

	/* Rotate into GCRS. */
	Trxpv(r, pvc, &pv)

	/* ICRS <-> GCRS parameters. */
	Apcs(date1, date2, pv, ebpv, ehp, astrom)

	/* Store the CIO based BPN matrix. */
	Cr(r, &astrom.Bpn)

}

/*
Apco13 Prepare for ICRS <-> observed, terrestrial

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed
coordinates.  The caller supplies UTC, site coordinates, ambient air
conditions and observing wavelength, and SOFA models are used to
obtain the Earth ephemeris, CIP/CIO and refraction constants.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the ICRS/CIRS transformations.

Given:

	utc1   float64     UTC as a 2-part...
	utc2   float64     ...quasi Julian Date (Notes 1,2)
	dut1   float64     UT1-UTC (seconds, Note 3)
	elong  float64     longitude (radians, east +ve, Note 4)
	phi    float64     latitude (geodetic, radians, Note 4)
	hm     float64     height above ellipsoid (m, geodetic, Notes 4,6)
	xp,yp  float64     polar motion coordinates (radians, Note 5)
	phpa   float64     pressure at the observer (hPa = mB, Note 6)
	tc     float64     ambient temperature at the observer (deg C)
	rh     float64     relative humidity at the observer (range 0-1)
	wl     float64     wavelength (micrometers, Note 7)

Returned:

	astrom ASTROM     star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       adjusted longitude (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)
	eo     float64       equation of the origins (ERA-GST)

Returned (function value):

	int        status: +1 = dubious year (Note 2)
	                    0 = OK
	                   -1 = unacceptable date

Notes:

 1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 2. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See Dat for further details.

 3. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 4. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 5. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto
    the local meridian.

 6. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

 7. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

 8. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

 9. In cases where the caller wishes to supply his own Earth
    ephemeris, Earth rotation information and refraction constants,
    the function Apco can be used instead of the present function.

    10)This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

    11)The context structure astrom produced by this function is used
    by Atioq, Atoiq, Atciq* and Aticq*.

Called:

	Utctai    UTC to TAI
	Taitt     TAI to TT
	Utcut1    UTC to UT1
	Epv00     Earth position and velocity
	Pnm06a    classical NPB matrix, IAU 2006/2000A
	Bpn2xy    extract CIP X,Y coordinates from NPB matrix
	S06       the CIO locator s, given X,Y, IAU 2006
	Era00     Earth rotation angle, IAU 2000
	Sp00      the TIO locator s', IERS 2000
	Refco     refraction constants for given ambient conditions
	Apco      astrometry parameters, ICRS-observed
	Eors      equation of the origins, given NPB matrix and s
*/
func Apco13(utc1, utc2, dut1 float64, elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64, astrom *ASTROM, eo *float64) int {
	var j int
	var tai1, tai2, tt1, tt2, ut11, ut12 float64
	var ehpv, ebpv [2][3]float64
	var r [3][3]float64
	var x, y, s, theta, sp, refa, refb float64

	/* UTC to other time scales. */
	j = Utctai(utc1, utc2, &tai1, &tai2)
	if j < 0 {
		return -1
	}

	Taitt(tai1, tai2, &tt1, &tt2)

	j = Utcut1(utc1, utc2, dut1, &ut11, &ut12)
	if j < 0 {
		return -1
	}

	/* Earth barycentric & heliocentric position/velocity (au, au/d). */
	Epv00(tt1, tt2, &ehpv, &ebpv)

	/* Form the equinox based BPN matrix, IAU 2006/2000A. */
	Pnm06a(tt1, tt2, &r)

	/* Extract CIP X,Y. */
	Bpn2xy(r, &x, &y)

	/* Obtain CIO locator s. */
	s = S06(tt1, tt2, x, y)

	/* Earth rotation angle. */
	theta = Era00(ut11, ut12)

	/* TIO locator s'. */
	sp = Sp00(tt1, tt2)

	/* Refraction constants A and B. */
	Refco(phpa, tc, rh, wl, &refa, &refb)

	/* Compute the star-independent astrometry parameters. */
	Apco(tt1, tt2, ebpv, ehpv[0], x, y, s, theta,
		elong, phi, hm, xp, yp, sp, refa, refb, astrom)

	/* Equation of the origins. */
	*eo = Eors(r, s)

	/* Return any warning status. */
	return j

}

/*
Apcs Prepare for ICRS <-> CIRS, space, special

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is supplied by the
caller.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

Given:

	date1  float64       TDB as a 2-part...
	date2  float64       ...Julian Date (Note 1)
	pv     [2][3]float64 observer's geocentric pos/vel (m, m/s)
	ebpv   [2][3]float64 Earth barycentric PV (au, au/day)
	ehp    [3]float64    Earth heliocentric P (au)

Returned:

	astrom ASTROM   star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       unchanged
	xpl    float64       unchanged
	ypl    float64       unchanged
	sphi   float64       unchanged
	cphi   float64       unchanged
	diurab float64       unchanged
	eral   float64       unchanged
	refa   float64       unchanged
	refb   float64       unchanged

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. Providing separate arguments for (i) the observer's geocentric
    position and velocity and (ii) the Earth ephemeris is done for
    convenience in the geocentric, terrestrial and Earth orbit cases.
    For deep space applications it maybe more convenient to specify
    zero geocentric position and velocity and to supply the
    observer's position and velocity information directly instead of
    with respect to the Earth.  However, note the different units:
    m and m/s for the geocentric vectors, au and au/day for the
    heliocentric and barycentric vectors.

 4. In cases where the caller does not wish to provide the Earth
    ephemeris, the function Apcs13 can be used instead of the
    present function.  This computes the Earth ephemeris using the
    SOFA function Epv00.

 5. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 6. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Cp        copy p-vector
	Pm        modulus of p-vector
	Pn        decompose p-vector into modulus and direction
	Ir        initialize r-matrix to identity
*/
func Apcs(date1, date2 float64, pv, ebpv [2][3]float64, ehp [3]float64, astrom *ASTROM) {
	/* au/d to m/s */
	const AUDMS = DAU / DAYSEC

	/* Light time for 1 au (day) */
	const CR = AULT / DAYSEC

	var i int
	var dp, dv float64
	var pb, vb, ph [3]float64
	var v2, w float64

	/* Time since reference epoch, years (for proper motion calculation). */
	astrom.Pmt = ((date1 - DJ00) + date2) / DJY

	/* Adjust Earth ephemeris to observer. */
	for i = 0; i < 3; i++ {
		dp = pv[0][i] / DAU
		dv = pv[1][i] / AUDMS
		pb[i] = ebpv[0][i] + dp
		vb[i] = ebpv[1][i] + dv
		ph[i] = ehp[i] + dp
	}

	/* Barycentric position of observer (au). */
	Cp(pb, &astrom.Eb)

	/* Heliocentric direction and distance (unit vector and au). */
	Pn(ph, &astrom.Em, &astrom.Eh)

	/* Barycentric vel. in units of c, and reciprocal of Lorenz factor. */
	v2 = 0.0
	for i = 0; i < 3; i++ {
		w = vb[i] * CR
		astrom.V[i] = w
		v2 += w * w
	}
	astrom.Bm1 = sqrt(1.0 - v2)

	/* Reset the NPB matrix. */
	Ir(&astrom.Bpn)
}

/*
Apcs13 Prepare for ICRS <-> CIRS, space

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is from SOFA models.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

Given:

	date1  float64       TDB as a 2-part...
	date2  float64       ...Julian Date (Note 1)
	pv     [2][3]float64 observer's geocentric pos/vel (Note 3)

Returned:

	astrom ASTROM   star-independent astrometry parameters:
	pmt    float64       PM time interval (SSB, Julian years)
	eb     [3]float64    SSB to observer (vector, au)
	eh     [3]float64    Sun to observer (unit vector)
	em     float64       distance from Sun to observer (au)
	v      [3]float64    barycentric observer velocity (vector, c)
	bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	bpn    [3][3]float64 bias-precession-nutation matrix
	along  float64       unchanged
	xpl    float64       unchanged
	ypl    float64       unchanged
	sphi   float64       unchanged
	cphi   float64       unchanged
	diurab float64       unchanged
	eral   float64       unchanged
	refa   float64       unchanged
	refb   float64       unchanged

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. All the vectors are with respect to BCRS axes.

 3. The observer's position and velocity pv are geocentric but with
    respect to BCRS axes, and in units of m and m/s.  No assumptions
    are made about proximity to the Earth, and the function can be
    used for deep space applications as well as Earth orbit and
    terrestrial.

 4. In cases where the caller wishes to supply his own Earth
    ephemeris, the function Apcs can be used instead of the present
    function.

 5. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 6. The context structure astrom produced by this function is used by
    Atciq* and Aticq*.

Called:

	Epv00     Earth position and velocity
	Apcs      astrometry parameters, ICRS-GCRS, space observer
*/
func Apcs13(date1, date2 float64, pv [2][3]float64, astrom *ASTROM) {
	var ehpv, ebpv [2][3]float64

	/* Earth barycentric & heliocentric position/velocity (au, au/d). */
	Epv00(date1, date2, &ehpv, &ebpv)

	/* Compute the star-independent astrometry parameters. */
	Apcs(date1, date2, pv, ebpv, ehpv[0], astrom)
}

/*
Aper Insert ERA into context

In the star-independent astrometry parameters, update only the
Earth rotation angle, supplied by the caller explicitly.

Given:

	theta   float64       Earth rotation angle (radians, Note 2)
	astrom  ASTROM        star-independent astrometry parameters:
	 pmt    float64       not used
	 eb     [3]float64    not used
	 eh     [3]float64    not used
	 em     float64       not used
	 v      [3]float64    not used
	 bm1    float64       not used
	 bpn    [3][3]float64 not used
	 along  float64       longitude + s' (radians)
	 xpl    float64       not used
	 ypl    float64       not used
	 sphi   float64       not used
	 cphi   float64       not used
	 diurab float64       not used
	 eral   float64       not used
	 refa   float64       not used
	 refb   float64       not used

Returned:

	astrom  ASTROM     star-independent astrometry parameters:
	 pmt    float64       unchanged
	 eb     [3]float64    unchanged
	 eh     [3]float64    unchanged
	 em     float64       unchanged
	 v      [3]float64    unchanged
	 bm1    float64       unchanged
	 bpn    [3][3]float64 unchanged
	 along  float64       unchanged
	 xpl    float64       unchanged
	 ypl    float64       unchanged
	 sphi   float64       unchanged
	 cphi   float64       unchanged
	 diurab float64       unchanged
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       unchanged
	 refb   float64       unchanged

Notes:

 1. This function exists to enable sidereal-tracking applications to
    avoid wasteful recomputation of the bulk of the astrometry
    parameters:  only the Earth rotation is updated.

 2. For targets expressed as equinox based positions, such as
    classical geocentric apparent (RA,Dec), the supplied theta can be
    Greenwich apparent sidereal time rather than Earth rotation
    angle.

 3. The function Aper13 can be used instead of the present
    function, and starts from UT1 rather than ERA itself.

 4. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.
*/
func Aper(theta float64, astrom *ASTROM) {
	astrom.Eral = theta + astrom.Along
}

/*
Aper13 Update context for Earth rotation

In the star-independent astrometry parameters, update only the
Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

Given:

	ut11    float64       UT1 as a 2-part...
	ut12    float64       ...Julian Date (Note 1)
	astrom  ASTROM        star-independent astrometry parameters:
	 pmt    float64       not used
	 eb     [3]float64    not used
	 eh     [3]float64    not used
	 em     float64       not used
	 v      [3]float64    not used
	 bm1    float64       not used
	 bpn    [3][3]float64 not used
	 along  float64       longitude + s' (radians)
	 xpl    float64       not used
	 ypl    float64       not used
	 sphi   float64       not used
	 cphi   float64       not used
	 diurab float64       not used
	 eral   float64       not used
	 refa   float64       not used
	 refb   float64       not used

Returned:

	astrom  ASTROM     star-independent astrometry parameters:
	 pmt    float64       unchanged
	 eb     [3]float64    unchanged
	 eh     [3]float64    unchanged
	 em     float64       unchanged
	 v      [3]float64    unchanged
	 bm1    float64       unchanged
	 bpn    [3][3]float64 unchanged
	 along  float64       unchanged
	 xpl    float64       unchanged
	 ypl    float64       unchanged
	 sphi   float64       unchanged
	 cphi   float64       unchanged
	 diurab float64       unchanged
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       unchanged
	 refb   float64       unchanged

Notes:

 1. The UT1 date (n.b. not UTC) ut11+ut12 is a Julian Date,
    apportioned in any convenient way between the arguments ut11 and
    ut12.  For example, JD(UT1)=2450123.7 could be expressed in any
    of these ways, among others:

    ut11           ut12

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 and MJD methods are good compromises
    between resolution and convenience.  The date & time method is
    best matched to the algorithm used:  maximum precision is
    delivered when the ut11 argument is for 0hrs UT1 on the day in
    question and the ut12 argument lies in the range 0 to 1, or vice
    versa.

 2. If the caller wishes to provide the Earth rotation angle itself,
    the function Aper can be used instead.  One use of this
    technique is to substitute Greenwich apparent sidereal time and
    thereby to support equinox based transformations directly.

 3. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

Called:

	Aper      astrometry parameters: update ERA
	Era00     Earth rotation angle, IAU 2000
*/
func Aper13(ut11, ut12 float64, astrom *ASTROM) {
	Aper(Era00(ut11, ut12), astrom)
}

/*
Apio Prepare for CIRS <-> observed, terrestrial, special

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed
coordinates.  The caller supplies the Earth orientation information
and the refraction constants as well as the site coordinates.

Given:

	sp     float64      the TIO locator s' (radians, Note 1)
	theta  float64      Earth rotation angle (radians)
	elong  float64      longitude (radians, east +ve, Note 2)
	phi    float64      geodetic latitude (radians, Note 2)
	hm     float64      height above ellipsoid (m, geodetic Note 2)
	xp,yp  float64      polar motion coordinates (radians, Note 3)
	refa   float64      refraction constant A (radians, Note 4)
	refb   float64      refraction constant B (radians, Note 4)

Returned:

	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       unchanged
	 eb     [3]float64    unchanged
	 eh     [3]float64    unchanged
	 em     float64       unchanged
	 v      [3]float64    unchanged
	 bm1    float64       unchanged
	 bpn    [3][3]float64 unchanged
	 along  float64       adjusted longitude (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Notes:

 1. sp, the TIO locator s', is a tiny quantity needed only by the
    most precise applications.  It can either be set to zero or
    predicted using the SOFA function Sp00.

 2. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 3. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many applications,
    xp and yp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto the
    local meridian.

 4. The refraction constants refa and refb are for use in a
    dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
    (i.e. refracted) zenith distance and dZ is the amount of
    refraction.

 5. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

 6. In cases where the caller does not wish to provide the Earth
    rotation information and refraction constants, the function
    Apio13 can be used instead of the present function.  This
    starts from UTC and weather readings etc. and computes suitable
    values using other SOFA functions.

 7. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

 8. The context structure astrom produced by this function is used by
    Atioq and Atoiq.

Called:

	Ir        initialize r-matrix to identity
	Rz        rotate around Z-axis
	Ry        rotate around Y-axis
	Rx        rotate around X-axis
	Anpm      normalize angle into range +/- pi
	Pvtob     position/velocity of terrestrial station
*/
func Apio(sp, theta float64, elong, phi, hm, xp, yp float64, refa, refb float64, astrom *ASTROM) {
	var r [3][3]float64
	var a, b, eral, c float64
	var pv [2][3]float64

	/* Form the rotation matrix, CIRS to apparent [HA,Dec]. */
	Ir(&r)
	Rz(theta+sp, &r)
	Ry(-xp, &r)
	Rx(-yp, &r)
	Rz(elong, &r)

	/* Solve for local Earth rotation angle. */
	a = r[0][0]
	b = r[0][1]
	// eral = ( a != 0.0 || b != 0.0 ) ?  atan2(b, a) : 0.0;
	if a != 0.0 || b != 0.0 {
		eral = atan2(b, a)
	} else {
		eral = 0.0
	}
	astrom.Eral = eral

	/* Solve for polar motion [X,Y] with respect to local meridian. */
	a = r[0][0]
	c = r[0][2]
	astrom.Xpl = atan2(c, sqrt(a*a+b*b))
	a = r[1][2]
	b = r[2][2]
	// astrom.Ypl = ( a != 0.0 || b != 0.0 ) ? -atan2(a, b) : 0.0;
	if a != 0.0 || b != 0.0 {
		astrom.Ypl = -atan2(a, b)
	} else {
		astrom.Ypl = 0.0
	}

	/* Adjusted longitude. */
	astrom.Along = Anpm(eral - theta)

	/* Functions of latitude. */
	astrom.Sphi = sin(phi)
	astrom.Cphi = cos(phi)

	/* Observer's geocentric position and velocity (m, m/s, CIRS). */
	Pvtob(elong, phi, hm, xp, yp, sp, theta, &pv)

	/* Magnitude of diurnal aberration vector. */
	astrom.Diurab = sqrt(pv[1][0]*pv[1][0]+pv[1][1]*pv[1][1]) / CMPS

	/* Refraction constants. */
	astrom.Refa = refa
	astrom.Refb = refb
}

/*
Apio13 Prepare for CIRS <-> observed, terrestrial

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed
coordinates.  The caller supplies UTC, site coordinates, ambient air
conditions and observing wavelength.

Given:

	utc1   float64      UTC as a 2-part...
	utc2   float64      ...quasi Julian Date (Notes 1,2)
	dut1   float64      UT1-UTC (seconds)
	elong  float64      longitude (radians, east +ve, Note 3)
	phi    float64      geodetic latitude (radians, Note 3)
	hm     float64      height above ellipsoid (m, geodetic Notes 4,6)
	xp,yp  float64      polar motion coordinates (radians, Note 5)
	phpa   float64      pressure at the observer (hPa = mB, Note 6)
	tc     float64      ambient temperature at the observer (deg C)
	rh     float64      relative humidity at the observer (range 0-1)
	wl     float64      wavelength (micrometers, Note 7)

Returned:

	astrom  ASTROM     star-independent astrometry parameters:
	 pmt    float64       unchanged
	 eb     [3]float64    unchanged
	 eh     [3]float64    unchanged
	 em     float64       unchanged
	 v      [3]float64    unchanged
	 bm1    float64       unchanged
	 bpn    [3][3]float64 unchanged
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned (function value):

	int         status: +1 = dubious year (Note 2)
	                     0 = OK
	                    -1 = unacceptable date

Notes:

 1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 2. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See Dat for further details.

 3. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 4. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 5. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many applications,
    xp and yp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto
    the local meridian.

 6. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to the
    pressure and that an accurate phpa value is important for
    precise work.

 7. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

 8. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

 9. In cases where the caller wishes to supply his own Earth
    rotation information and refraction constants, the function
    Apc can be used instead of the present function.

    10)This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

    functions         observer        transformation

    Apcg Apcg13      geocentric      ICRS <-> GCRS
    Apci Apci13      terrestrial     ICRS <-> CIRS
    Apco Apco13      terrestrial     ICRS <-> observed
    Apcs Apcs13      space           ICRS <-> GCRS
    Aper Aper13      terrestrial     update Earth rotation
    Apio Apio13      terrestrial     CIRS <-> observed

    Those with names ending in "13" use contemporary SOFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

    11)The context structure astrom produced by this function is used
    by Atioq and Atoiq.

Called:

	Utctai    UTC to TAI
	Taitt     TAI to TT
	Utcut1    UTC to UT1
	Sp00      the TIO locator s', IERS 2000
	Era00     Earth rotation angle, IAU 2000
	Refco     refraction constants for given ambient conditions
	Apio      astrometry parameters, CIRS-observed
*/
func Apio13(utc1, utc2, dut1 float64, elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64, astrom *ASTROM) int {
	var j int
	var tai1, tai2, tt1, tt2, ut11, ut12, sp, theta, refa, refb float64

	/* UTC to other time scales. */
	j = Utctai(utc1, utc2, &tai1, &tai2)
	if j < 0 {
		return -1
	}
	Taitt(tai1, tai2, &tt1, &tt2)
	j = Utcut1(utc1, utc2, dut1, &ut11, &ut12)
	if j < 0 {
		return -1
	}

	/* TIO locator s'. */
	sp = Sp00(tt1, tt2)

	/* Earth rotation angle. */
	theta = Era00(ut11, ut12)

	/* Refraction constants A and B. */
	Refco(phpa, tc, rh, wl, &refa, &refb)

	/* CIRS <-> observed astrometry parameters. */
	Apio(sp, theta, elong, phi, hm, xp, yp, refa, refb, astrom)

	/* Return any warning status. */
	return j
}

/*
Atcc13 J2000.0 catalog ICRS -> astrometric ICRS

Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
astrometric place.

Given:

	rc     float64   ICRS right ascension at J2000.0 (radians, Note 1)
	dc     float64   ICRS declination at J2000.0 (radians, Note 1)
	pr     float64   RA proper motion (radians/year, Note 2)
	pd     float64   Dec proper motion (radians/year)
	px     float64   parallax (arcsec)
	rv     float64   radial velocity (km/s, +ve if receding)
	date1  float64   TDB as a 2-part...
	date2  float64   ...Julian Date (Note 3)

Returned:

	ra,da  float64  ICRS astrometric RA,Dec (radians)

Notes:

 1. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require a
    preliminary call to Pmsafe before use.

 2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 3. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

Called:

	Apci13    astrometry parameters, ICRS-CIRS, 2013
	Atccq     quick catalog ICRS to astrometric
*/
func Atcc13(rc, dc float64, pr, pd, px, rv float64, date1, date2 float64, ra, da *float64) {
	/* Star-independent astrometry parameters */
	var astrom ASTROM

	var w float64

	/* The transformation parameters. */
	Apci13(date1, date2, &astrom, &w)

	/* Catalog ICRS (epoch J2000.0) to astrometric. */
	Atccq(rc, dc, pr, pd, px, rv, &astrom, ra, da)
}

/*
Atccq Astrometric ICRS -> J2000.0 catalog ICRS

Quick transformation of a star's ICRS catalog entry (epoch J2000.0)
into ICRS astrometric place, given precomputed star-independent
astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions Apci[13], Apcg[13], Apco[13] or Apcs[13].

If the parallax and proper motions are zero the transformation has
no effect.

Given:

	rc,dc  float64       ICRS RA,Dec at J2000.0 (radians)
	pr     float64       RA proper motion (radians/year, Note 3)
	pd     float64       Dec proper motion (radians/year)
	px     float64       parallax (arcsec)
	rv     float64       radial velocity (km/s, +ve if receding)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	ra,da  float64    ICRS astrometric RA,Dec (radians)

Notes:

 1. All the vectors are with respect to BCRS axes.

 2. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require a
    preliminary call to Pmsafe before use.

 3. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

Called:

	Pmpx      proper motion and parallax
	C2s       p-vector to spherical
	Anp       normalize angle into range 0 to 2pi
*/
func Atccq(rc, dc float64, pr, pd, px, rv float64, astrom *ASTROM, ra, da *float64) {
	var p [3]float64
	var w float64

	/* Proper motion and parallax, giving BCRS coordinate direction. */
	Pmpx(rc, dc, pr, pd, px, rv, astrom.Pmt, astrom.Eb, &p)

	/* ICRS astrometric RA,Dec. */
	C2s(p, &w, da)
	*ra = Anp(w)
}

/*
Atci13 Catalog -> CIRS

Transform ICRS star data, epoch J2000.0, to CIRS.

Given:

	rc     float64   ICRS right ascension at J2000.0 (radians, Note 1)
	dc     float64   ICRS declination at J2000.0 (radians, Note 1)
	pr     float64   RA proper motion (radians/year, Note 2)
	pd     float64   Dec proper motion (radians/year)
	px     float64   parallax (arcsec)
	rv     float64   radial velocity (km/s, +ve if receding)
	date1  float64   TDB as a 2-part...
	date2  float64   ...Julian Date (Note 3)

Returned:

	ri,di  float64  CIRS geocentric RA,Dec (radians)
	eo     float64  equation of the origins (ERA-GST, Note 5)

Notes:

 1. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require a
    preliminary call to Pmsafe before use.

 2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 3. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 4. The available accuracy is better than 1 milliarcsecond, limited
    mainly by the precession-nutation model that is used, namely
    2000A/2006.  Very close to solar system bodies, additional
    errors of up to several milliarcseconds can occur because of
    unmodeled light deflection;  however, the Sun's contribution is
    taken into account, to first order.  The accuracy limitations of
    the SOFA function Epv00 (used to compute Earth position and
    velocity) can contribute aberration errors of up to
    5 microarcseconds.  Light deflection at the Sun's limb is
    uncertain at the 0.4 mas level.

 5. Should the transformation to (equinox based) apparent place be
    required rather than (CIO based) intermediate place, subtract the
    equation of the origins from the returned right ascension:
    RA = RI - EO. (The Anp function can then be applied, as
    required, to keep the result in the conventional 0-2pi range.)

Called:

	Apci13    astrometry parameters, ICRS-CIRS, 2013
	Atciq     quick ICRS to CIRS
*/
func Atci13(rc, dc float64, pr, pd, px, rv float64, date1, date2 float64, ri, di, eo *float64) {
	/* Star-independent astrometry parameters */
	var astrom ASTROM

	/* The transformation parameters. */
	Apci13(date1, date2, &astrom, eo)

	/* ICRS (epoch J2000.0) to CIRS. */
	Atciq(rc, dc, pr, pd, px, rv, &astrom, ri, di)
}

/*
Atciq Quick ICRS -> CIRS

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions Apci[13], Apcg[13], Apco[13] or Apcs[13].

If the parallax and proper motions are zero the Atciqz function
can be used instead.

Given:

	rc,dc  float64       ICRS RA,Dec at J2000.0 (radians)
	pr     float64       RA proper motion (radians/year, Note 3)
	pd     float64       Dec proper motion (radians/year)
	px     float64       parallax (arcsec)
	rv     float64       radial velocity (km/s, +ve if receding)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	ri,di   float64    CIRS RA,Dec (radians)

Notes:

 1. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require a
    preliminary call to Pmsafe before use.

 2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

Called:

	Pmpx      proper motion and parallax
	Ldsun     light deflection by the Sun
	Ab        stellar aberration
	Rxp       product of r-matrix and pv-vector
	C2s       p-vector to spherical
	Anp       normalize angle into range 0 to 2pi
*/
func Atciq(rc, dc float64, pr, pd, px, rv float64, astrom *ASTROM, ri, di *float64) {
	var pco, pnat, ppr, pi [3]float64
	var w float64

	/* Proper motion and parallax, giving BCRS coordinate direction. */
	Pmpx(rc, dc, pr, pd, px, rv, astrom.Pmt, astrom.Eb, &pco)

	/* Light deflection by the Sun, giving BCRS natural direction. */
	Ldsun(pco, astrom.Eh, astrom.Em, &pnat)

	/* Aberration, giving GCRS proper direction. */
	Ab(pnat, astrom.V, astrom.Em, astrom.Bm1, &ppr)

	/* Bias-precession-nutation, giving CIRS proper direction. */
	Rxp(astrom.Bpn, ppr, &pi)

	/* CIRS RA,Dec. */
	C2s(pi, &w, di)
	*ri = Anp(w)
}

/*
Atciqn Quick ICRS -> CIRS, multiple deflections

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters plus a list of light-
deflecting bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions Apci[13], Apcg[13], Apco[13] or Apcs[13].

If the only light-deflecting body to be taken into account is the
Sun, the Atciq function can be used instead.  If in addition the
parallax and proper motions are zero, the Atciqz function can be
used.

Given:

	rc,dc  float64       ICRS RA,Dec at J2000.0 (radians)
	pr     float64       RA proper motion (radians/year, Note 3)
	pd     float64       Dec proper motion (radians/year)
	px     float64       parallax (arcsec)
	rv     float64       radial velocity (km/s, +ve if receding)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)
	n      int          number of bodies (Note 3)
	b      LDBODY[n]    data for each of the n bodies (Notes 3,4):
	 bm     float64       mass of the body (solar masses, Note 5)
	 dl     float64       deflection limiter (Note 6)
	 pv     [2][3]float64 barycentric PV of the body (au, au/day)

Returned:

	ri,di   float64    CIRS RA,Dec (radians)

Notes:

 1. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require a
    preliminary call to Pmsafe before use.

 2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 3. The struct b contains n entries, one for each body to be
    considered.  If n = 0, no gravitational light deflection will be
    applied, not even for the Sun.

 4. The struct b should include an entry for the Sun as well as for
    any planet or other body to be taken into account.  The entries
    should be in the order in which the light passes the body.

 5. In the entry in the b struct for body i, the mass parameter
    b[i].bm can, as required, be adjusted in order to allow for such
    effects as quadrupole field.

 6. The deflection limiter parameter b[i].dl is phi^2/2, where phi is
    the angular separation (in radians) between star and body at
    which limiting is applied.  As phi shrinks below the chosen
    threshold, the deflection is artificially reduced, reaching zero
    for phi = 0.   Example values suitable for a terrestrial
    observer, together with masses, are as follows:

    body i     b[i].bm        b[i].dl

    Sun        1.0            6e-6
    Jupiter    0.00095435     3e-9
    Saturn     0.00028574     3e-10

 7. For efficiency, validation of the contents of the b array is
    omitted.  The supplied masses must be greater than zero, the
    position and velocity vectors must be right, and the deflection
    limiter greater than zero.

Called:

	Pmpx      proper motion and parallax
	Ldn       light deflection by n bodies
	Ab        stellar aberration
	Rxp       product of r-matrix and pv-vector
	C2s       p-vector to spherical
	Anp       normalize angle into range 0 to 2pi
*/
func Atciqn(rc, dc float64, pr, pd, px, rv float64, astrom *ASTROM, n int, b []LDBODY, ri, di *float64) {
	var pco, pnat, ppr, pi [3]float64
	var w float64

	/* Proper motion and parallax, giving BCRS coordinate direction. */
	Pmpx(rc, dc, pr, pd, px, rv, astrom.Pmt, astrom.Eb, &pco)

	/* Light deflection, giving BCRS natural direction. */
	Ldn(n, b, astrom.Eb, pco, &pnat)

	/* Aberration, giving GCRS proper direction. */
	Ab(pnat, astrom.V, astrom.Em, astrom.Bm1, &ppr)

	/* Bias-precession-nutation, giving CIRS proper direction. */
	Rxp(astrom.Bpn, ppr, &pi)

	/* CIRS RA,Dec. */
	C2s(pi, &w, di)
	*ri = Anp(w)
}

/*
Atciqz Quick astrometric ICRS -> CIRS

Quick ICRS to CIRS transformation, given precomputed star-
independent astrometry parameters, and assuming zero parallax and
proper motion.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions Apci[13], Apcg[13], Apco[13] or Apcs[13].

The corresponding function for the case of non-zero parallax and
proper motion is Atciq.

Given:

	rc,dc  float64       ICRS astrometric RA,Dec (radians)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	ri,di  float64     CIRS RA,Dec (radians)

Note:

	All the vectors are with respect to BCRS axes.

References:

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013).

	Klioner, Sergei A., "A practical relativistic model for micro-
	arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

Called:

	S2c       spherical coordinates to unit vector
	Ldsun     light deflection due to Sun
	Ab        stellar aberration
	Rxp       product of r-matrix and p-vector
	C2s       p-vector to spherical
	Anp       normalize angle into range +/- pi
*/
func Atciqz(rc, dc float64, astrom *ASTROM, ri, di *float64) {
	var pco, pnat, ppr, pi [3]float64
	var w float64

	/* BCRS coordinate direction (unit vector). */
	S2c(rc, dc, &pco)

	/* Light deflection by the Sun, giving BCRS natural direction. */
	Ldsun(pco, astrom.Eh, astrom.Em, &pnat)

	/* Aberration, giving GCRS proper direction. */
	Ab(pnat, astrom.V, astrom.Em, astrom.Bm1, &ppr)

	/* Bias-precession-nutation, giving CIRS proper direction. */
	Rxp(astrom.Bpn, ppr, &pi)

	/* CIRS RA,Dec. */
	C2s(pi, &w, di)
	*ri = Anp(w)
}

/*
Atco13 ICRS -> observed

ICRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

SOFA models are used for the Earth ephemeris, bias-precession-
nutation, Earth orientation and refraction.

Given:

	rc,dc  float64   ICRS right ascension at J2000.0 (radians, Note 1)
	pr     float64   RA proper motion (radians/year, Note 2)
	pd     float64   Dec proper motion (radians/year)
	px     float64   parallax (arcsec)
	rv     float64   radial velocity (km/s, +ve if receding)
	utc1   float64   UTC as a 2-part...
	utc2   float64   ...quasi Julian Date (Notes 3-4)
	dut1   float64   UT1-UTC (seconds, Note 5)
	elong  float64   longitude (radians, east +ve, Note 6)
	phi    float64   latitude (geodetic, radians, Note 6)
	hm     float64   height above ellipsoid (m, geodetic, Notes 6,8)
	xp,yp  float64   polar motion coordinates (radians, Note 7)
	phpa   float64   pressure at the observer (hPa = mB, Note 8)
	tc     float64   ambient temperature at the observer (deg C)
	rh     float64   relative humidity at the observer (range 0-1)
	wl     float64   wavelength (micrometers, Note 9)

Returned:

	aob    float64  observed azimuth (radians: N=0,E=90)
	zob    float64  observed zenith distance (radians)
	hob    float64  observed hour angle (radians)
	dob    float64  observed declination (radians)
	rob    float64  observed right ascension (CIO-based, radians)
	eo     float64  equation of the origins (ERA-GST)

Returned (function value):

	int      status: +1 = dubious year (Note 4)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require
    a preliminary call to Pmsafe before use.

 2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 3. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 4. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See Dat for further details.

 5. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 6. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 7. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

 8. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB),
    is available, an adequate estimate of hm can be obtained from
    the expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

 9. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

    10)The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted observed
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions Atco13 and
    Atoc13 are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

    11)"Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

    12)It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

Called:

	Apco13    astrometry parameters, ICRS-observed, 2013
	Atciq     quick ICRS to CIRS
	Atioq     quick CIRS to observed
*/
func Atco13(rc, dc float64, pr, pd, px, rv float64, utc1, utc2, dut1 float64,
	elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64,
	aob, zob, hob *float64, dob, rob, eo *float64) int {
	var j int
	var astrom ASTROM
	var ri, di float64

	/* Star-independent astrometry parameters. */
	j = Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
		phpa, tc, rh, wl, &astrom, eo)

	/* Abort if bad UTC. */
	if j < 0 {
		return j
	}

	/* Transform ICRS to CIRS. */
	Atciq(rc, dc, pr, pd, px, rv, &astrom, &ri, &di)

	/* Transform CIRS to observed. */
	Atioq(ri, di, &astrom, aob, zob, hob, dob, rob)

	/* Return OK/warning status. */
	return j
}

/*
Atic13 CIRS -> ICRS

Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

Given:

	ri,di  float64  CIRS geocentric RA,Dec (radians)
	date1  float64  TDB as a 2-part...
	date2  float64  ...Julian Date (Note 1)

Returned:

	rc,dc  float64  ICRS astrometric RA,Dec (radians)
	eo     float64  equation of the origins (ERA-GST, Note 4)

Notes:

 1. The TDB date date1+date2 is a Julian Date, apportioned in any
    convenient way between the two arguments.  For example,
    JD(TDB)=2450123.7 could be expressed in any of these ways, among
    others:

    date1          date2

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.  For most
    applications of this function the choice will not be at all
    critical.

    TT can be used instead of TDB without any significant impact on
    accuracy.

 2. Iterative techniques are used for the aberration and light
    deflection corrections so that the functions Atic13 (or
    Aticq) and Atci13 (or Atciq) are accurate inverses;
    even at the edge of the Sun's disk the discrepancy is only about
    1 nanoarcsecond.

 3. The available accuracy is better than 1 milliarcsecond, limited
    mainly by the precession-nutation model that is used, namely
    2000A/2006.  Very close to solar system bodies, additional
    errors of up to several milliarcseconds can occur because of
    unmodeled light deflection;  however, the Sun's contribution is
    taken into account, to first order.  The accuracy limitations of
    the SOFA function Epv00 (used to compute Earth position and
    velocity) can contribute aberration errors of up to
    5 microarcseconds.  Light deflection at the Sun's limb is
    uncertain at the 0.4 mas level.

 4. Should the transformation to (equinox based) J2000.0 mean place
    be required rather than (CIO based) ICRS coordinates, subtract the
    equation of the origins from the returned right ascension:
    RA = RI - EO.  (The Anp function can then be applied, as
    required, to keep the result in the conventional 0-2pi range.)

Called:

	Apci13    astrometry parameters, ICRS-CIRS, 2013
	Aticq     quick CIRS to ICRS astrometric
*/
func Atic13(ri, di float64, date1, date2 float64, rc, dc, eo *float64) {
	/* Star-independent astrometry parameters */
	var astrom ASTROM

	/* Star-independent astrometry parameters. */
	Apci13(date1, date2, &astrom, eo)

	/* CIRS to ICRS astrometric. */
	Aticq(ri, di, &astrom, rc, dc)
}

/*
Aticq Quick CIRS -> ICRS

Quick CIRS RA,Dec to ICRS astrometric place, given the star-
independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling one of the functions Apci[13], Apcg[13], Apco[13]
or Apcs[13].

Given:

	ri,di  float64       CIRS RA,Dec (radians)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	rc,dc  float64     ICRS astrometric RA,Dec (radians)

Notes:

 1. Only the Sun is taken into account in the light deflection
    correction.

 2. Iterative techniques are used for the aberration and light
    deflection corrections so that the functions Atic13 (or
    Aticq) and Atci13 (or Atciq) are accurate inverses;
    even at the edge of the Sun's disk the discrepancy is only about
    1 nanoarcsecond.

Called:

	S2c       spherical coordinates to unit vector
	Trxp      product of transpose of r-matrix and p-vector
	Zp        zero p-vector
	Ab        stellar aberration
	Ldsun     light deflection by the Sun
	C2s       p-vector to spherical
	Anp       normalize angle into range +/- pi
*/
func Aticq(ri, di float64, astrom *ASTROM, rc, dc *float64) {
	var j, i int
	var pi, ppr, pnat, pco, d, before, after [3]float64
	var w, r2, r float64

	/* CIRS RA,Dec to Cartesian. */
	S2c(ri, di, &pi)

	/* Bias-precession-nutation, giving GCRS proper direction. */
	Trxp(astrom.Bpn, pi, &ppr)

	/* Aberration, giving GCRS natural direction. */
	Zp(&d)
	for j = 0; j < 2; j++ {
		r2 = 0.0
		for i = 0; i < 3; i++ {
			w = ppr[i] - d[i]
			before[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			before[i] /= r
		}
		Ab(before, astrom.V, astrom.Em, astrom.Bm1, &after)
		r2 = 0.0
		for i = 0; i < 3; i++ {
			d[i] = after[i] - before[i]
			w = ppr[i] - d[i]
			pnat[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			pnat[i] /= r
		}
	}

	/* Light deflection by the Sun, giving BCRS coordinate direction. */
	Zp(&d)
	for j = 0; j < 5; j++ {
		r2 = 0.0
		for i = 0; i < 3; i++ {
			w = pnat[i] - d[i]
			before[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			before[i] /= r
		}
		Ldsun(before, astrom.Eh, astrom.Em, &after)
		r2 = 0.0
		for i = 0; i < 3; i++ {
			d[i] = after[i] - before[i]
			w = pnat[i] - d[i]
			pco[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			pco[i] /= r
		}
	}

	/* ICRS astrometric RA,Dec. */
	C2s(pco, &w, dc)
	*rc = Anp(w)
}

/*
Aticqn Quick CIRS -> ICRS, multiple deflections

Quick CIRS to ICRS astrometric place transformation, given the star-
independent astrometry parameters plus a list of light-deflecting
bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling one of the functions Apci[13], Apcg[13], Apco[13]
or Apcs[13].

If the only light-deflecting body to be taken into account is the
Sun, the Aticq function can be used instead.

Given:

	ri,di  double        CIRS RA,Dec (radians)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)
	n      int            number of bodies (Note 3)
	b      LDBODY[n]     data for each of the n bodies (Notes 3,4):
	 bm     float64       mass of the body (solar masses, Note 5)
	 dl     float64       deflection limiter (Note 6)
	 pv     [2][3]float64 barycentric PV of the body (au, au/day)

Returned:

	rc,dc  float64     ICRS astrometric RA,Dec (radians)

Notes:

 1. Iterative techniques are used for the aberration and light
    deflection corrections so that the functions Aticqn and
    Atciqn are accurate inverses; even at the edge of the Sun's
    disk the discrepancy is only about 1 nanoarcsecond.

 2. If the only light-deflecting body to be taken into account is the
    Sun, the Aticq function can be used instead.

 3. The struct b contains n entries, one for each body to be
    considered.  If n = 0, no gravitational light deflection will be
    applied, not even for the Sun.

 4. The struct b should include an entry for the Sun as well as for
    any planet or other body to be taken into account.  The entries
    should be in the order in which the light passes the body.

 5. In the entry in the b struct for body i, the mass parameter
    b[i].bm can, as required, be adjusted in order to allow for such
    effects as quadrupole field.

 6. The deflection limiter parameter b[i].dl is phi^2/2, where phi is
    the angular separation (in radians) between star and body at
    which limiting is applied.  As phi shrinks below the chosen
    threshold, the deflection is artificially reduced, reaching zero
    for phi = 0.   Example values suitable for a terrestrial
    observer, together with masses, are as follows:

    body i     b[i].bm        b[i].dl

    Sun        1.0            6e-6
    Jupiter    0.00095435     3e-9
    Saturn     0.00028574     3e-10

 7. For efficiency, validation of the contents of the b array is
    omitted.  The supplied masses must be greater than zero, the
    position and velocity vectors must be right, and the deflection
    limiter greater than zero.

Called:

	S2c       spherical coordinates to unit vector
	Trxp      product of transpose of r-matrix and p-vector
	Zp        zero p-vector
	Ab        stellar aberration
	Ldn       light deflection by n bodies
	C2s       p-vector to spherical
	Anp       normalize angle into range +/- pi
*/
func Aticqn(ri, di float64, astrom *ASTROM, n int, b []LDBODY, rc, dc *float64) {
	var j, i int
	var pi, ppr, pnat, pco, d, before, after [3]float64
	var w, r2, r float64

	/* CIRS RA,Dec to Cartesian. */
	S2c(ri, di, &pi)

	/* Bias-precession-nutation, giving GCRS proper direction. */
	Trxp(astrom.Bpn, pi, &ppr)

	/* Aberration, giving GCRS natural direction. */
	Zp(&d)
	for j = 0; j < 2; j++ {
		r2 = 0.0
		for i = 0; i < 3; i++ {
			w = ppr[i] - d[i]
			before[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			before[i] /= r
		}
		Ab(before, astrom.V, astrom.Em, astrom.Bm1, &after)
		r2 = 0.0
		for i = 0; i < 3; i++ {
			d[i] = after[i] - before[i]
			w = ppr[i] - d[i]
			pnat[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			pnat[i] /= r
		}
	}

	/* Light deflection, giving BCRS coordinate direction. */
	Zp(&d)
	for j = 0; j < 5; j++ {
		r2 = 0.0
		for i = 0; i < 3; i++ {
			w = pnat[i] - d[i]
			before[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			before[i] /= r
		}
		Ldn(n, b, astrom.Eb, before, &after)
		r2 = 0.0
		for i = 0; i < 3; i++ {
			d[i] = after[i] - before[i]
			w = pnat[i] - d[i]
			pco[i] = w
			r2 += w * w
		}
		r = sqrt(r2)
		for i = 0; i < 3; i++ {
			pco[i] /= r
		}
	}

	/* ICRS astrometric RA,Dec. */
	C2s(pco, &w, dc)
	*rc = Anp(w)
}

/*
Atio13 CIRS -> observed

CIRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

Given:

	ri     float64   CIRS right ascension (CIO-based, radians)
	di     float64   CIRS declination (radians)
	utc1   float64   UTC as a 2-part...
	utc2   float64   ...quasi Julian Date (Notes 1,2)
	dut1   float64   UT1-UTC (seconds, Note 3)
	elong  float64   longitude (radians, east +ve, Note 4)
	phi    float64   geodetic latitude (radians, Note 4)
	hm     float64   height above ellipsoid (m, geodetic Notes 4,6)
	xp,yp  float64   polar motion coordinates (radians, Note 5)
	phpa   float64   pressure at the observer (hPa = mB, Note 6)
	tc     float64   ambient temperature at the observer (deg C)
	rh     float64   relative humidity at the observer (range 0-1)
	wl     float64   wavelength (micrometers, Note 7)

Returned:

	aob    float64  observed azimuth (radians: N=0,E=90)
	zob    float64  observed zenith distance (radians)
	hob    float64  observed hour angle (radians)
	dob    float64  observed declination (radians)
	rob    float64  observed right ascension (CIO-based, radians)

Returned (function value):

	int      status: +1 = dubious year (Note 2)
	                  0 = OK
	                 -1 = unacceptable date

Notes:

 1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 2. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See Dat for further details.

 3. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 4. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 5. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

 6. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

 7. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

 8. "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

 9. The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    10)The complementary functions Atio13 and Atoi13 are self-
    consistent to better than 1 microarcsecond all over the
    celestial sphere.

    11)It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

Called:

	Apio13    astrometry parameters, CIRS-observed, 2013
	Atioq     quick CIRS to observed
*/
func Atio13(ri, di float64, utc1, utc2, dut1 float64, elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64,
	aob, zob, hob *float64, dob, rob *float64) int {
	var j int
	var astrom ASTROM

	/* Star-independent astrometry parameters for CIRS->observed. */
	j = Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
		phpa, tc, rh, wl, &astrom)

	/* Abort if bad UTC. */
	if j < 0 {
		return j
	}

	/* Transform CIRS to observed. */
	Atioq(ri, di, &astrom, aob, zob, hob, dob, rob)

	/* Return OK/warning status. */
	return j
}

/*
Atioq Quick CIRS -> observed

Quick CIRS to observed place transformation.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling Apio[13] or Apco[13].

Given:

	ri     float64       CIRS right ascension
	di     float64       CIRS declination
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	aob    float64  observed azimuth (radians: N=0,E=90)
	zob    float64  observed zenith distance (radians)
	hob    float64  observed hour angle (radians)
	dob    float64  observed declination (radians)
	rob    float64  observed right ascension (CIO-based, radians)

Notes:

 1. This function returns zenith distance rather than altitude in
    order to reflect the fact that no allowance is made for
    depression of the horizon.

 2. The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted observed
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions Atioq and
    Atoiq are self-consistent to better than 1 microarcsecond all
    over the celestial sphere.  With refraction included, consistency
    falls off at high zenith distances, but is still better than
    0.05 arcsec at 85 degrees.

 3. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

 4. The CIRS RA,Dec is obtained from a star catalog mean place by
    allowing for space motion, parallax, the Sun's gravitational lens
    effect, annual aberration and precession-nutation.  For star
    positions in the ICRS, these effects can be applied by means of
    the Atci13 (etc.) functions.  Starting from classical "mean
    place" systems, additional transformations will be needed first.

 5. "Observed" Az,El means the position that would be seen by a
    perfect geodetically aligned theodolite.  This is obtained from
    the CIRS RA,Dec by allowing for Earth orientation and diurnal
    aberration, rotating from equator to horizon coordinates, and
    then adjusting for refraction.  The HA,Dec is obtained by
    rotating back into equatorial coordinates, and is the position
    that would be seen by a perfect equatorial with its polar axis
    aligned to the Earth's axis of rotation.  Finally, the RA is
    obtained by subtracting the HA from the local ERA.

 6. The star-independent CIRS-to-observed-place parameters in ASTROM
    may be computed with Apio[13] or Apco[13].  If nothing has
    changed significantly except the time, Aper[13] may be used to
    perform the requisite adjustment to the astrom structure.

Called:

	S2c       spherical coordinates to unit vector
	C2s       p-vector to spherical
	Anp       normalize angle into range 0 to 2pi
*/
func Atioq(ri, di float64, astrom *ASTROM, aob, zob *float64, hob, dob, rob *float64) {
	/* Minimum cos(alt) and sin(alt) for refraction purposes */
	const CELMIN = 1e-6
	const SELMIN = 0.05

	var v [3]float64
	var x, y, z, sx, cx, sy, cy, xhd, yhd, zhd, f,
		xhdt, yhdt, zhdt, xaet, yaet, zaet, azobs, r, tz, w, del,
		cosdel, xaeo, yaeo, zaeo, zdobs, hmobs, dcobs, raobs float64

	/* CIRS RA,Dec to Cartesian -HA,Dec. */
	S2c(ri-astrom.Eral, di, &v)
	x = v[0]
	y = v[1]
	z = v[2]

	/* Polar motion. */
	sx = sin(astrom.Xpl)
	cx = cos(astrom.Xpl)
	sy = sin(astrom.Ypl)
	cy = cos(astrom.Ypl)
	xhd = cx*x + sx*z
	yhd = sx*sy*x + cy*y - cx*sy*z
	zhd = -sx*cy*x + sy*y + cx*cy*z

	/* Diurnal aberration. */
	f = (1.0 - astrom.Diurab*yhd)
	xhdt = f * xhd
	yhdt = f * (yhd + astrom.Diurab)
	zhdt = f * zhd

	/* Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90). */
	xaet = astrom.Sphi*xhdt - astrom.Cphi*zhdt
	yaet = yhdt
	zaet = astrom.Cphi*xhdt + astrom.Sphi*zhdt

	/* Azimuth (N=0,E=90). */
	// azobs = ( xaet != 0.0 || yaet != 0.0 ) ? atan2(yaet,-xaet) : 0.0;
	if xaet != 0.0 || yaet != 0.0 {
		azobs = atan2(yaet, -xaet)
	} else {
		azobs = 0.0
	}

	/* ---------- */
	/* Refraction */
	/* ---------- */

	/* Cosine and sine of altitude, with precautions. */
	r = sqrt(xaet*xaet + yaet*yaet)
	// r = r > CELMIN ? r : CELMIN;
	if r <= CELMIN {
		r = CELMIN
	}
	// z = zaet > SELMIN ? zaet : SELMIN;
	if zaet > SELMIN {
		z = zaet
	} else {
		z = SELMIN
	}

	/* A*tan(z)+B*tan^3(z) model, with Newton-Raphson correction. */
	tz = r / z
	w = astrom.Refb * tz * tz
	del = (astrom.Refa + w) * tz /
		(1.0 + (astrom.Refa+3.0*w)/(z*z))

	/* Apply the change, giving observed vector. */
	cosdel = 1.0 - del*del/2.0
	f = cosdel - del*z/r
	xaeo = xaet * f
	yaeo = yaet * f
	zaeo = cosdel*zaet + del*r

	/* Observed ZD. */
	zdobs = atan2(sqrt(xaeo*xaeo+yaeo*yaeo), zaeo)

	/* Az/El vector to HA,Dec vector (both right-handed). */
	v[0] = astrom.Sphi*xaeo + astrom.Cphi*zaeo
	v[1] = yaeo
	v[2] = -astrom.Cphi*xaeo + astrom.Sphi*zaeo

	/* To spherical -HA,Dec. */
	C2s(v, &hmobs, &dcobs)

	/* Right ascension (with respect to CIO). */
	raobs = astrom.Eral + hmobs

	/* Return the results. */
	*aob = Anp(azobs)
	*zob = zdobs
	*hob = -hmobs
	*dob = dcobs
	*rob = Anp(raobs)
}

/*
Atoc13 Observed -> astrometric ICRS

Observed place at a groundbased site to to ICRS astrometric RA,Dec.
The caller supplies UTC, site coordinates, ambient air conditions
and observing wavelength.

Given:

	typ    string    type of coordinates - "R", "H" or "A" (Notes 1,2)
	ob1    float64   observed Az, HA or RA (radians; Az is N=0,E=90)
	ob2    float64   observed ZD or Dec (radians)
	utc1   float64   UTC as a 2-part...
	utc2   float64   ...quasi Julian Date (Notes 3,4)
	dut1   float64   UT1-UTC (seconds, Note 5)
	elong  float64   longitude (radians, east +ve, Note 6)
	phi    float64   geodetic latitude (radians, Note 6)
	hm     float64   height above ellipsoid (m, geodetic Notes 6,8)
	xp,yp  float64   polar motion coordinates (radians, Note 7)
	phpa   float64   pressure at the observer (hPa = mB, Note 8)
	tc     float64   ambient temperature at the observer (deg C)
	rh     float64   relative humidity at the observer (range 0-1)
	wl     float64   wavelength (micrometers, Note 9)

Returned:

	rc,dc  float64   ICRS astrometric RA,Dec (radians)

Returned (function value):

	int    status: +1 = dubious year (Note 4)
	                0 = OK
	               -1 = unacceptable date

Notes:

 1. "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

 2. Only the first character of the type argument is significant.
    "R" or "r" indicates that ob1 and ob2 are the observed right
    ascension and declination;  "H" or "h" indicates that they are
    hour angle (west +ve) and declination;  anything else ("A" or
    "a" is recommended) indicates that ob1 and ob2 are azimuth
    (north zero, east 90 deg) and zenith distance.

 3. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 4. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See Dat for further details.

 5. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 6. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 7. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

 8. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

 9. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

    10)The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions Atco13 and
    Atoc13 are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

    11)It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

Called:

	Apco13    astrometry parameters, ICRS-observed
	Atoiq     quick observed to CIRS
	Aticq     quick CIRS to ICRS
*/
func Atoc13(typ string, ob1, ob2 float64, utc1, utc2, dut1 float64,
	elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64,
	rc, dc *float64) int {
	var j int
	var astrom ASTROM
	var eo, ri, di float64

	/* Star-independent astrometry parameters. */
	j = Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
		phpa, tc, rh, wl, &astrom, &eo)

	/* Abort if bad UTC. */
	if j < 0 {
		return j
	}

	/* Transform observed to CIRS. */
	Atoiq(typ, ob1, ob2, &astrom, &ri, &di)

	/* Transform CIRS to ICRS. */
	Aticq(ri, di, &astrom, rc, dc)

	/* Return OK/warning status. */
	return j
}

/*
Atoi13 Observed -> CIRS

Observed place to CIRS.  The caller supplies UTC, site coordinates,
ambient air conditions and observing wavelength.

Given:

	typ    string   type of coordinates - "R", "H" or "A" (Notes 1,2)
	ob1    float64   observed Az, HA or RA (radians; Az is N=0,E=90)
	ob2    float64   observed ZD or Dec (radians)
	utc1   float64   UTC as a 2-part...
	utc2   float64   ...quasi Julian Date (Notes 3,4)
	dut1   float64   UT1-UTC (seconds, Note 5)
	elong  float64   longitude (radians, east +ve, Note 6)
	phi    float64   geodetic latitude (radians, Note 6)
	hm     float64   height above the ellipsoid (meters, Notes 6,8)
	xp,yp  float64   polar motion coordinates (radians, Note 7)
	phpa   float64   pressure at the observer (hPa = mB, Note 8)
	tc     float64   ambient temperature at the observer (deg C)
	rh     float64   relative humidity at the observer (range 0-1)
	wl     float64   wavelength (micrometers, Note 9)

Returned:

	ri     float64  CIRS right ascension (CIO-based, radians)
	di     float64  CIRS declination (radians)

Returned (function value):

	int    status: +1 = dubious year (Note 2)
	                0 = OK
	               -1 = unacceptable date

Notes:

 1. "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

 2. Only the first character of the type argument is significant.
    "R" or "r" indicates that ob1 and ob2 are the observed right
    ascension and declination;  "H" or "h" indicates that they are
    hour angle (west +ve) and declination;  anything else ("A" or
    "a" is recommended) indicates that ob1 and ob2 are azimuth
    (north zero, east 90 deg) and zenith distance.

 3. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function Dtf2d to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

 4. The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See Dat for further details.

 5. UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

 6. The geographical coordinates are with respect to the WGS84
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

 7. The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

 8. If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    hm = -29.3 * tsl * log ( phpa / 1013.25 );

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

 9. The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

    10)The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions Atio13 and
    Atoi13 are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

    12)It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

Called:

	Apio13    astrometry parameters, CIRS-observed, 2013
	Atoiq     quick observed to CIRS
*/
func Atoi13(typ string, ob1, ob2 float64, utc1, utc2, dut1 float64,
	elong, phi, hm, xp, yp float64, phpa, tc, rh, wl float64,
	ri, di *float64) int {
	var j int
	var astrom ASTROM

	/* Star-independent astrometry parameters for CIRS->observed. */
	j = Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
		phpa, tc, rh, wl, &astrom)

	/* Abort if bad UTC. */
	if j < 0 {
		return j
	}

	/* Transform observed to CIRS. */
	Atoiq(typ, ob1, ob2, &astrom, ri, di)

	/* Return OK/warning status. */
	return j
}

/*
Atoiq Quick observed -> CIRS

Quick observed place to CIRS, given the star-independent astrometry
parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling Apio[13] or Apco[13].

Given:

	typ    string     type of coordinates: "R", "H" or "A" (Note 1)
	ob1    float64     observed Az, HA or RA (radians; Az is N=0,E=90)
	ob2    float64     observed ZD or Dec (radians)
	astrom ASTROM        star-independent astrometry parameters:
	 pmt    float64       PM time interval (SSB, Julian years)
	 eb     [3]float64    SSB to observer (vector, au)
	 eh     [3]float64    Sun to observer (unit vector)
	 em     float64       distance from Sun to observer (au)
	 v      [3]float64    barycentric observer velocity (vector, c)
	 bm1    float64       sqrt(1-|v|^2): reciprocal of Lorenz factor
	 bpn    [3][3]float64 bias-precession-nutation matrix
	 along  float64       longitude + s' (radians)
	 xpl    float64       polar motion xp wrt local meridian (radians)
	 ypl    float64       polar motion yp wrt local meridian (radians)
	 sphi   float64       sine of geodetic latitude
	 cphi   float64       cosine of geodetic latitude
	 diurab float64       magnitude of diurnal aberration vector
	 eral   float64       "local" Earth rotation angle (radians)
	 refa   float64       refraction constant A (radians)
	 refb   float64       refraction constant B (radians)

Returned:

	ri     float64    CIRS right ascension (CIO-based, radians)
	di     float64    CIRS declination (radians)

Notes:

 1. "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  This is related to
    the observed HA,Dec via the standard rotation, using the geodetic
    latitude (corrected for polar motion), while the observed HA and
    RA are related simply through the Earth rotation angle and the
    site longitude.  "Observed" RA,Dec or HA,Dec thus means the
    position that would be seen by a perfect equatorial with its
    polar axis aligned to the Earth's axis of rotation.  By removing
    from the observed place the effects of atmospheric refraction and
    diurnal aberration, the CIRS RA,Dec is obtained.

 2. Only the first character of the type argument is significant.
    "R" or "r" indicates that ob1 and ob2 are the observed right
    ascension and declination;  "H" or "h" indicates that they are
    hour angle (west +ve) and declination;  anything else ("A" or
    "a" is recommended) indicates that ob1 and ob2 are azimuth (north
    zero, east 90 deg) and zenith distance.  (Zenith distance is used
    rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)

 3. The accuracy of the result is limited by the corrections for
    refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted intermediate
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better than
    20 arcmin (optical) or 25 arcmin (radio) at the horizon.

    Without refraction, the complementary functions Atioq and
    Atoiq are self-consistent to better than 1 microarcsecond all
    over the celestial sphere.  With refraction included, consistency
    falls off at high zenith distances, but is still better than
    0.05 arcsec at 85 degrees.

 4. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

Called:

	S2c       spherical coordinates to unit vector
	C2s       p-vector to spherical
	Anp       normalize angle into range 0 to 2pi
*/
func Atoiq(typ string, ob1, ob2 float64, astrom *ASTROM, ri, di *float64) {
	/* Minimum sin(alt) for refraction purposes */
	const SELMIN = 0.05

	var c byte
	var c1, c2, sphi, cphi, ce, xaeo, yaeo, zaeo,
		xmhdo, ymhdo, zmhdo, az, sz, zdo, refa, refb, tz, dref,
		zdt, xaet, yaet, zaet, xmhda, ymhda, zmhda,
		f, xhd, yhd, zhd, sx, cx, sy, cy, hma float64
	var v [3]float64

	/* Coordinate type. */
	c = typ[0]

	/* Coordinates. */
	c1 = ob1
	c2 = ob2

	/* Sin, cos of latitude. */
	sphi = astrom.Sphi
	cphi = astrom.Cphi

	/* Standardize coordinate type. */
	if c == 'r' || c == 'R' {
		c = 'R'
	} else if c == 'h' || c == 'H' {
		c = 'H'
	} else {
		c = 'A'
	}

	/* If Az,ZD, convert to Cartesian (S=0,E=90). */
	if c == 'A' {
		ce = sin(c2)
		xaeo = -cos(c1) * ce
		yaeo = sin(c1) * ce
		zaeo = cos(c2)

	} else {

		/* If RA,Dec, convert to HA,Dec. */
		if c == 'R' {
			c1 = astrom.Eral - c1
		}

		/* To Cartesian -HA,Dec. */
		S2c(-c1, c2, &v)
		xmhdo = v[0]
		ymhdo = v[1]
		zmhdo = v[2]

		/* To Cartesian Az,El (S=0,E=90). */
		xaeo = sphi*xmhdo - cphi*zmhdo
		yaeo = ymhdo
		zaeo = cphi*xmhdo + sphi*zmhdo
	}

	/* Azimuth (S=0,E=90). */
	// az = ( xaeo != 0.0 || yaeo != 0.0 ) ? atan2(yaeo,xaeo) : 0.0;
	if xaeo != 0.0 || yaeo != 0.0 {
		az = atan2(yaeo, xaeo)
	} else {
		az = 0.0
	}

	/* Sine of observed ZD, and observed ZD. */
	sz = sqrt(xaeo*xaeo + yaeo*yaeo)
	zdo = atan2(sz, zaeo)

	/*
	   Refraction
	   ----------
	*/

	/* Fast algorithm using two constant model. */
	refa = astrom.Refa
	refb = astrom.Refb
	// tz = sz / ( zaeo > SELMIN ? zaeo : SELMIN );
	if zaeo > SELMIN {
		tz = sz / zaeo
	} else {
		tz = sz / SELMIN
	}
	dref = (refa + refb*tz*tz) * tz
	zdt = zdo + dref

	/* To Cartesian Az,ZD. */
	ce = sin(zdt)
	xaet = cos(az) * ce
	yaet = sin(az) * ce
	zaet = cos(zdt)

	/* Cartesian Az,ZD to Cartesian -HA,Dec. */
	xmhda = sphi*xaet + cphi*zaet
	ymhda = yaet
	zmhda = -cphi*xaet + sphi*zaet

	/* Diurnal aberration. */
	f = (1.0 + astrom.Diurab*ymhda)
	xhd = f * xmhda
	yhd = f * (ymhda - astrom.Diurab)
	zhd = f * zmhda

	/* Polar motion. */
	sx = sin(astrom.Xpl)
	cx = cos(astrom.Xpl)
	sy = sin(astrom.Ypl)
	cy = cos(astrom.Ypl)
	v[0] = cx*xhd + sx*sy*yhd - sx*cy*zhd
	v[1] = cy*yhd + sy*zhd
	v[2] = sx*xhd - cx*sy*yhd + cx*cy*zhd

	/* To spherical -HA,Dec. */
	C2s(v, &hma, di)

	/* Right ascension. */
	*ri = Anp(astrom.Eral + hma)
}

/*
Ld Light deflection by a single solar-system body

Apply light deflection by a solar-system body, as part of
transforming coordinate direction into natural direction.

Given:

	bm     float64     mass of the gravitating body (solar masses)
	p      [3]float64  direction from observer to source (unit vector)
	q      [3]float64  direction from body to source (unit vector)
	e      [3]float64  direction from body to observer (unit vector)
	em     float64     distance from body to observer (au)
	dlim   float64     deflection limiter (Note 4)

Returned:

	p1     [3]float64  observer to deflected source (unit vector)

Notes:

 1. The algorithm is based on Expr. (70) in Klioner (2003) and
    Expr. (7.63) in the Explanatory Supplement (Urban & Seidelmann
    2013), with some rearrangement to minimize the effects of machine
    precision.

 2. The mass parameter bm can, as required, be adjusted in order to
    allow for such effects as quadrupole field.

 3. The barycentric position of the deflecting body should ideally
    correspond to the time of closest approach of the light ray to
    the body.

 4. The deflection limiter parameter dlim is phi^2/2, where phi is
    the angular separation (in radians) between source and body at
    which limiting is applied.  As phi shrinks below the chosen
    threshold, the deflection is artificially reduced, reaching zero
    for phi = 0.

 5. The returned vector p1 is not normalized, but the consequential
    departure from unit magnitude is always negligible.

 6. The arguments p and p1 can be the same array.

 7. To accumulate total light deflection taking into account the
    contributions from several bodies, call the present function for
    each body in succession, in decreasing order of distance from the
    observer.

 8. For efficiency, validation is omitted.  The supplied vectors must
    be of unit magnitude, and the deflection limiter non-zero and
    positive.

References:

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013).

	Klioner, Sergei A., "A practical relativistic model for micro-
	arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

Called:

	Pdp       scalar product of two p-vectors
	Pxp       vector product of two p-vectors
*/
func Ld(bm float64, p, q, e [3]float64, em, dlim float64, p1 *[3]float64) {

	var i int
	var qpe, eq, peq [3]float64
	var qdqpe, w float64

	/* q . (q + e). */
	for i = 0; i < 3; i++ {
		qpe[i] = q[i] + e[i]
	}
	qdqpe = Pdp(q, qpe)

	/* 2 x G x bm / ( em x c^2 x ( q . (q + e) ) ). */
	w = bm * SRS / em / gmax(qdqpe, dlim)

	/* p x (e x q). */
	Pxp(e, q, &eq)
	Pxp(p, eq, &peq)

	/* Apply the deflection. */
	for i = 0; i < 3; i++ {
		p1[i] = p[i] + w*peq[i]
	}
}

/*
Ldn Light deflection by multiple solar-system bodies

For a star, apply light deflection by multiple solar-system bodies,
as part of transforming coordinate direction into natural direction.

Given:

	n    int            number of bodies (note 1)
	b    LDBODY[n]      data for each of the n bodies (Notes 1,2):
	 bm   float64         mass of the body (solar masses, Note 3)
	 dl   float64         deflection limiter (Note 4)
	 pv   [2][3]float64   barycentric PV of the body (au, au/day)
	ob   [3]float64     barycentric position of the observer (au)
	sc   [3]float64     observer to star coord direction (unit vector)

Returned:

	sn    [3]float64      observer to deflected star (unit vector)

Notes:

 1. The array b contains n entries, one for each body to be
    considered.  If n = 0, no gravitational light deflection will be
    applied, not even for the Sun.

 2. The array b should include an entry for the Sun as well as for
    any planet or other body to be taken into account.  The entries
    should be in the order in which the light passes the body.

 3. In the entry in the b array for body i, the mass parameter
    b[i].bm can, as required, be adjusted in order to allow for such
    effects as quadrupole field.

 4. The deflection limiter parameter b[i].dl is phi^2/2, where phi is
    the angular separation (in radians) between star and body at
    which limiting is applied.  As phi shrinks below the chosen
    threshold, the deflection is artificially reduced, reaching zero
    for phi = 0.   Example values suitable for a terrestrial
    observer, together with masses, are as follows:

    body i     b[i].bm        b[i].dl

    Sun        1.0            6e-6
    Jupiter    0.00095435     3e-9
    Saturn     0.00028574     3e-10

 5. For cases where the starlight passes the body before reaching the
    observer, the body is placed back along its barycentric track by
    the light time from that point to the observer.  For cases where
    the body is "behind" the observer no such shift is applied.  If
    a different treatment is preferred, the user has the option of
    instead using the Ld function.  Similarly, Ld can be used
    for cases where the source is nearby, not a star.

 6. The returned vector sn is not normalized, but the consequential
    departure from unit magnitude is always negligible.

 7. The arguments sc and sn can be the same array.

 8. For efficiency, validation is omitted.  The supplied masses must
    be greater than zero, the position and velocity vectors must be
    right, and the deflection limiter greater than zero.

Reference:

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013), Section 7.2.4.

Called:

	Cp        copy p-vector
	Pdp       scalar product of two p-vectors
	Pmp       p-vector minus p-vector
	Ppsp      p-vector plus scaled p-vector
	Pn        decompose p-vector into modulus and direction
	Ld        light deflection by a solar-system body
*/
func Ldn(n int, b []LDBODY, ob, sc [3]float64, sn *[3]float64) {
	/* Light time for 1 au (days) */
	const CR = AULT / DAYSEC

	var i int
	var v, ev, e [3]float64
	var dt, em float64

	/* Star direction prior to deflection. */
	Cp(sc, sn)

	/* Body by body. */
	for i = 0; i < n; i++ {

		/* Body to observer vector at epoch of observation (au). */
		Pmp(ob, b[i].Pv[0], &v)

		/* Minus the time since the light passed the body (days). */
		dt = Pdp(*sn, v) * CR

		/* Neutralize if the star is "behind" the observer. */
		dt = gmin(dt, 0.0)

		/* Backtrack the body to the time the light was passing the body. */
		Ppsp(v, -dt, b[i].Pv[1], &ev)

		/* Body to observer vector as magnitude and direction. */
		Pn(ev, &em, &e)

		/* Apply light deflection for this body. */
		Ld(b[i].Bm, *sn, *sn, e, em, b[i].Dl, sn)

		/* Next body. */
	}
}

/*
Ldsun Light deflection by the Sun

Deflection of starlight by the Sun.

Given:

	p      [3]float64  direction from observer to star (unit vector)
	e      [3]float64  direction from Sun to observer (unit vector)
	em     float64     distance from Sun to observer (au)

Returned:

	p1     [3]float64  observer to deflected star (unit vector)

Notes:

 1. The source is presumed to be sufficiently distant that its
    directions seen from the Sun and the observer are essentially
    the same.

 2. The deflection is restrained when the angle between the star and
    the center of the Sun is less than a threshold value, falling to
    zero deflection for zero separation.  The chosen threshold value
    is within the solar limb for all solar-system applications, and
    is about 5 arcminutes for the case of a terrestrial observer.

 3. The arguments p and p1 can be the same array.

Called:

	Ld        light deflection by a solar-system body
*/
func Ldsun(p, e [3]float64, em float64, p1 *[3]float64) {
	var em2, dlim float64

	/* Deflection limiter (smaller for distant observers). */
	em2 = em * em
	if em2 < 1.0 {
		em2 = 1.0
	}
	// dlim = 1e-6 / (em2 > 1.0 ? em2 : 1.0);
	if em2 > 1.0 {
		dlim = 1e-6 / em2
	} else {
		dlim = 1e-6
	}

	/* Apply the deflection. */
	Ld(1.0, p, p, e, em, dlim, p1)
}

/*
Pmpx Apply proper motion and parallax

Given:

	rc,dc  float64     ICRS RA,Dec at catalog epoch (radians)
	pr     float64     RA proper motion (radians/year, Note 1)
	pd     float64     Dec proper motion (radians/year)
	px     float64     parallax (arcsec)
	rv     float64     radial velocity (km/s, +ve if receding)
	pmt    float64     proper motion time interval (SSB, Julian years)
	pob    [3]float64  SSB to observer vector (au)

Returned:

	pco    [3]float64  coordinate direction (BCRS unit vector)

Notes:

 1. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 2. The proper motion time interval is for when the starlight
    reaches the solar system barycenter.

 3. To avoid the need for iteration, the Roemer effect (i.e. the
    small annual modulation of the proper motion coming from the
    changing light time) is applied approximately, using the
    direction of the star at the catalog epoch.

References:

	1984 Astronomical Almanac, pp B39-B41.

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013), Section 7.2.

Called:

	Pdp       scalar product of two p-vectors
	Pn        decompose p-vector into modulus and direction
*/
func Pmpx(rc, dc float64, pr, pd, px, rv, pmt float64, pob [3]float64, pco *[3]float64) {
	/* Km/s to au/year */
	const VF = DAYSEC * DJM / DAU

	/* Light time for 1 au, Julian years */
	const AULTY = AULT / DAYSEC / DJY

	var i int
	var sr, cr, sd, cd, x, y, z, dt, pxr, w, pdz float64
	var p, pm [3]float64

	/* Spherical coordinates to unit vector (and useful functions). */
	sr = sin(rc)
	cr = cos(rc)
	sd = sin(dc)
	cd = cos(dc)
	x = cr * cd
	p[0] = x
	y = sr * cd
	p[1] = y
	z = sd
	p[2] = z

	/* Proper motion time interval (y) including Roemer effect. */
	dt = pmt + Pdp(p, pob)*AULTY

	/* Space motion (radians per year). */
	pxr = px * DAS2R
	w = VF * rv * pxr
	pdz = pd * z
	pm[0] = -pr*y - pdz*cr + w*x
	pm[1] = pr*x - pdz*sr + w*y
	pm[2] = pd*cd + w*z

	/* Coordinate direction of star (unit vector, BCRS). */
	for i = 0; i < 3; i++ {
		p[i] += dt*pm[i] - pxr*pob[i]
	}
	Pn(p, &w, pco)
}

/*
Pmsafe Apply proper motion, with zero-parallax precautions

Star proper motion:  update star catalog data for space motion, with
special handling to handle the zero parallax case.

Given:

	ra1    float64      right ascension (radians), before
	dec1   float64      declination (radians), before
	pmr1   float64      RA proper motion (radians/year), before
	pmd1   float64      Dec proper motion (radians/year), before
	px1    float64      parallax (arcseconds), before
	rv1    float64      radial velocity (km/s, +ve = receding), before
	ep1a   float64      "before" epoch, part A (Note 1)
	ep1b   float64      "before" epoch, part B (Note 1)
	ep2a   float64      "after" epoch, part A (Note 1)
	ep2b   float64      "after" epoch, part B (Note 1)

Returned:

	ra2    float64      right ascension (radians), after
	dec2   float64      declination (radians), after
	pmr2   float64      RA proper motion (radians/year), after
	pmd2   float64      Dec proper motion (radians/year), after
	px2    float64      parallax (arcseconds), after
	rv2    float64      radial velocity (km/s, +ve = receding), after

Returned (function value):

	int    status:
	         -1 = system error (should not occur)
	          0 = no warnings or errors
	          1 = distance overridden (Note 6)
	          2 = excessive velocity (Note 7)
	          4 = solution didn't converge (Note 8)
	       else = binary logical OR of the above warnings

Notes:

 1. The starting and ending TDB epochs ep1a+ep1b and ep2a+ep2b are
    Julian Dates, apportioned in any convenient way between the two
    parts (A and B).  For example, JD(TDB)=2450123.7 could be
    expressed in any of these ways, among others:

    epNa            epNb

    2450123.7           0.0       (JD method)
    2451545.0       -1421.3       (J2000 method)
    2400000.5       50123.2       (MJD method)
    2450123.5           0.2       (date & time method)

    The JD method is the most natural and convenient to use in cases
    where the loss of several decimal digits of resolution is
    acceptable.  The J2000 method is best matched to the way the
    argument is handled internally and will deliver the optimum
    resolution.  The MJD method and the date & time methods are both
    good compromises between resolution and convenience.

 2. In accordance with normal star-catalog conventions, the object's
    right ascension and declination are freed from the effects of
    secular aberration.  The frame, which is aligned to the catalog
    equator and equinox, is Lorentzian and centered on the SSB.

    The proper motions are the rate of change of the right ascension
    and declination at the catalog epoch and are in radians per TDB
    Julian year.

    The parallax and radial velocity are in the same frame.

 3. Care is needed with units.  The star coordinates are in radians
    and the proper motions in radians per Julian year, but the
    parallax is in arcseconds.

 4. The RA proper motion is in terms of coordinate angle, not true
    angle.  If the catalog uses arcseconds for both RA and Dec proper
    motions, the RA proper motion will need to be divided by cos(Dec)
    before use.

 5. Straight-line motion at constant speed, in the inertial frame, is
    assumed.

 6. An extremely small (or zero or negative) parallax is overridden
    to ensure that the object is at a finite but very large distance,
    but not so large that the proper motion is equivalent to a large
    but safe speed (about 0.1c using the chosen constant).  A warning
    status of 1 is added to the status if this action has been taken.

 7. If the space velocity is a significant fraction of c (see the
    constant VMAX in the function Starpv), it is arbitrarily set
    to zero.  When this action occurs, 2 is added to the status.

 8. The relativistic adjustment carried out in the Starpv function
    involves an iterative calculation.  If the process fails to
    converge within a set number of iterations, 4 is added to the
    status.

Called:

	Seps      angle between two points
	Starpm    update star catalog data for space motion
*/
func Pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1 float64,
	ep1a, ep1b, ep2a, ep2b float64,
	ra2, dec2, pmr2, pmd2, px2, rv2 *float64) int {

	/* Minimum allowed parallax (arcsec) */
	const PXMIN = 5e-7

	/* Factor giving maximum allowed transverse speed of about 1% c */
	const F = 326.0

	var jpx, j int
	var pm, px1a float64

	/* Proper motion in one year (radians). */
	pm = Seps(ra1, dec1, ra1+pmr1, dec1+pmd1)

	/* Override the parallax to reduce the chances of a warning status. */
	jpx = 0
	px1a = px1
	pm *= F
	if px1a < pm {
		jpx = 1
		px1a = pm
	}
	if px1a < PXMIN {
		jpx = 1
		px1a = PXMIN
	}

	/* Carry out the transformation using the modified parallax. */
	j = Starpm(ra1, dec1, pmr1, pmd1, px1a, rv1,
		ep1a, ep1b, ep2a, ep2b,
		ra2, dec2, pmr2, pmd2, px2, rv2)

	/* Revise and return the status. */
	if (j % 2) != 1 {
		j += jpx
	}
	return j
}

/*
Pvstar Star position+velocity vector to catalog coordinates

Convert star position+velocity vector to catalog coordinates.

Given (Note 1):

	pv     [2][3]float64   pv-vector (au, au/day)

Returned (Note 2):

	ra     float64         right ascension (radians)
	dec    float64         declination (radians)
	pmr    float64         RA proper motion (radians/year)
	pmd    float64         Dec proper motion (radians/year)
	px     float64         parallax (arcsec)
	rv     float64         radial velocity (km/s, positive = receding)

Returned (function value):

	int     status:
	          0 = OK
	         -1 = superluminal speed (Note 5)
	         -2 = null position vector

Notes:

 1. The specified pv-vector is the coordinate direction (and its rate
    of change) for the date at which the light leaving the star
    reached the solar-system barycenter.

 2. The star data returned by this function are "observables" for an
    imaginary observer at the solar-system barycenter.  Proper motion
    and radial velocity are, strictly, in terms of barycentric
    coordinate time, TCB.  For most practical applications, it is
    permissible to neglect the distinction between TCB and ordinary
    "proper" time on Earth (TT/TAI).  The result will, as a rule, be
    limited by the intrinsic accuracy of the proper-motion and
    radial-velocity data;  moreover, the supplied pv-vector is likely
    to be merely an intermediate result (for example generated by the
    function Starpv), so that a change of time unit will cancel
    out overall.

    In accordance with normal star-catalog conventions, the object's
    right ascension and declination are freed from the effects of
    secular aberration.  The frame, which is aligned to the catalog
    equator and equinox, is Lorentzian and centered on the SSB.

    Summarizing, the specified pv-vector is for most stars almost
    identical to the result of applying the standard geometrical
    "space motion" transformation to the catalog data.  The
    differences, which are the subject of the Stumpff paper cited
    below, are:

    (i) In stars with significant radial velocity and proper motion,
    the constantly changing light-time distorts the apparent proper
    motion.  Note that this is a classical, not a relativistic,
    effect.

    (ii) The transformation complies with special relativity.

 3. Care is needed with units.  The star coordinates are in radians
    and the proper motions in radians per Julian year, but the
    parallax is in arcseconds; the radial velocity is in km/s, but
    the pv-vector result is in au and au/day.

 4. The proper motions are the rate of change of the right ascension
    and declination at the catalog epoch and are in radians per Julian
    year.  The RA proper motion is in terms of coordinate angle, not
    true angle, and will thus be numerically larger at high
    declinations.

 5. Straight-line motion at constant speed in the inertial frame is
    assumed.  If the speed is greater than or equal to the speed of
    light, the function aborts with an error status.

 6. The inverse transformation is performed by the function Starpv.

Called:

	Pn        decompose p-vector into modulus and direction
	Pdp       scalar product of two p-vectors
	Sxp       multiply p-vector by scalar
	Pmp       p-vector minus p-vector
	Pm        modulus of p-vector
	Ppp       p-vector plus p-vector
	Pv2s      pv-vector to spherical
	Anp       normalize angle into range 0 to 2pi

Reference:

	Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
*/
func Pvstar(pv [2][3]float64, ra, dec, pmr, pmd, px, rv *float64) int {

	var r, vr, vt, bett, betr, d, w, del, a, rad, decd, rd float64

	var pu, ur, ut, usr, ust [3]float64

	/* Isolate the radial component of the velocity (au/day, inertial). */
	Pn(pv[0], &r, &pu)
	vr = Pdp(pu, pv[1])
	Sxp(vr, pu, &ur)

	/* Isolate the transverse component of the velocity (au/day, inertial). */
	Pmp(pv[1], ur, &ut)
	vt = Pm(ut)

	/* Special-relativity dimensionless parameters. */
	bett = vt / DC
	betr = vr / DC

	/* The observed-to-inertial correction terms. */
	d = 1.0 + betr
	w = betr*betr + bett*bett
	if d == 0.0 || w > 1.0 {
		return -1
	}
	del = -w / (sqrt(1.0-w) + 1.0)

	/* Scale inertial tangential velocity vector into observed (au/d). */
	Sxp(1.0/d, ut, &ust)

	/* Compute observed radial velocity vector (au/d). */
	Sxp(DC*(betr-del)/d, pu, &usr)

	/* Combine the two to obtain the observed velocity vector (au/day). */
	Ppp(usr, ust, &pv[1])

	/* Cartesian to spherical. */
	Pv2s(pv, &a, dec, &r, &rad, &decd, &rd)
	if r == 0.0 {
		return -2
	}

	/* Return RA in range 0 to 2pi. */
	*ra = Anp(a)

	/* Return proper motions in radians per year. */
	*pmr = rad * DJY
	*pmd = decd * DJY

	/* Return parallax in arcsec. */
	*px = DR2AS / r

	/* Return radial velocity in km/s. */
	*rv = 1e-3 * rd * DAU / DAYSEC

	/* Success. */
	return 0
}

/*
Pvtob Observatory position and velocity

Position and velocity of a terrestrial observing station.

Given:

	elong   float64       longitude (radians, east +ve, Note 1)
	phi     float64       latitude (geodetic, radians, Note 1)
	hm      float64       height above ref. ellipsoid (geodetic, m)
	xp,yp   float64       coordinates of the pole (radians, Note 2)
	sp      float64       the TIO locator s' (radians, Note 2)
	theta   float64       Earth rotation angle (radians, Note 3)

Returned:

	pv      [2][3]float64 position/velocity vector (m, m/s, CIRS)

Notes:

 1. The terrestrial coordinates are with respect to the WGS84
    reference ellipsoid.

 2. xp and yp are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions), measured along the
    meridians 0 and 90 deg west respectively.  sp is the TIO locator
    s', in radians, which positions the Terrestrial Intermediate
    Origin on the equator.  For many applications, xp, yp and
    (especially) sp can be set to zero.

 3. If theta is Greenwich apparent sidereal time instead of Earth
    rotation angle, the result is with respect to the true equator
    and equinox of date, i.e. with the x-axis at the equinox rather
    than the celestial intermediate origin.

 4. The velocity units are meters per UT1 second, not per SI second.
    This is unlikely to have any practical consequences in the modern
    era.

 5. No validation is performed on the arguments.  Error cases that
    could lead to arithmetic exceptions are trapped by the Gd2gc
    function, and the result set to zeros.

References:

	McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
	IERS Technical Note No. 32, BKG (2004)

	Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
	the Astronomical Almanac, 3rd ed., University Science Books
	(2013), Section 7.4.3.3.

Called:

	Gd2gc     geodetic to geocentric transformation
	Pom00     polar motion matrix
	Trxp      product of transpose of r-matrix and p-vector
*/
func Pvtob(elong, phi, hm float64, xp, yp, sp, theta float64, pv *[2][3]float64) {
	/* Earth rotation rate in radians per UT1 second */
	const OM = 1.00273781191135448 * D2PI / DAYSEC

	var xyz, xyzm [3]float64
	var rpm [3][3]float64
	var x, y, z, s, c float64

	/* Geodetic to geocentric transformation (WGS84). */
	Gd2gc(1, elong, phi, hm, &xyzm)

	/* Polar motion and TIO position. */
	Pom00(xp, yp, sp, &rpm)
	Trxp(rpm, xyzm, &xyz)
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]

	/* Functions of ERA. */
	s = sin(theta)
	c = cos(theta)

	/* Position. */
	pv[0][0] = c*x - s*y
	pv[0][1] = s*x + c*y
	pv[0][2] = z

	/* Velocity. */
	pv[1][0] = OM * (-s*x - c*y)
	pv[1][1] = OM * (c*x - s*y)
	pv[1][2] = 0.0
}

/*
Refco Refraction constants

Determine the constants A and B in the atmospheric refraction model
dZ = A tan Z + B tan^3 Z.

Z is the "observed" zenith distance (i.e. affected by refraction)
and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
zenith distance.

Given:

	phpa   float64    pressure at the observer (hPa = millibar)
	tc     float64    ambient temperature at the observer (deg C)
	rh     float64    relative humidity at the observer (range 0-1)
	wl     float64    wavelength (micrometers)

Returned:

	refa   float64   tan Z coefficient (radians)
	refb   float64   tan^3 Z coefficient (radians)

Notes:

 1. The model balances speed and accuracy to give good results in
    applications where performance at low altitudes is not paramount.
    Performance is maintained across a range of conditions, and
    applies to both optical/IR and radio.

 2. The model omits the effects of (i) height above sea level (apart
    from the reduced pressure itself), (ii) latitude (i.e. the
    flattening of the Earth), (iii) variations in tropospheric lapse
    rate and (iv) dispersive effects in the radio.

    The model was tested using the following range of conditions:

    lapse rates 0.0055, 0.0065, 0.0075 deg/meter
    latitudes 0, 25, 50, 75 degrees
    heights 0, 2500, 5000 meters ASL
    pressures mean for height -10% to +5% in steps of 5%
    temperatures -10 deg to +20 deg with respect to 280 deg at SL
    relative humidity 0, 0.5, 1
    wavelengths 0.4, 0.6, ... 2 micron, + radio
    zenith distances 15, 45, 75 degrees

    The accuracy with respect to raytracing through a model
    atmosphere was as follows:

    worst         RMS

    optical/IR           62 mas       8 mas
    radio               319 mas      49 mas

    For this particular set of conditions:

    lapse rate 0.0065 K/meter
    latitude 50 degrees
    sea level
    pressure 1005 mb
    temperature 280.15 K
    humidity 80%
    wavelength 5740 Angstroms

    the results were as follows:

    ZD       raytrace       Refco   Saastamoinen

    10         10.27        10.27        10.27
    20         21.19        21.20        21.19
    30         33.61        33.61        33.60
    40         48.82        48.83        48.81
    45         58.16        58.18        58.16
    50         69.28        69.30        69.27
    55         82.97        82.99        82.95
    60        100.51       100.54       100.50
    65        124.23       124.26       124.20
    70        158.63       158.68       158.61
    72        177.32       177.37       177.31
    74        200.35       200.38       200.32
    76        229.45       229.43       229.42
    78        267.44       267.29       267.41
    80        319.13       318.55       319.10

    deg        arcsec       arcsec       arcsec

    The values for Saastamoinen's formula (which includes terms
    up to tan^5) are taken from Hohenkerk and Sinclair (1985).

 3. A wl value in the range 0-100 selects the optical/IR case and is
    wavelength in micrometers.  Any value outside this range selects
    the radio case.

 4. Outlandish input parameters are silently limited to
    mathematically safe values.  Zero pressure is permissible, and
    causes zeroes to be returned.

 5. The algorithm draws on several sources, as follows:

    a) The formula for the saturation vapour pressure of water as
    a function of temperature and temperature is taken from
    Equations (A4.5-A4.7) of Gill (1982).

    b) The formula for the water vapour pressure, given the
    saturation pressure and the relative humidity, is from
    Crane (1976), Equation (2.5.5).

    c) The refractivity of air is a function of temperature,
    total pressure, water-vapour pressure and, in the case
    of optical/IR, wavelength.  The formulae for the two cases are
    developed from Hohenkerk & Sinclair (1985) and Rueger (2002).
    The IAG (1999) optical refractivity for dry air is used.

    d) The formula for beta, the ratio of the scale height of the
    atmosphere to the geocentric distance of the observer, is
    an adaption of Equation (9) from Stone (1996).  The
    adaptations, arrived at empirically, consist of (i) a small
    adjustment to the coefficient and (ii) a humidity term for the
    radio case only.

    e) The formulae for the refraction constants as a function of
    n-1 and beta are from Green (1987), Equation (4.31).

References:

	Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral
	Atmosphere", Methods of Experimental Physics: Astrophysics 12B,
	Academic Press, 1976.

	Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press,
	1982.

	Green, R.M., "Spherical Astronomy", Cambridge University Press,
	1987.

	Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63,
	1985.

	IAG Resolutions adopted at the XXIIth General Assembly in
	Birmingham, 1999, Resolution 3.

	Rueger, J.M., "Refractive Index Formulae for Electronic Distance
	Measurement with Radio and Millimetre Waves", in Unisurv Report
	S-68, School of Surveying and Spatial Information Systems,
	University of New South Wales, Sydney, Australia, 2002.

	Stone, Ronald C., P.A.S.P. 108, 1051-1058, 1996.
*/
func Refco(phpa, tc, rh, wl float64, refa, refb *float64) {
	var optic bool
	var p, t, r, w, ps, pw, tk, wlsq, gamma, beta float64

	/* Decide whether optical/IR or radio case:  switch at 100 microns. */
	optic = (wl <= 100.0)

	/* Restrict parameters to safe values. */
	t = gmax(tc, -150.0)
	t = gmin(t, 200.0)
	p = gmax(phpa, 0.0)
	p = gmin(p, 10000.0)
	r = gmax(rh, 0.0)
	r = gmin(r, 1.0)
	w = gmax(wl, 0.1)
	w = gmin(w, 1e6)

	/* Water vapour pressure at the observer. */
	if p > 0.0 {
		ps = pow(10.0, (0.7859+0.03477*t)/
			(1.0+0.00412*t)) *
			(1.0 + p*(4.5e-6+6e-10*t*t))
		pw = r * ps / (1.0 - (1.0-r)*ps/p)
	} else {
		pw = 0.0
	}

	/* Refractive index minus 1 at the observer. */
	tk = t + 273.15
	if optic {
		wlsq = w * w
		gamma = ((77.53484e-6+
			(4.39108e-7+3.666e-9/wlsq)/wlsq)*p -
			11.2684e-6*pw) / tk
	} else {
		gamma = (77.6890e-6*p - (6.3938e-6-0.375463/tk)*pw) / tk
	}

	/* Formula for beta from Stone, with empirical adjustments. */
	beta = 4.4474e-6 * tk
	if !optic {
		beta -= 0.0074 * pw * beta
	}

	/* Refraction constants from Green. */
	*refa = gamma * (1.0 - beta)
	*refb = -gamma * (beta - gamma/2.0)
}

/*
Starpm Proper motion between two epochs

Star proper motion:  update star catalog data for space motion.

Given:

	ra1    float64     right ascension (radians), before
	dec1   float64     declination (radians), before
	pmr1   float64     RA proper motion (radians/year), before
	pmd1   float64     Dec proper motion (radians/year), before
	px1    float64     parallax (arcseconds), before
	rv1    float64     radial velocity (km/s, +ve = receding), before
	ep1a   float64     "before" epoch, part A (Note 1)
	ep1b   float64     "before" epoch, part B (Note 1)
	ep2a   float64     "after" epoch, part A (Note 1)
	ep2b   float64     "after" epoch, part B (Note 1)

Returned:

	ra2    float64     right ascension (radians), after
	dec2   float64     declination (radians), after
	pmr2   float64     RA proper motion (radians/year), after
	pmd2   float64     Dec proper motion (radians/year), after
	px2    float64     parallax (arcseconds), after
	rv2    float64     radial velocity (km/s, +ve = receding), after

Returned (function value):

	int     status:
	          -1 = system error (should not occur)
	           0 = no warnings or errors
	           1 = distance overridden (Note 6)
	           2 = excessive velocity (Note 7)
	           4 = solution didn't converge (Note 8)
	        else = binary logical OR of the above warnings

Notes:

 1. The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
    Julian Dates, apportioned in any convenient way between the two
    parts (A and B).  For example, JD(TDB)=2450123.7 could be
    expressed in any of these ways, among others:

    epna          epnb

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

 2. In accordance with normal star-catalog conventions, the object's
    right ascension and declination are freed from the effects of
    secular aberration.  The frame, which is aligned to the catalog
    equator and equinox, is Lorentzian and centered on the SSB.

    The proper motions are the rate of change of the right ascension
    and declination at the catalog epoch and are in radians per TDB
    Julian year.

    The parallax and radial velocity are in the same frame.

 3. Care is needed with units.  The star coordinates are in radians
    and the proper motions in radians per Julian year, but the
    parallax is in arcseconds.

 4. The RA proper motion is in terms of coordinate angle, not true
    angle.  If the catalog uses arcseconds for both RA and Dec proper
    motions, the RA proper motion will need to be divided by cos(Dec)
    before use.

 5. Straight-line motion at constant speed, in the inertial frame,
    is assumed.

 6. An extremely small (or zero or negative) parallax is interpreted
    to mean that the object is on the "celestial sphere", the radius
    of which is an arbitrary (large) value (see the Starpv
    function for the value used).  When the distance is overridden in
    this way, the status, initially zero, has 1 added to it.

 7. If the space velocity is a significant fraction of c (see the
    constant VMAX in the function Starpv), it is arbitrarily set
    to zero.  When this action occurs, 2 is added to the status.

 8. The relativistic adjustment carried out in the Starpv function
    involves an iterative calculation.  If the process fails to
    converge within a set number of iterations, 4 is added to the
    status.

Called:

	Starpv    star catalog data to space motion pv-vector
	Pvu       update a pv-vector
	Pdp       scalar product of two p-vectors
	Pvstar    space motion pv-vector to star catalog data
*/
func Starpm(ra1, dec1, pmr1, pmd1, px1, rv1 float64, ep1a, ep1b, ep2a, ep2b float64,
	ra2, dec2, pmr2, pmd2, px2, rv2 *float64) int {

	var pv1, pv, pv2 [2][3]float64
	var tl1, dt, r2, rdv, v2, c2mv2, tl2 float64

	var j1, j2 int

	/* RA,Dec etc. at the "before" epoch to space motion pv-vector. */
	j1 = Starpv(ra1, dec1, pmr1, pmd1, px1, rv1, &pv1)

	/* Light time when observed (days). */
	tl1 = Pm(pv1[0]) / DC

	/* Time interval, "before" to "after" (days). */
	dt = (ep2a - ep1a) + (ep2b - ep1b)

	/* Move star along track from the "before" observed position to the */
	/* "after" geometric position. */
	Pvu(dt+tl1, pv1, &pv)

	/* From this geometric position, deduce the observed light time (days) */
	/* at the "after" epoch (with theoretically unneccessary error check). */
	r2 = Pdp(pv[0], pv[0])
	rdv = Pdp(pv[0], pv[1])
	v2 = Pdp(pv[1], pv[1])
	c2mv2 = DC*DC - v2
	if c2mv2 <= 0 {
		return -1
	}
	tl2 = (-rdv + sqrt(rdv*rdv+c2mv2*r2)) / c2mv2

	/* Move the position along track from the observed place at the */
	/* "before" epoch to the observed place at the "after" epoch. */
	Pvu(dt+(tl1-tl2), pv1, &pv2)

	/* Space motion pv-vector to RA,Dec etc. at the "after" epoch. */
	j2 = Pvstar(pv2, ra2, dec2, pmr2, pmd2, px2, rv2)

	/* Final status. */

	if j2 == 0 {
		return j1
	}
	return -1
}

/*
Starpv Star catalog coordinates to position+velocity vector

Convert star catalog coordinates to position+velocity vector.

Given (Note 1):

	ra     float64        right ascension (radians)
	dec    float64        declination (radians)
	pmr    float64        RA proper motion (radians/year)
	pmd    float64        Dec proper motion (radians/year)
	px     float64        parallax (arcseconds)
	rv     float64        radial velocity (km/s, positive = receding)

Returned (Note 2):

	pv     [2][3]float64  pv-vector (au, au/day)

Returned (function value):

	int     status:
	           0 = no warnings
	           1 = distance overridden (Note 6)
	           2 = excessive speed (Note 7)
	           4 = solution didn't converge (Note 8)
	        else = binary logical OR of the above

Notes:

 1. The star data accepted by this function are "observables" for an
    imaginary observer at the solar-system barycenter.  Proper motion
    and radial velocity are, strictly, in terms of barycentric
    coordinate time, TCB.  For most practical applications, it is
    permissible to neglect the distinction between TCB and ordinary
    "proper" time on Earth (TT/TAI).  The result will, as a rule, be
    limited by the intrinsic accuracy of the proper-motion and
    radial-velocity data;  moreover, the pv-vector is likely to be
    merely an intermediate result, so that a change of time unit
    would cancel out overall.

    In accordance with normal star-catalog conventions, the object's
    right ascension and declination are freed from the effects of
    secular aberration.  The frame, which is aligned to the catalog
    equator and equinox, is Lorentzian and centered on the SSB.

 2. The resulting position and velocity pv-vector is with respect to
    the same frame and, like the catalog coordinates, is freed from
    the effects of secular aberration.  Should the "coordinate
    direction", where the object was located at the catalog epoch, be
    required, it may be obtained by calculating the magnitude of the
    position vector pv[0][0-2] dividing by the speed of light in
    au/day to give the light-time, and then multiplying the space
    velocity pv[1][0-2] by this light-time and adding the result to
    pv[0][0-2].

    Summarizing, the pv-vector returned is for most stars almost
    identical to the result of applying the standard geometrical
    "space motion" transformation.  The differences, which are the
    subject of the Stumpff paper referenced below, are:

    (i) In stars with significant radial velocity and proper motion,
    the constantly changing light-time distorts the apparent proper
    motion.  Note that this is a classical, not a relativistic,
    effect.

    (ii) The transformation complies with special relativity.

 3. Care is needed with units.  The star coordinates are in radians
    and the proper motions in radians per Julian year, but the
    parallax is in arcseconds; the radial velocity is in km/s, but
    the pv-vector result is in au and au/day.

 4. The RA proper motion is in terms of coordinate angle, not true
    angle.  If the catalog uses arcseconds for both RA and Dec proper
    motions, the RA proper motion will need to be divided by cos(Dec)
    before use.

 5. Straight-line motion at constant speed, in the inertial frame,
    is assumed.

 6. An extremely small (or zero or negative) parallax is interpreted
    to mean that the object is on the "celestial sphere", the radius
    of which is an arbitrary (large) value (see the constant PXMIN).
    When the distance is overridden in this way, the status,
    initially zero, has 1 added to it.

 7. If the space velocity is a significant fraction of c (see the
    constant VMAX), it is arbitrarily set to zero.  When this action
    occurs, 2 is added to the status.

 8. The relativistic adjustment involves an iterative calculation.
    If the process fails to converge within a set number (IMAX) of
    iterations, 4 is added to the status.

 9. The inverse transformation is performed by the function
    Pvstar.

Called:

	S2pv      spherical coordinates to pv-vector
	Pm        modulus of p-vector
	Zp        zero p-vector
	Pn        decompose p-vector into modulus and direction
	Pdp       scalar product of two p-vectors
	Sxp       multiply p-vector by scalar
	Pmp       p-vector minus p-vector
	Ppp       p-vector plus p-vector
*/
func Starpv(ra, dec, pmr, pmd, px, rv float64, pv *[2][3]float64) int {
	/* Smallest allowed parallax */
	const PXMIN = 1e-7

	/* Largest allowed speed (fraction of c) */
	const VMAX = 0.5

	/* Maximum number of iterations for relativistic solution */
	const IMAX = 100

	var i, iwarn int
	var w, r, rd, rad, decd, v,
		vsr, vst, betst, betsr, bett, betr,
		dd, ddel,
		d, del,
		odd, oddel,
		od, odel float64
	d, del = 0.0, 0.0     /* to prevent */
	odd, oddel = 0.0, 0.0 /* compiler   */
	od, odel = 0.0, 0.0   /* warnings   */
	var pu, usr, ust, ur, ut [3]float64

	/* Distance (au). */
	if px >= PXMIN {
		w = px
		iwarn = 0
	} else {
		w = PXMIN
		iwarn = 1
	}
	r = DR2AS / w

	/* Radial speed (au/day). */
	rd = DAYSEC * rv * 1e3 / DAU

	/* Proper motion (radian/day). */
	rad = pmr / DJY
	decd = pmd / DJY

	/* To pv-vector (au,au/day). */
	S2pv(ra, dec, r, rad, decd, rd, pv)

	/* If excessive velocity, arbitrarily set it to zero. */
	v = Pm(pv[1])
	if (v / DC) > VMAX {
		Zp(&pv[1])
		iwarn += 2
	}

	/* Isolate the radial component of the velocity (au/day). */
	Pn(pv[0], &w, &pu)
	vsr = Pdp(pu, pv[1])
	Sxp(vsr, pu, &usr)

	/* Isolate the transverse component of the velocity (au/day). */
	Pmp(pv[1], usr, &ust)
	vst = Pm(ust)

	/* Special-relativity dimensionless parameters. */
	betsr = vsr / DC
	betst = vst / DC

	/* Determine the observed-to-inertial relativistic correction terms. */
	bett = betst
	betr = betsr
	for i = 0; i < IMAX; i++ {
		d = 1.0 + betr
		w = betr*betr + bett*bett
		del = -w / (sqrt(1.0-w) + 1.0)
		betr = d*betsr + del
		bett = d * betst
		if i > 0 {
			dd = fabs(d - od)
			ddel = fabs(del - odel)
			if (i > 1) && (dd >= odd) && (ddel >= oddel) {
				break
			}
			odd = dd
			oddel = ddel
		}
		od = d
		odel = del
	}
	if i >= IMAX {
		iwarn += 4
	}

	// Scale observed tangential velocity vector into inertial (au/d). */
	Sxp(d, ust, &ut)

	/* Compute inertial radial velocity vector (au/d). */
	Sxp(DC*(d*betsr+del), pu, &ur)

	/* Combine the two to obtain the inertial space velocity vector. */
	Ppp(ur, ut, &pv[1])

	/* Return the status. */
	return iwarn
}
