// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Coord

// Galactic Coordinates

/*
Icrs2g Transformation from ICRS to Galactic Coordinates.

Given:
    dr     float64      ICRS right ascension (radians)
    dd     float64      ICRS declination (radians)

Returned:
    dl     float64      galactic longitude (radians)
    db     float64      galactic latitude (radians)

Notes:

 1) The IAU 1958 system of Galactic coordinates was defined with
    respect to the now obsolete reference system FK4 B1950.0.  When
    interpreting the system in a modern context, several factors have
    to be taken into account:

     . The inclusion in FK4 positions of the E-terms of aberration.
     . The distortion of the FK4 proper motion system by differential
       Galactic rotation.
     . The use of the B1950.0 equinox rather than the now-standard
       J2000.0.
     . The frame bias between ICRS and the J2000.0 mean place system.

    The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
    matrix that transforms directly between ICRS and Galactic
    coordinates with the above factors taken into account.  The
    matrix is derived from three angles, namely the ICRS coordinates
    of the Galactic pole and the longitude of the ascending node of
    the galactic equator on the ICRS equator.  They are given in
    degrees to five decimal places and for canonical purposes are
    regarded as exact.  In the Hipparcos Catalogue the matrix
    elements are given to 10 decimal places (about 20 microarcsec).
    In the present SOFA function the matrix elements have been
    recomputed from the canonical three angles and are given to 30
    decimal places.

 2) The inverse transformation is performed by the function G2icrs.

Called:
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi
    S2c       spherical coordinates to unit vector
    Rxp       product of r-matrix and p-vector
    C2s       p-vector to spherical

Reference:
    Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
    catalogues.  Astrometric and photometric star catalogues
    derived from the ESA Hipparcos Space Astrometry Mission.  ESA
    Publications Division, Noordwijk, Netherlands.
*/
func Icrs2g(dr, dd float64, dl, db *float64) {
	var v1, v2 [3]float64

	/*
	   L2,B2 system of galactic coordinates in the form presented in the
	   Hipparcos Catalogue.  In degrees:

	   P = 192.85948    right ascension of the Galactic north pole in ICRS
	   Q =  27.12825    declination of the Galactic north pole in ICRS
	   R =  32.93192    Galactic longitude of the ascending node of
	                    the Galactic equator on the ICRS equator

	   ICRS to galactic rotation matrix, obtained by computing
	   R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
	*/
	r := [3][3]float64{
		{-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116},
		{+0.494109427875583673525222371358, -0.444829629960011178146614061616, +0.746982244497218890527388004556},
		{-0.867666149019004701181616534570, -0.198076373431201528180486091412, +0.455983776175066922272100478348},
	}

	/* Spherical to Cartesian. */
	S2c(dr, dd, &v1)

	/* ICRS to Galactic. */
	Rxp(r, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, dl, db)

	/* Express in conventional ranges. */
	*dl = Anp(*dl)
	*db = Anpm(*db)
}

/*
G2icrsTransformation from Galactic Coordinates to ICRS.

Given:
    dl     float64      galactic longitude (radians)
    db     float64      galactic latitude (radians)

Returned:
    dr     float64      ICRS right ascension (radians)
    dd     float64      ICRS declination (radians)

Notes:

 1) The IAU 1958 system of Galactic coordinates was defined with
    respect to the now obsolete reference system FK4 B1950.0.  When
    interpreting the system in a modern context, several factors have
    to be taken into account:

     . The inclusion in FK4 positions of the E-terms of aberration.
     . The distortion of the FK4 proper motion system by differential
       Galactic rotation.
     . The use of the B1950.0 equinox rather than the now-standard
       J2000.0.
     . The frame bias between ICRS and the J2000.0 mean place system.

    The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
    matrix that transforms directly between ICRS and Galactic
    coordinates with the above factors taken into account.  The
    matrix is derived from three angles, namely the ICRS coordinates
    of the Galactic pole and the longitude of the ascending node of
    the galactic equator on the ICRS equator.  They are given in
    degrees to five decimal places and for canonical purposes are
    regarded as exact.  In the Hipparcos Catalogue the matrix
    elements are given to 10 decimal places (about 20 microarcsec).
    In the present SOFA function the matrix elements have been
    recomputed from the canonical three angles and are given to 30
    decimal places.

 2) The inverse transformation is performed by the function Icrs2g.

Called:
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi
    S2c       spherical coordinates to unit vector
    Trxp      product of transpose of r-matrix and p-vector
    C2s       p-vector to spherical

Reference:
    Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
    catalogues.  Astrometric and photometric star catalogues
    derived from the ESA Hipparcos Space Astrometry Mission.  ESA
    Publications Division, Noordwijk, Netherlands.
*/
func G2icrs(dl, db float64, dr, dd *float64) {
	var v1, v2 [3]float64

	/*
	   L2,B2 system of galactic coordinates in the form presented in the
	   Hipparcos Catalogue.  In degrees:

	   P = 192.85948    right ascension of the Galactic north pole in ICRS
	   Q =  27.12825    declination of the Galactic north pole in ICRS
	   R =  32.93192    Galactic longitude of the ascending node of
	                    the Galactic equator on the ICRS equator

	   ICRS to galactic rotation matrix, obtained by computing
	   R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
	*/
	r := [3][3]float64{
		{-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116},
		{+0.494109427875583673525222371358, -0.444829629960011178146614061616, +0.746982244497218890527388004556},
		{-0.867666149019004701181616534570, -0.198076373431201528180486091412, +0.455983776175066922272100478348},
	}

	/* Spherical to Cartesian. */
	S2c(dl, db, &v1)

	/* Galactic to ICRS. */
	Trxp(r, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, dr, dd)

	/* Express in conventional ranges. */
	*dr = Anp(*dr)
	*dd = Anpm(*dd)
}

/*
Ae2hd Horizon to equatorial coordinates, transform azimuth and altitude to hour angle and declination.

Given:
    az       float64       azimuth
    el       float64       altitude (informally, elevation)
    phi      float64       site latitude

Returned:
    ha       float64       hour angle (local)
    dec      float64       declination

Notes:

 1) All the arguments are angles in radians.

 2) The sign convention for azimuth is north zero, east +pi/2.

 3) HA is returned in the range +/-pi.  Declination is returned in
    the range +/-pi/2.

 4) The latitude phi is pi/2 minus the angle between the Earth's
    rotation axis and the adopted zenith.  In many applications it
    will be sufficient to use the published geodetic latitude of the
    site.  In very precise (sub-arcsecond) applications, phi can be
    corrected for polar motion.

 5) The azimuth az must be with respect to the rotational north pole,
    as opposed to the ITRS pole, and an azimuth with respect to north
    on a map of the Earth's surface will need to be adjusted for
    polar motion if sub-arcsecond accuracy is required.

 6) Should the user wish to work with respect to the astronomical
    zenith rather than the geodetic zenith, phi will need to be
    adjusted for deflection of the vertical (often tens of
    arcseconds), and the zero point of ha will also be affected.

 7) The transformation is the same as Ve = Ry(phi-pi/2)*Rz(pi)*Vh,
    where Ve and Vh are lefthanded unit vectors in the (ha,dec) and
    (az,el) systems respectively and Rz and Ry are rotations about
    first the z-axis and then the y-axis.  (n.b. Rz(pi) simply
    reverses the signs of the x and y components.)  For efficiency,
    the algorithm is written out rather than calling other utility
    functions.  For applications that require even greater
    efficiency, additional savings are possible if constant terms
    such as functions of latitude are computed once and for all.

 8) Again for efficiency, no range checking of arguments is carried
    out.
*/
func Ae2hd(az, el, phi float64, ha, dec *float64) {
	var sa, ca, se, ce, sp, cp, x, y, z, r float64

	/* Useful trig functions. */
	sa = sin(az)
	ca = cos(az)
	se = sin(el)
	ce = cos(el)
	sp = sin(phi)
	cp = cos(phi)

	/* HA,Dec unit vector. */
	x = -ca*ce*sp + se*cp
	y = -sa * ce
	z = ca*ce*cp + se*sp

	/* To spherical. */
	r = sqrt(x*x + y*y)
	if r != 0.0 {
		*ha = atan2(y, x)
	} else {
		*ha = 0.0
	}
	*dec = atan2(z, r)
}

/*
Hd2ae Equatorial to horizon coordinates: transform hour angle and declination to azimuth and altitude.

Given:
    ha       float64       hour angle (local)
    dec      float64       declination
    phi      float64       site latitude

Returned:
    az      float64       azimuth
    el      float64       altitude (informally, elevation)

Notes:

 1) All the arguments are angles in radians.

 2) Azimuth is returned in the range 0-2pi;  north is zero, and east
    is +pi/2.  Altitude is returned in the range +/- pi/2.

 3) The latitude phi is pi/2 minus the angle between the Earth's
    rotation axis and the adopted zenith.  In many applications it
    will be sufficient to use the published geodetic latitude of the
    site.  In very precise (sub-arcsecond) applications, phi can be
    corrected for polar motion.

 4) The returned azimuth az is with respect to the rotational north
    pole, as opposed to the ITRS pole, and for sub-arcsecond
    accuracy will need to be adjusted for polar motion if it is to
    be with respect to north on a map of the Earth's surface.

 5) Should the user wish to work with respect to the astronomical
    zenith rather than the geodetic zenith, phi will need to be
    adjusted for deflection of the vertical (often tens of
    arcseconds), and the zero point of the hour angle ha will also
    be affected.

 6) The transformation is the same as Vh = Rz(pi)*Ry(pi/2-phi)*Ve,
    where Vh and Ve are lefthanded unit vectors in the (az,el) and
    (ha,dec) systems respectively and Ry and Rz are rotations about
    first the y-axis and then the z-axis.  (n.b. Rz(pi) simply
    reverses the signs of the x and y components.)  For efficiency,
    the algorithm is written out rather than calling other utility
    functions.  For applications that require even greater
    efficiency, additional savings are possible if constant terms
    such as functions of latitude are computed once and for all.

 7) Again for efficiency, no range checking of arguments is carried
    out.
*/
func Hd2ae(ha, dec, phi float64, az, el *float64) {
	var sh, ch, sd, cd, sp, cp, x, y, z, r, a float64

	/* Useful trig functions. */
	sh = sin(ha)
	ch = cos(ha)
	sd = sin(dec)
	cd = cos(dec)
	sp = sin(phi)
	cp = cos(phi)

	/* Az,Alt unit vector. */
	x = -ch*cd*sp + sd*cp
	y = -sh * cd
	z = ch*cd*cp + sd*sp

	/* To spherical. */
	r = sqrt(x*x + y*y)
	if r != 0.0 {
		a = atan2(y, x)
	} else {
		a = 0.0
	}
	if a < 0.0 {
		*az = a + D2PI
	} else {
		*az = a
	}

	*el = atan2(z, r)
}

/*
Hd2pa Parallactic angle for a given hour angle and declination.

Given:
    ha     float64     hour angle
    dec    float64     declination
    phi    float64     site latitude

Returned (function value):
    float64     parallactic angle

Notes:

 1) All the arguments are angles in radians.

 2) The parallactic angle at a point in the sky is the position
    angle of the vertical, i.e. the angle between the directions to
    the north celestial pole and to the zenith respectively.

 3) The result is returned in the range -pi to +pi.

 4) At the pole itself a zero result is returned.

 5) The latitude phi is pi/2 minus the angle between the Earth's
    rotation axis and the adopted zenith.  In many applications it
    will be sufficient to use the published geodetic latitude of the
    site.  In very precise (sub-arcsecond) applications, phi can be
    corrected for polar motion.

 6) Should the user wish to work with respect to the astronomical
    zenith rather than the geodetic zenith, phi will need to be
    adjusted for deflection of the vertical (often tens of
    arcseconds), and the zero point of the hour angle ha will also
    be affected.

Reference:
    Smart, W.M., "Spherical Astronomy", Cambridge University Press,
    6th edition (Green, 1977), p49.
*/
func Hd2pa(ha, dec, phi float64) float64 {
	var cp, cqsz, sqsz float64

	cp = cos(phi)
	sqsz = cp * sin(ha)
	cqsz = sin(phi)*cos(dec) - cp*sin(dec)*cos(ha)
	if sqsz != 0.0 || cqsz != 0.0 {
		return atan2(sqsz, cqsz)
	}
	return 0.0
}

/*
Ecm06 ICRS equatorial to ecliptic rotation matrix, IAU 2006.

Given:
    date1,date2  float64         TT as a 2-part Julian date (Note 1)

Returned:
    rm           [3][3]float64   ICRS to ecliptic rotation matrix

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

 2) The matrix is in the sense

    	E_ep = rm x P_ICRS,

    where P_ICRS is a vector with respect to ICRS right ascension
    and declination axes and E_ep is the same vector with respect to
    the (inertial) ecliptic and equinox of date.

 3) P_ICRS is a free vector, merely a direction, typically of unit
    magnitude, and not bound to any particular spatial origin, such
    as the Earth, Sun or SSB.  No assumptions are made about whether
    it represents starlight and embodies astrometric effects such as
    parallax or aberration.  The transformation is approximately that
    between mean J2000.0 right ascension and declination and ecliptic
    longitude and latitude, with only frame bias (always less than
    25 mas) to disturb this classical picture.

Called:
    Obl06     mean obliquity, IAU 2006
    Pmat06    PB matrix, IAU 2006
    Ir        initialize r-matrix to identity
    Rx        rotate around X-axis
    Rxr       product of two r-matrices
*/
func Ecm06(date1, date2 float64, rm *[3][3]float64) {
	var ob float64
	var bp, e [3][3]float64

	/* Obliquity, IAU 2006. */
	ob = Obl06(date1, date2)

	/* Precession-bias matrix, IAU 2006. */
	Pmat06(date1, date2, &bp)

	/* Equatorial of date to ecliptic matrix. */
	Ir(&e)
	Rx(ob, &e)

	/* ICRS to ecliptic coordinates rotation matrix, IAU 2006. */
	Rxr(e, bp, rm)

}

/*
Eqec06 Transformation from ICRS equatorial coordinates to ecliptic coordinates
(mean equinox and ecliptic of date) using IAU 2006 precession model.

Given:
    date1,date2 float64 TT as a 2-part Julian date (Note 1)
    dr,dd       float64 ICRS right ascension and declination (radians)

Returned:
    dl,db       float64 ecliptic longitude and latitude (radians)

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

 2) No assumptions are made about whether the coordinates represent
    starlight and embody astrometric effects such as parallax or
    aberration.

 3) The transformation is approximately that from mean J2000.0 right
    ascension and declination to ecliptic longitude and latitude
    (mean equinox and ecliptic of date), with only frame bias (always
    less than 25 mas) to disturb this classical picture.

Called:
    S2c       spherical coordinates to unit vector
    Ecm06     J2000.0 to ecliptic rotation matrix, IAU 2006
    Rxp       product of r-matrix and p-vector
    C2s       unit vector to spherical coordinates
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi
*/
func Eqec06(date1, date2 float64, dr, dd float64, dl, db *float64) {
	var rm [3][3]float64
	var v1, v2 [3]float64
	var a, b float64

	/* Spherical to Cartesian. */
	S2c(dr, dd, &v1)

	/* Rotation matrix, ICRS equatorial to ecliptic. */
	Ecm06(date1, date2, &rm)

	/* The transformation from ICRS to ecliptic. */
	Rxp(rm, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, &a, &b)

	/* Express in conventional ranges. */
	*dl = Anp(a)
	*db = Anpm(b)
}

/*
Eceq06 Transformation from ecliptic coordinates (mean equinox and ecliptic
of date) to ICRS RA,Dec, using the IAU 2006 precession model.

Given:
    date1,date2 float64 TT as a 2-part Julian date (Note 1)
    dl,db       float64 ecliptic longitude and latitude (radians)

Returned:
    dr,dd       float64 ICRS right ascension and declination (radians)

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

 2) No assumptions are made about whether the coordinates represent
    starlight and embody astrometric effects such as parallax or
    aberration.

 3) The transformation is approximately that from ecliptic longitude
    and latitude (mean equinox and ecliptic of date) to mean J2000.0
    right ascension and declination, with only frame bias (always
    less than 25 mas) to disturb this classical picture.

Called:
    S2c       spherical coordinates to unit vector
    Ecm06     J2000.0 to ecliptic rotation matrix, IAU 2006
    Trxp      product of transpose of r-matrix and p-vector
    C2s       unit vector to spherical coordinates
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi
*/
func Eceq06(date1, date2 float64, dl, db float64, dr, dd *float64) {
	var rm [3][3]float64
	var v1, v2 [3]float64
	var a, b float64

	/* Spherical to Cartesian. */
	S2c(dl, db, &v1)

	/* Rotation matrix, ICRS equatorial to ecliptic. */
	Ecm06(date1, date2, &rm)

	/* The transformation from ecliptic to ICRS. */
	Trxp(rm, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, &a, &b)

	/* Express in conventional ranges. */
	*dr = Anp(a)
	*dd = Anpm(b)
}

/*
Ltecm ICRS equatorial to ecliptic rotation matrix, long-term.

Given:
    epj     float64         Julian epoch (TT)

Returned:
    rm      [3][3]float64   ICRS to ecliptic rotation matrix

Notes:

 1) The matrix is in the sense

    	E_ep = rm x P_ICRS,

    where P_ICRS is a vector with respect to ICRS right ascension
    and declination axes and E_ep is the same vector with respect to
    the (inertial) ecliptic and equinox of epoch epj.

 2) P_ICRS is a free vector, merely a direction, typically of unit
    magnitude, and not bound to any particular spatial origin, such
    as the Earth, Sun or SSB.  No assumptions are made about whether
    it represents starlight and embodies astrometric effects such as
    parallax or aberration.  The transformation is approximately that
    between mean J2000.0 right ascension and declination and ecliptic
    longitude and latitude, with only frame bias (always less than
    25 mas) to disturb this classical picture.

 3) The Vondrak et al. (2011, 2012) 400 millennia precession model
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
func Ltecm(epj float64, rm *[3][3]float64) {
	/* Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33) */
	dx := -0.016617 * DAS2R
	de := -0.0068192 * DAS2R
	dr := -0.0146 * DAS2R

	var p, z, w, x, y [3]float64
	var s float64

	/* Equator pole. */
	Ltpequ(epj, &p)

	/* Ecliptic pole (bottom row of equatorial to ecliptic matrix). */
	Ltpecl(epj, &z)

	/* Equinox (top row of matrix). */
	Pxp(p, z, &w)
	Pn(w, &s, &x)

	/* Middle row of matrix. */
	Pxp(z, x, &y)

	/* Combine with frame bias. */
	rm[0][0] = x[0] - x[1]*dr + x[2]*dx
	rm[0][1] = x[0]*dr + x[1] + x[2]*de
	rm[0][2] = -x[0]*dx - x[1]*de + x[2]
	rm[1][0] = y[0] - y[1]*dr + y[2]*dx
	rm[1][1] = y[0]*dr + y[1] + y[2]*de
	rm[1][2] = -y[0]*dx - y[1]*de + y[2]
	rm[2][0] = z[0] - z[1]*dr + z[2]*dx
	rm[2][1] = z[0]*dr + z[1] + z[2]*de
	rm[2][2] = -z[0]*dx - z[1]*de + z[2]
}

/*
Lteqec Transformation from ICRS equatorial coordinates to ecliptic coordinates
(mean equinox and ecliptic of date) using a long-term precession model.

Given:
    epj     float64     Julian epoch (TT)
    dr,dd   float64     ICRS right ascension and declination (radians)

Returned:
    dl,db   float64     ecliptic longitude and latitude (radians)

Notes:
 1) No assumptions are made about whether the coordinates represent
    starlight and embody astrometric effects such as parallax or
    aberration.

 2) The transformation is approximately that from mean J2000.0 right
    ascension and declination to ecliptic longitude and latitude
    (mean equinox and ecliptic of date), with only frame bias (always
    less than 25 mas) to disturb this classical picture.

 3) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

Called:
    S2c       spherical coordinates to unit vector
    Ltecm     J2000.0 to ecliptic rotation matrix, long term
    Rxp       product of r-matrix and p-vector
    C2s       unit vector to spherical coordinates
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Lteqec(epj float64, dr, dd float64, dl, db *float64) {
	var rm [3][3]float64
	var v1, v2 [3]float64
	var a, b float64

	/* Spherical to Cartesian. */
	S2c(dr, dd, &v1)

	/* Rotation matrix, ICRS equatorial to ecliptic. */
	Ltecm(epj, &rm)

	/* The transformation from ICRS to ecliptic. */
	Rxp(rm, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, &a, &b)

	/* Express in conventional ranges. */
	*dl = Anp(a)
	*db = Anpm(b)
}

/*
Lteceq Transformation from ecliptic coordinates (mean equinox and ecliptic
of date) to ICRS RA,Dec, using a long-term precession model.

Given:
    epj     float64     Julian epoch (TT)
    dl,db   float64     ecliptic longitude and latitude (radians)

Returned:
    dr,dd   float64     ICRS right ascension and declination (radians)

Notes:
 1) No assumptions are made about whether the coordinates represent
    starlight and embody astrometric effects such as parallax or
    aberration.

 2) The transformation is approximately that from ecliptic longitude
    and latitude (mean equinox and ecliptic of date) to mean J2000.0
    right ascension and declination, with only frame bias (always
    less than 25 mas) to disturb this classical picture.

 3) The Vondrak et al. (2011, 2012) 400 millennia precession model
    agrees with the IAU 2006 precession at J2000.0 and stays within
    100 microarcseconds during the 20th and 21st centuries.  It is
    accurate to a few arcseconds throughout the historical period,
    worsening to a few tenths of a degree at the end of the
    +/- 200,000 year time span.

Called:
    S2c       spherical coordinates to unit vector
    Ltecm     J2000.0 to ecliptic rotation matrix, long term
    Trxp      product of transpose of r-matrix and p-vector
    C2s       unit vector to spherical coordinates
    Anp       normalize angle into range 0 to 2pi
    Anpm      normalize angle into range +/- pi

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1
*/
func Lteceq(epj float64, dl, db float64, dr, dd *float64) {
	var rm [3][3]float64
	var v1, v2 [3]float64
	var a, b float64

	/* Spherical to Cartesian. */
	S2c(dl, db, &v1)

	/* Rotation matrix, ICRS equatorial to ecliptic. */
	Ltecm(epj, &rm)

	/* The transformation from ecliptic to ICRS. */
	Trxp(rm, v1, &v2)

	/* Cartesian to spherical. */
	C2s(v2, &a, &b)

	/* Express in conventional ranges. */
	*dr = Anp(a)
	*dd = Anpm(b)
}

/*
Eform a,f for a nominated Earth reference ellipsoid

Given:
    n    int         ellipsoid identifier (Note 1)

Returned:
    a    float64      equatorial radius (meters, Note 2)
    f    float64      flattening (Note 2)

Returned (function value):
    int   status:  0 = OK
                  -1 = illegal identifier (Note 3)

Notes:

 1) The identifier n is a number that specifies the choice of
    reference ellipsoid.  The following are supported:

    n    ellipsoid

    1     WGS84
    2     GRS80
    3     WGS72

    The n value has no significance outside the SOFA software.  For
    convenience, symbols WGS84 etc. are defined in sofam.h.

 2) The ellipsoid parameters are returned in the form of equatorial
    radius in meters (a) and flattening (f).  The latter is a number
    around 0.00335, i.e. around 1/298.

 3) For the case where an unsupported n value is supplied, zero a and
    f are returned, as well as error status.

References:

    Department of Defense World Geodetic System 1984, National
    Imagery and Mapping Agency Technical Report 8350.2, Third
    Edition, p3-2.

    Moritz, H., Bull. Geodesique 66-2, 187 (1992).

    The Department of Defense World Geodetic System 1972, World
    Geodetic System Committee, May 1974.

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    p220.
*/
func Eform(n int, a, f *float64) int {
	/* Look up a and f for the specified reference ellipsoid. */
	switch n {

	case WGS84:
		*a = 6378137.0
		*f = 1.0 / 298.257223563
		break

	case GRS80:
		*a = 6378137.0
		*f = 1.0 / 298.257222101
		break

	case WGS72:
		*a = 6378135.0
		*f = 1.0 / 298.26
		break

	default:

		/* Invalid identifier. */
		*a = 0.0
		*f = 0.0
		return -1

	}

	/* OK status. */
	return 0
}

/*
Gc2gd Transform geocentric coordinates to geodetic using the specified
reference ellipsoid.

Given:
    n       int         ellipsoid identifier (Note 1)
    xyz     [3]float64  geocentric vector (Note 2)

Returned:
    elong   float64     longitude (radians, east +ve, Note 3)
    phi     float64     latitude (geodetic, radians, Note 3)
    height  float64     height above ellipsoid (geodetic, Notes 2,3)

Returned (function value):
    int     status:  0 = OK
                    -1 = illegal identifier (Note 3)
                    -2 = internal error (Note 3)

Notes:

 1) The identifier n is a number that specifies the choice of
    reference ellipsoid.  The following are supported:

        n     ellipsoid

        1     WGS84
        2     GRS80
        3     WGS72

    The n value has no significance outside the SOFA software.  For
    convenience, symbols WGS84 etc. are defined in sofam.h.

 2) The geocentric vector (xyz, given) and height (height, returned)
    are in meters.

 3) An error status -1 means that the identifier n is illegal.  An
    error status -2 is theoretically impossible.  In all error cases,
    all three results are set to -1e9.

 4) The inverse transformation is performed in the function iauGd2gc.

Called:
    Eform     Earth reference ellipsoids
    Gc2gde    geocentric to geodetic transformation, general
*/
func Gc2gd(n int, xyz [3]float64, elong, phi, height *float64) int {
	var j int
	var a, f float64

	/* Obtain reference ellipsoid parameters. */
	j = Eform(n, &a, &f)

	/* If OK, transform x,y,z to longitude, geodetic latitude, height. */
	if j == 0 {
		j = Gc2gde(a, f, xyz, elong, phi, height)
		if j < 0 {
			j = -2
		}
	}

	/* Deal with any errors. */
	if j < 0 {
		*elong = -1e9
		*phi = -1e9
		*height = -1e9
	}

	/* Return the status. */
	return j
}

/*
Gc2gde Transform geocentric coordinates to geodetic for a reference
ellipsoid of specified form.

Given:
    a       float64     equatorial radius (Notes 2,4)
    f       float64     flattening (Note 3)
    xyz     [3]float64  geocentric vector (Note 4)

Returned:
    elong   float64     longitude (radians, east +ve)
    phi     float64     latitude (geodetic, radians)
    height  float64     height above ellipsoid (geodetic, Note 4)

Returned (function value):
    int     status:  0 = OK
                    -1 = illegal f
                    -2 = illegal a

Notes:

 1) This function is based on the GCONV2H Fortran subroutine by
    Toshio Fukushima (see reference).

 2) The equatorial radius, a, can be in any units, but meters is
    the conventional choice.

 3) The flattening, f, is (for the Earth) a value around 0.00335,
    i.e. around 1/298.

 4) The equatorial radius, a, and the geocentric vector, xyz,
    must be given in the same units, and determine the units of
    the returned height, height.

 5) If an error occurs (status < 0), elong, phi and height are
    unchanged.

 6) The inverse transformation is performed in the function
    iauGd2gce.

 7) The transformation for a standard ellipsoid (such as WGS84) can
    more conveniently be performed by calling iauGc2gd, which uses a
    numerical code to identify the required A and F values.

Reference:

    Fukushima, T., "Transformation from Cartesian to geodetic
    coordinates accelerated by Halley's method", J.Geodesy (2006)
    79: 689-693
*/
func Gc2gde(a, f float64, xyz [3]float64, elong, phi, height *float64) int {
	var aeps2, e2, e4t, ec2, ec, b, x, y, z, p2, absz, p, s0, pn, zc,
		c0, c02, c03, s02, s03, a02, a0, a03, d0, f0, b0, s1,
		cc, s12, cc2 float64

	/* ------------- */
	/* Preliminaries */
	/* ------------- */

	/* Validate ellipsoid parameters. */
	if f < 0.0 || f >= 1.0 {
		return -1
	}
	if a <= 0.0 {
		return -2
	}

	/* Functions of ellipsoid parameters (with further validation of f). */
	aeps2 = a * a * 1e-32
	e2 = (2.0 - f) * f
	e4t = e2 * e2 * 1.5
	ec2 = 1.0 - e2
	if ec2 <= 0.0 {
		return -1
	}
	ec = sqrt(ec2)
	b = a * ec

	/* Cartesian components. */
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]

	/* Distance from polar axis squared. */
	p2 = x*x + y*y

	/* Longitude. */
	// *elong = p2 > 0.0 ? atan2(y, x) : 0.0;
	if p2 > 0.0 {
		*elong = atan2(y, x)
	} else {
		*elong = 0.0
	}

	/* Unsigned z-coordinate. */
	absz = fabs(z)

	/* Proceed unless polar case. */
	if p2 > aeps2 {

		/* Distance from polar axis. */
		p = sqrt(p2)

		/* Normalization. */
		s0 = absz / a
		pn = p / a
		zc = ec * s0

		/* Prepare Newton correction factors. */
		c0 = ec * pn
		c02 = c0 * c0
		c03 = c02 * c0
		s02 = s0 * s0
		s03 = s02 * s0
		a02 = c02 + s02
		a0 = sqrt(a02)
		a03 = a02 * a0
		d0 = zc*a03 + e2*s03
		f0 = pn*a03 - e2*c03

		/* Prepare Halley correction factor. */
		b0 = e4t * s02 * c02 * pn * (a0 - ec)
		s1 = d0*f0 - b0*s0
		cc = ec * (f0*f0 - b0*c0)

		/* Evaluate latitude and height. */
		*phi = atan(s1 / cc)
		s12 = s1 * s1
		cc2 = cc * cc
		*height = (p*cc + absz*s1 - a*sqrt(ec2*s12+cc2)) / sqrt(s12+cc2)
	} else {

		/* Exception: pole. */
		*phi = DPI / 2.0
		*height = absz - b
	}

	/* Restore sign of latitude. */
	if z < 0 {
		*phi = -*phi
	}

	/* OK status. */
	return 0
}

/*
Gd2gc Transform geodetic coordinates to geocentric using the specified
reference ellipsoid.

Given:
    n       int         ellipsoid identifier (Note 1)
    elong   float64     longitude (radians, east +ve)
    phi     float64     latitude (geodetic, radians, Note 3)
    height  float64     height above ellipsoid (geodetic, Notes 2,3)

Returned:
    xyz     [3]float64  geocentric vector (Note 2)

Returned (function value):
    int     status:  0 = OK
                    -1 = illegal identifier (Note 3)
                    -2 = illegal case (Note 3)

Notes:

 1) The identifier n is a number that specifies the choice of
    reference ellipsoid.  The following are supported:

        n     ellipsoid

        1     WGS84
        2     GRS80
        3     WGS72

    The n value has no significance outside the SOFA software.  For
    convenience, symbols WGS84 etc. are defined in sofam.h.

 2) The height (height, given) and the geocentric vector (xyz,
    returned) are in meters.

 3) No validation is performed on the arguments elong, phi and
    height.  An error status -1 means that the identifier n is
    illegal.  An error status -2 protects against cases that would
    lead to arithmetic exceptions.  In all error cases, xyz is set
    to zeros.

 4) The inverse transformation is performed in the function iauGc2gd.

Called:
    Eform     Earth reference ellipsoids
    Gd2gce    geodetic to geocentric transformation, general
    Zp        zero p-vector
*/
func Gd2gc(n int, elong, phi, height float64, xyz *[3]float64) int {
	var j int
	var a, f float64

	/* Obtain reference ellipsoid parameters. */
	j = Eform(n, &a, &f)

	/* If OK, transform longitude, geodetic latitude, height to x,y,z. */
	if j == 0 {
		j = Gd2gce(a, f, elong, phi, height, xyz)
		if j != 0 {
			j = -2
		}
	}

	/* Deal with any errors. */
	if j != 0 {
		Zp(xyz)
	}

	/* Return the status. */
	return j
}

/*
Gd2gce Transform geodetic coordinates to geocentric for a reference
ellipsoid of specified form.

Given:
    a       float64     equatorial radius (Notes 1,4)
    f       float64     flattening (Notes 2,4)
    elong   float64     longitude (radians, east +ve)
    phi     float64     latitude (geodetic, radians, Note 4)
    height  float64     height above ellipsoid (geodetic, Notes 3,4)

Returned:
    xyz     [3]float64  geocentric vector (Note 3)

Returned (function value):
    int     status:  0 = OK
                    -1 = illegal case (Note 4)
Notes:

 1) The equatorial radius, a, can be in any units, but meters is
    the conventional choice.

 2) The flattening, f, is (for the Earth) a value around 0.00335,
    i.e. around 1/298.

 3) The equatorial radius, a, and the height, height, must be
    given in the same units, and determine the units of the
    returned geocentric vector, xyz.

 4) No validation is performed on individual arguments.  The error
    status -1 protects against (unrealistic) cases that would lead
    to arithmetic exceptions.  If an error occurs, xyz is unchanged.

 5) The inverse transformation is performed in the function
    iauGc2gde.

 6) The transformation for a standard ellipsoid (such as WGS84) can
    more conveniently be performed by calling iauGd2gc,  which uses a
    numerical code to identify the required a and f values.

References:

    Green, R.M., Spherical Astronomy, Cambridge University Press,
    (1985) Section 4.5, p96.

    Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 4.22, p202.
*/
func Gd2gce(a, f float64, elong, phi, height float64, xyz *[3]float64) int {
	var sp, cp, w, d, ac, as, r float64

	/* Functions of geodetic latitude. */
	sp = sin(phi)
	cp = cos(phi)
	w = 1.0 - f
	w = w * w
	d = cp*cp + w*sp*sp
	if d <= 0.0 {
		return -1
	}
	ac = a / sqrt(d)
	as = w * ac

	/* Geocentric vector. */
	r = (ac + height) * cp
	xyz[0] = r * cos(elong)
	xyz[1] = r * sin(elong)
	xyz[2] = (as + height) * sp

	/* Success. */
	return 0
}
