// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

/*
Fk425 Convert B1950.0 FK4 star catalog data to J2000.0 FK5

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system.

Given: (all B1950.0, FK4)
    r1950,d1950    float64   B1950.0 RA,Dec (rad)
    dr1950,dd1950  float64   B1950.0 proper motions (rad/trop.yr)
    p1950          float64   parallax (arcsec)
    v1950          float64   radial velocity (km/s, +ve = moving away)

Returned: (all J2000.0, FK5)
    r2000,d2000    float64   J2000.0 RA,Dec (rad)
    dr2000,dd2000  float64   J2000.0 proper motions (rad/Jul.yr)
    p2000          float64   parallax (arcsec)
    v2000          float64   radial velocity (km/s, +ve = moving away)

Notes:

 1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
    and are per year rather than per century.

 2) The conversion is somewhat complicated, for several reasons:

    . Change of standard epoch from B1950.0 to J2000.0.

    . An intermediate transition date of 1984 January 1.0 TT.

    . A change of precession model.

    . Change of time unit for proper motion (tropical to Julian).

    . FK4 positions include the E-terms of aberration, to simplify
      the hand computation of annual aberration.  FK5 positions
      assume a rigorous aberration computation based on the Earth's
      barycentric velocity.

    . The E-terms also affect proper motions, and in particular cause
      objects at large distances to exhibit fictitious proper
      motions.

    The algorithm is based on Smith et al. (1989) and Yallop et al.
    (1989), which presented a matrix method due to Standish (1982) as
    developed by Aoki et al. (1983), using Kinoshita's development of
    Andoyer's post-Newcomb precession.  The numerical constants from
    Seidelmann (1992) are used canonically.

 3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
    Conversions for different epochs and equinoxes would require
    additional treatment for precession, proper motion and E-terms.

 4) In the FK4 catalog the proper motions of stars within 10 degrees
    of the poles do not embody differential E-terms effects and
    should, strictly speaking, be handled in a different manner from
    stars outside these regions.  However, given the general lack of
    homogeneity of the star data available for routine astrometry,
    the difficulties of handling positions that may have been
    determined from astrometric fields spanning the polar and non-
    polar regions, the likelihood that the differential E-terms
    effect was not taken into account when allowing for proper motion
    in past astrometry, and the undesirability of a discontinuity in
    the algorithm, the decision has been made in this SOFA algorithm
    to include the effects of differential E-terms on the proper
    motions for all stars, whether polar or not.  At epoch J2000.0,
    and measuring "on the sky" rather than in terms of RA change, the
    errors resulting from this simplification are less than
    1 milliarcsecond in position and 1 milliarcsecond per century in
    proper motion.

Called:
    Anp       normalize angle into range 0 to 2pi
    Pv2s      pv-vector to spherical coordinates
    Pdp       scalar product of two p-vectors
    Pvmpv     pv-vector minus pv_vector
    Pvppv     pv-vector plus pv_vector
    S2pv      spherical coordinates to pv-vector
    Sxp       multiply p-vector by scalar

References:

    Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
    FK4-based positions of stars to epoch J2000.0 positions in
    accordance with the new IAU resolutions".  Astron.Astrophys.
    128, 263-267.

    Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
    Astronomical Almanac", ISBN 0-935702-68-7.

    Smith, C.A. et al., 1989, "The transformation of astrometric
    catalog systems to the equinox J2000.0".  Astron.J. 97, 265.

    Standish, E.M., 1982, "Conversion of positions and proper motions
    from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
    115, 1, 20-22.

    Yallop, B.D. et al., 1989, "Transformation of mean star places
    from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
    Astron.J. 97, 274.
*/
func Fk425(r1950, d1950 float64, dr1950, dd1950 float64, p1950, v1950 float64,
	r2000, d2000 *float64, dr2000, dd2000 *float64, p2000, v2000 *float64) {
	/* Radians per year to arcsec per century */
	const PMF = 100.0 * DR2AS

	/* Small number to avoid arithmetic problems */
	const TINY = 1e-30

	/* Miscellaneous */
	var r, d, ur, ud, px, rv, pxvf, w, rd float64
	var i, j, k, l int

	/* Pv-vectors */
	var r0, pv1, pv2 [2][3]float64

	/*
	 CANONICAL CONSTANTS (Seidelmann 1992)
	*/

	/* Km per sec to AU per tropical century */
	/* = 86400 * 36524.2198782 / 149597870.7 */
	const VF = 21.095

	/* Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot) */
	a := [2][3]float64{
		{-1.62557e-6, -0.31919e-6, -0.13843e-6},
		{+1.245e-3, -1.580e-3, -0.659e-3},
	}

	/* 3x2 matrix of pv-vectors (cf. Seidelmann 3.591-4, matrix M) */
	em := [2][3][2][3]float64{

		{{{+0.9999256782, -0.0111820611, -0.0048579477},
			{+0.00000242395018, -0.00000002710663, -0.00000001177656}},

			{{+0.0111820610, +0.9999374784, -0.0000271765},
				{+0.00000002710663, +0.00000242397878, -0.00000000006587}},

			{{+0.0048579479, -0.0000271474, +0.9999881997},
				{+0.00000001177656, -0.00000000006582, +0.00000242410173}}},

		{{{-0.000551, -0.238565, +0.435739},
			{+0.99994704, -0.01118251, -0.00485767}},

			{{+0.238514, -0.002667, -0.008541},
				{+0.01118251, +0.99995883, -0.00002718}},

			{{-0.435623, +0.012254, +0.002117},
				{+0.00485767, -0.00002714, +1.00000956}}},
	}

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* The FK4 data (units radians and arcsec per tropical century). */
	r = r1950
	d = d1950
	ur = dr1950 * PMF
	ud = dd1950 * PMF
	px = p1950
	rv = v1950

	/* Express as a pv-vector. */
	pxvf = px * VF
	w = rv * pxvf
	S2pv(r, d, 1.0, ur, ud, w, &r0)

	/* Allow for E-terms (cf. Seidelmann 3.591-2). */
	Pvmpv(r0, a, &pv1)
	Sxp(Pdp(r0[0], a[0]), r0[0], &pv2[0])
	Sxp(Pdp(r0[0], a[1]), r0[0], &pv2[1])
	Pvppv(pv1, pv2, &pv1)

	/* Convert pv-vector to Fricke system (cf. Seidelmann 3.591-3). */
	for i = 0; i < 2; i++ {
		for j = 0; j < 3; j++ {
			w = 0.0
			for k = 0; k < 2; k++ {
				for l = 0; l < 3; l++ {
					w += em[i][j][k][l] * pv1[k][l]
				}
			}
			pv2[i][j] = w
		}
	}

	/* Revert to catalog form. */
	Pv2s(pv2, &r, &d, &w, &ur, &ud, &rd)
	if px > TINY {
		rv = rd / pxvf
		px = px / w
	}

	/* Return the results. */
	*r2000 = Anp(r)
	*d2000 = d
	*dr2000 = ur / PMF
	*dd2000 = ud / PMF
	*v2000 = rv
	*p2000 = px
}

/*
Fk45z Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion in the FK5 system

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system,
in such a way that the FK5 proper motion is zero.  Because such a
star has, in general, a non-zero proper motion in the FK4 system,
the function requires the epoch at which the position in the FK4
system was determined.

Given:
    r1950,d1950    float64   B1950.0 FK4 RA,Dec at epoch (rad)
    bepoch         float64   Besselian epoch (e.g. 1979.3)

Returned:
    r2000,d2000    float64   J2000.0 FK5 RA,Dec (rad)

Notes:

 1) The epoch bepoch is strictly speaking Besselian, but if a
    Julian epoch is supplied the result will be affected only to a
    negligible extent.

 2) The method is from Appendix 2 of Aoki et al. (1983), but using
    the constants of Seidelmann (1992).  See the function Fk425
    for a general introduction to the FK4 to FK5 conversion.

 3) Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only
    is provided for.  Conversions for different starting and/or
    ending epochs would require additional treatment for precession,
    proper motion and E-terms.

 4) In the FK4 catalog the proper motions of stars within 10 degrees
    of the poles do not embody differential E-terms effects and
    should, strictly speaking, be handled in a different manner from
    stars outside these regions.  However, given the general lack of
    homogeneity of the star data available for routine astrometry,
    the difficulties of handling positions that may have been
    determined from astrometric fields spanning the polar and non-
    polar regions, the likelihood that the differential E-terms
    effect was not taken into account when allowing for proper motion
    in past astrometry, and the undesirability of a discontinuity in
    the algorithm, the decision has been made in this SOFA algorithm
    to include the effects of differential E-terms on the proper
    motions for all stars, whether polar or not.  At epoch 2000.0,
    and measuring "on the sky" rather than in terms of RA change, the
    errors resulting from this simplification are less than
    1 milliarcsecond in position and 1 milliarcsecond per century in
    proper motion.

References:

    Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
    FK4-based positions of stars to epoch J2000.0 positions in
    accordance with the new IAU resolutions".  Astron.Astrophys.
    128, 263-267.

    Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
    Astronomical Almanac", ISBN 0-935702-68-7.

Called:
    Anp       normalize angle into range 0 to 2pi
    C2s       p-vector to spherical
    Epb2jd    Besselian epoch to Julian date
    Epj       Julian date to Julian epoch
    Pdp       scalar product of two p-vectors
    Pmp       p-vector minus p-vector
    Ppsp      p-vector plus scaled p-vector
    Pvu       update a pv-vector
    S2c       spherical to p-vector
*/
func Fk45z(r1950, d1950, bepoch float64, r2000, d2000 *float64) {
	/* Radians per year to arcsec per century */
	const PMF = 100.0 * DR2AS

	/* Position and position+velocity vectors */
	var r0, p [3]float64
	var pv [2][3]float64

	/* Miscellaneous */
	var w, djm0, djm float64
	var i, j, k int

	/*
		CANONICAL CONSTANTS (Seidelmann 1992)
	*/

	/* Vectors A and Adot (Seidelmann 3.591-2) */
	a := [3]float64{-1.62557e-6, -0.31919e-6, -0.13843e-6}
	ad := [3]float64{+1.245e-3, -1.580e-3, -0.659e-3}

	/* 3x2 matrix of p-vectors (cf. Seidelmann 3.591-4, matrix M) */
	em := [2][3][3]float64{
		{{+0.9999256782, -0.0111820611, -0.0048579477},
			{+0.0111820610, +0.9999374784, -0.0000271765},
			{+0.0048579479, -0.0000271474, +0.9999881997}},
		{{-0.000551, -0.238565, +0.435739},
			{+0.238514, -0.002667, -0.008541},
			{-0.435623, +0.012254, +0.002117}},
	}

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* Spherical coordinates to p-vector. */
	S2c(r1950, d1950, &r0)

	/* Adjust p-vector A to give zero proper motion in FK5. */
	w = (bepoch - 1950) / PMF
	Ppsp(a, w, ad, &p)

	/* Remove E-terms. */
	Ppsp(p, -Pdp(r0, p), r0, &p)
	Pmp(r0, p, &p)

	/* Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3). */
	for i = 0; i < 2; i++ {
		for j = 0; j < 3; j++ {
			w = 0.0
			for k = 0; k < 3; k++ {
				w += em[i][j][k] * p[k]
			}
			pv[i][j] = w
		}
	}

	/* Allow for fictitious proper motion. */
	Epb2jd(bepoch, &djm0, &djm)
	w = (Epj(djm0, djm) - 2000.0) / PMF
	Pvu(w, pv, &pv)

	/* Revert to spherical coordinates. */
	C2s(pv[0], &w, d2000)
	*r2000 = Anp(w)
}

/*
Fk524 Convert J2000.0 FK5 star catalog data to B1950.0 FK4

Given: (all J2000.0, FK5)
    r2000,d2000    float64   J2000.0 RA,Dec (rad)
    dr2000,dd2000  float64   J2000.0 proper motions (rad/Jul.yr)
    p2000          float64   parallax (arcsec)
    v2000          float64   radial velocity (km/s, +ve = moving away)

Returned: (all B1950.0, FK4)
    r1950,d1950    float64   B1950.0 RA,Dec (rad)
    dr1950,dd1950  float64   B1950.0 proper motions (rad/trop.yr)
    p1950          float64   parallax (arcsec)
    v1950          float64   radial velocity (km/s, +ve = moving away)

Notes:

 1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
    and are per year rather than per century.

 2) The conversion is somewhat complicated, for several reasons:

    . Change of standard epoch from J2000.0 to B1950.0.

    . An intermediate transition date of 1984 January 1.0 TT.

    . A change of precession model.

    . Change of time unit for proper motion (Julian to tropical).

    . FK4 positions include the E-terms of aberration, to simplify
      the hand computation of annual aberration.  FK5 positions
      assume a rigorous aberration computation based on the Earth's
      barycentric velocity.

    . The E-terms also affect proper motions, and in particular cause
      objects at large distances to exhibit fictitious proper
      motions.

 3) The algorithm is based on Smith et al. (1989) and Yallop et al.
    (1989), which presented a matrix method due to Standish (1982) as
    developed by Aoki et al. (1983), using Kinoshita's development of
    Andoyer's post-Newcomb precession.  The numerical constants from
    Seidelmann (1992) are used canonically.

 4) In the FK4 catalog the proper motions of stars within 10 degrees
    of the poles do not embody differential E-terms effects and
    should, strictly speaking, be handled in a different manner from
    stars outside these regions.  However, given the general lack of
    homogeneity of the star data available for routine astrometry,
    the difficulties of handling positions that may have been
    determined from astrometric fields spanning the polar and non-
    polar regions, the likelihood that the differential E-terms
    effect was not taken into account when allowing for proper motion
    in past astrometry, and the undesirability of a discontinuity in
    the algorithm, the decision has been made in this SOFA algorithm
    to include the effects of differential E-terms on the proper
    motions for all stars, whether polar or not.  At epoch J2000.0,
    and measuring "on the sky" rather than in terms of RA change, the
    errors resulting from this simplification are less than
    1 milliarcsecond in position and 1 milliarcsecond per century in
    proper motion.

Called:
    Anp       normalize angle into range 0 to 2pi
    Pdp       scalar product of two p-vectors
    Pm        modulus of p-vector
    Pmp       p-vector minus p-vector
    Ppp       p-vector pluus p-vector
    Pv2s      pv-vector to spherical coordinates
    S2pv      spherical coordinates to pv-vector
    Sxp       multiply p-vector by scalar

References:

    Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
    FK4-based positions of stars to epoch J2000.0 positions in
    accordance with the new IAU resolutions".  Astron.Astrophys.
    128, 263-267.

    Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
    Astronomical Almanac", ISBN 0-935702-68-7.

    Smith, C.A. et al., 1989, "The transformation of astrometric
    catalog systems to the equinox J2000.0".  Astron.J. 97, 265.

    Standish, E.M., 1982, "Conversion of positions and proper motions
    from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
    115, 1, 20-22.

    Yallop, B.D. et al., 1989, "Transformation of mean star places
    from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
    Astron.J. 97, 274.
*/
func Fk524(r2000, d2000 float64, dr2000, dd2000 float64, p2000, v2000 float64,
	r1950, d1950 *float64, dr1950, dd1950 *float64, p1950, v1950 *float64) {
	/* Radians per year to arcsec per century */
	const PMF = 100.0 * DR2AS

	/* Small number to avoid arithmetic problems */
	const TINY = 1e-30

	/* Miscellaneous */
	var r, d, ur, ud, px, rv, pxvf, w, rd float64
	var i, j, k, l int

	/* Vectors, p and pv */
	var r0, r1, pv [2][3]float64
	var p1, p2 [3]float64

	/*
	 CANONICAL CONSTANTS (Seidelmann 1992)
	*/

	/* Km per sec to AU per tropical century */
	/* = 86400 * 36524.2198782 / 149597870.7 */
	const VF = 21.095

	/* Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot) */
	a := [2][3]float64{
		{-1.62557e-6, -0.31919e-6, -0.13843e-6},
		{+1.245e-3, -1.580e-3, -0.659e-3},
	}

	/* 3x2 matrix of pv-vectors (cf. Seidelmann 3.592-1, matrix M^-1) */
	em := [2][3][2][3]float64{

		{{
			{+0.9999256795, +0.0111814828, +0.0048590039},
			{-0.00000242389840, -0.00000002710544, -0.00000001177742},
		}, {
			{-0.0111814828, +0.9999374849, -0.0000271771},
			{+0.00000002710544, -0.00000242392702, +0.00000000006585},
		}, {
			{-0.0048590040, -0.0000271557, +0.9999881946},
			{+0.00000001177742, +0.00000000006585, -0.00000242404995},
		},
		},

		{{
			{-0.000551, +0.238509, -0.435614},
			{+0.99990432, +0.01118145, +0.00485852},
		}, {
			{-0.238560, -0.002667, +0.012254},
			{-0.01118145, +0.99991613, -0.00002717},
		}, {
			{+0.435730, -0.008541, +0.002117},
			{-0.00485852, -0.00002716, +0.99996684},
		}},
	}

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* The FK5 data (units radians and arcsec per Julian century). */
	r = r2000
	d = d2000
	ur = dr2000 * PMF
	ud = dd2000 * PMF
	px = p2000
	rv = v2000

	/* Express as a pv-vector. */
	pxvf = px * VF
	w = rv * pxvf
	S2pv(r, d, 1.0, ur, ud, w, &r0)

	/* Convert pv-vector to Bessel-Newcomb system (cf. Seidelmann 3.592-1). */
	for i = 0; i < 2; i++ {
		for j = 0; j < 3; j++ {
			w = 0.0
			for k = 0; k < 2; k++ {
				for l = 0; l < 3; l++ {
					w += em[i][j][k][l] * r0[k][l]
				}
			}
			r1[i][j] = w
		}
	}

	/* Apply E-terms (equivalent to Seidelmann 3.592-3, one iteration). */

	/* Direction. */
	w = Pm(r1[0])
	Sxp(Pdp(r1[0], a[0]), r1[0], &p1)
	Sxp(w, a[0], &p2)
	Pmp(p2, p1, &p1)
	Ppp(r1[0], p1, &p1)

	/* Recompute length. */
	w = Pm(p1)

	/* Direction. */
	Sxp(Pdp(r1[0], a[0]), r1[0], &p1)
	Sxp(w, a[0], &p2)
	Pmp(p2, p1, &p1)
	Ppp(r1[0], p1, &pv[0])

	/* Derivative. */
	Sxp(Pdp(r1[0], a[1]), pv[0], &p1)
	Sxp(w, a[1], &p2)
	Pmp(p2, p1, &p1)
	Ppp(r1[1], p1, &pv[1])

	/* Revert to catalog form. */
	Pv2s(pv, &r, &d, &w, &ur, &ud, &rd)
	if px > TINY {
		rv = rd / pxvf
		px = px / w
	}

	/* Return the results. */
	*r1950 = Anp(r)
	*d1950 = d
	*dr1950 = ur / PMF
	*dd1950 = ud / PMF
	*p1950 = px
	*v1950 = rv
}

/*
Fk52h Transform FK5  (J2000.0) star data into the Hipparcos frame

Given (all FK5, equinox J2000.0, epoch J2000.0):
    r5      float64    RA (radians)
    d5      float64    Dec (radians)
    dr5     float64    proper motion in RA (dRA/dt, rad/Jyear)
    dd5     float64    proper motion in Dec (dDec/dt, rad/Jyear)
    px5     float64    parallax (arcsec)
    rv5     float64    radial velocity (km/s, positive = receding)

Returned (all Hipparcos, epoch J2000.0):
    rh      float64    RA (radians)
    dh      float64    Dec (radians)
    drh     float64    proper motion in RA (dRA/dt, rad/Jyear)
    ddh     float64    proper motion in Dec (dDec/dt, rad/Jyear)
    pxh     float64    parallax (arcsec)
    rvh     float64    radial velocity (km/s, positive = receding)

Notes:

 1) This function transforms FK5 star positions and proper motions
    into the system of the Hipparcos catalog.

 2) The proper motions in RA are dRA/dt rather than
    cos(Dec)*dRA/dt, and are per year rather than per century.

 3) The FK5 to Hipparcos transformation is modeled as a pure
    rotation and spin;  zonal errors in the FK5 catalog are not
    taken into account.

 4) See also H2fk5, Fk5hz, Hfk5z.

Called:
    Starpv    star catalog data to space motion pv-vector
    Fk5hip    FK5 to Hipparcos rotation and spin
    Rxp       product of r-matrix and p-vector
    Pxp       vector product of two p-vectors
    Ppp       p-vector plus p-vector
    Pvstar    space motion pv-vector to star catalog data

Reference:

    F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
*/
func Fk52h(r5, d5 float64, dr5, dd5, px5, rv5 float64,
	rh, dh *float64, drh, ddh, pxh, rvh *float64) {

	var i int
	var pv5, pvh [2][3]float64
	var r5h [3][3]float64
	var s5h, wxp, vv [3]float64

	/* FK5 barycentric position/velocity pv-vector (normalized). */
	Starpv(r5, d5, dr5, dd5, px5, rv5, &pv5)

	/* FK5 to Hipparcos orientation matrix and spin vector. */
	Fk5hip(&r5h, &s5h)

	/* Make spin units per day instead of per year. */
	for i = 0; i < 3; i++ {
		s5h[i] /= 365.25
	}

	/* Orient the FK5 position into the Hipparcos system. */
	Rxp(r5h, pv5[0], &pvh[0])

	/* Apply spin to the position giving an extra space motion component. */
	Pxp(pv5[0], s5h, &wxp)

	/* Add this component to the FK5 space motion. */
	Ppp(wxp, pv5[1], &vv)

	/* Orient the FK5 space motion into the Hipparcos system. */
	Rxp(r5h, vv, &pvh[1])

	/* Hipparcos pv-vector to spherical. */
	Pvstar(pvh, rh, dh, drh, ddh, pxh, rvh)
}

/*
Fk54z Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero proper motion in FK5 system and zero parallax

Given:
    r2000,d2000    float64   J2000.0 FK5 RA,Dec (rad)
    bepoch         float64   Besselian epoch (e.g. 1950.0)

Returned:
    r1950,d1950    float64   B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH
    dr1950,dd1950  float64   B1950.0 FK4 proper motions (rad/trop.yr)

Notes:

 1) In contrast to the Fk524 function, here the FK5 proper
    motions, the parallax and the radial velocity are presumed zero.

 2) This function converts a star position from the IAU 1976 FK5
   (Fricke) system to the former FK4 (Bessel-Newcomb) system, for
    cases such as distant radio sources where it is presumed there is
    zero parallax and no proper motion.  Because of the E-terms of
    aberration, such objects have (in general) non-zero proper motion
    in FK4, and the present function returns those fictitious proper
    motions.

 3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
    Conversions involving other equinoxes would require additional
    treatment for precession.

 4) The position returned by this function is in the B1950.0 FK4
    reference system but at Besselian epoch BEPOCH.  For comparison
    with catalogs the BEPOCH argument will frequently be 1950.0. (In
    this context the distinction between Besselian and Julian epoch
    is insignificant.)

 5) The RA component of the returned (fictitious) proper motion is
    dRA/dt rather than cos(Dec)*dRA/dt.

Called:
    Anp       normalize angle into range 0 to 2pi
    C2s       p-vector to spherical
    Fk524     FK4 to FK5
    S2c       spherical to p-vector
*/
func Fk54z(r2000, d2000, bepoch float64, r1950, d1950 *float64, dr1950, dd1950 *float64) {
	var r, d, pr, pd, px, rv, w float64
	var p, v [3]float64
	var i int

	/* FK5 equinox J2000.0 to FK4 equinox B1950.0. */
	Fk524(r2000, d2000, 0.0, 0.0, 0.0, 0.0,
		&r, &d, &pr, &pd, &px, &rv)

	/* Spherical to Cartesian. */
	S2c(r, d, &p)

	/* Fictitious proper motion (radians per year). */
	v[0] = -pr*p[1] - pd*cos(r)*sin(d)
	v[1] = pr*p[0] - pd*sin(r)*sin(d)
	v[2] = pd * cos(d)

	/* Apply the motion. */
	w = bepoch - 1950.0
	for i = 0; i < 3; i++ {
		p[i] += w * v[i]
	}

	/* Cartesian to spherical. */
	C2s(p, &w, d1950)
	*r1950 = Anp(w)

	/* Fictitious proper motion. */
	*dr1950 = pr
	*dd1950 = pd
}

/*
Fk5hip FK5 orientation and spin with respect to Hipparcos

Returned:
    r5h   [3][3]float64  r-matrix: FK5 rotation wrt Hipparcos (Note 2)
    s5h   [3]float64     r-vector: FK5 spin wrt Hipparcos (Note 3)

Notes:

 1) This function models the FK5 to Hipparcos transformation as a
    pure rotation and spin;  zonal errors in the FK5 catalogue are
    not taken into account.

 2) The r-matrix r5h operates in the sense:

          P_Hipparcos = r5h x P_FK5

    where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
    the equivalent Hipparcos p-vector.

 3) The r-vector s5h represents the time derivative of the FK5 to
    Hipparcos rotation.  The units are radians per year (Julian,
    TDB).

Called:
    Rv2m      r-vector to r-matrix

Reference:

    F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
*/
func Fk5hip(r5h *[3][3]float64, s5h *[3]float64) {
	var v [3]float64

	/* FK5 wrt Hipparcos orientation and spin (radians, radians/year) */
	var epx, epy, epz float64
	var omx, omy, omz float64

	epx = -19.9e-3 * DAS2R
	epy = -9.1e-3 * DAS2R
	epz = 22.9e-3 * DAS2R

	omx = -0.30e-3 * DAS2R
	omy = 0.60e-3 * DAS2R
	omz = 0.70e-3 * DAS2R

	/* FK5 to Hipparcos orientation expressed as an r-vector. */
	v[0] = epx
	v[1] = epy
	v[2] = epz

	/* Re-express as an r-matrix. */
	Rv2m(v, r5h)

	/* Hipparcos wrt FK5 spin expressed as an r-vector. */
	s5h[0] = omx
	s5h[1] = omy
	s5h[2] = omz
}

/*
Fk5hz FK5 to Hipparcos assuming zero Hipparcos proper motion

Transform an FK5 (J2000.0) star position into the system of the
Hipparcos catalogue, assuming zero Hipparcos proper motion.

Given:
    r5           float64   FK5 RA (radians), equinox J2000.0, at date
    d5           float64   FK5 Dec (radians), equinox J2000.0, at date
    date1,date2  float64   TDB date (Notes 1,2)

Returned:
    rh           float64   Hipparcos RA (radians)
    dh           float64   Hipparcos Dec (radians)

Notes:

 1) This function converts a star position from the FK5 system to
    the Hipparcos system, in such a way that the Hipparcos proper
    motion is zero.  Because such a star has, in general, a non-zero
    proper motion in the FK5 system, the function requires the date
    at which the position in the FK5 system was determined.

 2) The TT date date1+date2 is a Julian Date, apportioned in any
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

 3) The FK5 to Hipparcos transformation is modeled as a pure
    rotation and spin;  zonal errors in the FK5 catalogue are not
    taken into account.

 4) The position returned by this function is in the Hipparcos
    reference system but at date date1+date2.

 5) See also Fk52h, H2fk5, Hfk5z.

Called:
    S2c       spherical coordinates to unit vector
    Fk5hip    FK5 to Hipparcos rotation and spin
    Sxp       multiply p-vector by scalar
    Rv2m      r-vector to r-matrix
    Trxp      product of transpose of r-matrix and p-vector
    Pxp       vector product of two p-vectors
    C2s       p-vector to spherical
    Anp       normalize angle into range 0 to 2pi

Reference:

    F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
*/
func Fk5hz(r5, d5 float64, date1, date2 float64, rh, dh *float64) {
	var t, w float64
	var p5e, s5h, vst, p5, ph [3]float64
	var r5h, rst [3][3]float64

	/* Interval from given date to fundamental epoch J2000.0 (JY). */
	t = -((date1 - DJ00) + date2) / DJY

	/* FK5 barycentric position vector. */
	S2c(r5, d5, &p5e)

	/* FK5 to Hipparcos orientation matrix and spin vector. */
	Fk5hip(&r5h, &s5h)

	/* Accumulated Hipparcos wrt FK5 spin over that interval. */
	Sxp(t, s5h, &vst)

	/* Express the accumulated spin as a rotation matrix. */
	Rv2m(vst, &rst)

	/* Derotate the vector's FK5 axes back to date. */
	Trxp(rst, p5e, &p5)

	/* Rotate the vector into the Hipparcos system. */
	Rxp(r5h, p5, &ph)

	/* Hipparcos vector to spherical. */
	C2s(ph, &w, dh)
	*rh = Anp(w)
}

/*
H2fk5 Transform Hipparcos star data into the FK5 (J2000.0) frame

Given (all Hipparcos, epoch J2000.0):
    rh      float64    RA (radians)
    dh      float64    Dec (radians)
    drh     float64    proper motion in RA (dRA/dt, rad/Jyear)
    ddh     float64    proper motion in Dec (dDec/dt, rad/Jyear)
    pxh     float64    parallax (arcsec)
    rvh     float64    radial velocity (km/s, positive = receding)

Returned (all FK5, equinox J2000.0, epoch J2000.0):
    r5      float64    RA (radians)
    d5      float64    Dec (radians)
    dr5     float64    proper motion in RA (dRA/dt, rad/Jyear)
    dd5     float64    proper motion in Dec (dDec/dt, rad/Jyear)
    px5     float64    parallax (arcsec)
    rv5     float64    radial velocity (km/s, positive = receding)

Notes:

 1) This function transforms Hipparcos star positions and proper
    motions into FK5 J2000.0.

 2) The proper motions in RA are dRA/dt rather than
    cos(Dec)*dRA/dt, and are per year rather than per century.

 3) The FK5 to Hipparcos transformation is modeled as a pure
    rotation and spin;  zonal errors in the FK5 catalog are not
    taken into account.

 4) See also Fk52h, Fk5hz, Hfk5z.

Called:
    Starpv    star catalog data to space motion pv-vector
    Fk5hip    FK5 to Hipparcos rotation and spin
    Rv2m      r-vector to r-matrix
    Rxp       product of r-matrix and p-vector
    Trxp      product of transpose of r-matrix and p-vector
    Pxp       vector product of two p-vectors
    Pmp       p-vector minus p-vector
    Pvstar    space motion pv-vector to star catalog data

Reference:

    F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
*/
func H2fk5(rh, dh float64, drh, ddh, pxh, rvh float64,
	r5, d5 *float64, dr5, dd5, px5, rv5 *float64) {

	var i int
	var pvh, pv5 [2][3]float64
	var r5h [3][3]float64
	var s5h, sh, wxp, vv [3]float64

	/* Hipparcos barycentric position/velocity pv-vector (normalized). */
	Starpv(rh, dh, drh, ddh, pxh, rvh, &pvh)

	/* FK5 to Hipparcos orientation matrix and spin vector. */
	Fk5hip(&r5h, &s5h)

	/* Make spin units per day instead of per year. */
	for i = 0; i < 3; i++ {
		s5h[i] /= 365.25
	}

	/* Orient the spin into the Hipparcos system. */
	Rxp(r5h, s5h, &sh)

	/* De-orient the Hipparcos position into the FK5 system. */
	Trxp(r5h, pvh[0], &pv5[0])

	/* Apply spin to the position giving an extra space motion component. */
	Pxp(pvh[0], sh, &wxp)

	/* Subtract this component from the Hipparcos space motion. */
	Pmp(pvh[1], wxp, &vv)

	/* De-orient the Hipparcos space motion into the FK5 system. */
	Trxp(r5h, vv, &pv5[1])

	/* FK5 pv-vector to spherical. */
	Pvstar(pv5, r5, d5, dr5, dd5, px5, rv5)
}

/*
Hfk5z Hipparcos to FK5 assuming zero Hipparcos proper motion

Transform a Hipparcos star position into FK5 J2000.0, assuming
zero Hipparcos proper motion.

Given:
    rh            float64    Hipparcos RA (radians)
    dh            float64    Hipparcos Dec (radians)
    date1,date2   float64    TDB date (Note 1)

Returned (all FK5, equinox J2000.0, date date1+date2):
    r5            float64    RA (radians)
    d5            float64    Dec (radians)
    dr5           float64    FK5 RA proper motion (rad/year, Note 4)
    dd5           float64    Dec proper motion (rad/year, Note 4)

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

 2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

 3) The FK5 to Hipparcos transformation is modeled as a pure rotation
    and spin;  zonal errors in the FK5 catalogue are not taken into
    account.

 4) It was the intention that Hipparcos should be a close
    approximation to an inertial frame, so that distant objects have
    zero proper motion;  such objects have (in general) non-zero
    proper motion in FK5, and this function returns those fictitious
    proper motions.

 5) The position returned by this function is in the FK5 J2000.0
    reference system but at date date1+date2.

 6) See also Fk52h, H2fk5, Fk5zhz.

Called:
    S2c       spherical coordinates to unit vector
    Fk5hip    FK5 to Hipparcos rotation and spin
    Rxp       product of r-matrix and p-vector
    Sxp       multiply p-vector by scalar
    Rxr       product of two r-matrices
    Trxp      product of transpose of r-matrix and p-vector
    Pxp       vector product of two p-vectors
    Pv2s      pv-vector to spherical
    Anp       normalize angle into range 0 to 2pi

Reference:

    F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
*/
func Hfk5z(rh, dh float64, date1, date2 float64, r5, d5, dr5, dd5 *float64) {
	var t, w, r, v float64
	var ph, s5h, sh, vst, vv [3]float64
	var pv5e [2][3]float64
	var r5h, rst, r5ht [3][3]float64

	/* Time interval from fundamental epoch J2000.0 to given date (JY). */
	t = ((date1 - DJ00) + date2) / DJY

	/* Hipparcos barycentric position vector (normalized). */
	S2c(rh, dh, &ph)

	/* FK5 to Hipparcos orientation matrix and spin vector. */
	Fk5hip(&r5h, &s5h)

	/* Rotate the spin into the Hipparcos system. */
	Rxp(r5h, s5h, &sh)

	/* Accumulated Hipparcos wrt FK5 spin over that interval. */
	Sxp(t, s5h, &vst)

	/* Express the accumulated spin as a rotation matrix. */
	Rv2m(vst, &rst)

	/* Rotation matrix:  accumulated spin, then FK5 to Hipparcos. */
	Rxr(r5h, rst, &r5ht)

	/* De-orient & de-spin the Hipparcos position into FK5 J2000.0. */
	Trxp(r5ht, ph, &pv5e[0])

	/* Apply spin to the position giving a space motion. */
	Pxp(sh, ph, &vv)

	/* De-orient & de-spin the Hipparcos space motion into FK5 J2000.0. */
	Trxp(r5ht, vv, &pv5e[1])

	/* FK5 position/velocity pv-vector to spherical. */
	Pv2s(pv5e, &w, d5, &r, dr5, dd5, &v)
	*r5 = Anp(w)
}
