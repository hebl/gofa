// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

/*
Tpxes Project celestial to tangent plane, spherical

In the tangent plane projection, given celestial spherical
coordinates for a star and the tangent point, solve for the star's
rectangular coordinates in the tangent plane.

Given:
    a,b       float64  star's spherical coordinates
    a0,b0     float64  tangent point's spherical coordinates

Returned:
    xi,eta    float64  rectangular coordinates of star image (Note 2)

Returned (function value):
    int     status: 0 = OK
                    1 = star too far from axis
                    2 = antistar on tangent plane
                    3 = antistar too far from axis

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the spherical coordinates are observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  For right-handed spherical coordinates,
    (xi,eta) are also right-handed.  The units of (xi,eta) are,
    effectively, radians at the tangent point.

 3) All angular arguments are in radians.

 4) This function is a member of the following set:

      spherical      vector         solve for

    > iauTpxes <    iauTpxev         xi,eta
      iauTpsts      iauTpstv          star
      iauTpors      iauTporv         origin

References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
    1987, Chapter 13.
*/
func Tpxes(a, b, a0, b0 float64, xi, eta *float64) int {
	const TINY = 1e-6

	var j int
	var sb0, sb, cb0, cb, da, sda, cda, d float64

	/* Functions of the spherical coordinates. */
	sb0 = sin(b0)
	sb = sin(b)
	cb0 = cos(b0)
	cb = cos(b)
	da = a - a0
	sda = sin(da)
	cda = cos(da)

	/* Reciprocal of star vector length to tangent plane. */
	d = sb*sb0 + cb*cb0*cda

	/* Check for error cases. */
	if d > TINY {
		j = 0
	} else if d >= 0.0 {
		j = 1
		d = TINY
	} else if d > -TINY {
		j = 2
		d = -TINY
	} else {
		j = 3
	}

	/* Return the tangent plane coordinates (even in dubious cases). */
	*xi = cb * sda / d
	*eta = (sb*cb0 - cb*sb0*cda) / d

	/* Return the status. */
	return j
}

/*
Tpxev Project celestial to tangent plane, vector

In the tangent plane projection, given celestial direction cosines
for a star and the tangent point, solve for the star's rectangular
coordinates in the tangent plane.

Given:
    v         [3]float64  direction cosines of star (Note 4)
    v0        [3]float64  direction cosines of tangent point (Note 4)

Returned:
    xi,eta  float64     tangent plane coordinates of star

Returned (function value):
    int     status: 0 = OK
                    1 = star too far from axis
                    2 = antistar on tangent plane
                    3 = antistar too far from axis

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the direction cosines represent observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  If the direction cosines are with
    respect to a right-handed triad, (xi,eta) are also right-handed.
    The units of (xi,eta) are, effectively, radians at the tangent
    point.

 3) The method used is to extend the star vector to the tangent
    plane and then rotate the triad so that (x,y) becomes (xi,eta).
    Writing (a,b) for the celestial spherical coordinates of the
    star, the sequence of rotations is (a+pi/2) around the z-axis
    followed by (pi/2-b) around the x-axis.

 4) If vector v0 is not of unit length, or if vector v is of zero
    length, the results will be wrong.

 5) If v0 points at a pole, the returned (xi,eta) will be based on
    the arbitrary assumption that the longitude coordinate of the
    tangent point is zero.

 6) This function is a member of the following set:

    spherical      vector         solve for

    iauTpxes    > iauTpxev <       xi,eta
    iauTpsts      iauTpstv          star
    iauTpors      iauTporv         origin

References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
        1987, Chapter 13.
*/
func Tpxev(v, v0 [3]float64, xi, eta *float64) int {
	const TINY = 1e-6
	var j int
	var x, y, z, x0, y0, z0, r2, r, w, d float64

	/* Star and tangent point. */
	x = v[0]
	y = v[1]
	z = v[2]
	x0 = v0[0]
	y0 = v0[1]
	z0 = v0[2]

	/* Deal with polar case. */
	r2 = x0*x0 + y0*y0
	r = sqrt(r2)
	if r == 0.0 {
		r = 1e-20
		x0 = r
	}

	/* Reciprocal of star vector length to tangent plane. */
	w = x*x0 + y*y0
	d = w + z*z0

	/* Check for error cases. */
	if d > TINY {
		j = 0
	} else if d >= 0.0 {
		j = 1
		d = TINY
	} else if d > -TINY {
		j = 2
		d = -TINY
	} else {
		j = 3
	}

	/* Return the tangent plane coordinates (even in dubious cases). */
	d *= r
	*xi = (y*x0 - x*y0) / d
	*eta = (z*r2 - z0*w) / d

	/* Return the status. */
	return j
}

/*
Tpsts Project tangent plane to celestial, spherical

In the tangent plane projection, given the star's rectangular
coordinates and the spherical coordinates of the tangent point,
solve for the spherical coordinates of the star.

Given:
    xi,eta    float64  rectangular coordinates of star image (Note 2)
    a0,b0     float64  tangent point's spherical coordinates

Returned:
    a,b     float64  star's spherical coordinates

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the spherical coordinates are observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  If the spherical coordinates are with
    respect to a right-handed triad, (xi,eta) are also right-handed.
    The units of (xi,eta) are, effectively, radians at the tangent
    point.

 3) All angular arguments are in radians.

 4) This function is a member of the following set:

      spherical      vector         solve for

      iauTpxes      iauTpxev         xi,eta
    > iauTpsts <    iauTpstv          star
      iauTpors      iauTporv         origin

Called:
    Anp       normalize angle into range 0 to 2pi

References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
    1987, Chapter 13.
*/
func Tpsts(xi, eta, a0, b0 float64, a, b *float64) {
	var sb0, cb0, d float64

	sb0 = sin(b0)
	cb0 = cos(b0)
	d = cb0 - eta*sb0
	*a = Anp(atan2(xi, d) + a0)
	*b = atan2(sb0+eta*cb0, sqrt(xi*xi+d*d))
}

/*
Tpstv Project tangent plane to celestial, vector

In the tangent plane projection, given the star's rectangular
coordinates and the direction cosines of the tangent point, solve
for the direction cosines of the star.

Given:
    xi,eta  float64     rectangular coordinates of star image (Note 2)
    v0      [3]float64  tangent point's direction cosines

Returned:
    v       [3]float64  star's direction cosines

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the direction cosines represent observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  If the direction cosines are with
    respect to a right-handed triad, (xi,eta) are also right-handed.
    The units of (xi,eta) are, effectively, radians at the tangent
    point.

 3) The method used is to complete the star vector in the (xi,eta)
    based triad and normalize it, then rotate the triad to put the
    tangent point at the pole with the x-axis aligned to zero
    longitude.  Writing (a0,b0) for the celestial spherical
    coordinates of the tangent point, the sequence of rotations is
    (b-pi/2) around the x-axis followed by (-a-pi/2) around the
    z-axis.

 4) If vector v0 is not of unit length, the returned vector v will
    be wrong.

 5) If vector v0 points at a pole, the returned vector v will be
    based on the arbitrary assumption that the longitude coordinate
    of the tangent point is zero.

 6) This function is a member of the following set:

        spherical      vector         solve for

        iauTpxes      iauTpxev         xi,eta
        iauTpsts    > iauTpstv <        star
        iauTpors      iauTporv         origin

 References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
    1987, Chapter 13.
*/
func Tpstv(xi, eta float64, v0 [3]float64, v *[3]float64) {
	var x, y, z, f, r float64

	/* Tangent point. */
	x = v0[0]
	y = v0[1]
	z = v0[2]

	/* Deal with polar case. */
	r = sqrt(x*x + y*y)
	if r == 0.0 {
		r = 1e-20
		x = r
	}

	/* Star vector length to tangent plane. */
	f = sqrt(1.0 + xi*xi + eta*eta)

	/* Apply the transformation and normalize. */
	v[0] = (x - (xi*y+eta*x*z)/r) / f
	v[1] = (y + (xi*x-eta*y*z)/r) / f
	v[2] = (z + eta*r) / f
}

/*
Tpors Solve for tangent point, spherical

In the tangent plane projection, given the rectangular coordinates
of a star and its spherical coordinates, determine the spherical
coordinates of the tangent point.

Given:
    xi,eta     float64  rectangular coordinates of star image (Note 2)
    a,b        float64  star's spherical coordinates (Note 3)

Returned:
    a01,b01  float64  tangent point's spherical coordinates, Soln. 1
    a02,b02  float64  tangent point's spherical coordinates, Soln. 2

Returned (function value):
    int     number of solutions:
            0 = no solutions returned (Note 5)
            1 = only the first solution is useful (Note 6)
            2 = both solutions are useful (Note 6)

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the spherical coordinates are observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  If the spherical coordinates are with
    respect to a right-handed triad, (xi,eta) are also right-handed.
    The units of (xi,eta) are, effectively, radians at the tangent
    point.

 3) All angular arguments are in radians.

 4) The angles a01 and a02 are returned in the range 0-2pi.  The
    angles b01 and b02 are returned in the range +/-pi, but in the
    usual, non-pole-crossing, case, the range is +/-pi/2.

 5) Cases where there is no solution can arise only near the poles.
    For example, it is clearly impossible for a star at the pole
    itself to have a non-zero xi value, and hence it is meaningless
    to ask where the tangent point would have to be to bring about
    this combination of xi and dec.

 6) Also near the poles, cases can arise where there are two useful
    solutions.  The return value indicates whether the second of the
    two solutions returned is useful;  1 indicates only one useful
    solution, the usual case.

 7) The basis of the algorithm is to solve the spherical triangle PSC,
    where P is the north celestial pole, S is the star and C is the
    tangent point.  The spherical coordinates of the tangent point are
    [a0,b0];  writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), side c
    is then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
    found) is (pi/2-b0).  Angle C is given by sin(C) = xi/rho and
    cos(C) = eta/rho.  Angle P (to be found) is the longitude
    difference between star and tangent point (a-a0).

 8) This function is a member of the following set:

        spherical      vector         solve for

        iauTpxes      iauTpxev         xi,eta
        iauTpsts      iauTpstv          star
      > iauTpors <    iauTporv         origin

 Called:
    Anp       normalize angle into range 0 to 2pi

 References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
    1987, Chapter 13.
*/
func Tpors(xi, eta, a, b float64, a01, b01, a02, b02 *float64) int {
	var xi2, r, sb, cb, rsb, rcb, w2, w, s, c float64

	xi2 = xi * xi
	r = sqrt(1.0 + xi2 + eta*eta)
	sb = sin(b)
	cb = cos(b)
	rsb = r * sb
	rcb = r * cb
	w2 = rcb*rcb - xi2
	if w2 >= 0.0 {
		w = sqrt(w2)
		s = rsb - eta*w
		c = rsb*eta + w
		if xi == 0.0 && w == 0.0 {
			w = 1.0
		}
		*a01 = Anp(a - atan2(xi, w))
		*b01 = atan2(s, c)
		w = -w
		s = rsb - eta*w
		c = rsb*eta + w
		*a02 = Anp(a - atan2(xi, w))
		*b02 = atan2(s, c)
		if fabs(rsb) < 1.0 {
			return 1
		} else {
			return 2
		}
	} else {
		return 0
	}
}

/*
Tporv Solve for tangent point, vector

In the tangent plane projection, given the rectangular coordinates
of a star and its direction cosines, determine the direction
cosines of the tangent point.

Given:
    xi,eta   float64    rectangular coordinates of star image (Note 2)
    v        [3]float64 star's direction cosines (Note 3)

Returned:
    v01      [3]float64 tangent point's direction cosines, Solution 1
    v02      [3]float64 tangent point's direction cosines, Solution 2

Returned (function value):
    int     number of solutions:
            0 = no solutions returned (Note 4)
            1 = only the first solution is useful (Note 5)
            2 = both solutions are useful (Note 5)

Notes:

 1) The tangent plane projection is also called the "gnomonic
    projection" and the "central projection".

 2) The eta axis points due north in the adopted coordinate system.
    If the direction cosines represent observed (RA,Dec), the tangent
    plane coordinates (xi,eta) are conventionally called the
    "standard coordinates".  If the direction cosines are with
    respect to a right-handed triad, (xi,eta) are also right-handed.
    The units of (xi,eta) are, effectively, radians at the tangent
    point.

 3) The vector v must be of unit length or the result will be wrong.

 4) Cases where there is no solution can arise only near the poles.
    For example, it is clearly impossible for a star at the pole
    itself to have a non-zero xi value, and hence it is meaningless
    to ask where the tangent point would have to be.

 5) Also near the poles, cases can arise where there are two useful
    solutions.  The return value indicates whether the second of the
    two solutions returned is useful;  1 indicates only one useful
    solution, the usual case.

 6) The basis of the algorithm is to solve the spherical triangle
    PSC, where P is the north celestial pole, S is the star and C is
    the tangent point.  Calling the celestial spherical coordinates
    of the star and tangent point (a,b) and (a0,b0) respectively, and
    writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), and
    transforming the vector v into (a,b) in the normal way, side c is
    then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
    found) is (pi/2-b0), while angle C is given by sin(C) = xi/rho
    and cos(C) = eta/rho;  angle P (to be found) is (a-a0).  After
    solving the spherical triangle, the result (a0,b0) can be
    expressed in vector form as v0.

 7) This function is a member of the following set:

        spherical      vector         solve for

        iauTpxes      iauTpxev         xi,eta
        iauTpsts      iauTpstv          star
        iauTpors    > iauTporv <       origin

 References:

    Calabretta M.R. & Greisen, E.W., 2002, "Representations of
    celestial coordinates in FITS", Astron.Astrophys. 395, 1077

    Green, R.M., "Spherical Astronomy", Cambridge University Press,
    1987, Chapter 13.
*/
func Tporv(xi, eta float64, v [3]float64, v01, v02 *[3]float64) int {
	var x, y, z, rxy2, xi2, eta2p1, r, rsb, rcb, w2, w, c float64

	x = v[0]
	y = v[1]
	z = v[2]
	rxy2 = x*x + y*y
	xi2 = xi * xi
	eta2p1 = eta*eta + 1.0
	r = sqrt(xi2 + eta2p1)
	rsb = r * z
	rcb = r * sqrt(x*x+y*y)
	w2 = rcb*rcb - xi2
	if w2 > 0.0 {
		w = sqrt(w2)
		c = (rsb*eta + w) / (eta2p1 * sqrt(rxy2*(w2+xi2)))
		v01[0] = c * (x*w + y*xi)
		v01[1] = c * (y*w - x*xi)
		v01[2] = (rsb - eta*w) / eta2p1
		w = -w
		c = (rsb*eta + w) / (eta2p1 * sqrt(rxy2*(w2+xi2)))
		v02[0] = c * (x*w + y*xi)
		v02[1] = c * (y*w - x*xi)
		v02[2] = (rsb - eta*w) / eta2p1
		if fabs(rsb) < 1.0 {
			return 1
		} else {
			return 2
		}
	} else {
		return 0
	}
}
