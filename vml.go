// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Vector/Matrix Library

// 1. Initialization (4)

// Operations involving p-vectors and r-matrices (3)

/*
Zp Zero a p-vector.

Returned:
    p        [3]float64      zero p-vector
*/
func Zp(p *[3]float64) {
	p[0] = 0.0
	p[1] = 0.0
	p[2] = 0.0
}

/*
Zr Initialize an r-matrix to the null matrix.

Returned:
    r        [3][3]float64    r-matrix
*/
func Zr(r *[3][3]float64) {
	r[0][0] = 0.0
	r[0][1] = 0.0
	r[0][2] = 0.0
	r[1][0] = 0.0
	r[1][1] = 0.0
	r[1][2] = 0.0
	r[2][0] = 0.0
	r[2][1] = 0.0
	r[2][2] = 0.0
}

/*
Ir Initialize an r-matrix to the identity matrix.

Returned:
    r    [3][3]float64    r-matrix
*/
func Ir(r *[3][3]float64) {
	r[0][0] = 1.0
	r[0][1] = 0.0
	r[0][2] = 0.0
	r[1][0] = 0.0
	r[1][1] = 1.0
	r[1][2] = 0.0
	r[2][0] = 0.0
	r[2][1] = 0.0
	r[2][2] = 1.0
}

/*
Zpv  Zero a pv-vector.

Returned:
    pv       [2][3]float64      zero pv-vector

Called:
    Zp        zero p-vector
*/
func Zpv(pv *[2][3]float64) {
	Zp(&pv[0])
	Zp(&pv[1])
}

// 2. Copy/Extend/Extract (5)

// Operations involving p-vectors and r-matrices (2)

/*
Cp Copy a p-vector.

Given:
    p        [3]float64     p-vector to be copied

Returned:
    c        [3]float64     copy
*/
func Cp(p [3]float64, c *[3]float64) {
	c[0] = p[0]
	c[1] = p[1]
	c[2] = p[2]
}

/*
Cr Copy an r-matrix.

Given:
    r        [3][3]float64    r-matrix to be copied

Returned:
    c        [3][3]float64    copy

Called:
    Cp        copy p-vector
*/
func Cr(r [3][3]float64, c *[3][3]float64) {
	Cp(r[0], &c[0])
	Cp(r[1], &c[1])
	Cp(r[2], &c[2])
}

// Operations involving pv-vectors (3)

/*
Cpv Copy a position/velocity vector.

Given:
    pv     [2][3]float64    position/velocity vector to be copied

Returned:
    c      [2][3]float64    copy
*/
func Cpv(pv [2][3]float64, c *[2][3]float64) {
	Cp(pv[0], &c[0])
	Cp(pv[1], &c[1])
}

/*
P2pv Extend a p-vector to a pv-vector by appending a zero velocity.

Given:
    p        [3]float64       p-vector

Returned:
    pv       [2][3]float64    pv-vector

Called:
    Cp        copy p-vector
    Zp        zero p-vector
*/
func P2pv(p [3]float64, pv *[2][3]float64) {
	Cp(p, &pv[0])
	Zp(&pv[1])
}

/*
Pv2p Discard velocity component of a pv-vector.

Given:
    pv      [2][3]float64     pv-vector

Returned:
    p       [3]float64        p-vector

Called:
    Cp       copy p-vector
*/
func Pv2p(pv [2][3]float64, p *[3]float64) {
	Cp(pv[0], p)
}

// 3. Build Rotations (3)

/*
Rx Rotate an r-matrix about the x-axis.

Given:
    phi    float64          angle (radians)

Given and returned:
    r      [3][3]float64    r-matrix, rotated

Notes:

 1) Calling this function with positive phi incorporates in the
    supplied r-matrix r an additional rotation, about the x-axis,
    anticlockwise as seen looking towards the origin from positive x.

 2) The additional rotation can be represented by this matrix:

        (  1        0            0      )
        (                               )
        (  0   + cos(phi)   + sin(phi)  )
        (                               )
        (  0   - sin(phi)   + cos(phi)  )
*/
func Rx(phi float64, r *[3][3]float64) {
	var s, c, a10, a11, a12, a20, a21, a22 float64

	s = sin(phi)
	c = cos(phi)

	a10 = c*r[1][0] + s*r[2][0]
	a11 = c*r[1][1] + s*r[2][1]
	a12 = c*r[1][2] + s*r[2][2]
	a20 = -s*r[1][0] + c*r[2][0]
	a21 = -s*r[1][1] + c*r[2][1]
	a22 = -s*r[1][2] + c*r[2][2]

	r[1][0] = a10
	r[1][1] = a11
	r[1][2] = a12
	r[2][0] = a20
	r[2][1] = a21
	r[2][2] = a22
}

/*
Ry Rotate an r-matrix about the y-axis.

Given:
    theta  float64          angle (radians)

Given and returned:
    r      [3][3]float64    r-matrix, rotated

Notes:

 1) Calling this function with positive theta incorporates in the
    supplied r-matrix r an additional rotation, about the y-axis,
    anticlockwise as seen looking towards the origin from positive y.

 2) The additional rotation can be represented by this matrix:

        (  + cos(theta)     0      - sin(theta)  )
        (                                        )
        (       0           1           0        )
        (                                        )
        (  + sin(theta)     0      + cos(theta)  )
*/
func Ry(theta float64, r *[3][3]float64) {
	var s, c, a00, a01, a02, a20, a21, a22 float64

	s = sin(theta)
	c = cos(theta)

	a00 = c*r[0][0] - s*r[2][0]
	a01 = c*r[0][1] - s*r[2][1]
	a02 = c*r[0][2] - s*r[2][2]
	a20 = s*r[0][0] + c*r[2][0]
	a21 = s*r[0][1] + c*r[2][1]
	a22 = s*r[0][2] + c*r[2][2]

	r[0][0] = a00
	r[0][1] = a01
	r[0][2] = a02
	r[2][0] = a20
	r[2][1] = a21
	r[2][2] = a22
}

/*
Rz Rotate an r-matrix about the z-axis.

Given:
    psi    float64          angle (radians)

Given and returned:
    r      [3][3]float64    r-matrix, rotated

Notes:

 1) Calling this function with positive psi incorporates in the
    supplied r-matrix r an additional rotation, about the z-axis,
    anticlockwise as seen looking towards the origin from positive z.

 2) The additional rotation can be represented by this matrix:

        (  + cos(psi)   + sin(psi)     0  )
        (                                 )
        (  - sin(psi)   + cos(psi)     0  )
        (                                 )
        (       0            0         1  )
*/
func Rz(psi float64, r *[3][3]float64) {
	var s, c, a00, a01, a02, a10, a11, a12 float64

	s = sin(psi)
	c = cos(psi)

	a00 = c*r[0][0] + s*r[1][0]
	a01 = c*r[0][1] + s*r[1][1]
	a02 = c*r[0][2] + s*r[1][2]
	a10 = -s*r[0][0] + c*r[1][0]
	a11 = -s*r[0][1] + c*r[1][1]
	a12 = -s*r[0][2] + c*r[1][2]

	r[0][0] = a00
	r[0][1] = a01
	r[0][2] = a02
	r[1][0] = a10
	r[1][1] = a11
	r[1][2] = a12

}

// 4. Spherical/Cartesian Conversions (6)

// Operations involving p-vectors and r-matrices (4)

/*
S2c Convert spherical coordinates to Cartesian.

Given:
    theta    float64       longitude angle (radians)
    phi      float64       latitude angle (radians)

Returned:
    c        [3]float64    direction cosines
*/
func S2c(theta, phi float64, c *[3]float64) {
	var cp float64

	cp = cos(phi)
	c[0] = cos(theta) * cp
	c[1] = sin(theta) * cp
	c[2] = sin(phi)
}

/*
C2s P-vector to spherical coordinates.

Given:
    p      [3]float64    p-vector

Returned:
    theta  float64       longitude angle (radians)
    phi    float64       latitude angle (radians)

Notes:

 1) The vector p can have any magnitude; only its direction is used.

 2) If p is null, zero theta and phi are returned.

 3) At either pole, zero theta is returned.
*/
func C2s(p [3]float64, theta *float64, phi *float64) {
	var x, y, z, d2 float64

	x = p[0]
	y = p[1]
	z = p[2]
	d2 = x*x + y*y

	if d2 == 0.0 {
		*theta = 0.0
	} else {
		*theta = atan2(y, x)
	}

	if z == 0.0 {
		*phi = 0.0
	} else {
		*phi = atan2(z, sqrt(d2))
	}
}

/*
S2p Convert spherical polar coordinates to p-vector.

Given:
    theta   float64       longitude angle (radians)
    phi     float64       latitude angle (radians)
    r       float64       radial distance

Returned:
    p       [3]float64    Cartesian coordinates

Called:
    S2c       spherical coordinates to unit vector
    Sxp       multiply p-vector by scalar
*/
func S2p(theta, phi, r float64, p *[3]float64) {
	u := [3]float64{}

	S2c(theta, phi, &u)
	Sxp(r, u, p)
}

/*
P2s P-vector to spherical polar coordinates.

Given:
    p        [3]float64    p-vector

Returned:
    theta    float64       longitude angle (radians)
    phi      float64       latitude angle (radians)
    r        float64       radial distance

Notes:

 1) If P is null, zero theta, phi and r are returned.

 2) At either pole, zero theta is returned.

Called:
    C2s       p-vector to spherical
    Pm        modulus of p-vector
*/
func P2s(p [3]float64, theta *float64, phi *float64, r *float64) {
	C2s(p, theta, phi)
	*r = Pm(p)
}

// Operations involving pv-vectors (2)

/*
S2pv Convert position/velocity from spherical to Cartesian coordinates.

Given:
    theta    float64          longitude angle (radians)
    phi      float64          latitude angle (radians)
    r        float64          radial distance
    td       float64          rate of change of theta
    pd       float64          rate of change of phi
    rd       float64          rate of change of r

Returned:
    pv       [2][3]float64    pv-vector
*/
func S2pv(theta, phi, r float64, td, pd, rd float64, pv *[2][3]float64) {
	var st, ct, sp, cp, rcp, x, y, rpd, w float64

	st = sin(theta)
	ct = cos(theta)
	sp = sin(phi)
	cp = cos(phi)
	rcp = r * cp
	x = rcp * ct
	y = rcp * st
	rpd = r * pd
	w = rpd*sp - cp*rd

	pv[0][0] = x
	pv[0][1] = y
	pv[0][2] = r * sp
	pv[1][0] = -y*td - w*ct
	pv[1][1] = x*td - w*st
	pv[1][2] = rpd*cp + sp*rd
}

/*
Pv2s Convert position/velocity from Cartesian to spherical coordinates.

Given:
    pv       [2][3]float64  pv-vector

Returned:
    theta    float64        longitude angle (radians)
    phi      float64        latitude angle (radians)
    r        float64        radial distance
    td       float64        rate of change of theta
    pd       float64        rate of change of phi
    rd       float64        rate of change of r

Notes:

 1) If the position part of pv is null, theta, phi, td and pd
    are indeterminate.  This is handled by extrapolating the
    position through unit time by using the velocity part of
    pv.  This moves the origin without changing the direction
    of the velocity component.  If the position and velocity
    components of pv are both null, zeroes are returned for all
    six results.

 2) If the position is a pole, theta, td and pd are indeterminate.
    In such cases zeroes are returned for all three.
*/
func Pv2s(pv [2][3]float64, theta, phi, r *float64, td, pd, rd *float64) {
	var x, y, z, xd, yd, zd, rxy2, rxy, r2, rtrue, rw, xyp float64

	/* Components of position/velocity vector. */
	x = pv[0][0]
	y = pv[0][1]
	z = pv[0][2]
	xd = pv[1][0]
	yd = pv[1][1]
	zd = pv[1][2]

	/* Component of r in XY plane squared. */
	rxy2 = x*x + y*y

	/* Modulus squared. */
	r2 = rxy2 + z*z

	/* Modulus. */
	rtrue = sqrt(r2)

	/* If null vector, move the origin along the direction of movement. */
	rw = rtrue
	if rtrue == 0.0 {
		x = xd
		y = yd
		z = zd
		rxy2 = x*x + y*y
		r2 = rxy2 + z*z
		rw = sqrt(r2)
	}

	/* Position and velocity in spherical coordinates. */
	rxy = sqrt(rxy2)
	xyp = x*xd + y*yd
	if rxy2 != 0.0 {
		*theta = atan2(y, x)
		*phi = atan2(z, rxy)
		*td = (x*yd - y*xd) / rxy2
		*pd = (zd*rxy2 - z*xyp) / (r2 * rxy)
	} else {
		*theta = 0.0
		if z != 0.0 {
			*phi = atan2(z, rxy)
		} else {
			*phi = 0.0
		}

		*td = 0.0
		*pd = 0.0
	}
	*r = rtrue

	if rw != 0.0 {
		*rd = (xyp + z*zd) / rw
	} else {
		*rd = 0.0
	}

}

// 5. Operations on Vectors (17)

// Operations involving p-vectors and r-matrices (8)

/*
Ppp P-vector addition.

Given:
    a        [3]float64      first p-vector
    b        [3]float64      second p-vector

Returned:
    apb      [3]float64      a + b

Note:
    It is permissible to re-use the same array for any of the
    arguments.
*/
func Ppp(a, b [3]float64, apb *[3]float64) {
	apb[0] = a[0] + b[0]
	apb[1] = a[1] + b[1]
	apb[2] = a[2] + b[2]
}

/*
Pmp P-vector subtraction.

Given:
    a        [3]float64      first p-vector
    b        [3]float64      second p-vector

Returned:
    amb      [3]float64      a - b

Note:
    It is permissible to re-use the same array for any of the
    arguments.
*/
func Pmp(a, b [3]float64, amb *[3]float64) {
	amb[0] = a[0] - b[0]
	amb[1] = a[1] - b[1]
	amb[2] = a[2] - b[2]
}

/*
Ppsp P-vector plus scaled p-vector.

Given:
    a      [3]float64     first p-vector
    s      float64        scalar (multiplier for b)
    b      [3]float64     second p-vector

Returned:
    apsb   [3]float64     a + s*b

Note:
    It is permissible for any of a, b and apsb to be the same array.

Called:
    Sxp       multiply p-vector by scalar
    Ppp       p-vector plus p-vector
*/
func Ppsp(a [3]float64, s float64, b [3]float64, apsb *[3]float64) {
	var sb [3]float64

	/* s*b. */
	Sxp(s, b, &sb)

	/* a + s*b. */
	Ppp(a, sb, apsb)

}

/*
Pdp p-vector inner (=scalar=dot) product.

Given:
    a      [3]float64     first p-vector
    b      [3]float64     second p-vector

Returned (function value):
    float64        a . b
*/
func Pdp(a, b [3]float64) float64 {
	w := a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
	return w
}

/*
Pxp p-vector outer (=vector=cross) product.

Given:
    a        [3]float64      first p-vector
    b        [3]float64      second p-vector

Returned:
    axb      [3]float64      a x b

Note:
    It is permissible to re-use the same array for any of the
    arguments.
*/
func Pxp(a, b [3]float64, axb *[3]float64) {
	var xa, ya, za, xb, yb, zb float64

	xa = a[0]
	ya = a[1]
	za = a[2]
	xb = b[0]
	yb = b[1]
	zb = b[2]
	axb[0] = ya*zb - za*yb
	axb[1] = za*xb - xa*zb
	axb[2] = xa*yb - ya*xb
}

/*
Pm Modulus of p-vector.

Given:
    p      [3]float64     p-vector

Returned (function value):
    float64        modulus
*/
func Pm(p [3]float64) float64 {
	return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
}

/*
Pn Convert a p-vector into modulus and unit vector.

Given:
    p        [3]float64      p-vector

Returned:
    r        float64         modulus
    u        [3]float64      unit vector

Notes:

 1) If p is null, the result is null.  Otherwise the result is a unit
    vector.

 2) It is permissible to re-use the same array for any of the
    arguments.

Called:
    Pm        modulus of p-vector
    Zp        zero p-vector
    Sxp       multiply p-vector by scalar
*/
func Pn(p [3]float64, r *float64, u *[3]float64) {
	var w float64

	/* Obtain the modulus and test for zero. */
	w = Pm(p)
	if w == 0.0 {
		/* Null vector. */
		Zp(u)
	} else {
		/* Unit vector. */
		Sxp(1.0/w, p, u)
	}

	/* Return the modulus. */
	*r = w
}

/*
Sxp Multiply a p-vector by a scalar.

Given:
    s      float64        scalar
    p      [3]float64     p-vector

Returned:
    sp     [3]float64     s * p

Note:
    It is permissible for p and sp to be the same array.
*/
func Sxp(s float64, p [3]float64, sp *[3]float64) {
	sp[0] = s * p[0]
	sp[1] = s * p[1]
	sp[2] = s * p[2]
}

// Operations involving pv-vectors (9)

/*
Pvppv Add one pv-vector to another.

Given:
    a        [2][3]float64      first pv-vector
    b        [2][3]float64      second pv-vector

Returned:
    apb      [2][3]float64      a + b

Note:
    It is permissible to re-use the same array for any of the
    arguments.

Called:
    Ppp       p-vector plus p-vector
*/
func Pvppv(a, b [2][3]float64, apb *[2][3]float64) {
	Ppp(a[0], b[0], &apb[0])
	Ppp(a[1], b[1], &apb[1])
}

/*
Pvmpv Subtract one pv-vector from another.

Given:
    a       [2][3]float64      first pv-vector
    b       [2][3]float64      second pv-vector

Returned:
    amb     [2][3]float64      a - b

Note:
    It is permissible to re-use the same array for any of the
    arguments.

Called:
    Pmp       p-vector minus p-vector
*/
func Pvmpv(a, b [2][3]float64, amb *[2][3]float64) {
	Pmp(a[0], b[0], &amb[0])
	Pmp(a[1], b[1], &amb[1])
}

/*
Pvdpv Inner (=scalar=dot) product of two pv-vectors.

Given:
    a        [2][3]float64      first pv-vector
    b        [2][3]float64      second pv-vector

Returned:
    adb      [2]float64         a . b (see note)

Note:

    If the position and velocity components of the two pv-vectors are
    ( ap, av ) and ( bp, bv ), the result, a . b, is the pair of
    numbers ( ap . bp , ap . bv + av . bp ).  The two numbers are the
    dot-product of the two p-vectors and its derivative.

Called:
    Pdp       scalar product of two p-vectors
*/
func Pvdpv(a, b [2][3]float64, adb *[2]float64) {
	var adbd, addb float64

	/* a . b = constant part of result. */
	adb[0] = Pdp(a[0], b[0])

	/* a . bdot */
	adbd = Pdp(a[0], b[1])

	/* adot . b */
	addb = Pdp(a[1], b[0])

	/* Velocity part of result. */
	adb[1] = adbd + addb
}

/*
Pvxpv Outer (=vector=cross) product of two pv-vectors.

Given:
    a        [2][3]float64      first pv-vector
    b        [2][3]float64      second pv-vector

Returned:
    axb      [2][3]float64      a x b

Notes:

 1) If the position and velocity components of the two pv-vectors are
    ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
    vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
    cross-product of the two p-vectors and its derivative.

 2) It is permissible to re-use the same array for any of the
    arguments.

Called:
    Cpv       copy pv-vector
    Pxp       vector product of two p-vectors
    Ppp       p-vector plus p-vector
*/
func Pvxpv(a, b [2][3]float64, axb *[2][3]float64) {
	var wa, wb [2][3]float64
	var axbd, adxb [3]float64

	/* Make copies of the inputs. */
	Cpv(a, &wa)
	Cpv(b, &wb)

	/* a x b = position part of result. */
	Pxp(wa[0], wb[0], &axb[0])

	/* a x bdot + adot x b = velocity part of result. */
	Pxp(wa[0], wb[1], &axbd)
	Pxp(wa[1], wb[0], &adxb)
	Ppp(axbd, adxb, &axb[1])
}

/*
Pvm Modulus of pv-vector.

Given:
	pv     [2][3]float64   pv-vector

Returned:
	r      float64         modulus of position component
	s      float64         modulus of velocity component

Called:
	Pm        modulus of p-vector
*/
func Pvm(pv [2][3]float64, r, s *float64) {
	/* Distance. */
	*r = Pm(pv[0])

	/* Speed. */
	*s = Pm(pv[1])
}

/*
Sxpv Multiply a pv-vector by a scalar.

Given:
	s       float64          scalar
	pv      [2][3]float64    pv-vector

Returned:
	spv     [2][3]float64    s * pv

Note:
	It is permissible for pv and spv to be the same array.

Called:
	S2xpv     multiply pv-vector by two scalars
*/
func Sxpv(s float64, pv [2][3]float64, spv *[2][3]float64) {
	S2xpv(s, s, pv, spv)
}

/*
S2xpv Multiply a pv-vector by two scalars.

Given:
	s1     float64         scalar to multiply position component by
	s2     float64         scalar to multiply velocity component by
	pv     [2][3]float64   pv-vector

Returned:
	spv    [2][3]float64   pv-vector: p scaled by s1, v scaled by s2

Note:
	It is permissible for pv and spv to be the same array.

Called:
	Sxp       multiply p-vector by scalar
*/
func S2xpv(s1, s2 float64, pv [2][3]float64, spv *[2][3]float64) {
	Sxp(s1, pv[0], &spv[0])
	Sxp(s2, pv[1], &spv[1])
}

/*
Pvu Update a pv-vector.

Given:
	dt       float64           time interval
	pv       [2][3]float64     pv-vector

Returned:
	upv      [2][3]float64     p updated, v unchanged

Notes:

 1) "Update" means "refer the position component of the vector
    to a new date dt time units from the existing date".

 2) The time units of dt must match those of the velocity.

 3) It is permissible for pv and upv to be the same array.

Called:
	Ppsp      p-vector plus scaled p-vector
	Cp        copy p-vector
*/
func Pvu(dt float64, pv [2][3]float64, upv *[2][3]float64) {
	Ppsp(pv[0], dt, pv[1], &upv[0])
	Cp(pv[1], &upv[1])
}

/*
Pvup Update a pv-vector, discarding the velocity component.

Status:  vector/matrix support function.

Given:
	dt       float64            time interval
	pv       [2][3]float64      pv-vector

Returned:
	p        [3]float64         p-vector

Notes:

 1) "Update" means "refer the position component of the vector to a
    new date dt time units from the existing date".

 2) The time units of dt must match those of the velocity.
*/
func Pvup(dt float64, pv [2][3]float64, p *[3]float64) {
	p[0] = pv[0][0] + dt*pv[1][0]
	p[1] = pv[0][1] + dt*pv[1][1]
	p[2] = pv[0][2] + dt*pv[1][2]
}

// 6. Operations on matrices (2)

/*
Rxr Multiply two r-matrices.

Given:
	a        [3][3]float64    first r-matrix
	b        [3][3]float64    second r-matrix

Returned:
	atb      [3][3]float64    a * b

Note:
	It is permissible to re-use the same array for any of the
	arguments.

Called:
	Cr        copy r-matrix
*/
func Rxr(a, b [3][3]float64, atb *[3][3]float64) {

	var w float64
	var wm [3][3]float64

	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			w = 0.0
			for k := 0; k < 3; k++ {
				w += a[i][k] * b[k][j]
			}
			wm[i][j] = w
		}
	}
	Cr(wm, atb)
}

/*
Tr Transpose an r-matrix.

Given:
	r        [3][3]float64    r-matrix

Returned:
	rt       [3][3]float64    transpose

Note:
	It is permissible for r and rt to be the same array.

Called:
	Cr        copy r-matrix
*/
func Tr(r [3][3]float64, rt *[3][3]float64) {
	var wm [3][3]float64

	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			wm[i][j] = r[j][i]
		}
	}
	Cr(wm, rt)
}

// 7. Matrix-vector products (4)

// Operations involving p-vectors and r-matrices (2)

/*
Rxp Multiply a p-vector by an r-matrix.

Given:
	r        [3][3]float64    r-matrix
	p        [3]float64       p-vector

Returned:
	rp       [3]float64       r * p

Note:
	It is permissible for p and rp to be the same array.

Called:
	Cp        copy p-vector
*/
func Rxp(r [3][3]float64, p [3]float64, rp *[3]float64) {
	var w float64
	var wrp [3]float64
	//var i, j int

	/* Matrix r * vector p. */
	for j := 0; j < 3; j++ {
		w = 0.0
		for i := 0; i < 3; i++ {
			w += r[j][i] * p[i]
		}
		wrp[j] = w
	}

	/* Return the result. */
	Cp(wrp, rp)
}

/*
Trxp Multiply a p-vector by the transpose of an r-matrix.

Given:
	r        [3][3]float64   r-matrix
	p        [3]float64      p-vector

Returned:
	trp      [3]float64      r^T * p

Note:
	It is permissible for p and trp to be the same array.

Called:
	Tr        transpose r-matrix
	Rxp       product of r-matrix and p-vector
*/
func Trxp(r [3][3]float64, p [3]float64, trp *[3]float64) {
	var tr [3][3]float64

	/* Transpose of matrix r. */
	Tr(r, &tr)

	/* Matrix tr * vector p -> vector trp. */
	Rxp(tr, p, trp)
}

// Operations involving pv-vectors (2)

/*
Rxpv Multiply a pv-vector by an r-matrix.

Given:
	r        [3][3]float64    r-matrix
	pv       [2][3]float64    pv-vector

Returned:
	rpv      [2][3]float64    r * pv

Notes:

 1) The algorithm is for the simple case where the r-matrix r is not
    a function of time.  The case where r is a function of time leads
    to an additional velocity component equal to the product of the
    derivative of r and the position vector.

 2) It is permissible for pv and rpv to be the same array.

Called:
	Rxp       product of r-matrix and p-vector
*/
func Rxpv(r [3][3]float64, pv [2][3]float64, rpv *[2][3]float64) {
	Rxp(r, pv[0], &rpv[0])
	Rxp(r, pv[1], &rpv[1])
}

/*
Trxpv Multiply a pv-vector by the transpose of an r-matrix.

Given:
	r        [3][3]float64    r-matrix
	pv       [2][3]float64    pv-vector

Returned:
	trpv     [2][3]float64    r^T * pv

Notes:

 1) The algorithm is for the simple case where the r-matrix r is not
    a function of time.  The case where r is a function of time leads
    to an additional velocity component equal to the product of the
    derivative of the transpose of r and the position vector.

 2) It is permissible for pv and rpv to be the same array.

Called:
	Tr        transpose r-matrix
	Rxpv      product of r-matrix and pv-vector
*/
func Trxpv(r [3][3]float64, pv [2][3]float64, trpv *[2][3]float64) {
	var tr [3][3]float64

	/* Transpose of matrix r. */
	Tr(r, &tr)

	/* Matrix tr * vector pv -> vector trpv. */
	Rxpv(tr, pv, trpv)
}

// 9. Rotation vectors (2)

/*
Rv2m Form the r-matrix corresponding to a given r-vector.

Given:
	w        [3]float64      rotation vector (Note 1)

Returned:
	r        [3][3]float64    rotation matrix

Notes:

 1) A rotation matrix describes a rotation through some angle about
    some arbitrary axis called the Euler axis.  The "rotation vector"
    supplied to This function has the same direction as the Euler
    axis, and its magnitude is the angle in radians.

 2) If w is null, the identity matrix is returned.

 3) The reference frame rotates clockwise as seen looking along the
    rotation vector from the origin.
*/
func Rv2m(w [3]float64, r *[3][3]float64) {
	var x, y, z, phi, s, c, f float64

	/* Euler angle (magnitude of rotation vector) and functions. */
	x = w[0]
	y = w[1]
	z = w[2]
	phi = sqrt(x*x + y*y + z*z)
	s = sin(phi)
	c = cos(phi)
	f = 1.0 - c

	/* Euler axis (direction of rotation vector), perhaps null. */
	if phi > 0.0 {
		x /= phi
		y /= phi
		z /= phi
	}

	/* Form the rotation matrix. */
	r[0][0] = x*x*f + c
	r[0][1] = x*y*f + z*s
	r[0][2] = x*z*f - y*s
	r[1][0] = y*x*f - z*s
	r[1][1] = y*y*f + c
	r[1][2] = y*z*f + x*s
	r[2][0] = z*x*f + y*s
	r[2][1] = z*y*f - x*s
	r[2][2] = z*z*f + c

}

/*
Rm2v Express an r-matrix as an r-vector.

Given:
	r        [3][3]float64    rotation matrix

Returned:
	w        [3]float64       rotation vector (Note 1)

Notes:

 1) A rotation matrix describes a rotation through some angle about
    some arbitrary axis called the Euler axis.  The "rotation vector"
    returned by this function has the same direction as the Euler axis,
    and its magnitude is the angle in radians.  (The magnitude and
    direction can be separated by means of the function iauPn.)

 2) If r is null, so is the result.  If r is not a rotation matrix
    the result is undefined;  r must be proper (i.e. have a positive
    determinant) and real orthogonal (inverse = transpose).

 3) The reference frame rotates clockwise as seen looking along
    the rotation vector from the origin.
*/
func Rm2v(r [3][3]float64, w *[3]float64) {
	var x, y, z, s2, c2, phi, f float64

	x = r[1][2] - r[2][1]
	y = r[2][0] - r[0][2]
	z = r[0][1] - r[1][0]
	s2 = sqrt(x*x + y*y + z*z)
	if s2 > 0 {
		c2 = r[0][0] + r[1][1] + r[2][2] - 1.0
		phi = atan2(s2, c2)
		f = phi / s2
		w[0] = x * f
		w[1] = y * f
		w[2] = z * f
	} else {
		w[0] = 0.0
		w[1] = 0.0
		w[2] = 0.0
	}
}
