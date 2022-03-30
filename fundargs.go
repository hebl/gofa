// Copyright 2022 HE Boliang
// All rights reserved.

package gofa

// Fundamental Arguments (14)

/*
Fad03  mean elongation of the Moon from the Sun.

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    D, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Fad03(t float64) float64 {
	var a float64

	/* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
	a = fmod(1072260.703692+
		t*(1602961601.2090+
			t*(-6.3706+
				t*(0.006593+
					t*(-0.00003169)))), TURNAS) * DAS2R

	return a
}

/*
Fae03  mean longitude of Earth.

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Earth, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Fae03(t float64) float64 {
	var a float64

	/* Mean longitude of Earth (IERS Conventions 2003). */
	a = fmod(1.753470314+628.3075849991*t, D2PI)

	return a
}

/*
Faf03 mean longitude of the Moon minus mean longitude of the ascending node.

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    F, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Faf03(t float64) float64 {
	var a float64

	/* Mean longitude of the Moon minus that of the ascending node */
	/* (IERS Conventions 2003).                                    */
	a = fmod(335779.526232+
		t*(1739527262.8478+
			t*(-12.7512+
				t*(-0.001037+
					t*(0.00000417)))), TURNAS) * DAS2R

	return a
}

/*
Faju03 Mean longitude of Jupiter

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Jupiter, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Faju03(t float64) float64 {
	var a float64

	/* Mean longitude of Jupiter (IERS Conventions 2003). */
	a = fmod(0.599546497+52.9690962641*t, D2PI)

	return a
}

/*
Fal03    Mean anomaly of the Moon

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    l, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Fal03(t float64) float64 {
	var a float64

	/* Mean anomaly of the Moon (IERS Conventions 2003). */
	a = fmod(485868.249036+
		t*(1717915923.2178+
			t*(31.8792+
				t*(0.051635+
					t*(-0.00024470)))), TURNAS) * DAS2R

	return a
}

/*
Falp03 Mean anomaly of the Sun

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    l', radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Falp03(t float64) float64 {
	var a float64

	/* Mean anomaly of the Sun (IERS Conventions 2003). */
	a = fmod(1287104.793048+
		t*(129596581.0481+
			t*(-0.5532+
				t*(0.000136+
					t*(-0.00001149)))), TURNAS) * DAS2R

	return a
}

/*
Fama03 Mean longitude of Mars

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Mars, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Fama03(t float64) float64 {
	var a float64

	/* Mean longitude of Mars (IERS Conventions 2003). */
	a = fmod(6.203480913+334.0612426700*t, D2PI)

	return a
}

/*
Fame03 Mean longitude of Mercury

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Mercury, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Fame03(t float64) float64 {
	var a float64

	/* Mean longitude of Mercury (IERS Conventions 2003). */
	a = fmod(4.402608842+2608.7903141574*t, D2PI)

	return a
}

/*
Fane03 Mean longitude of Neptune

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Neptune, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is adapted from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Fane03(t float64) float64 {
	var a float64

	/* Mean longitude of Neptune (IERS Conventions 2003). */
	a = fmod(5.311886287+3.8133035638*t, D2PI)

	return a
}

/*
Faom03 Mean longitude of the Moon's ascending node

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    Omega, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683.
*/
func Faom03(t float64) float64 {
	var a float64

	/* Mean longitude of the Moon's ascending node */
	/* (IERS Conventions 2003).                    */
	a = fmod(450160.398036+
		t*(-6962890.5431+
			t*(7.4722+
				t*(0.007702+
					t*(-0.00005939)))), TURNAS) * DAS2R

	return a
}

/*
Fapa03 General accumulated precession in longitude

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    general precession in longitude, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003).  It
    is taken from Kinoshita & Souchay (1990) and comes originally
    from Lieske et al. (1977).

References:

    Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
    48, 187

    Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    Astron.Astrophys. 58, 1-16

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)
*/
func Fapa03(t float64) float64 {
	var a float64

	/* General accumulated precession in longitude. */
	a = (0.024381750 + 0.00000538691*t) * t

	return a
}

/*
Fasa03 Mean longitude of Saturn

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Saturn, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Fasa03(t float64) float64 {
	var a float64

	/* Mean longitude of Saturn (IERS Conventions 2003). */
	a = fmod(0.874016757+21.3299104960*t, D2PI)

	return a
}

/*
Faur03 Mean longitude of Uranus

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned  (function value):
    float64    mean longitude of Uranus, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    is adapted from Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*/
func Faur03(t float64) float64 {
	var a float64

	/* Mean longitude of Uranus (IERS Conventions 2003). */
	a = fmod(5.481293872+7.4781598567*t, D2PI)

	return a
}

/*
Fave03 Mean longitude of Venus

Fundamental argument, IERS Conventions (2003)

Given:
    t     float64    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
    float64    mean longitude of Venus, radians (Note 2)

Notes:

 1) Though t is strictly TDB, it is usually more convenient to use
    TT, which makes no significant difference.

 2) The expression used is as adopted in IERS Conventions (2003) and
    comes from Souchay et al. (1999) after Simon et al. (1994).

References:

    McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

    Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

    Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111
*/
func Fave03(t float64) float64 {
	var a float64

	/* Mean longitude of Venus (IERS Conventions 2003). */
	a = fmod(3.176146697+1021.3285546211*t, D2PI)

	return a
}
