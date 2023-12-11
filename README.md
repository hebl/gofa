<a href="https://pkg.go.dev/github.com/hebl/gofa"><img src="https://pkg.go.dev/badge/github.com/hebl/gofa.svg" alt="Go Reference"></a>

# GoFA

GoFA is a pure Golang derived from the International Astronomical
Union's (*IAU*) "`Standards Of Fundamental Astronomy (SOFA)`"  library
official C release [http://iausofa.org](http://iausofa.org).

GoFA documents are in: <https://pkg.go.dev/github.com/hebl/gofa>

GoFA functions are just translated from SOFA ANSI C functions.

Current version is 1.19. It is based on SOFA version `19` ([2023-10-11](http://iausofa.org/2023_1011_C/))

Examples in manuals are in the subdirectory: [examples](examples)

All function are tested. The [test](test) functions are derived from `t_sofa_c.c` .

## Installation

```shell
go get -u github.com/hebl/gofa
```

or

```shell
go install github.com/hebl/gofa@latest
```

## Functions

### Vector/Matrix Library [vml.go](vml.go)

- Initialization (4)
- Copy/Extend/Extract (5)
- Build rotations (3)
- Operations on vectors (17)
- Operations on matrices (2)
- Matrix-vector products (4)
- Rotation vectors (2)

### Angle [angle.go](angle.go)

- Spherical/Cartesian conversions (6)
- Separation and position-angle (4)
- Operations on angles (8)

### Astrometry (38) [astrometry.go](astrometry.go)

### Calendars (7) [jd.go](jd.go)

### Time Scales (20) [ts.go](ts.go)

### Coordinates [coord.go](coord.go)

- Ecliptic Coordinates (6)
- Galactic Coordinates (2)
- Horizon/Equatorial Coordinates (3)
- Geocentric/Geodetic Transformations (5)

### Earth Rotation and Sidereal Time (15) [erast.go](erast.go)

### Ephemerides (3) [ephem.go](ephem.go)

### Fundamental Arguments (14) [fundargs.go](fundargs.go)

### Gnomonic Projections (6) [projection.go](projection.go)

### Precession/Nutation/Polar Motion (64) [pn.go](pn.go)

### Star Catalog Conversions (9) [catalog.go](catalog.go)

## Versions

### v1.19

Version 1.19 include 192 routines for astronomy library, 55 routines for vector/matrix library.

## License

- [SOFA License](sofa_copyr.txt)

2023-10-23
