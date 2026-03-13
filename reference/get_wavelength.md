# Computes the body wavelength based on the phase at each point and the arc length

Computes the spatial derivative in phase relative to arc length using
several different possible methods:

- 'deriv' Computes the derivative directly using
  [`deriv()`](https://tytell.github.io/fishmechr/reference/deriv.md)

- 'slope' Fits a line to the phase relative to arc length and uses the
  slope of that line as an estimate of the wavelength

- 'cycle' Looks for pairs of points along the body where the phase
  differs by a full cycle. The arc length between those points is one
  wavelength.

- 'halfcycle' Looks for pairs of points along the body where the phase
  differs by one half cycle. The arc length between those points is half
  of a wavelength.

## Usage

``` r
get_wavelength(
  arclen,
  ph,
  unwrap = TRUE,
  method = "deriv",
  ignore_arclen_vals = NULL,
  sort_arclen = FALSE,
  check_reasonableness = TRUE,
  mod = 2 * pi,
  traveling_wave_dir = -1
)
```

## Arguments

- arclen:

  Arc length

- ph:

  Phase

- unwrap:

  (TRUE or FALSE) Unwrap the phase along the body

- method:

  ('deriv', 'slope', 'cycle', 'halfcycle') as explained above

- ignore_arclen_vals:

  NULL or a function that returns TRUE or FALSE for certain values of
  arclength where the phase estimate is not reliable. This is often used
  to exclude points near the head (e.g.,
  `ignore_arclen_vals = \(s) s < 0.3`)

- sort_arclen:

  (TRUE or FALSE) Sort the phase values by arc length before computing
  the wavelength. Useful if arc lengths arrive out of order. (default =
  FALSE)

- check_reasonableness:

  (TRUE or FALSE) Check that the phase decreases along the body as
  expected for a head-to-tail traveling wave, and warn if it does not.
  (default = TRUE)

- mod:

  Modulus for the phase variable

- traveling_wave_dir:

  (1 or -1) Defines the direction of the traveling wave. For a normal
  head-to-tail traveling wave, use -1 (default)

## Value

The wavelength as a vector the same size as the phase variable

## Examples

``` r
s <- seq(0, 1, by=0.1)
# artificial data with a wavelength of exactly 0.6
ph <- 2*pi*((1 - s) / 0.6)

get_wavelength(s, ph, ignore_arclen_vals = \(s) s < 0.3)
#>  [1]  NA  NA  NA 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6
```
