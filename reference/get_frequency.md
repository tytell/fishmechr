# Estimates the cycle frequency based on time and phase

Estimates frequency based on time and phase by one of two methods:

- 'deriv' Takes the derivative of phase vs time to get frequency

- 'slope' Fits a line to phase vs time and uses the slope as an estimate
  of frequency.

## Usage

``` r
get_frequency(
  t,
  ph,
  unwrap = FALSE,
  method = "deriv",
  mod = 2 * pi,
  check_reasonableness = TRUE
)
```

## Arguments

- t:

  Time

- ph:

  Phase

- unwrap:

  (TRUE or FALSE) Unwrap the phase. Note that the function will not work
  unless the phase is unwrapped, so you should only set unwrap to FALSE
  if the phase has been unwrapped earlier.

- method:

  ('deriv' or 'slope') Method to estimate the frequency. 'deriv' takes
  the derivative of phase vs time; 'slope' fits a line to phase vs time
  and uses the slope.

- mod:

  Modulus for the phase. Default is `2*pi`.

- check_reasonableness:

  Runs some simple checks to make sure the data make sense. Checks to
  make sure phase is mostly increasing and warns if it isn't.

## Value

The frequency estimate, either as a vector the same size as `ph` (for
the 'deriv' algorithm) or as a single value (for the 'slope' algorithm)

## Examples

``` r
t <- seq(0, 3, by=0.1)
# example phase that has a frequency of exactly 2.3 Hz
ph <- 2*pi*(t * 2.3)

get_frequency(t, ph)
#>  [1] 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3
#> [20] 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3 2.3
```
