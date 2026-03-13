# Gets cycle numbers from a phase variable

Given a phase variable that increases, each time the phase passes
through 2 k pi, a new cycle starts. This function unwraps the phase, so
that it increases steadily, then takes the floor of the phase divided by
2 pi (or another modulus), so that it gets an integer for each cycle.

## Usage

``` r
get_cycle_numbers(
  ph,
  unwrap = FALSE,
  mod = 2 * pi,
  exclude_partial_cycles = TRUE
)
```

## Arguments

- ph:

  Phase variable

- unwrap:

  (TRUE or FALSE) Unwrap the phase variable. Note that the function will
  not work unless the phase is unwrapped, so you should only set unwrap
  to FALSE if the phase has been unwrapped earlier.

- mod:

  Modulus for the phase. Default is `2*pi`.

- exclude_partial_cycles:

  (TRUE or FALSE) Exclude cycles in which the phase does not advance
  from close to 0 to close to 2pi.

## Value

Integer cycle numbers with the same size as `ph`

## Details

Optionally, it will try to exclude partial cycles, setting the cycle
number to NA for cycles that do not start at a phase close to 0 and
progress to a phase close to 2 pi. "Close" here is defined based on the
average change in phase.

## See also

[`get_body_cycle_numbers_df()`](https://tytell.github.io/fishmechr/reference/get_body_cycle_numbers_df.md)

## Examples

``` r
# example phase that advances by slightly more than 3 cycles, modulo 2pi
ph <- seq(0, 20, by = pi/10) %% (2*pi)
get_cycle_numbers(ph, unwrap=TRUE)
#>  [1]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1
#> [26]  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2
#> [51]  2  2  2  2  2  2  2  2  2  2 NA NA NA NA
```
