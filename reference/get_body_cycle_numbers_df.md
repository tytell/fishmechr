# Gets oscillation cycle numbers for a midline data set

Extracts the phase for a specific point along the body (defined by
`pointval`) and uses
[`get_cycle_numbers()`](https://tytell.github.io/fishmechr/reference/get_cycle_numbers.md)
to get the cycle number. Joins the cycle number data with the original
data set so that cycle number is defined for each point along the body

## Usage

``` r
get_body_cycle_numbers_df(
  .data,
  ph,
  pointval,
  .frame = frame,
  .point = point,
  .out = NULL,
  overwrite = TRUE,
  ...
)
```

## Arguments

- .data:

  Data frame containing the midlines

- ph:

  Phase variable

- pointval:

  Specific point on the body to use to define the phase and and cycle
  number

- .frame, .point:

  Columns identifying frames and points (defaults are `frame` and
  `point`)

- .out:

  Name of the output column. Needs to have one element specifying the
  name for the cycle variable, or as a named list with an element named
  `cycle`). Default is 'cycle'

- overwrite:

  TRUE or FALSE to overwrite existing columns, if present.

- ...:

  Extra parameters to supply to
  [`get_cycle_numbers()`](https://tytell.github.io/fishmechr/reference/get_cycle_numbers.md)

## Value

A data frame containing a new column with the cycle numbers named
'cycle' or as specified in .out.

## Examples

``` r
library(dplyr)
# get phase at each frame for one point, then label each complete tail-beat cycle
# this is a minimal example; check the vignette for more details. In particular,
# you generally need to do much more smoothing on the data set for it to give
# good output
lampreydata |>
   group_by(frame) |>
   mutate(arclen = arclength(mxmm, mymm),
          curve_ang = curvature(arclen, mxmm, mymm)) |>
   group_by(point) |>
   mutate(ph_c = hilbert_phase(curve_ang)) |>
   ungroup() |>
   get_body_cycle_numbers_df(ph_c, pointval = 18)
#> Warning: There were 21 warnings in `mutate()`.
#> The first warning was:
#> ℹ In argument: `ph_c = hilbert_phase(curve_ang)`.
#> ℹ In group 2: `point = 2`.
#> Caused by warning in `hilbert_phase()`:
#> ! Phase seems to go backwards a lot, which may indicate an overly noisy signal
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 20 remaining warnings.
#> # A tibble: 1,600 × 9
#>        t frame point  mxmm  mymm arclen curve_ang  ph_c cycle
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl>     <dbl> <dbl> <dbl>
#>  1  0.02     1     1    NA    NA     NA        NA    NA    NA
#>  2  0.02     1     2    NA    NA     NA        NA    NA    NA
#>  3  0.02     1     3    NA    NA     NA        NA    NA    NA
#>  4  0.02     1     4    NA    NA     NA        NA    NA    NA
#>  5  0.02     1     5    NA    NA     NA        NA    NA    NA
#>  6  0.02     1     6    NA    NA     NA        NA    NA    NA
#>  7  0.02     1     7    NA    NA     NA        NA    NA    NA
#>  8  0.02     1     8    NA    NA     NA        NA    NA    NA
#>  9  0.02     1     9    NA    NA     NA        NA    NA    NA
#> 10  0.02     1    10    NA    NA     NA        NA    NA    NA
#> # ℹ 1,590 more rows
```
