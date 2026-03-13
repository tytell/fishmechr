# Interpolates and smooths a 2D curve at new arc length

For a 2D curve with (x,y) coordinates parameterized by the arc length,
interpolate new (x,y) coordinates at new arc lengths. Smooth the input
data with a smoothing spline.

## Usage

``` r
interpolate_points_df(
  .data,
  arclen,
  x,
  y,
  arclen_out = NULL,
  spar = 0.8,
  tailmethod = "extrapolate",
  fill_gaps = 0,
  .suffix = "_s",
  .out = NULL,
  overwrite = TRUE,
  .frame = frame,
  .point = point
)
```

## Arguments

- .data:

  Data frame

- arclen:

  Name of the input arc length column in `.data`

- x, y:

  Name of the columns that contain the coordinates of points on the
  curve

- arclen_out:

  Vector containing the new arc length

- spar:

  Smoothing parameter (ranges from 0 for no smoothing to 1 for high
  smoothing; see
  [`smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html) for
  more details.)

- tailmethod:

  ('keep', 'extrapolate', or 'NA') Methods to estimate the position of
  the tail tip if the last value of `arclength_out` is longer than the
  maximum arc length in the current frame.

  - 'keep' to keep the existing tail point, even if it is not at the
    requested arc length

  - 'extrapolate' to extrapolate a tail tip position, assuming that the
    curve continues straight

  - 'NA' to use replace the tail point with NA in this case.

- fill_gaps:

  Fill internal missing points of this size or smaller. (0, default,
  means no filling; 1 means to fill single missing points)

- .suffix:

  (default = '\_s') Suffix to append to the names of the arclen, x, and
  y columns after smoothing and interpolation.

- .out:

  Names of the output columns. Defaults are (arclen = 'arclen_s', xs =
  'xs', ys = 'ys'). Overrides the `.suffix` parameter if it is included.

- overwrite:

  TRUE or FALSE to overwrite existing columns

- .frame:

  Name of the frame variable in the data frame

- .point:

  Name of the point variable in the data frame

## Value

A data frame with updated columns for the smoothed and iterpolated arc
length, x and y coordinates.

## Details

Operates on each frame (as defined in the `.frame` parameter)
individually.

## Examples

``` r
library(dplyr)
# compute arc length, then interpolate all midlines to 20 evenly-spaced points
lampreydata |>
  group_by(frame) |>
  mutate(arclen = arclength(mxmm, mymm)) |>
  ungroup() |>
  interpolate_points_df(arclen, mxmm, mymm)
#> # A tibble: 1,600 × 9
#> # Groups:   frame [80]
#>        t frame point  mxmm  mymm arclen arclen_s mxmm_s mymm_s
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl>    <dbl>  <dbl>  <dbl>
#>  1  0.02     1     1    NA    NA     NA     0        NA     NA
#>  2  0.02     1     2    NA    NA     NA     8.13     NA     NA
#>  3  0.02     1     3    NA    NA     NA    16.3      NA     NA
#>  4  0.02     1     4    NA    NA     NA    24.4      NA     NA
#>  5  0.02     1     5    NA    NA     NA    32.5      NA     NA
#>  6  0.02     1     6    NA    NA     NA    40.7      NA     NA
#>  7  0.02     1     7    NA    NA     NA    48.8      NA     NA
#>  8  0.02     1     8    NA    NA     NA    56.9      NA     NA
#>  9  0.02     1     9    NA    NA     NA    65.1      NA     NA
#> 10  0.02     1    10    NA    NA     NA    73.2      NA     NA
#> # ℹ 1,590 more rows
```
