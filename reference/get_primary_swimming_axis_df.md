# Gets the primary swimming axis for many midlines

Processes midlines from many frames of a video

## Usage

``` r
get_primary_swimming_axis_df(
  .data,
  tm,
  x,
  y,
  cutoff = NULL,
  overwrite = TRUE,
  .out = NULL,
  .frame = frame,
  .point = point,
  check_reasonableness = TRUE,
  na.rm = FALSE
)
```

## Arguments

- .data:

  Data frame containing the midline data

- tm:

  Column containing the time data. If a cutoff frequency is passed in,
  then this variable will be used to get the sampling frequency.

- x, y:

  Columns containing the x and y coordinates of each point along the
  midline.

- cutoff:

  (optional) If this parameter is included, smooth the swimming axis
  data with a low-pass filter with a cutoff at this frequency.

- overwrite:

  (default TRUE). Overwrite output columns if they exist

- .out:

  Names of the output columns. Needs to have four elements specifying
  the names for the x and y coordinates of the swim axis and the
  parallel and perpendicular components of the excursion, in that order.
  Or it can be a named list containing at least some of the elements
  `swimaxis_x`, `swimaxis_y`, `exc_x`, `exc`, in any order. If the
  return elements aren't in the named list, the defaults are
  'swimaxis_x', 'swimaxis_y', 'exc_x', and 'exc')

- .frame, .point:

  Columns identifying frames and points (defaults are `frame` and
  `point`)

- check_reasonableness:

  (default TRUE). Run some checks that the data are reasonable before
  processing.

- na.rm:

  (default FALSE) Remove NA points before computing the SVD

## Value

A data frame containing the original variables along with

- XX_ctr,YY_ctr: The center of each midline at each time, where XX and
  YY are the original names of the x and y coordinates.

- exc,exc_x: The new midlines centered and projected on to the swimming
  direction and the perpendicular axis. `exc` is useful as the lateral
  excursion of the swimming undulation.

## Details

Uses
[`get_primary_swimming_axis()`](https://tytell.github.io/fishmechr/reference/get_primary_swimming_axis.md)
to compute the swimming axis for a midline. Then optionally smooths the
axis using a Butterworth filter, and then projects the midlines on to
the new time-varying axes.

## Examples

``` r
library(dplyr)
# subtract the center of mass, then estimate the primary swimming axis
lampreydata |>
  group_by(frame) |>
  mutate(arclen = arclength(mxmm, mymm)) |>
  ungroup() |>
  get_midline_center_df(arclen, mxmm, mymm) |>
  mutate(mxmm_ctr = mxmm - xcom, mymm_ctr = mymm - ycom) |>
  get_primary_swimming_axis_df(t, mxmm_ctr, mymm_ctr)
#> ℹ Estimating center of mass as the centroid of x and y
#> # A tibble: 1,600 × 15
#>        t frame point  mxmm  mymm arclen  xcom  ycom  nsum mxmm_ctr mymm_ctr
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl> <dbl> <dbl> <int>    <dbl>    <dbl>
#>  1  0.02     1     1    NA    NA     NA   NaN   NaN     0       NA       NA
#>  2  0.02     1     2    NA    NA     NA   NaN   NaN     0       NA       NA
#>  3  0.02     1     3    NA    NA     NA   NaN   NaN     0       NA       NA
#>  4  0.02     1     4    NA    NA     NA   NaN   NaN     0       NA       NA
#>  5  0.02     1     5    NA    NA     NA   NaN   NaN     0       NA       NA
#>  6  0.02     1     6    NA    NA     NA   NaN   NaN     0       NA       NA
#>  7  0.02     1     7    NA    NA     NA   NaN   NaN     0       NA       NA
#>  8  0.02     1     8    NA    NA     NA   NaN   NaN     0       NA       NA
#>  9  0.02     1     9    NA    NA     NA   NaN   NaN     0       NA       NA
#> 10  0.02     1    10    NA    NA     NA   NaN   NaN     0       NA       NA
#> # ℹ 1,590 more rows
#> # ℹ 4 more variables: swimaxis_x <dbl>, swimaxis_y <dbl>, exc_x <dbl>,
#> #   exc <dbl>
```
