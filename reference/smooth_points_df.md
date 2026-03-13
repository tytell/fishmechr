# Smooths locations of points over time

Smooths columns specified in `cols` for each individual point over time,
using a smoothing spline. Optionally fills in gaps of less than
`fillgaps` frames.

## Usage

``` r
smooth_points_df(
  .data,
  cols,
  spar,
  .out = NULL,
  .frame = frame,
  .point = point,
  fillgaps = 0
)
```

## Arguments

- .data:

  Data frame containing the midlines.

- cols:

  Columns containing the components to be smoothed. Often these will be
  the x and y coordinates of the midline.

- spar:

  Smoothing parameter passed to
  [`smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html).
  Values range from 0 (no smoothing) to 1 (heavy smoothing).

- .out:

  Names of the output columns. Should either be a list with the same
  number of elements as cols, or a glue specification as in `across` for
  the `.names` parameter. The default (NULL) means that the output
  columns will have the same name as the original column with an 's'
  appended at the end.

- .frame, .point:

  Columns identifying frames and points (defaults are `frame` and
  `point`)

- fillgaps:

  Longest gap to interpolate over. default is 0, which means not to fill
  gaps

## Value

A data frame containing the smoothed points

## Examples

``` r
library(dplyr)
# smooth x and y coordinates of the lamprey midline over time
lampreydata |>
  smooth_points_df(c(mxmm, mymm), spar = 0.6)
#> # A tibble: 1,600 × 8
#> # Groups:   point [20]
#>        t frame point  mxmm  mymm gaplen mxmms mymms
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl> <dbl> <dbl>
#>  1  0.02     1     1    NA    NA      0  425.  137.
#>  2  0.02     1     2    NA    NA      0  433.  136.
#>  3  0.02     1     3    NA    NA      0  441.  137.
#>  4  0.02     1     4    NA    NA      0  449.  137.
#>  5  0.02     1     5    NA    NA      0  457.  138.
#>  6  0.02     1     6    NA    NA      0  464.  141.
#>  7  0.02     1     7    NA    NA      0  471.  146.
#>  8  0.02     1     8    NA    NA      0  477.  152.
#>  9  0.02     1     9    NA    NA      0  485.  156.
#> 10  0.02     1    10    NA    NA      0  493.  157.
#> # ℹ 1,590 more rows
```
