# Calculate arc length along a 2D curve

Computes the straight-line distance between points on the curve and then
adds them up to get the arc length.

## Usage

``` r
arclength(x, y, na.skip = FALSE)
```

## Arguments

- x, y:

  Coordinates of the curve.

- na.skip:

  (TRUE or FALSE) If TRUE, skip NAs and compute arc length for the
  non-NA points, returning NA for the NA positions. If FALSE (default),
  return all NAs if any input point is NA.

## Value

Arc length along the curve.

## Examples

``` r
# compute arc length in each frame from the lamprey data set
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
lampreydata |>
  group_by(frame) |>
  mutate(arclen = arclength(mxmm, mymm))
#> # A tibble: 1,600 × 6
#> # Groups:   frame [80]
#>        t frame point  mxmm  mymm arclen
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl>
#>  1  0.02     1     1    NA    NA     NA
#>  2  0.02     1     2    NA    NA     NA
#>  3  0.02     1     3    NA    NA     NA
#>  4  0.02     1     4    NA    NA     NA
#>  5  0.02     1     5    NA    NA     NA
#>  6  0.02     1     6    NA    NA     NA
#>  7  0.02     1     7    NA    NA     NA
#>  8  0.02     1     8    NA    NA     NA
#>  9  0.02     1     9    NA    NA     NA
#> 10  0.02     1    10    NA    NA     NA
#> # ℹ 1,590 more rows
```
