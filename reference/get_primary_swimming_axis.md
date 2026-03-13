# Gets the main swimming axis from a midline

Computes the main swimming axis of a midline as a unit vector, using the
singular value decomposition
([`svd()`](https://rdrr.io/r/base/svd.html)). This only works well if
the midlines are centered around zero, so it optionally subtracts off
the mean of x and y. For more sophisticated centering algorithms, see
[`get_midline_center_df()`](https://tytell.github.io/fishmechr/reference/get_midline_center_df.md).

## Usage

``` r
get_primary_swimming_axis(x, y, center = TRUE, na.rm = FALSE)
```

## Arguments

- x, y:

  Coordinates of the midline

- center:

  (TRUE or FALSE) Subtract the mean from the x and y coordinates

- na.rm:

  (default FALSE) Remove NA points before computing the SVD

## Value

A data frame with the following columns:

- `swimaxis_x`, `swimaxis_y` x and y components of the swimming axis
  vector

- `swimaxis_xctr`, `swimaxis_yctr` Mean x and y values that were
  subtracted before running the SVD

## See also

[`get_midline_center_df()`](https://tytell.github.io/fishmechr/reference/get_midline_center_df.md)

## Examples

``` r
library(dplyr)
library(tidyr)
# run the algorithm across multiple midlines at different times
lampreydata |>
  group_by(t) |>
  summarize(swimaxis = get_primary_swimming_axis(mxmm, mymm)) |>
  unnest(swimaxis)
#> # A tibble: 80 × 5
#>        t swimaxis_x swimaxis_y swimaxis_xctr swimaxis_yctr
#>    <dbl>      <dbl>      <dbl>         <dbl>         <dbl>
#>  1  0.02     NA         NA               NA            NA 
#>  2  0.04     NA         NA               NA            NA 
#>  3  0.06      0.989      0.150          487.          147.
#>  4  0.08      0.990      0.143          482.          146.
#>  5  0.1       0.991      0.136          477.          145.
#>  6  0.12      0.991      0.130          473.          143.
#>  7  0.14      0.993      0.121          467.          142.
#>  8  0.16      0.993      0.118          463.          141.
#>  9  0.18      0.993      0.116          458.          140.
#> 10  0.2       0.993      0.117          452.          139.
#> # ℹ 70 more rows
```
