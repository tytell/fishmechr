# Gets the center of a midline for many midlines in a data frame

Estimates the center of a midline based on mass distribution, volume
distribution, or body width.

## Usage

``` r
get_midline_center_df(
  .data,
  arclen,
  x,
  y,
  mass,
  width,
  height,
  .out = NULL,
  .frame = frame,
  .point = point,
  excludepoints = c(),
  method = "mutate",
  overwrite = TRUE
)
```

## Arguments

- .data:

  Data frame containing the midlines.

- arclen:

  The column containing the arc length. See
  [`arclength()`](https://tytell.github.io/fishmechr/reference/arclength.md)

- x, y:

  Columns containing the x and y coordinates of the midline. There
  should be N points.

- mass:

  (optional) Column containing the mass of each segment, with N-1
  segments.

- width:

  (optional) Column containing the horizontal plane width of the body at
  each midline point (N points)

- height:

  (optional) Column containing the dorso-ventral height of the body at
  each midline point (N points)

- .out:

  Names of the output columns. Needs to have two elements specifying the
  names for the x and y coordinates of center position. Or it can be a
  named list containing at least some of the elements `xctr` and `yctr`.
  If the return elements aren't in the named list, the defaults are
  'xcom' and 'ycom')

- .frame, .point:

  Columns identifying frames and points (defaults are `frame` and
  `point`)

- excludepoints:

  Exclude these points when estimating center. Some points (like the tip
  of the tail) have relatively little mass and are hard to track, so can
  introduce errors.

- method:

  'mutate' or 'summarize'. If summarize, returns one center position for
  each frame. If mutate, returns a same center position repeated for
  each point in a frame.

- overwrite:

  TRUE or FALSE to overwrite existing columns, if present.

## Value

A data frame containing the original variables along with xcom, ycom (or
names as specified in `.out`). The center of each midline in each frame.

## Details

Given a mass distribution, it produces an estimates of the true center
of mass. If given the body width and height, it assumes that the body
has an oval cross section with varying width and height, and it
estimates the volume distribution. This method will give a good estimate
of the center of mass if the body has close to uniform density. If given
just the width, it uses the width to estimate a weight average centroid
position.

## Examples

``` r
library(dplyr)
lampreydata |>
  group_by(frame) |>
  mutate(arclen = arclength(mxmm, mymm),
         width = interpolate_width(fishwidth$s, fishwidth$ammowidth, arclen)) |>
  ungroup() |>
  get_midline_center_df(arclen, mxmm, mymm, width = width)
#> ℹ Estimating center of mass based on width
#> # A tibble: 1,600 × 10
#>        t frame point  mxmm  mymm arclen width  xcom  ycom  nsum
#>    <dbl> <int> <int> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
#>  1  0.02     1     1    NA    NA     NA    NA    NA    NA     0
#>  2  0.02     1     2    NA    NA     NA    NA    NA    NA     0
#>  3  0.02     1     3    NA    NA     NA    NA    NA    NA     0
#>  4  0.02     1     4    NA    NA     NA    NA    NA    NA     0
#>  5  0.02     1     5    NA    NA     NA    NA    NA    NA     0
#>  6  0.02     1     6    NA    NA     NA    NA    NA    NA     0
#>  7  0.02     1     7    NA    NA     NA    NA    NA    NA     0
#>  8  0.02     1     8    NA    NA     NA    NA    NA    NA     0
#>  9  0.02     1     9    NA    NA     NA    NA    NA    NA     0
#> 10  0.02     1    10    NA    NA     NA    NA    NA    NA     0
#> # ℹ 1,590 more rows
```
