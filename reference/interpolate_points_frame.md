# Interpolates x and y points on a curve to different arc lengths

For a 2D curve with (x,y) coordinates parameterized by the arc length,
interpolate new (x,y) coordinates at new arc lengths. Smooth the input
data with a smoothing spline.

## Usage

``` r
interpolate_points_frame(arclen, x, y, arclen_out, spar = 0.1, fill_gaps = 0)
```

## Arguments

- arclen:

  Input arc length

- x, y:

  Coordinates of points on the curve

- arclen_out:

  New arc length

- spar:

  Smoothing parameter (ranges from 0 for no smoothing to 1 for high
  smoothing; see
  [`smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html) for
  more details.)

- fill_gaps:

  Fill internal missing points of this size or smaller. (0, default,
  means no filling; 1 means to fill single missing points)

## Value

A tibble containing the new interpolated and smoothed x and y
coordinates as columns `xs` and `ys`

## Examples

``` r
# get one frame of lamprey midline data
df1 <- lampreydata[lampreydata$frame == 3, ]
# compute arc length along the midline
s <- with(df1, arclength(mxmm, mymm))
# define 20 evenly-spaced output arc lengths from head to tail
s_new <- seq(0, max(s, na.rm = TRUE), length.out = 20)
# interpolate and smooth the midline at the new arc lengths
with(df1, interpolate_points_frame(s, mxmm, mymm, s_new))
#> # A tibble: 20 × 2
#>       xs    ys
#>    <dbl> <dbl>
#>  1  415.  136.
#>  2  423.  136.
#>  3  431.  137.
#>  4  439.  137.
#>  5  446.  139.
#>  6  454.  142.
#>  7  461.  146.
#>  8  468.  149.
#>  9  476.  151.
#> 10  484.  150.
#> 11  491.  148.
#> 12  499.  146.
#> 13  507.  145.
#> 14  515.  145.
#> 15  522.  147.
#> 16  529.  151.
#> 17  536.  155.
#> 18  543.  158.
#> 19  551.  160.
#> 20  559.  159.
```
