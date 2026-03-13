# Interpolate the location of a peak based on three points

Uses a parabolic approximation to determine the location of a peak from
3 points.

## Usage

``` r
interpolate_peak_location(y, x = c(-1, 0, 1))
```

## Arguments

- y:

  3 points on the curve, where the peak/trough should be at `y[2]`, so
  that `y[1]` and `y[3]` are lower/higher, respectively, than `y[2]`

- x:

  x coordinates. Defaults to `x = c(-1, 0, 1)`, so that the output is an
  offset of the peak location from what was originally detected. It
  could also be the true `x` coordinates, which is important if they're
  unevenly spaced.

## Value

The location of the peak, either as an offset (with default `x`) or as a
true `x` coordinate

## Examples

``` r
y <- c(1, 2, 0.5)
interpolate_peak_location(y)
#> [1] -0.1
```
