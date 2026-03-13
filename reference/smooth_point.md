# Applies a smoothing spline to a data series, potentially with gaps

Builds a smoothing spline for y(x), where `y` may contain gaps (NAs).
Ignores the NAs. The uses the spline to interpolate values at the same
`x` coordinates, potentially filling in the gaps.

## Usage

``` r
smooth_point(x, y, spar, goodout = NULL)
```

## Arguments

- x:

  x coordinate

- y:

  y coordinate, potentially with NAs

- spar:

  Smoothing parameter (see
  [`smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html))

- goodout:

  Logical vector of where in the x coordinate to interpolate

## Value

The smoothed values

## Examples

``` r
# smooth a noisy sine wave with two missing points
x <- 1:20
y <- sin(2 * pi * x / 10) + rnorm(20, sd = 0.1)
y[c(9, 10)] <- NA
smooth_point(x, y, spar = 0.5)
#>  [1]  0.84327603  0.73127994  0.53785165  0.22529583 -0.16112910 -0.50470913
#>  [7] -0.68509595 -0.62930570          NA          NA  0.39144429  0.61552250
#> [13]  0.62610588  0.41431753  0.05795349 -0.32174107 -0.59994154 -0.67803047
#> [19] -0.54409923 -0.29123282
```
