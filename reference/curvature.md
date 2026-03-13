# Estimates curvature for a single curve

Estimates curvature either directly through derivatives of the x and y
coordinates relative to arc length, or as the derivative of segment
angle relative to arc length.

## Usage

``` r
curvature(s, x, y, method = "angle")
```

## Arguments

- s:

  Arc length along the body.

- x, y:

  Coordinates of each point along the body

- method:

  ("xy" or "angle") for the direct formula or for the angle derivative

## Value

Curvature.

## Details

Assumes that points are in order along the curve from head to tail.
