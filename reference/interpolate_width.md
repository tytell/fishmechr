# Interpolates and scales fish body width

Interpolates the width for a new arc length and scales it based on body
length. Assumes that the input width and arc length have the same units
(they could be in fractions of body length, cm, or pixels, as long as
they are the same). Once the width is estimated at the new arc length,
scales it based on the new maximum length.

## Usage

``` r
interpolate_width(arclen0, width0, arclen, scale_to_body_length = TRUE)
```

## Arguments

- arclen0:

  Arc length for the width measurement. The first value should be at the
  head and the last value should be at the tail tip.

- width0:

  Width measurement. Should have the same units as `arclen0`

- arclen:

  New arc length

- scale_to_body_length:

  TRUE or FALSE to scale the interpolated width by multiplying by body
  length. This only works if `arclen` is in real units (like cm) so that
  the last value in `arclen` is equal to the total length of the fish.

## Value

Width at the new values of arc length, scaled for the new length

## Details

Width here is defined as the distance from one side of the body to the
other (like a diameter), not from the center to a side (like a radius).
