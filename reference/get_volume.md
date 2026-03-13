# Gets the volume of segments of a cylindrical body with elliptical cross section

Used for estimating the center of mass of a fish. If we know the width
and height profile, and we assume that the cross section is elliptical,
then we can estimate the volume of each segment as the volume of a
truncated elliptical cone.

## Usage

``` r
get_volume(arclen, width, height)
```

## Arguments

- arclen, width, height:

  Arc length, width and height. Should have the same units. N points

## Value

Volume of each segment (N-1 values). Last value will be NA

## Details

The formula for such a cone is V = pi s (w h + 1/2 dw h + 1/2 dh w + 1/3
dw dh) where s is the length of the cone, w and h are the half width and
height, and dw and dh are the difference in width or height from one end
to the other (e.g., dw = w(s) - w(0) if w is a function of s)

## Examples

``` r
# volume of lamprey body segments using measured width and estimated height
h <- seq(0.05, 0.03, length.out = nrow(fishwidth))  # height tapers head to tail
get_volume(fishwidth$s, fishwidth$ammowidth, h)
#>  [1] 1.345691e-04 1.387190e-04 1.408014e-04 1.398325e-04 1.376938e-04
#>  [6] 1.355068e-04 1.310126e-04 1.239205e-04 1.165924e-04 1.084350e-04
#> [11] 9.917619e-05 8.992416e-05 8.073924e-05 7.074856e-05 6.033214e-05
#> [16] 5.017277e-05 3.835812e-05 2.531397e-05 9.638956e-06           NA
```
