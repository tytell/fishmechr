# Compute phase of an oscillation by locating peaks and zero crossings.

For an oscillatory signal, we can define a positive peak as having phase
0 and a negative peak as having phase pi, with the zero crossings or
intermediate points having phase pi/2 and 3pi/2. Then, if we find those
peaks and zero crossings, we can interpolate to find the phase at any
point.

## Usage

``` r
peak_phase(
  x,
  unwrap = TRUE,
  check_reasonableness = TRUE,
  check_ordering = TRUE,
  interpolate_peaks = TRUE,
  interpolate_zeros = TRUE,
  zero_mode = "midpoint",
  ...
)
```

## Arguments

- x:

  Oscillatory signal

- unwrap:

  TRUE or FALSE to unwrap the phase signal so that it increases smoothly
  over time (default = TRUE)

- check_reasonableness:

  TRUE or FALSE to run some checks on the input data to make sure that
  the output will be reasonable (default = TRUE)

- check_ordering:

  TRUE or FALSE to check the order of peaks and zeros. A good signal
  should have a positive peak followed by a downward zero, then a
  negative peak followed by an upward zero. (default = TRUE)

- interpolate_peaks:

  TRUE or FALSE to interpolate the locations of peaks using a 3-point
  parabola (see
  [`interpolate_peak_location()`](https://tytell.github.io/fishmechr/reference/interpolate_peak_location.md)).
  (default = TRUE)

- interpolate_zeros:

  TRUE or FALSE to interpolate the locations of zero crossings linearly
  (default = TRUE)

- zero_mode:

  'midpoint', 'zero', or 'none' or NA. Define zero crossings as the
  point halfway between a positive and a negative peak ('midpoint'), or
  as a genuine zero crossing ('zero'). If 'none' or NA, do not detect
  zeros.

- ...:

  Other parameters to supply to
  [`pracma::findpeaks()`](https://rdrr.io/pkg/pracma/man/findpeaks.html).
  'minpeakdistance', the minimum number of index values between peaks,
  is often the most useful.

## Value

The phase of the oscillatory signal.

## Details

Uses
[`pracma::findpeaks()`](https://rdrr.io/pkg/pracma/man/findpeaks.html)
to locate peaks.

## Examples

``` r
t <- seq(0, 4, by=0.1)
# signal that increases in frequency over time
x <- cos(2*pi*t*(2.2 + 0.2*t))
ph <- peak_phase(x)
```
