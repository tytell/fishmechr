# Compute phase of an oscillation using the Hilbert transform

Given a value that oscillates (`x`), first computes the analytic signal
using the Hilbert transform, then the phase of that signal.

## Usage

``` r
hilbert_phase(x, na.skip = TRUE, unwrap = TRUE, check_reasonableness = TRUE)
```

## Arguments

- x:

  Oscillating signal

- na.skip:

  TRUE or FALSE to skip NAs. See
  [`skip_na()`](https://tytell.github.io/fishmechr/reference/skip_na.md)
  for differences with na.omit (default = TRUE)

- unwrap:

  TRUE or FALSE to unwrap the phase signal so that it increases smoothly
  over time (default = TRUE)

- check_reasonableness:

  TRUE or FALSE to run some checks on the input data to make sure that
  the output will be reasonable (default = TRUE)

## Value

Phase (mod 2pi)

## Examples

``` r
t <- seq(0, 4, by=0.1)
# signal that increases in frequency over time
x <- cos(2*pi*t*(2.2 + 0.2*t))
ph <- hilbert_phase(x)
#> Warning: For a good phase estimate, the signal should oscillate around 0
```
