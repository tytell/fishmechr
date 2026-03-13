# Constructs a smoothing filter

Wrapper function for
[`gsignal::butter`](https://rdrr.io/pkg/gsignal/man/butter.html) to
construct a bandpass (or low or high pass) filter at a particular
frequency. Uses the Sos form for the filter, which works better
numerically particularly for very low frequencies.

## Usage

``` r
build_filter(lo = NULL, hi = NULL, sampfreq, order = 13)
```

## Arguments

- lo:

  Lower frequency cutoff in Hz. If only `lo` is specified, then
  `build_filter` constructs a high pass filter

- hi:

  Upper frequency cutoff in Hz. If only `hi` is specified, then
  `build_filter` constructs a high pass filter

- sampfreq:

  Sampling frequency for the data in Hz.

- order:

  (optional) Order for the filter

## Value

Filter parameters in Sos form

## Examples

``` r
# Low pass filter with a cutoff at 0.5Hz for data sampled at 100Hz
lopass <- build_filter(lo=0.5, sampfreq=100)

# Band pass filter that passes frequencies between 0.5 and 10Hz
bandpass <- build_filter(lo=0.5, hi=10, sampfreq=100)
```
