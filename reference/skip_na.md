# Skip NAs when running a function on a vector

`skip_na()` is a helper function related to
[`na.omit()`](https://rdrr.io/r/stats/na.fail.html). It runs a function
`f` on a vector that may contain NAs or NaNs, skipping all the NAs, and
returns the results as a vector of the same length as `x` with the NAs
in the same places.

## Usage

``` r
skip_na(x, f, min.len = 1, ...)
```

## Arguments

- x:

  A vector that may have NAs

- f:

  The function to run on the vector

- min.len:

  Minimum number of non-NA values required to run `f` (default 1)

- ...:

  Other parameters to supply to the function

## Value

A vector with the same length as `x` with NAs in the same places

## Examples

``` r
x <- c(1,2,3,NA,4,5,6,NA,7,8)
skip_na(x, cumsum)
#>  [1]  1  3  6 NA 10 15 21 NA 28 36
# should return a vector the same length as x with NAs in position 4 and 8
```
