# Apply a filter constructed with [build_filter](https://tytell.github.io/fishmechr/reference/build_filter.md)

Wrapper function for
[`gsignal::filtfilt`](https://rdrr.io/pkg/gsignal/man/filtfilt.html) to
apply a filter to a dataset. Potentially skips NAs in the data set (see
[`skip_na()`](https://tytell.github.io/fishmechr/reference/skip_na.md)
for details).

## Usage

``` r
apply_filter(filt, x, na.skip = TRUE)
```

## Arguments

- filt:

  Filter as returned from
  [`build_filter()`](https://tytell.github.io/fishmechr/reference/build_filter.md)
  or [`gsignal::butter()`](https://rdrr.io/pkg/gsignal/man/butter.html)

- x:

  Vector of data to filter

- na.skip:

  (TRUE or FALSE) to skip NAs in the data set.

## Value

Filtered data set

## Examples

``` r
filt <- build_filter(hi = 5, sampfreq = 100)
# 2 Hz signal sampled at 100 Hz with high-frequency noise and two NA gaps
x <- sin(2 * pi * (1:100) / 100 * 2) + rnorm(100, sd = 0.3)
x[c(30, 31, 70)] <- NA
apply_filter(filt, x)
#>   [1] -0.29517703 -0.10756466  0.07573123  0.25044203  0.41265865  0.55899098
#>   [7]  0.68669141  0.79373439  0.87884816  0.94149856  0.98182903  1.00056481
#>  [13]  0.99889164  0.97832137  0.94055738  0.88737184  0.82050518  0.74159538
#>  [19]  0.65214116  0.55349953  0.44691468  0.33357170  0.21466646  0.09148133
#>  [25] -0.03454406 -0.16175689 -0.28825857 -0.41190136 -0.53031691          NA
#>  [31]          NA -0.64097645 -0.74127890 -0.82866039 -0.90071681 -0.95532930
#>  [37] -0.99078244 -1.00586557 -0.99994878 -0.97302769 -0.92573340 -0.85930711
#>  [43] -0.77554154 -0.67669407 -0.56537826 -0.44444199 -0.31684111 -0.18551718
#>  [49] -0.05328724  0.07724791  0.20377241  0.32429850  0.43718375  0.54112428
#>  [55]  0.63512726  0.71846771  0.79063522  0.85127669  0.90014105  0.93703123
#>  [61]  0.96176760  0.97416615  0.97403266  0.96117336  0.93542039  0.89666986
#>  [67]  0.84492863  0.78036561  0.70336286          NA  0.61456173  0.51489948
#>  [73]  0.40563258  0.28834363  0.16492992  0.03757312 -0.09130969 -0.21913148
#>  [79] -0.34321444 -0.46087896 -0.56953763 -0.66678925 -0.75050777 -0.81892118
#>  [85] -0.87067605 -0.90488427 -0.92114969 -0.91957336 -0.90073766 -0.86567063
#>  [91] -0.81579287 -0.75285057 -0.67883843 -0.59591702 -0.50632882 -0.41231703
#>  [97] -0.31605091 -0.21956050 -0.12468278 -0.03302087
```
