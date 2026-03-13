# Find gaps in a data series

Looks for NAs (gaps) in a data series and counts the number of NAs.
Useful for `smooth_points_df`, which can fill gaps up to a certain
length. 0 means no gap, 1 means a single NA, and so forth.

## Usage

``` r
find_gaps_df(.data, cols, .frame = frame, .out = c("gaplen"))
```

## Arguments

- .data:

  Data frame containing the midlines.

- cols:

  Columns containing the components to be smoothed. Often these will be
  the x and y coordinates of the midline.

- .frame:

  Column identifying frames (default is `frame`)

- .out:

  Name of the column to contain the length of the gap.

## Value

The data frame with a new column.

## Examples

``` r
# create a data frame with two NA gaps of different lengths
df <- data.frame(frame = 1:10,
                 x = c(1, 2, NA, NA, 5, 6, 7, NA, 9, 10))
find_gaps_df(df, x)
#>    frame  x gaplen
#> 1      1  1      0
#> 2      2  2      0
#> 3      3 NA      2
#> 4      4 NA      2
#> 5      5  5      0
#> 6      6  6      0
#> 7      7  7      0
#> 8      8 NA      1
#> 9      9  9      0
#> 10    10 10      0
```
