# Check if a data frame has columns that we might overwrite

Helper function for other functions in the package that create new
columns in a data frame, to check if the columns are already present in
the data frame.

## Usage

``` r
check_if_overwrite_columns(df, newcols, overwrite)
```

## Arguments

- df:

  Data frame

- newcols:

  Names of the columns, as strings

- overwrite:

  TRUE or FALSE to overwrite the columns

## Details

Produces a warning if the columns are present but `overwrite` is true,
and an error if `overwrite` is false.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(a=c(1,2,3), b=c(1,2,3))

# this should give a warning
check_if_overwrite_columns(df, c('a', 'd'), overwrite=TRUE)
} # }
```
