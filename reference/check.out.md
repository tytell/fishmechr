# Helper function to check the `.out` and `.out_default` variables used in this package

- If `.out` is NULL, returns `.out_default`.

- If `.out` is a vector containing strings, checks to make sure that
  there are the same number of strings in `.out` as items in
  `.out_default`

- If `.out` is a list with named elements, checks to make sure the names
  are present in `.out_default`. Ignores names not present in
  `.out_default`, but gives a warning. If any elements in `.out_default`
  are not in `.out`, uses the values from `.out_default`. Sorts the
  items in the same order as in `.out_default`.

## Usage

``` r
check.out(.data, .out, .out_default, overwrite)
```

## Arguments

- .data:

  Data frame

- .out:

  Name for output columns

- .out_default:

  Named list containing default names for output columns

- overwrite:

  TRUE or FALSE to overwrite columns.

## Value

The updated `.out` list.

## Details

Checks if there are columns with the same names in the data frame
`.data`. Gives a warning if there are such columns and `overwrite` is
TRUE, and an error if `overwrite` is FALSE.
