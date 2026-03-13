# Pivots a kinematics dataset into long format

Converts a wide-format data frame where each point's variables are
separate columns (e.g. `head.x`, `head.y`, `tail.x`, `tail.y`) into a
long format with one row per point per frame. The point column is
returned as a factor with levels in the order given by `pointnames`.

## Usage

``` r
pivot_kinematics_longer(df, pointnames, point_to = "point", sep = "\\.")
```

## Arguments

- df:

  The data frame

- pointnames:

  The names of the points, in order from head to tail

- point_to:

  The name of the column to put the point names in

- sep:

  The separator between the point name and variable name, as a regular
  expression (default = `"\\."`)

## Value

The long data set

## Examples

``` r
df <- data.frame(
  frame = 1:3,
  head.x = c(0, 0, 0), head.y = c(0, 1, 2),
  tail.x = c(5, 5, 5), tail.y = c(0, 1, 2)
)
pivot_kinematics_longer(df, c("head", "tail"))
#> # A tibble: 6 × 4
#>   frame point     x     y
#>   <int> <fct> <dbl> <dbl>
#> 1     1 head      0     0
#> 2     1 tail      5     0
#> 3     2 head      0     1
#> 4     2 tail      5     1
#> 5     3 head      0     2
#> 6     3 tail      5     2
```
