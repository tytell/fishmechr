# Estimate first or second derivatives for dy/dx.

Uses central differencing where possible.

## Usage

``` r
deriv(x, y, ord = 1, method = "direct", ends = "forwardback")
```

## Arguments

- x:

  x variable. Does not need to be evenly spaced.

- y:

  y variable.

- ord:

  Order of the derivative (1 or 2).

- method:

  Method for taking second derivatives. Either

  - 'direct' (default) Uses a direct formula, based on a central
    difference of forward and backward differences, from
    <https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/>

  - 'repeat' Repeat two first derivatives.

- ends:

  ('forwardback', 'NA', or 'drop') How to handle the endpoints where
  central differencing is not possible. `'forwardback'` (default) uses
  forward differencing at the first point and backward differencing at
  the last. `'NA'` sets endpoints to NA. `'drop'` removes the endpoints.

## Value

Derivative of y relative to x.
