arclength <- function(x, y)
{
  assertthat::assert_that(length(x) == length(y))

  s <- rep(0, length(x))
  ds <- sqrt(diff(x)^2 + diff(y)^2)
  s[2:length(s)] <- cumsum(ds)

  s
}

#' Estimate first or second derivatives for dy/dx.
#'
#' Uses central differencing where possible.
#'
#' @param x x variable. Does not need to be evenly spaced.
#' @param y y variable.
#' @param ord Order of the derivative (1 or 2).
#' @param method Method for taking second derivatives. Either
#'   * 'direct' (default) Uses a direct formula, based on a central difference of
#'     forward and backward differences, from [https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/]
#'   * 'repeat' Repeat two first derivatives.
deriv <- function(x, y, ord = 1, method = 'direct', ends = 'forwardback')
{
  assertthat::assert_that(ord %in% c(1, 2),
                          method %in% c('direct', 'repeat'),
                          ends %in% c('forwardback', 'NA', NA, 'drop'))

  if (ord == 1) {
    # standard central difference formula for first derivative
    D <- (lead(y) - lag(y)) / (lead(x) - lag(x))

    if (ends == "forwardback") {
      # forward difference for the first point
      D[1] <- (y[2] - y[1]) / (x[2] - x[1])

      # backward difference for the last
      n <- length(x)
      D[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])
    } else if (ends == "drop") {
      D <- D[2:n-1]
    }
  } else if (ord == 2) {
    if (method == 'direct') {
      # direct formula for the second derivative, given uneven spacing in x
      # see https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/
      D <- 2*lag(y) / ((x - lag(x))*(lead(x) - lag(x))) -
        2*y / ((lead(x) - x)*(x - lag(x))) +
        2*lead(y) / ((lead(x) - x)*(lead(x) - lag(x)))

      if (ends == "drop") {
        D <- D[2:n-1]
      }
    } else if (method == 'repeat') {
      # second derivative by repeating first derivatives
      dydx <- deriv(x, y, ord = 1)
      D <- deriv(x, dydx, ord = 1)
    }
  }
  D
}

#' Estimates curvature for a single curve
#'
#' Assumes that points are in order
curvature <- function(s, x, y, method="angle")
{
  assertthat::assert_that(method %in% c("angle", "xy"))

  if (method == "xy") {
    dx <- deriv(s, x)
    dy <- deriv(s, y)
    ddx <- deriv(s, x, 2, method='direct')
    ddy <- deriv(s, y, 2, method='direct')

    (dx*ddy - dy*ddx) / ((dx^2 + dy^2)^1.5)
  } else if (method == "angle") {
    dx <- diff(x)
    dy <- diff(y)

    ang <- atan2(dy, dx)
    ds <- diff(s, lag = 2) / 2

    c(NA, diff(ang) / ds, NA)
  }
}
