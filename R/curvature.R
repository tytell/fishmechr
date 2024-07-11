#' Calculate arc length along a 2D curve
#'
#' Computes the straight-line distance between points on the curve and
#' then adds them up to get the arc length.
#'
#' @param x,y Coordinates of the curve.
#'
#' @returns Arc length along the curve.
#' @export
arclength <- function(x, y)
{
  assertthat::assert_that(length(x) == length(y))

  if (any(!is.finite(x)) | any(!is.finite(y))) {
    s <- rep(NA, length(x))
  } else {
    s <- rep(0, length(x))
    ds <- sqrt(diff(x)^2 + diff(y)^2)
    s[2:length(s)] <- cumsum(ds)
  }
  s
}

interpolate_points_frame <- function(arclen, x,y, arclen_out,
                                     spar = 0.1)
{
  assertthat::assert_that(length(x) == length(y))

  xs <- numeric(length(x))
  ys <- numeric(length(y))

  good <- is.finite(x) & is.finite(y)

  if (any(good)) {
    spx <- smooth.spline(arclen[good], x[good], spar = spar)
    spy <- smooth.spline(arclen[good], y[good], spar = spar)

    xs[good] = predict(spx, x = arclen_out)$y
    ys[good] = predict(spy, x = arclen_out)$y
  }

  xs[!good] <- NA
  ys[!good] <- NA

  tibble(xs = xs, ys = ys)
}

interpolate_points_df <- function(.data, arclen, x,y,
                                  arclen_out = NULL,
                                  spar = 0.8,
                                  tailmethod = "extrapolate",
                                  .suffix = "_s",
                                  .out = NULL,
                              .frame = frame, .point = point)
{
  assertthat::assert_that(tailmethod %in% c("keep", "extrapolate", "NA"))

    if (missing(.frame)) {
    .frame <- enquo(.frame)
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.frame)),
                            msg = "Default column 'frame' not present. Use .frame to specify the name of the frame column")
  }
  if (missing(.point)) {
    .point <- enquo(.point)
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.point)),
                            msg = "Default column 'point' not present. Use .point to specify the name of the point column")
  }

  if (is.null(arclen_out)) {
    arclen_out <-
      .data |>
      group_by({{.point}}) |>
      summarize(s = median({{arclen}}, na.rm = TRUE))
    arclen_out <- arclen_out$s
  }
  assertthat::assert_that(all(is.finite(arclen_out)))

  if (is.null(.out)) {
    .out = c(paste0(rlang::as_name(enquo(arclen)), .suffix),
             paste0(rlang::as_name(enquo(x)), .suffix),
             paste0(rlang::as_name(enquo(y)), .suffix))
  } else {
    assertthat::assert_that(length(.out) == 3,
                            msg = ".out must contain three names, for the arc length, x, and y positions")
  }

  df <- .data |>
    group_by({{.frame}}, .add = TRUE) |>
    mutate(xys = list(interpolate_points_frame({{arclen}}, {{x}},{{y}}, arclen_out, spar)),
           "{.out[1]}" := arclen_out,
           "{.out[2]}" := xys[[1]]$xs,
           "{.out[3]}" := xys[[1]]$ys) |>
    select(-xys)

  if (tailmethod == "keep") {
    df <- df |>
      mutate(
        "{.out[1]}" := if_else(row_number({{arclen}}) == n(), {{arclen}}, .data[[.out[1]]]),
        "{.out[2]}" := if_else(row_number({{x}}) == n(), {{x}}, .data[[.out[2]]]),
        "{.out[3]}" := if_else(row_number({{y}}) == n(), {{y}}, .data[[.out[3]]]))
  } else if (tailmethod == "NA") {
    df <- df |>
      mutate(
        "{.out[2]}" := if_else(row_number({{x}}) == n() &&
                                 .data[[.out[1]]] > {{arclen}}, NA, .data[[.out[2]]]),
        "{.out[3]}" := if_else(row_number({{y}}) == n() &&
                                 .data[[.out[1]]] > {{arclen}}, NA, .data[[.out[3]]]))
  }

  df
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
#'
#' @returns Derivative of y relative to x.
#' @export
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
#' Estimates curvature either directly through derivatives of the x and y
#' coordinates relative to arc length, or as the derivative of segment angle
#' relative to arc length.
#'
#' Assumes that points are in order along the curve from head to tail.
#'
#' @param s Arc length along the body.
#' @param x,y Coordinates of each point along the body
#' @param method ("xy" or "angle") for the direct formula or for the angle
#'     derivative
#' @returns Curvature.
#' @export
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
