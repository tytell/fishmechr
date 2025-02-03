#' Calculate arc length along a 2D curve
#'
#' Computes the straight-line distance between points on the curve and
#' then adds them up to get the arc length.
#'
#' @param x,y Coordinates of the curve.
#'
#' @returns Arc length along the curve.
#' @export
#'
#' @concept pipeline
#'
#' @examples
#' # compute arc length in each frame from the lamprey data set
#' lampreydata |>
#'   group_by(frame) |>
#'   mutate(arclen = arclength(mxmm, mymm))
#'
arclength <- function(x, y,
                      na.skip = FALSE)
{
  assertthat::assert_that(length(x) == length(y))

  good <- is.finite(x) & is.finite(y)
  if (all(good)) {
    s <- rep(0, length(x))
    ds <- sqrt(diff(x)^2 + diff(y)^2)
    s[2:length(s)] <- cumsum(ds)
  }
  else if (any(!good) & !na.skip) {
    s <- rep(NA, length(x))
  }
  else {
    x1 <- x[good]
    y1 <- y[good]

    s1 <- rep(0, length(x1))

    ds1 <- sqrt(diff(x1)^2 + diff(y1)^2)
    s1[2:length(s1)] <- cumsum(ds1)

    s <- rep(NA, length(x))
    s[good] <- s1
  }

  s
}

#' Interpolates x and y points on a curve to different arc lengths
#'
#' For a 2D curve with (x,y) coordinates parameterized by the arc length,
#' interpolate new (x,y) coordinates at new arc lengths. Smooth the input
#' data with a smoothing spline.
#'
#' @param arclen Input arc length
#' @param x,y Coordinates of points on the curve
#' @param arclen_out New arc length
#' @param spar Smoothing parameter (ranges from 0 for no smoothing to 1 for
#'   high smoothing; see [smooth.spline()] for more details.)
#' @param fill_gaps Fill internal missing points of this size or smaller.
#'   (0, default, means no filling; 1 means to fill single missing points)
#'
#' @returns A tibble containing the new interpolated and smoothed x and y
#'   coordinates as columns `xs` and `ys`
#' @export
#'
#' @examples
interpolate_points_frame <- function(arclen, x,y, arclen_out,
                                     spar = 0.1, fill_gaps = 0)
{
  assertthat::assert_that(length(x) == length(y))

  xs <- numeric(length(x))
  ys <- numeric(length(y))

  good <- is.finite(arclen) & is.finite(x) & is.finite(y)
  fillpts <- good

  if (sum(good) >= 4) {
    if (any(!good) & (fill_gaps > 0)) {
      gapn <- rep_len(0, length.out = length(arclen))

      for (i in seq(2, length(good))) {
        if (!good[i]) {
          gapn[i] <- gapn[i-1] + 1
        }
      }

      # don't fill gaps that go all the way to the end
      i <- length(good)
      while (!good[i] & (i > 1)) {
        fillpts[i] <- FALSE
        i <- i-1
      }
      i <- i-1

      gapsz <- 0
      while (i > 1) {
        if ((gapn[i] > 0) & (gapn[i+1] == 0)) {
          gapsz <- gapn[i]
        } else if (gapn[i] == 0) {
          gapsz <- 0
        }
        gapn[i] <- gapsz
        i <- i-1
      }

      nfill <- sum((gapn > 0) & (gapn <= fill_gaps), na.rm = TRUE)
      if (nfill > 0) {
        fillpts[gapn <= fill_gaps] <- TRUE
      }
    }

    spx <- smooth.spline(arclen[good], x[good], spar = spar)
    spy <- smooth.spline(arclen[good], y[good], spar = spar)

    xs[fillpts] = predict(spx, x = arclen_out[fillpts])$y
    ys[fillpts] = predict(spy, x = arclen_out[fillpts])$y
  }

  xs[!fillpts] <- NA
  ys[!fillpts] <- NA

  tibble(xs = xs, ys = ys)
}

#' Interpolates and smooths a 2D curve at new arc length
#'
#' For a 2D curve with (x,y) coordinates parameterized by the arc length,
#' interpolate new (x,y) coordinates at new arc lengths. Smooth the input
#' data with a smoothing spline.
#'
#' Operates on each frame (as defined in the `.frame` parameter) individually.
#'
#' @param .data Data frame
#' @param arclen Name of the input arc length column in `.data`
#' @param x,y Name of the columns that contain the coordinates of points on the
#'   curve
#' @param arclen_out Vector containing the new arc length
#' @param spar Smoothing parameter (ranges from 0 for no smoothing to 1 for
#'   high smoothing; see [smooth.spline()] for more details.)
#' @param tailmethod ('keep', 'extrapolate', or 'NA') Methods to estimate the
#'   position of the tail tip if the last value of `arclength_out` is longer
#'   than the maximum arc length in the current frame.
#'   * 'keep' to keep the existing tail point, even if it is not at the
#'     requested arc length
#'   * 'extrapolate' to extrapolate a tail tip position, assuming that the
#'     curve continues straight
#'   * 'NA' to use replace the tail point with NA in this case.
#' @param fill_gaps Fill internal missing points of this size or smaller.
#'   (0, default, means no filling; 1 means to fill single missing points)
#' @param .suffix (default = '_s') Suffix to append to the names of the
#'   arclen, x, and y columns after smoothing and interpolation.
#' @param .out Names of the output columns. Defaults are
#'   (arclen = 'arclen_s', xs = 'xs', ys = 'ys'). Overrides the `.suffix`
#'   parameter if it is included.
#' @param overwrite TRUE or FALSE to overwrite existing columns
#' @param .frame Name of the frame variable in the data frame
#' @param .point Name of the point variable in the data frame
#'
#' @returns A data frame with updated columns for the smoothed and iterpolated
#'   arc length, x and y coordinates.
#' @export
#'
#' @concept pipeline
#' @examples
interpolate_points_df <- function(.data, arclen, x,y,
                                  arclen_out = NULL,
                                  spar = 0.8,
                                  tailmethod = "extrapolate",
                                  fill_gaps = 0,
                                  .suffix = "_s",
                                  .out = NULL,
                                  overwrite = TRUE,
                              .frame = frame, .point = point)
{
  assertthat::assert_that(tailmethod %in% c("keep", "extrapolate", "NA"))

  if (is.null(.out)) {
    .out = c(paste0(rlang::as_name(enquo(arclen)), .suffix),
             paste0(rlang::as_name(enquo(x)), .suffix),
             paste0(rlang::as_name(enquo(y)), .suffix))
  } else {
    .out <- check.out(.data, .out,
                      .out_default = c(arclen='arclen_s', xs='xs', ys='ys'),
                      overwrite=overwrite)
  }

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

  df <- .data |>
    group_by({{.frame}}, .add = TRUE) |>
    mutate(xys = list(interpolate_points_frame({{arclen}}, {{x}},{{y}},
                                               arclen_out, spar, fill_gaps = fill_gaps)),
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
    D <- (dplyr::lead(y) - dplyr::lag(y)) / (dplyr::lead(x) - dplyr::lag(x))

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
      D <- 2*dplyr::lag(y) / ((x - dplyr::lag(x))*(dplyr::lead(x) - dplyr::lag(x))) -
        2*y / ((dplyr::lead(x) - x)*(x - dplyr::lag(x))) +
        2*dplyr::lead(y) / ((dplyr::lead(x) - x)*(dplyr::lead(x) - dplyr::lag(x)))

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
#'
#' @concept pipeline
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
