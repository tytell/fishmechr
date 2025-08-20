#' Gets the center of a midline for many midlines in a data frame
#'
#' Estimates the center of a midline based on mass distribution, volume distribution,
#' or body width.
#'
#' Given a mass distribution, it produces an estimates of the true center of mass.
#' If given the body width and height, it assumes that the body has an oval cross
#' section with varying width and height, and it estimates the volume distribution.
#' This method will give a good estimate of the center of mass if the body has
#' close to uniform density. If given just the width, it uses the width to
#' estimate a weight average centroid position.
#'
#' @param .data Data frame containing the midlines.
#' @param arclen The column containing the arc length. See [arclength()]
#' @param x,y Columns containing the x and y coordinates of the midline. There
#'   should be N points.
#' @param mass (optional) Column containing the mass of each segment, with N-1
#'   segments.
#' @param width (optional) Column containing the horizontal plane width of the
#'   body at each midline point (N points)
#' @param height (optional) Column containing the dorso-ventral height of the
#'   body at each midline point (N points)
#' @param .out Names of the output columns. Needs to have two elements specifying
#'   the names for the x and y coordinates of center position. Or it can be a named
#'   list containing at least some of the elements `xctr` and `yctr`. If the
#'   return elements aren't in the named list, the defaults are 'xcom' and
#'   'ycom')
#' @param .frame,.point Columns identifying frames and points (defaults are `frame`
#'   and `point`)
#' @param excludepoints Exclude these points when estimating center. Some points
#'   (like the tip of the tail) have relatively little mass and are hard to track,
#'   so can introduce errors.
#' @param cutoff (optional) If this parameter is included, smooth the swimming
#'   axis data with a low-pass filter with a cutoff at this frequency.
#' @param method 'mutate' or 'summarize'. If summarize, returns one center position
#'   for each frame. If mutate, returns a same center position repeated for
#'   each point in a frame.
#' @param overwrite TRUE or FALSE to overwrite existing columns, if present.
#'
#' @returns A data frame containing the original variables along with
#'   xcom, ycom (or names as specified in `.out`). The center of each midline
#'   in each frame.
#' @export
#'
#' @concept pipeline
#' @examples
#' lampreydata |>
#'   group_by(frame) |>
#'   mutate(arclen = arclength(mxmm, mymm),
#'          width = interpolate_width(fishwidth$s, fishwidth$ammowidth, arclen)) |>
#'   get_midline_center_df(arclen, mxmm,mymm, width=width)
get_midline_center_df <- function(
  .data,
  arclen,
  x,
  y,
  mass,
  width,
  height,
  .out = NULL,
  .frame = frame,
  .point = point,
  excludepoints = c(),
  method = "mutate",
  overwrite = TRUE
) {
  assertthat::assert_that(method %in% c("mutate", "summarize", "summarise"))

  assertthat::assert_that(
    !is_grouped_df(.data),
    msg = "`get_midline_center_df` does not work on grouped data frames. Consider wrapping it in a call to `group_modify` to operate on groups separately"
  )

  .out <- check.out(
    .data,
    .out,
    .out_default = c(xctr = 'xcom', yctr = 'ycom'),
    overwrite = overwrite
  )

  if (missing(.frame)) {
    .frame <- enquo(.frame)
    assertthat::assert_that(
      assertthat::has_name(.data, rlang::as_name(.frame)),
      msg = "Default column 'frame' not present. Use .frame to specify the name of the frame column"
    )
  } else {
    .frame <- enquo(.frame)
  }

  if (length(excludepoints) > 0) {
    if (missing(.point)) {
      .point <- enquo(.point)
      assertthat::assert_that(
        assertthat::has_name(.data, rlang::as_name(.point)),
        msg = "Default column 'point' not present. Use .point to specify the name of the point column"
      )
    } else {
      .point <- enquo(.point)
    }

    if (any(!(excludepoints %in% rlang::eval_tidy(.point, .data)))) {
      cli::cli_alert_warning(
        'Some excluded points are not present in the data set'
      )
    }

    excludefcn <- function(df) filter(df, !(!!.point %in% excludepoints))
  } else {
    excludefcn <- function(df) df
  }

  if (method == "mutate") {
    fcn <- dplyr::mutate
  } else if (method %in% c("summarize", "summarise")) {
    fcn <- dplyr::summarise
  }

  if (!missing(mass)) {
    cli::cli_alert_info(
      "Estimating true center of mass based on mass distribution"
    )
    mass <- rlang::enquo(mass)

    com <- .data |>
      dplyr::group_by(!!.frame) |>
      excludefcn() |>
      summarize(
        M = sum(!!mass, na.rm = TRUE),
        data.table::`:=`(
          "{.out[1]}",
          sum(!!mass * ({{ x }} + dplyr::lead({{ x }})), na.rm = TRUE) / (2 * M)
        ),
        data.table::`:=`(
          "{.out[2]}",
          sum(!!mass * ({{ y }} + dplyr::lead({{ y }})), na.rm = TRUE) / (2 * M)
        ),
        nsum = sum(!is.na(!!mass) & !is.na({{ x }})),
        .groups = 'drop'
      ) |>
      select(-c(M))
  } else if (!missing(height) & !missing(width)) {
    cli::cli_alert_info("Estimating center of mass based on width and height")
    width <- rlang::enquo(width)
    height <- rlang::enquo(height)

    com <- .data |>
      excludefcn() |>
      group_by(!!.frame) |>
      mutate(V = get_volume({{ arclen }}, !!width, !!height)) |>
      summarize(
        sumV = sum(V, na.rm = TRUE),
        "{.out[1]}" := sum(V * ({{ x }} + dplyr::lead({{ x }})), na.rm = TRUE) /
          (2 * sumV),
        "{.out[2]}" := sum(V * ({{ y }} + dplyr::lead({{ y }})), na.rm = TRUE) /
          (2 * sumV),
        nsum = sum(!is.na(V) & !is.na({{ x }})),
        .groups = 'drop'
      ) |>
      select(-c(sumV))
  } else if (!missing(width) & missing(height)) {
    cli::cli_alert_info("Estimating center of mass based on width")
    width <- rlang::enquo(width)

    com <- .data |>
      group_by(!!.frame) |>
      excludefcn() |>
      summarize(
        sumw = sum(!!width, na.rm = TRUE),
        "{.out[1]}" := sum({{ x }} * !!width) / sumw,
        "{.out[2]}" := sum({{ y }} * !!width) / sumw,
        nsum = sum(!is.na({{ x }})),
        .groups = 'drop'
      ) |>
      select(-sumw)
  } else {
    cli::cli_alert_info("Estimating center of mass as the centroid of x and y")
    com <- .data |>
      group_by(!!.frame) |>
      excludefcn() |>
      summarize(
        "{.out[1]}" := mean({{ x }}, na.rm = TRUE),
        "{.out[2]}" := mean({{ y }}, na.rm = TRUE),
        nsum = sum(!is.na({{ x }})),
        .groups = 'drop'
      )
  }

  nmax <- max(com$nsum, na.rm = TRUE)
  if (any((com$nsum > 0) & (com$nsum < nmax), na.rm = TRUE)) {
    cli::cli_alert_warning(
      "Some frames have missing points. Dropping COM estimates for those frames"
    )
    com <- com |>
      mutate(across(any_of(.out), \(x) if_else(nsum == nmax, x, NA)))
  }
  if (method == "mutate") {
    .data <- .data |>
      ungroup() |>
      select(-any_of(.out)) |>
      left_join(com, by = c(rlang::quo_name(.frame)))
  } else {
    .data <- com
  }

  .data
}

#' Gets the volume of segments of a cylindrical body with elliptical cross section
#'
#' Used for estimating the center of mass of a fish. If we know the width and
#' height profile, and we assume that the cross section is elliptical, then
#' we can estimate the volume of each segment as the volume of a truncated
#' elliptical cone.
#'
#' The formula for such a cone is
#' V = pi s (w h + 1/2 dw h + 1/2 dh w + 1/3 dw dh)
#' where s is the length of the cone, w and h are the half width and height,
#' and dw and dh are the difference in width or height from
#' one end to the other (e.g., dw = w(s) - w(0) if w is a function of
#' s)
#'
#' @param arclen,width,height Arc length, width and height. Should have the same
#'   units. N points
#'
#' @returns Volume of each segment (N-1 values). Last value will be NA
#' @export
#'
#' @examples
get_volume <- function(arclen, width, height) {
  ds <- dplyr::lead(arclen) - arclen

  dw <- dplyr::lead(width) - width
  dh <- dplyr::lead(height) - height

  pi *
    ds *
    (width * height + 0.5 * dw * height + 0.5 * dh * width + 1 / 3 * dw * dh)
}

#' Interpolates and scales fish body width
#'
#' Interpolates the width for a new arc length and scales it based on body length.
#' Assumes that the input width and arc length have the same units (they could be
#' in fractions of body length, cm, or pixels, as long as they are the same).
#' Once the width is estimated at the new arc length, scales it based on the new
#' maximum length.
#'
#' Width here is defined as the distance from one side of the body to the other
#' (like a diameter), not from the center to a side (like a radius).
#'
#' @param arclen0 Arc length for the width measurement. The first value should
#'   be at the head and the last value should be at the tail tip.
#' @param width0 Width measurement. Should have the same units as `arclen0`
#' @param arclen New arc length
#' @param scale_to_body_length TRUE or FALSE to scale the interpolated width
#'   by multiplying by body length. This only works if `arclen` is in real units
#'   (like cm) so that the last value in `arclen` is equal to the total length
#'   of the fish.
#'
#' @concept pipeline
#' @returns Width at the new values of arc length, scaled for the new length
#' @export
interpolate_width <- function(
  arclen0,
  width0,
  arclen,
  scale_to_body_length = TRUE
) {
  if (sum(!is.na(arclen)) > 2) {
    len0 <- max(arclen0, na.rm = TRUE)
    len1 <- max(arclen, na.rm = TRUE)

    xy <- approx(arclen0 / len0 * len1, width0, xout = arclen)
    if (scale_to_body_length) {
      xy$y = xy$y / len0 * len1
    }
    xy$y
  } else {
    rep(NA, length(arclen))
  }
}
