#' Smooths locations of points over time
#'
#' Smooths columns specified in `cols` for each individual point over time, using
#' a smoothing spline. Optionally fills in gaps of less than `fillgaps` frames.
#'
#' @param .data Data frame containing the midlines.
#' @param cols Columns containing the components to be smoothed. Often these will
#'   be the x and y coordinates of the midline.
#' @param .out Names of the output columns. Currently ignores this parameter
#' @param .frame,.point Columns identifying frames and points (defaults are `frame`
#'   and `point`)
#' @param fillgaps Longest gap to interpolate over. default is 0, which means not
#'   to fill gaps
#'
#' @returns A data frame containing the smoothed points
#' @export
#'
#' @examples
smooth_points_df <- function(
  .data,
  cols,
  spar,
  .out = NULL,
  .frame = frame,
  .point = point,
  fillgaps = 0
) {
  assertthat::assert_that(
    !is_grouped_df(.data),
    msg = "`get_midline_center_df` does not work on grouped data frames. Consider wrapping it in a call to `group_modify` to operate on groups separately"
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
  if (missing(.point)) {
    .point <- enquo(.point)
    assertthat::assert_that(
      assertthat::has_name(.data, rlang::as_name(.point)),
      msg = "Default column 'point' not present. Use .point to specify the name of the point column"
    )
  } else {
    .point <- enquo(.point)
  }

  if (fillgaps > 0) {
    gapdata <- .data |>
      find_gaps_df({{ cols }}, .frame = {{ .frame }})
  } else {
    gapdata <- .data |>
      mutate(gaplen = 0)
  }

  gapdata |>
    group_by({{ .point }}) |>
    mutate(across(
      {{ cols }},
      \(y) smooth_point({{ .frame }}, y, gaplen <= fillgaps, spar = spar),
      .names = '{.col}s'
    ))
}

find_gaps_df <- function(.data, cols, .frame = frame) {
  gaps <- .data |>
    mutate(across({{ cols }}, is.na)) |>
    mutate(good = (rowSums(pick({{ cols }})) == 0))

  gaps <- gaps |>
    mutate(
      gapstart = if_else(!good & dplyr::lag(good), {{ .frame }}, NA),
      gapend = if_else(!good & dplyr::lead(good), {{ .frame }}, NA)
    ) |>
    fill(gapend, .direction = 'up') |>
    fill(gapstart, .direction = 'down') |>
    mutate(
      gaplen = if_else(!good, gapend - gapstart + 1, 0),
      gaplen = if_else(is.na(gaplen), 0, gaplen)
    )

  .data$gaplen = gaps$gaplen

  return(.data)
}

smooth_point <- function(x, y, goodout, spar) {
  ys <- numeric(length(y))
  good <- !is.na(y)

  if (missing(goodout)) {
    goodout <- good
  }

  ys <- rep(NA, length(y))
  if (any(good)) {
    sp <- smooth.spline(x[good], y[good], spar = spar)

    if (any(goodout)) {
      ys[goodout] = predict(sp, x = x[goodout])$y
    }
  }

  ys[!goodout] <- NA

  return(ys)
}
