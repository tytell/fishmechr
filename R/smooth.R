#' Smooths locations of points over time
#'
#' Smooths columns specified in `cols` for each individual point over time, using
#' a smoothing spline. Optionally fills in gaps of less than `fillgaps` frames.
#'
#' @param .data Data frame containing the midlines.
#' @param cols Columns containing the components to be smoothed. Often these will
#'   be the x and y coordinates of the midline.
#' @param .out Names of the output columns. Should either be a list with the same
#'   number of elements as cols, or a glue specification as in `dplyr::across` for
#'   the `.names` parameter. The default (NULL) means that the output columns will
#'   have the same name as the original column with an 's' appended at the end.
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

  # only look for gaps if they want us to fill them
  if (fillgaps > 0) {
    gapdata <- .data |>
      find_gaps_df({{ cols }}, .frame = {{ .frame }})
  } else {
    gapdata <- .data |>
      mutate(gaplen = 0)
  }

  if (is.null(.out)) {
    nms <- '{.col}s'
  } else {
    assertthat::assert_that(
      (length(.out) == 1) ||
        length(.out) == length(cols),
      msg = "`.out` must either have the same length as `cols` or must have a name spec suitable for use in `dplyr::across` in the `.names` parameter"
    )
    nms <- .out
  }

  gapdata |>
    group_by({{ .point }}) |>
    mutate(across(
      {{ cols }},
      \(y) smooth_point({{ .frame }}, y, gaplen <= fillgaps, spar = spar),
      .names = nms
    ))
}


#' Find gaps in a data series
#'
#' Looks for NAs (gaps) in a data series and counts the number of NAs.
#' Useful for `smooth_points_df`, which can fill gaps up to a certain length.
#' 0 means no gap, 1 means a single NA, and so forth.
#'
#' @param .data Data frame containing the midlines.
#' @param cols Columns containing the components to be smoothed. Often these will
#'   be the x and y coordinates of the midline.
#' @param .frame Column identifying frames (default is `frame`)
#' @param .out Name of the column to contain the length of the gap.
#'
#' @returns The data frame with a new column.
#' @export
#'
#' @examples
find_gaps_df <- function(.data, cols, .frame = frame, .out = c('gaplen')) {
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

#' Applies a smoothing spline to a data series, potentially with gaps
#'
#' Builds a smoothing spline for y(x), where `y` may contain gaps (NAs).
#' Ignores the NAs. The uses the spline to interpolate values at the
#' same `x` coordinates, potentially filling in the gaps.
#'
#' @param x x coordinate
#' @param y y coordinate, potentially with NAs
#' @param goodout Logical vector of where in the x coordinate to interpolate
#' @param spar Smoothing parameter (see `smooth.spline`)
#'
#' @returns The smoothed values
#' @export
#'
#' @examples
smooth_point <- function(x, y, goodout = NULL, spar) {
  ys <- numeric(length(y))
  good <- !is.na(y)

  if (missing(goodout) || is.null(goodout)) {
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
