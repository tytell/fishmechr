#' Gets the main swimming axis from a midline
#'
#' Computes the main swimming axis of a midline as a unit vector,
#' using the singular value
#' decomposition ([svd()]). This only works well if the midlines are centered
#' around zero, so it optionally subtracts off the mean of x and y. For more
#' sophisticated centering algorithms, see [get_midline_center_df()].
#'
#' @param x,y Coordinates of the midline
#' @param center (TRUE or FALSE) Subtract the mean from the x and y coordinates
#'
#' @returns A data frame with the following columns:
#'   * `swimaxis_x`, `swimaxis_y` x and y components of the swimming axis vector
#'   * `swimaxis_xctr`, `swimaxis_yctr` Mean x and y values that were subtracted
#'     before running the SVD
#'
#' @seealso [get_midline_center_df()]
#'
#' @export
#'
#' @examples
#' # run the algorithm across multiple midlines at different times
#' lampreydata |>
#'   group_by(t) |>
#'   summarize(swimaxis = get_primary_swimming_axis(mxmm, mymm)) |>
#'   unnest(swimaxis)
get_primary_swimming_axis <- function(x, y, center=TRUE)
{
  if (all(is.finite(x)) && all(is.finite(y))) {
    if (center) {
      x_ctr <- mean(x)
      y_ctr <- mean(y)

      x <- x - x_ctr
      y <- y - y_ctr
    }

    XY <- matrix(c(x, y), ncol=2)
    S <- svd(XY)

    ab <- data.frame(swimaxis_x = S$v[1,1], swimaxis_y = S$v[2,1])
    if (center) {
      ab <- dplyr::bind_cols(ab,
         data.frame(swimaxis_xctr = x_ctr, swimaxis_yctr = y_ctr))
    }
  } else {
    ab <- data.frame(swimaxis_x = NA, swimaxis_y = NA)
    if (center) {
      ab <- dplyr::bind_cols(ab,
                             data.frame(swimaxis_xctr = NA,
                                        swimaxis_yctr = NA))
    }
  }

  ab
}

#' Gets the primary swimming axis for many midlines
#'
#' Processes midlines from many frames of a video
#'
#' Uses [get_primary_swimming_axis()] to compute the swimming axis for a midline.
#' Then optionally smooths the axis using a Butterworth filter, and then
#' projects the midlines on to the new time-varying axes.
#'
#' @param .data Data frame containing the midline data
#' @param t Column containing the time data. If a cutoff frequency is passed in,
#'   then this variable will be used to get the sampling frequency.
#' @param x,y Columns containing the x and y coordinates of each point along
#'   the midline.
#' @param .out Names of the output columns. Needs to have four elements specifying
#'   the names for the x and y coordinates of the swim axis and the parallel and
#'   perpendicular components of the excursion, in that order. Or it can be a named
#'   list containing at least some of the elements `swimaxis_x`, `swimaxis_y`, `exc_x`,
#'   `exc`, in any order. If the return elements aren't in the named list, the
#'   defaults are 'swimaxis_x', 'swimaxis_y', 'exc_x', and 'exc_y')
#' @param .frame,.point Columns identifying frames and points (defaults are `frame`
#'   and `point`)
#' @param cutoff (optional) If this parameter is included, smooth the swimming
#'   axis data with a low-pass filter with a cutoff at this frequency.
#'
#' @returns A data frame containing the original variables along with
#' * XX_ctr,YY_ctr: The center of each midline at each time, where XX and YY are
#'   the original names of the x and y coordinates.
#' * a,b: The new midlines centered and projected on to the swimming direction
#'   and the perpendicular axis. `b` is useful as the lateral excursion of the
#'   swimming undulation.
#'
#' @concept pipeline
#' @export
get_primary_swimming_axis_df <- function(.data, t, x,y,
                                         .out = NULL,
                                         .frame=frame, .point=point,
                                         cutoff=NULL, overwrite=TRUE,
                                         check_reasonableness=TRUE)
{
  .out <- check.out(.data, .out,
                    .out_default = c(swimaxis_x='swimaxis_x', swimaxis_y='swimaxis_y',
                                     exc_x='exc_x', exc='exc'),
                    overwrite=overwrite)

  .frame <- enquo(.frame)
  if (missing(.frame)) {
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.frame)),
                            msg = "Default column 'frame' not present. Use .frame to specify the name of the frame column")
  }
  .point <- enquo(.point)
  if (missing(.point)) {
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.point)),
                            msg = "Default column 'point' not present. Use .point to specify the name of the point column")
  }

  if (check_reasonableness) {
    centering <-
      .data |>
      group_by(!!.frame) |>
      summarize(
        xctr = mean({{x}}, na.rm = TRUE),
        xsd = sd({{x}}, na.rm = TRUE),
        yctr = mean({{y}}, na.rm = TRUE),
        ysd = sd({{y}}, na.rm = TRUE)
        ) |>
      summarize(notcenterx = sum(xctr > xsd, na.rm = TRUE) / n(),
                notcentery = sum(yctr > ysd, na.rm = TRUE) / n())

    if (centering$notcenterx > 0.1 ||
        centering$notcentery > 0.1) {
      warning("Many frames seem not to be centered around zero. Did you remember to subtract the center of mass?")
    }
  }
  # run across all of the frames
  swimaxis <-
    .data |>
    group_by(!!.frame, .add = TRUE) |>
    summarize(swimaxis = get_primary_swimming_axis({{x}},{{y}}),
              t = first({{t}})) |>
    tidyr::unnest(swimaxis)

  # get the sampling rate
  dt <- swimaxis$t[2] - swimaxis$t[1]

  # filter the x and y components of the swimming axis
  if (!is.null(cutoff)) {
    # run the filter
    filt <- build_filter(hi = cutoff, 1/dt)

    swimaxis <- swimaxis |>
      mutate(swimaxis_x0 = swimaxis_x,
             swimaxis_y0 = swimaxis_y,
             across(c(swimaxis_x, swimaxis_y), \(x) apply_filter(filt, x)))

    # be careful to re-normalize, because the smoothed data doesn't necessarily
    # represent a unit vector
    swimaxis <- swimaxis |>
      mutate(swimaxis_mag = sqrt(swimaxis_x^2 + swimaxis_y^2),
             across(c(swimaxis_x, swimaxis_y), \(x) x / swimaxis_mag)) |>
      select(-swimaxis_mag)
  } else {
    swimaxis <- swimaxis |>
      mutate(swimaxis_x0 = swimaxis_x,
             swimaxis_y0 = swimaxis_y)
  }
  swimaxis <- swimaxis |>
    rename("{.out[1]}" := swimaxis_x,
           "{.out[2]}" := swimaxis_y)

  # the swimming axis is defined once per time, but we need to repeat it for
  # each point along the body, so we use a left_join to repeat the values
  ab <-
    left_join(
    .data |>
      ungroup(!!.frame) |>
      select(!!.frame, {{.point}}, {{x}}, {{y}}),
    swimaxis,
    by = rlang::as_name(.frame))

  # then this centers each midline and projects them on to the swimming axis
  # and the perpendicular axis
  ab <-
    ab |>
    mutate("{.out[3]}" := {{x}} * .data[[.out[1]]] + {{y}} * .data[[.out[2]]],
           "{.out[4]}" := -{{y}} * .data[[.out[1]]] + {{x}} * .data[[.out[2]]]) |>
    select(any_of(.out), !!.frame, {{.point}})

  left_join(
    .data |> select(-any_of(.out)),
    ab,
    by = c(rlang::as_name(.frame), rlang::as_name(.point))
    )
}
