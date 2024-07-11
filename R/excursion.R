#' Gets the main swimming axis from a midline
#'
#' Computes the main swimming axis of a midline as a unit vector,
#' using the singular value
#' decomposition ([svd()]). This only works well if the midlines are centered
#' around zero, so it optionally subtracts off the mean of x and y.
#'
#' @param x,y Coordinates of the midline
#' @param center (TRUE or FALSE) Subtract the mean from the x and y coordinates
#'
#' @returns A data frame with the following columns:
#'   * `swimaxis_x`, `swimaxis_y` x and y components of the swimming axis vector
#'   * `swimaxis_xctr`, `swimaxis_yctr` Mean x and y values that were subtracted
#'     before running the SVD
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
#' @param df Data frame containing the midline data
#' @param t Column containing the time data. If a cutoff frequency is passed in,
#'   then this variable will be used to get the sampling frequency.
#' @param point Column containing the point identification.
#' @param x,y Columns containing the x and y coordinates of each point along
#'   the midline.
#' @param cutoff (optional) If this parameter is included, smooth the swimming
#'   axis data with a low-pass filter with a cutoff at this frequency.
#'
#' @returns A data frame containing the original variables along with
#' * XX_ctr,YY_ctr: The center of each midline at each time, where XX and YY are
#'   the original names of the x and y coordinates.
#' * a,b: The new midlines centered and projected on to the swimming direction
#'   and the perpendicular axis. `b` is useful as the lateral excursion of the
#'   swimming undulation.
#' @export
get_primary_swimming_axis_df <- function(df, x,y, .frame=frame, .point=point,
                                         cutoff=NULL, overwrite=TRUE)
{

  newcols <- c("swimaxis_x", "swimaxis_y", "a", "b")
  if (any(newcols %in% colnames(df))) {
    dfname <- "Data frame"
    if (overwrite) {
      warning(dfname, " has columns that are assigned in 'get_primary_swimming_axis_df'. Overwriting")
    } else {
      stop(dfname, " has columns that are assigned in 'get_primary_swimming_axis_df'. Stopping")
    }
  }
  df <- df |>
    rename(t = {{t}},
           pt = {{point}})

  # run across all of the frames
  swimaxis <-
    df |>
    group_by(t) |>
    summarize(swimaxis = get_primary_swimming_axis({{x}},{{y}})) |>
    unnest(swimaxis)

  # get the sampling rate
  dt <- swimaxis$t[2] - swimaxis$t[1]

  # filter the x and y components of the swimming axis
  if (!is.null(cutoff)) {
    # run the filter
    filt <- build_filter(hi = cutoff, 1/dt)

    swimaxis <- swimaxis |>
      ungroup() |>
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
      ungroup() |>
      mutate(swimaxis_x0 = swimaxis_x,
             swimaxis_y0 = swimaxis_y)
  }

  # the swimming axis is defined once per time, but we need to repeat it for
  # each point along the body, so we use a left_join to repeat the values
  ab <-
    left_join(
    df |>
      ungroup() |>
      select(t, pt, {{x}}, {{y}}),
    swimaxis,
    by = "t")

  # then this centers each midline and projects them on to the swimming axis
  # and the perpendicular axis
  ab <-
    ab |>
    mutate(x_ctr = {{x}} - swimaxis_xctr,
           y_ctr = {{y}} - swimaxis_yctr,
           a = x_ctr * swimaxis_x + y_ctr * swimaxis_y,
           b = -y_ctr * swimaxis_x + x_ctr * swimaxis_y) |>
    select(t, pt, a,b, swimaxis_x, swimaxis_y,
           swimaxis_xctr, swimaxis_yctr)

  ctrnames <- c(paste0(rlang::as_name(enquo(x)), '_ctr'),
                paste0(rlang::as_name(enquo(y)), '_ctr'))

  # and join it back up with the original data frame, renaming variables back
  # to their original names
  left_join(
    df |> select(-any_of(c("a", "b", "swimaxis_x", "swimaxis_y",
                                   "swimaxis_xctr", "swimaxis_yctr"))),
    ab,
    by = c("t", "pt")
    ) |>
    select(-any_of(ctrnames)) |>
    rename("{{t}}" := t,
           "{{point}}" := pt,
           "{{x}}_ctr" := swimaxis_xctr,
           "{{y}}_ctr" := swimaxis_yctr)
}
