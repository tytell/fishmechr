#' Estimates the cycle frequency based on time and phase
#'
#' Estimates frequency based on time and phase by one of two methods:
#'   * 'deriv' Takes the derivative of phase vs time to get frequency
#'   * 'slope' Fits a line to phase vs time and uses the slope as an estimate of
#'     frequency.
#'
#' @param t Time
#' @param ph Phase
#' @param unwrap (TRUE or FALSE) Unwrap the phase. Note that the
#'   function will not work unless the phase is unwrapped, so you should only
#'   set unwrap to FALSE if the phase has been unwrapped earlier.
#' @param method ('deriv' or 'slope') Estimate the frequency based on
#' @param mod Modulus for the phase. Default is `2*pi`.
#' @param check_reasonableness Runs some simple checks to make sure the data
#'   make sense. Checks to make sure phase is mostly increasing and warns if
#'   it isn't.
#'
#' @returns The frequency estimate, either as a vector the same size as `ph`
#'   (for the 'deriv' algorithm) or as a single value (for the 'slope' algorithm)
#' @export
#' @concept pipeline
#' @examples
#' t <- seq(0, 3, by=0.1)
#' # example phase that has a frequency of exactly 2.3 Hz
#' ph <- 2*pi*(t * 2.3)
#'
#' get_frequency(t, ph)
get_frequency <- function(t, ph, unwrap=FALSE, method='deriv', mod=2*pi,
                          check_reasonableness=TRUE)
{
  if (unwrap) {
    ph <- skip_na(ph, gsignal::unwrap)
  }

  if (method == "deriv") {
    freq <- deriv(t, ph) / mod

    if (check_reasonableness) {
      fracneg <- sum(freq < 0, na.rm=TRUE) / length(freq)

      if (fracneg > 0.05)
        warning("Phase seems to go backwards a lot, which may indicate an overly noisy signal")
    }
  } else if (method == "slope") {
    if (any(is.finite(ph)) && any(is.finite(t))) {
      m <- stats::lm(ph ~ t, na.action = na.omit)
      freq <- coefficients(m)[2] / mod
    }
    else {
      freq <- NA
    }
  }
  freq
}

#' Gets cycle numbers from a phase variable
#'
#' Given a phase variable that increases, each time the phase passes through
#' 2 k pi, a new cycle starts. This function unwraps the phase, so that it
#' increases steadily, then takes the floor of the phase divided by 2 pi (or
#' another modulus), so that it gets an integer for each cycle.
#'
#' Optionally, it will try to exclude partial cycles, setting the cycle number
#' to NA for cycles that do not start at a phase close to 0 and progress to a
#' phase close to 2 pi. "Close" here is defined based on the average change in
#' phase.
#'
#' @param ph Phase variable
#' @param unwrap (TRUE or FALSE) Unwrap the phase variable. Note that the
#'   function will not work unless the phase is unwrapped, so you should only
#'   set unwrap to FALSE if the phase has been unwrapped earlier.
#' @param mod Modulus for the phase. Default is `2*pi`.
#' @param exclude_partial_cycles (TRUE or FALSE) Exclude cycles in which the
#'   phase does not advance from close to 0 to close to 2pi.
#'
#' @returns Integer cycle numbers with the same size as `ph`
#' @export
#'
#' @examples
#' # example phase that advances by slightly more than 3 cycles, modulo 2pi
#' ph <- seq(0, 20, by = pi/10) %% (2*pi)
#' get_cycle_numbers(ph, unwrap=TRUE)
#'
#' @seealso [get_body_cycle_numbers_df()]
get_cycle_numbers <- function(ph, unwrap=FALSE, mod=2*pi,
                              exclude_partial_cycles=TRUE)
{
  if (unwrap) {
    ph <- skip_na(ph, gsignal::unwrap)
  }

  ph <- ph / mod
  dph = mean(diff(ph), na.rm = TRUE)

  if (exclude_partial_cycles) {
    d <- tibble::tibble(ph, cyc = floor(ph)) |>
      dplyr::group_by(cyc) |>
      dplyr::mutate(lo = min(ph) - cyc,
                    hi = max(ph) - cyc,
                    good = (lo <= 1.5*dph) & (hi >= (1 - 1.5*dph)),
                    cyc = dplyr::if_else(good, cyc, NA))

    d$cyc
  } else {
    floor(ph)
  }
}

#' Gets oscillation cycle numbers for a midline data set
#'
#' Extracts the phase for a specific point along the body (defined by `pointval`)
#' and uses [get_cycle_numbers()] to get the cycle number. Joins the cycle
#' number data with the original data set so that cycle number is defined for
#' each point along the body
#'
#' @param .data Data frame containing the midlines
#' @param ph Phase variable
#' @param pointval Specific point on the body to use to define the phase and
#'   and cycle number
#' @param .frame,.point Columns identifying frames and points (defaults are `frame`
#'   and `point`)
#' @param overwrite TRUE or FALSE to overwrite existing columns, if present.
#' @param .out Name of the output column. Needs to have one element specifying
#'   the name for the cycle variable, or as a named list with an element named
#'   `cycle`). Default is 'cycle'
#' @param ... Extra parameters to supply to [get_cycle_numbers()]
#' @concept pipeline
#' @returns A data frame containing a new column with the cycle numbers named
#'   'cycle' or as specified in .out.
#' @export
get_body_cycle_numbers_df <- function(.data, ph, pointval,
                                      .frame=frame, .point=point,
                                      .out=NULL,
                                      overwrite=TRUE, ...)
{
  .out <- check.out(.data, .out, .out_default = c(cycle = 'cycle'),
                    overwrite = overwrite)

  .frame <- rlang::enquo(.frame)
  if (missing(.frame)) {
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.frame)),
                            msg = "Default column 'frame' not present. Use .frame to specify the name of the frame column")
  }
  .point <- rlang::enquo(.point)
  if (missing(.point)) {
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.point)),
                            msg = "Default column 'point' not present. Use .point to specify the name of the point column")
  }

  by <- rlang::as_name(.frame)

  cycdf <-
    .data |>
    dplyr::ungroup() |>
    dplyr::filter(!!.point == pointval) |>
    dplyr::mutate("{.out[1]}" := get_cycle_numbers({{ph}}, ...))

  dplyr::left_join(.data |>
                     dplyr::select(-dplyr::any_of(.out)),
                   cycdf |>
                     dplyr::select(dplyr::any_of(by), cycle),
                   by = by)
}
