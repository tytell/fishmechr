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

#' Computes the body wavelength based on the phase at each point and the arc length
#'
#' Computes the spatial derivative in phase relative to arc length using several
#' different possible methods:
#'   * 'deriv' Computes the derivative directly using [deriv()]
#'   * 'slope' Fits a line to the phase relative to arc length and uses the
#'     slope of that line as an estimate of the wavelength
#'   * 'cycle' Looks for pairs of points along the body where the phase differs by
#'     a full cycle. The arc length between those points is one wavelength.
#'   * 'halfcycle' Looks for pairs of points along the body where the phase
#'     differs by one half cycle. The arc length between those points is
#'     half of a wavelength.
#'
#' @param arclen Arc length
#' @param ph Phase
#' @param unwrap (TRUE or FALSE) Unwrap the phase along the body
#' @param method ('deriv', 'slope', 'cycle', 'halfcycle') as explained above
#' @param ignore_arclen_vals NULL or a function that returns TRUE or FALSE for
#'   certain values of arclength where the phase estimate is not reliable. This
#'   is often used to exclude points near the head (e.g., `ignore_arclen_vals =
#'   \(s) s < 0.3`)
#' @param mod Modulus for the phase variable
#' @param traveling_wave_dir (1 or -1) Defines the direction of the traveling
#'   wave. For a normal head-to-tail traveling wave, use -1 (default)
#'
#' @returns The wavelength as a vector the same size as the phase variable
#' @export
#' @concept pipeline
#' @examples
#' s <- seq(0, 1, by=0.1)
#' # artificial data with a wavelength of exactly 0.6
#' ph <- 2*pi*((1 - s) / 0.6)
#'
#' get_wavelength(s, ph, ignore_arclen_vals = \(s) s < 0.3)
get_wavelength <- function(arclen, ph, unwrap=TRUE, method='deriv',
                           ignore_arclen_vals=NULL,
                           sort_arclen=FALSE,
                           check_reasonableness=TRUE,
                           mod=2*pi, traveling_wave_dir=-1)
{
  assertthat::assert_that(method %in% c('deriv', 'slope', 'cycle', 'halfcycle'))

  if (unwrap) {
    ph <- skip_na(ph, gsignal::unwrap)
  }

  if (!is.null(ignore_arclen_vals)) {
    bad <- ignore_arclen_vals(arclen)
  } else {
    bad <- rep(FALSE, length(arclen))
  }
  good <- !bad

  if (sort_arclen) {
    ph <- sort_by(ph, arclen)
    arclen <- sort(arclen, na.last = TRUE)
  }

  if (method == "deriv") {
    k <- deriv(arclen, ph) / mod
    wavelen <- traveling_wave_dir / k      # defined at all points
    wavelen[bad] <- NA
  }
  else if (method == "slope") {
    m <- lm(ph[good] ~ arclen[good])
    k <- coefficients(m)[2] / mod
    wavelen <- traveling_wave_dir / k     # single value
  }
  else if (method %in% c('cycle', 'halfcycle')) {
    if (method == "halfcycle") {
      m <- mod/2
      c <- 0.5
    } else {
      m <- mod
      c <- 1
    }

    ph[bad] <- NA

    ph <- ph * traveling_wave_dir

    if (check_reasonableness) {
      phdir <- ph[2:length(ph)] - ph[1:(length(ph)-1)]
      if (sum(phdir > 0, na.rm = TRUE) < 0.8*sum(!is.na(phdir))) {
        if (traveling_wave_dir < 0) {
          warning("Phase does not seem to be decreasing along the body. Wavelength may be strange")
        } else {
          warning("Phase does not seem to be decreasing along the body. Wavelength may be strange")
        }
      }
    }
    k <- which.max(ph)
    ph <- ph - ph[k]
    ph <- ph %% m
    ph[k] <- m

    i <- which(diff(ph) < 0)

    wavelen <- rep(NA_real_, length(arclen))
    if (length(i) >= 1) {
      # interpolate to get the exact location of the wrap

      swrap1 <- arclen[i]
      swrap2 <- arclen[i+1]
      phwrap1 <- ph[i]
      phwrap2 <- ph[i+1] + m

      swrap <- swrap1 + (swrap2 - swrap1) / (phwrap2 - phwrap1) * (m - phwrap1)

      swrap <- c(swrap, arclen[length(arclen)])
      i <- c(i, length(arclen))

      wavelen1 <- diff(swrap) / c
      j <- round((i[1:length(i)-1] + i[2:length(i)])/2)

      wavelen[j] <- wavelen1   # defined in between the wrap points
    }
  }
  wavelen
}
