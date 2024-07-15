#' Skip NAs when running a function on a vector
#'
#' `skip_na()` is a helper function related to [na.omit()]. It runs a function
#' `f` on a vector that may contain NAs or NaNs, skipping all the NAs, and
#' returns the results as a vector of the same length as `x` with the NAs
#' in the same places.
#'
#' @param x A vector that may have NAs
#' @param f The function to run on the vector
#' @param ... Other parameters to supply to the function
#'
#' @return A vector with the same length as `x` with NAs in the same places
#' @export
skip_na <- function(x, f, min.len = 1, ...)
{
  good <- is.finite(x)

  if (sum(good) >= min.len) {
    xf0 <- f(x[good], ...)
    xf <- vector(typeof(xf0), length = length(x))

    xf[good] <- xf0
    xf[!good] <- NA
  } else {
    xf <- vector(typeof(x), length = length(x))
    xf[!good] <- NA
  }

  xf
}

#' Constructs a smoothing filter
#'
#' Wrapper function for `gsignal::butter` to construct a bandpass (or
#' low or high pass) filter at a particular frequency. Uses the Sos form for
#' the filter, which works better numerically particularly for very low
#' frequencies.
#'
#' @param lo Lower frequency cutoff in Hz. If only `lo` is specified, then
#'   `build_filter` constructs a high pass filter
#' @param hi Upper frequency cutoff in Hz. If only `hi` is specified, then
#'   `build_filter` constructs a high pass filter
#' @param sampfreq Sampling frequency for the data in Hz.
#' @param order (optional) Order for the filter
#'
#' @return Filter parameters in Sos form
#' @export
build_filter <- function(lo=NULL, hi=NULL, sampfreq, order=13)
{
  assertthat::assert_that(!(is.null(lo) & is.null(hi)),
                          msg = "You must provide at least one of lo or hi")

  # Nyquist frequency
  nyqfreq <- sampfreq / 2

  if (is.null(lo)) {
    cutoff <- hi
    type <- "low"
  }
  else if (is.null(hi)) {
    cutoff <- lo
    type <- "high"
  } else {
    cutoff <- c(lo, hi)
    type <- "pass"
  }

  gsignal::butter(order, cutoff / nyqfreq, type=type, output="Sos")
}

#' Apply a filter constructed with [build_filter]
#'
#' Wrapper function for `gsignal::filtfilt` to apply a filter to a dataset.
#' Potentially skips NAs in the data set (see [skip_na()] for details).
#'
#' @param filt Filter as returned from [build_filter()] or [gsignal::butter()]
#' @param x Vector of data to filter
#' @param na.skip (TRUE or FALSE) to skip NAs in the data set.
#'
#' @return Filtered data set
#' @export
apply_filter <- function(filt, x, na.skip=TRUE)
{
  if (na.skip) {
    skip_na(x, \(y) gsignal::filtfilt(filt, y))
  } else {
    gsignal::filtfilt(filt, x)
  }
}

check_if_overwrite_columns <- function(df, newcols, overwrite)
{
  if (any(newcols %in% colnames(df))) {
    dfname <- "Data frame"
    callingfcn <- as.list(sys.call(-1))[[1]]

    if (overwrite) {
      warning(dfname, " has columns that are assigned in '", callingfcn, "'. Overwriting")
    } else {
      stop(dfname, " has columns that are assigned in '", callingfcn, "'. Stopping")
    }
  }
}

check.out <- function(.data, .out, .out_default, overwrite)
{
  if (is.null(.out)) {
    .out <- .out_default
  } else if (length(.out) < length(.out_default)) {
    if (is.null(names(.out))) {
      stop(".out must have ", length(.out_default),
           "elements or be a named list for specific output columns")
    }
    if (!all(names(.out) %in% names(.out_default))) {
      nout <- names(.out)
      extra = !(nout %in% names(.out_default))

      nout_str = paste(nout[extra], collapse = ",")
      warning('Some names in .out (', nout_str, ') are not used')
    }
    for (n in names(.out_default)) {
      if (!(n %in% names(.out))) {
        .out[[n]] <- .out_default[[n]]
      }
    }
  }

  if (!is.null(names(.out))) {
    .out <- .out[names(.out_default)]
  }

  if (any(.out %in% colnames(.data))) {
    overwritenames <- .out[.out %in% colnames(.data)]
    overwritestr <- paste(overwritenames, collapse = ',')
    if (overwrite) {
      warning("Columns (", overwritestr, ") will be overwritten")
    } else {
      warning("Columns (", overwritestr, ") already present in data frame. Stopping.")
    }
  }

  .out
}
