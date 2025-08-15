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
#'
#' @examples
#' # Low pass filter with a cutoff at 0.5Hz for data sampled at 100Hz
#' lopass <- build_filter(lo=0.5, sampfreq=100)
#'
#' # Band pass filter that passes frequencies between 0.5 and 10Hz
#' bandpass <- build_filter(lo=0.5, hi=10, sampfreq=100)
build_filter <- function(lo = NULL, hi = NULL, sampfreq, order = 13) {
  assertthat::assert_that(
    !(is.null(lo) & is.null(hi)),
    msg = "You must provide at least one of lo or hi"
  )

  # Nyquist frequency
  nyqfreq <- sampfreq / 2

  if (is.null(lo)) {
    cutoff <- hi
    type <- "low"
  } else if (is.null(hi)) {
    cutoff <- lo
    type <- "high"
  } else {
    cutoff <- c(lo, hi)
    type <- "pass"
  }

  gsignal::butter(order, cutoff / nyqfreq, type = type, output = "Sos")
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
apply_filter <- function(filt, x, na.skip = TRUE) {
  if (na.skip) {
    skip_na(x, \(y) gsignal::filtfilt(filt, y))
  } else {
    gsignal::filtfilt(filt, x)
  }
}
