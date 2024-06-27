#' Title
#'
#' @param lo
#' @param hi
#' @param sampfreq
#' @param order
#'
#' @return
#' @export
#'
#' @examples
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

apply_filter <- function(filt, x, na.skip=TRUE)
{
  if (na.skip) {
    skip_na(x, \(y) gsignal::filtfilt(filt, y))
  } else {
    gsignal::filtfilt(filt, x)
  }
}

hilbert_phase <- function(x, na.skip=TRUE, unwrap=TRUE,
                          check_reasonableness=TRUE)
{
  if (check_reasonableness) {
    m <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)

    assertthat::assert_that(abs(m) < s,
                            msg="For a good phase estimate, the signal should oscillate around 0")
  }

  # get the analytic signal
  if (na.skip) {
    A <- skip_na(x, gsignal::hilbert)
  }
  else {
    A <- gsignal::hilbert(x)
  }

  ph <- Arg(A)
  ph <- skip_na(A, gsignal::unwrap)

  if (check_reasonableness) {
    d <- diff(ph)
    fracneg <- sum(d < 0, na.rm=TRUE) / length(d)

    assertthat::assert_that(fracneg < 0.05,
                            msg="Phase seems to go backwards a lot, which may indicate an overly noisy signal")
  }

  ph
}
