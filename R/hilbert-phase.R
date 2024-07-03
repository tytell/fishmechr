
hilbert_phase <- function(x, na.skip=TRUE, unwrap=TRUE,
                          check_reasonableness=TRUE)
{
  if (check_reasonableness) {
    finitefrac <- sum(is.finite(x)) / length(x)
    if (finitefrac > 0 && finitefrac < 0.9)
      warning("A large fraction of values are NA. Phase estimate may work poorly")

    if (finitefrac > 0) {
      m <- mean(x, na.rm=TRUE)
      s <- sd(x, na.rm=TRUE)

      if (abs(m) < s)
        warning("For a good phase estimate, the signal should oscillate around 0")
    }
  }

  # get the analytic signal
  if (na.skip) {
    A <- skip_na(x, gsignal::hilbert, min.len = 4)
  }
  else {
    A <- gsignal::hilbert(x)
  }

  # for some odd reason, R's Arg returns the phase angle as a complex number,
  # even though all the imaginary components are 0
  ph <- Re(Arg(A))
  ph <- skip_na(ph, gsignal::unwrap)

  if (check_reasonableness) {
    d <- diff(ph)
    fracneg <- sum(d < 0, na.rm=TRUE) / length(d)

    if (fracneg > 0.05)
      warning("Phase seems to go backwards a lot, which may indicate an overly noisy signal")
  }

  ph
}
