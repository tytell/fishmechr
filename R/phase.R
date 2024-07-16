
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

intepolate_peak_location <- function(y, x = c(-1, 0, 1))
{
  # formula from https://math.stackexchange.com/questions/889569/finding-a-parabola-from-three-points-algebraically
  den <- (x[1] - x[2])*(x[1] - x[3])*(x[2] - x[3])
  A <- (x[3]*(y[2] - y[1]) + x[2]*(y[1] - y[3]) + x[1]*(y[3] - y[2])) / den
  B <- (x[1]^2*(y[2] - y[3]) + x[3]^2*(y[1] - y[2]) + x[2]^2*(y[3] - y[1])) / den
  # C <- (x[1]^2*(x[2]*y[3] - x[3]*y[2]) +
  #         x[2]^2*(x[3]*y[1] - x[1]*y[3]) +
  #         x[3]^2*(x[1]*y[2] - x[2]*y[1])) / den

  xv <- -B / (2*A)

  xv
}

peak_phase <- function(x, na.skip=TRUE, unwrap=TRUE,
                       check_reasonableness=TRUE,
                       interpolate_peaks=TRUE, interpolate_zeros=TRUE,
                       zero_mode='midpoint',
                       ...)
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

  pkhi <- pracma::findpeaks(x, ...)
  pklo <- pracma::findpeaks(-x, ...)

  pk <- c(pkhi[,1], pklo[,1])
  ipk <- c(pkhi[,2], pklo[,2])
  pk <- sort_by(pk, ipk)
  ipk <- sort(ipk)

  ipkoff <- rep_len(0, length(ipk))

  for (j in seq_along(ipk)) {
    ipk1 <- ipk[j]
    ipkoff[j] <- interpolate_peak_location(x[(ipk1-1):(ipk1+1)])
  }

  if (zero_mode == 'midpoint') {
    izero <- rep_len(0, length(ipk)-1)
    izerooff <- rep_len(0, length(ipk)-1)
    for (j in seq(1, length(ipk)-1)) {
      half <- (pk[j] + pk[j+1])/2

      ind <- seq(ipk[j], ipk[j+1]-1)
      if (pk[j] > pk[j+1]) {
        iz1 <- which((x[ind] >= half) & (x[ind+1] < half))
      } else {
        iz1 <- which((x[ind] <= half) & (x[ind+1] > half))
      }
      iz1 <- iz1 - 1 + ipk[j]
      izero[j] <- iz1

      off1 <- 1 / (x[iz1+1] - x[iz1]) * (half - x[iz1])
      izerooff[j] <- off1
    }
  }
  else if (zero_mode == 'zero') {
    izerodown <- which((x[1:length(x)-1] >= 0) & (x[2:length(x)] < 0))
    izeroup <- which((x[1:length(x)-1] <= 0) & (x[2:length(x)] > 0))

    izero <- c(izerodown, izeroup)
    izero <- sort(izero)

    izerooff <- rep_len(0, length(izero))
    for (j in seq_along(izero)) {
      iz1 <- izero[j]
      off1 <- 1 / (x[iz1+1] - x[iz1]) * (0 - x[iz1])
      izerooff[j] <- off1
    }
  }

  ## CONTINUE here
  # Interpolate phase based on peaks and zeros
}
