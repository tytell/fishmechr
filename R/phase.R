
#' Compute phase of an oscillation using the Hilbert transform
#'
#' @param x
#' @param na.skip
#' @param unwrap
#' @param check_reasonableness
#'
#' @returns Phase
#' @export
#'
#' @examples
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

#' Interpolate the location of a peak based on three points
#'
#' Uses a parabolic approximation to determine the location of a
#' peak from 3 points.
#'
#' @param y
#' @param x
#'
#' @return
#' @export
#'
#' @examples
interpolate_peak_location <- function(y, x = c(-1, 0, 1))
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

#' Compute phase of an oscillation by locating peaks and zero crossings
#'
#' @param x
#' @param unwrap
#' @param check_reasonableness
#' @param check_ordering
#' @param interpolate_peaks
#' @param interpolate_zeros
#' @param zero_mode
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
peak_phase <- function(x, unwrap=TRUE,
                       check_reasonableness=TRUE,
                       check_ordering=TRUE,
                       interpolate_peaks=TRUE, interpolate_zeros=TRUE,
                       zero_mode='midpoint',
                       ...)
{
  if (check_reasonableness) {
    finitefrac <- sum(is.finite(x)) / length(x)
    if (finitefrac > 0 && finitefrac < 0.9)
      warning("A large fraction of values are NA. Phase estimate may work poorly")

    if ((finitefrac > 0) && (zero_mode == 'zero')) {
      m <- mean(x, na.rm=TRUE)
      s <- sd(x, na.rm=TRUE)

      if (abs(m) < s)
        warning("For a good phase estimate, the signal should oscillate around 0")
    }
  }

  pkhi <- pracma::findpeaks(x, ...)
  pklo <- pracma::findpeaks(-x, ...)

  signlevels <- c('hi', 'down', 'lo', 'up')

  pk <- c(pkhi[,1], -pklo[,1])
  ipk <- c(pkhi[,2], pklo[,2])
  pksign <- factor(c(rep_len('hi', nrow(pkhi)), rep_len('lo', nrow(pklo))),
                   levels = signlevels)
  pk <- sort_by(pk, ipk)
  pksign <- sort_by(pksign, ipk)
  ipk <- sort(ipk)

  ipkoff <- rep_len(0, length(ipk))

  if (interpolate_peaks) {
    for (j in seq_along(ipk)) {
      ipk1 <- ipk[j]
      ipkoff[j] <- interpolate_peak_location(x[(ipk1-1):(ipk1+1)])
    }
    good <- (ipkoff >= -1) & (ipkoff <= 1)

    pk <- pk[good]
    ipk <- ipk[good]
    pksign <- pksign[good]
    ipkoff <- ipkoff[good]
  }

  if (is.na(zero_mode) || (zero_mode == 'none')) {
    izero <- rep_len(NA, length(ipk)-1)
    izerooff <- rep_len(0, length(ipk)-1)
    zerosign <- factor(rep_len(NA, length(ipk)-1), levels = signlevels)
  }
  else if (zero_mode == 'midpoint') {
    izero <- rep_len(NA, length(ipk)-1)
    izerooff <- rep_len(0, length(ipk)-1)
    zerosign <- factor(rep_len(NA, length(ipk)-1), levels = signlevels)
    for (j in seq(1, length(ipk)-1)) {
      if (pksign[j] != pksign[j+1]) {
        half <- (pk[j] + pk[j+1])/2

        ind <- seq(ipk[j], ipk[j+1])
        if (pk[j] > pk[j+1]) {
          iz1 <- which((x[ind] >= half) & (x[ind+1] < half))
          zs1 <- 'down'
        } else {
          iz1 <- which((x[ind] <= half) & (x[ind+1] > half))
          zs1 <- 'up'
        }

        if (length(iz1) == 1) {
          # we shouldn't ever get iz1 being anything other than length 1,
          # but just to be careful, we check...
          iz1 <- iz1 - 1 + ipk[j]
          izero[j] <- iz1

          if (interpolate_zeros) {
            off1 <- 1 / (x[iz1+1] - x[iz1]) * (half - x[iz1])
            izerooff[j] <- off1
          }
          zerosign[j] <- zs1
        }
      }
    }
  }
  else if (zero_mode == 'zero') {
    izerodown <- which((x[1:length(x)-1] >= 0) & (x[2:length(x)] < 0))
    izeroup <- which((x[1:length(x)-1] <= 0) & (x[2:length(x)] > 0))

    izero <- c(izerodown, izeroup)
    zerosign <- factor(c(rep_len('down', length(izerodown)), rep_len('up', length(izeroup))),
                        levels = signlevels)
    zerosign <- sort_by(izerosign, izero)
    izero <- sort(izero)

    izerooff <- rep_len(0, length(izero))
    if (interpolate_zeros) {
      for (j in seq_along(izero)) {
        iz1 <- izero[j]
        off1 <- 1 / (x[iz1+1] - x[iz1]) * (0 - x[iz1])
        izerooff[j] <- off1
      }
    }
  }

  # Interpolate phase based on peaks and zeros
  ind <- c(ipk + ipkoff, izero + izerooff)
  sgn <- c(pksign, zerosign)

  sgn <- sgn[!is.na(ind)]
  ind <- ind[!is.na(ind)]

  sgn <- sort_by(sgn, ind)
  ind <- sort(ind)

  ph0 <- dplyr::case_when(sgn == 'hi'  ~  0,
                        sgn == 'down'  ~  pi/2,
                        sgn == 'lo'  ~  pi,
                        sgn == 'up'  ~  3*pi/2,
                        .default = NA)

  if (check_ordering) {
    good <- dplyr::case_when(sgn == 'hi'  ~ lag(sgn) == 'up' & lead(sgn) == 'down',
                             sgn == 'down'  ~  lag(sgn) == 'hi' & lead(sgn) == 'lo',
                             sgn == 'lo'  ~  lag(sgn) == 'down' & lead(sgn) == 'up',
                             sgn == 'up'  ~  lag(sgn) == 'lo' & lead(sgn) == 'hi',
                             .default = FALSE)
    good[is.na(good)] <- TRUE
    ph0[!good] <- NA
  }
  else {
    good <- is.finite(ph0) & is.finite(ind)
    ph0 <- ph0[good]
    ind <- ind[good]
  }

  ph0 <- skip_na(ph0, gsignal::unwrap)

  span <- seq(ceiling(ind[1]), floor(ind[length(ind)]))

  ph <- rep_len(NA, length(x))
  apx <- approx(ind, ph0, xout = span)

  ph[span] <- apx$y

  if (!unwrap) {
    ph <- ph %% (2*pi)
  }
  ph
}
