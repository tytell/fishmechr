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
      m <- lm(ph ~ t, na.action = na.omit)
      freq <- coefficients(m)[2] / mod
    }
    else {
      freq <- NA
    }
  }
  freq
}

get_cycle_numbers <- function(ph, unwrap=FALSE, mod=2*pi)
{
  if (unwrap) {
    ph <- skip_na(ph, gsignal::unwrap)
  }

  cyc <- floor(ph / mod)
}

get_wavelength <- function(s, ph, unwrap=TRUE, method='deriv',
                           ignore_s=NULL,
                           mod=2*pi, traveling_wave_dir=-1)
{
  assertthat::assert_that(method %in% c('deriv', 'slope', 'cycle', 'halfcycle'))

  if (unwrap) {
    ph <- skip_na(ph, gsignal::unwrap)
  }

  if (!is.null(ignore_s)) {
    bad <- ignore_s(s)
  } else {
    bad <- rep(FALSE, length(s))
  }
  good <- !bad

  if (method == "deriv") {
    k <- deriv(s, ph) / mod
    wavelen <- traveling_wave_dir / k      # defined at all points
    wavelen[bad] <- NA
  }
  else if (method == "slope") {
    m <- lm(ph[good] ~ s[good])
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
    k <- which.max(ph)
    ph <- ph - ph[k]
    ph <- ph %% m
    ph[k] <- m

    i <- which(diff(ph) < 0)

    wavelen <- rep(NA_real_, length(s))
    if (length(i) >= 1) {
      # interpolate to get the exact location of the wrap

      swrap1 <- s[i]
      swrap2 <- s[i+1]
      phwrap1 <- ph[i]
      phwrap2 <- ph[i+1] + m

      swrap <- swrap1 + (swrap2 - swrap1) / (phwrap2 - phwrap1) * (m - phwrap1)

      swrap <- c(swrap, s[length(s)])
      i <- c(i, length(s))

      wavelen1 <- diff(swrap) / c
      j <- round((i[1:length(i)-1] + i[2:length(i)])/2)

      wavelen[j] <- wavelen1   # defined in between the wrap points
    }
  }
  wavelen
}
