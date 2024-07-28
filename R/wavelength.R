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
