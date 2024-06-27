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
#'
#' @examples
skip_na <- function(x, f, ...)
{
  good <- !is.na(x) & is.finite(x)
  xf <- purrr::rep_along(NA, x)

  xf[good] <- f(x[good], ...)
  xf
}
