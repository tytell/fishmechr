% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cycle_frequency.R
\name{get_cycle_numbers}
\alias{get_cycle_numbers}
\title{Gets cycle numbers from a phase variable}
\usage{
get_cycle_numbers(
  ph,
  unwrap = FALSE,
  mod = 2 * pi,
  exclude_partial_cycles = TRUE
)
}
\arguments{
\item{ph}{Phase variable}

\item{unwrap}{(TRUE or FALSE) Unwrap the phase variable. Note that the
function will not work unless the phase is unwrapped, so you should only
set unwrap to FALSE if the phase has been unwrapped earlier.}

\item{mod}{Modulus for the phase. Default is \code{2*pi}.}

\item{exclude_partial_cycles}{(TRUE or FALSE) Exclude cycles in which the
phase does not advance from close to 0 to close to 2pi.}
}
\value{
Integer cycle numbers with the same size as \code{ph}
}
\description{
Given a phase variable that increases, each time the phase passes through
2 k pi, a new cycle starts. This function unwraps the phase, so that it
increases steadily, then takes the floor of the phase divided by 2 pi (or
another modulus), so that it gets an integer for each cycle.
}
\details{
Optionally, it will try to exclude partial cycles, setting the cycle number
to NA for cycles that do not start at a phase close to 0 and progress to a
phase close to 2 pi. "Close" here is defined based on the average change in
phase.
}
\examples{
# example phase that advances by slightly more than 3 cycles, modulo 2pi
ph <- seq(0, 20, by = pi/10) \%\% (2*pi)
get_cycle_numbers(ph, unwrap=TRUE)

}
\seealso{
\code{\link[=get_body_cycle_numbers_df]{get_body_cycle_numbers_df()}}
}
