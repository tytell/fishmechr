#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across any_of bind_cols case_when filter first group_by
#'   if_else is_grouped_df lag lead left_join mutate n pick rename row_number
#'   select summarise summarize ungroup
#' @importFrom rlang := .data as_name enquo eval_tidy quo_name
#' @importFrom stats approx coefficients lm median na.omit predict sd smooth.spline
#' @importFrom tidyr fill unnest
#' @importFrom tidyselect contains
## usethis namespace: end

utils::globalVariables(c(
  # smooth.R: dplyr NSE column names
  "frame", "point", "gapstart", "gapend", "gaplen", "good",
  # cycle_frequency.R
  "cyc", "lo", "hi", "cycle",
  # center_of_mass.R
  "M", "V", "sumV", "sumw", "nsum",
  # excursion.R
  "xctr", "xsd", "yctr", "ysd",
  "swimaxis_x", "swimaxis_y", "swimaxis_mag",
  # excursion.R / curvature.R
  "xys"
))

NULL
