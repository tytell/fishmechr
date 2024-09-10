#' Pivots a kinematics dataset into long format
#'
#' @param df The data frame
#' @param pointnames The names of the points, in order from head to tail
#' @param sep The separator, as a regular expression (default = "\\.")
#' @param point_to The name of the column to put the point names in
#'
#' @returns The long data set
#' @export
#'
#' @examples
pivot_kinematics_longer <- function(df, pointnames,
                                    point_to = "point",
                                    sep = '\\.')
{
  tidyr::pivot_longer(df, cols = contains(pointnames),
               names_to = c(point_to, ".value"),
               names_sep = sep) |>
    dplyr::mutate("{point_to}" := factor(.data[[point_to]], levels=pointnames))
}
