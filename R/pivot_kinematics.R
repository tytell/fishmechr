#' Pivots a kinematics dataset into long format
#'
#' Converts a wide-format data frame where each point's variables are separate
#' columns (e.g. `head.x`, `head.y`, `tail.x`, `tail.y`) into a long format
#' with one row per point per frame. The point column is returned as a factor
#' with levels in the order given by `pointnames`.
#'
#' @param df The data frame
#' @param pointnames The names of the points, in order from head to tail
#' @param sep The separator between the point name and variable name, as a
#'   regular expression (default = `"\\."`)
#' @param point_to The name of the column to put the point names in
#'
#' @returns The long data set
#' @export
#'
#' @examples
#' df <- data.frame(
#'   frame = 1:3,
#'   head.x = c(0, 0, 0), head.y = c(0, 1, 2),
#'   tail.x = c(5, 5, 5), tail.y = c(0, 1, 2)
#' )
#' pivot_kinematics_longer(df, c("head", "tail"))
pivot_kinematics_longer <- function(df, pointnames,
                                    point_to = "point",
                                    sep = '\\.')
{
  tidyr::pivot_longer(df, cols = contains(pointnames),
               names_to = c(point_to, ".value"),
               names_sep = sep) |>
    mutate("{point_to}" := factor(.data[[point_to]], levels=pointnames))
}
