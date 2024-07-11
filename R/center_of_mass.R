get_center_df <- function(.data, arclen, x, y, mass,width,height,
                          .out = c('xcom', 'ycom'),
                          .frame = frame, .point = point,
                          method = "mutate")
{
  assertthat::assert_that(method %in% c("mutate", "summarize", "summarise"))
  assertthat::assert_that(length(.out) == 2,
                          msg = ".out must contain two names, for the x and y center positions")

  if (missing(.frame)) {
    .frame <- enquo(.frame)
    assertthat::assert_that(assertthat::has_name(.data, rlang::as_name(.frame)),
                            msg = "Default column 'frame' not present. Use .frame to specify the name of the frame column")
  }
  if (missing(.point)) {
    .point <- enquo(.point)
  }

  if (method == "mutate")
    fcn <- dplyr::mutate
  else if (method %in% c("summarize", "summarise"))
    fcn <- dplyr::summarise

  if (!missing(mass)) {
    message("Estimating true center of mass based on mass distribution")
    mass <- enquo(mass)

    .data <- .data |>
      group_by({{.frame}}, .add = TRUE) |>
      fcn(M = sum(!!mass, na.rm = TRUE),
          "{.out[1]}" := sum(!!mass * ({{x}} + dplyr::lead({{x}})), na.rm=TRUE) / (2*M),
          "{.out[2]}" := sum(!!mass * ({{y}} + dplyr::lead({{y}})), na.rm=TRUE) / (2*M)) |>
      select(-c(M))
  }
  else if (!missing(height) & !missing(width)) {
    message("Estimating center of mass based on width and height")
    width <- enquo(width)
    height <- enquo(height)

    .data <- .data |>
      group_by({{.frame}}, .add = TRUE) |>
      fcn(ds = lead({{arclen}}) - {{arclen}},
          dw = lead(!!width) - !!width,
          dh = lead(!!height) - !!height,
          V = pi * ds * (!!width * !!height + 0.5*dw*!!height + 0.5*dh*!!width + 0.333333*dw*dh),
          sumV = sum(V, na.rm = TRUE),

          "{.out[1]}" := sum({{V}} * ({{x}} + dplyr::lead({{x}})), na.rm=TRUE) / (2*sumV),
          "{.out[2]}" := sum({{V}} * ({{y}} + dplyr::lead({{y}})), na.rm=TRUE) / (2*sumV)) |>
      select(-c(ds,dw,dh,V,sumV))
  }
  else if (!missing(width) & missing(height)) {
    message("Estimating center of mass based on width")
    width <- enquo(width)

    .data <- .data |>
      group_by({{.frame}}, .add = TRUE) |>
      fcn(sumw = sum(!!width, na.rm=TRUE),
          "{.out[1]}" := sum({{x}} * !!width) / sumw,
          "{.out[2]}" := sum({{y}} * !!width) / sumw) |>
      select(-sumw)
  } else {
    message("Estimating center of mass as the centroid of x and y")
    .data <- .data |>
      group_by({{.frame}}) |>
      fcn("{.out[1]}" := mean({{x}}, na.rm = TRUE),
          "{.out[2]}" := mean({{y}}, na.rm = TRUE))
  }

  if (method == "mutate") {
    .data <- .data |>
      ungroup({{.frame}})
  }

  .data
}

interpolate_width <- function(arclen0, width0, arclen)
{

  if (sum(!is.na(arclen)) > 2) {
    len0 <- max(arclen0, na.rm = TRUE)
    len1 <- max(arclen, na.rm = TRUE)

    xy <- approx(arclen0/len0*len1, width0/len0*len1, xout=arclen)
    xy$y
  } else {
    rep(NA, length(arclen))
  }
}
