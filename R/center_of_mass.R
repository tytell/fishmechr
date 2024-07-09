get_center_df <- function(df, t, s, x, y, m = NULL, w = NULL, h = NULL,
                       method = "mutate")
{
  assertthat::assert_that(method %in% c("mutate", "summarize", "summarise"))

  if (method == "mutate")
    fcn <- dplyr::mutate
  else if (method %in% c("summarize", "summarise"))
    fcn <- dplyr::summarise

  if (!is.null(m)) {
    message("Estimating true center of mass based on mass distribution")

    df |>
      group_by({{t}}) |>
      fcn(M = sum({{m}}, na.rm = TRUE),
          xcom = sum({{m}} * ({{x}} + dplyr::lead({{x}})), na.rm=TRUE) / (2*M),
          ycom = sum({{m}} * ({{y}} + dplyr::lead({{y}})), na.rm=TRUE) / (2*M))
  }
  else if (!is.null(h) & !is.null(w)) {
    message("Estimating center of mass based on width and height")

    df |>
      group_by({{t}}) |>
      fcn(ds = lead({{s}}) - {{s}},
          dw = lead({{w}}) - {{w}},
          dh = lead({{h}}) - {{h}},
          V = pi * ds * ({{w}} * {{h}} + 0.5*dw*{{h}} + 0.5*dh*{{w}} + 0.333333*dw*dh),
          sumV = sum(V, na.rm = TRUE),

          xcom = sum({{V}} * ({{x}} + dplyr::lead({{x}})), na.rm=TRUE) / (2*sumV),
          ycom = sum({{V}} * ({{y}} + dplyr::lead({{y}})), na.rm=TRUE) / (2*sumV))
  }
  else if (!is.null(w) & is.null(h)) {
    message("Estimating center of mass based on width")

    df |>
      group_by({{t}}) |>
      fcn(ds = lead({{s}}) - {{s}},
          xcom = 0.25 * sum(ds * ({{x}} + dplyr::lead({{x}})) * ({{w}} + dplyr::lead({{w}})), na.rm=TRUE),
          ycom = 0.25 * sum(ds * ({{y}} + dplyr::lead({{y}})) * ({{w}} + dplyr::lead({{w}})), na.rm=TRUE))
  } else {
    message("Estimating center of mass as the centroid of x and y")
    df |>
      group_by({{t}}) |>
      fcn(xcom = mean({{x}}, na.rm = TRUE),
          ycom = mean({{y}}, na.rm = TRUE))
  }
}

interpolate_width <- function(arclen0, width0, arclen)
{
  len0 <- max(arclen0, na.rm = TRUE)
  len1 <- max(arclen, na.rm = TRUE)

  if (sum(!is.na(arclen)) > 2) {
    xy <- approx(arclen0/len0*len1, width0/len0*len1, xout=arclen)
    xy$y
  } else {
    rep(NA, length(arclen))
  }
}
