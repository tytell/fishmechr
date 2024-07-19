test_that("check oval cylinder", {
  # oval cross section cylinder of length 2
  cyl <- tibble::tribble(
    ~arclen, ~width, ~height,
    0, 1, 0.5,
    1, 1, 0.5,
    2, 1, 0.5
  )

  expect_equal(
    get_volume(cyl, arclen, width, height) |>
      dplyr::summarize(V = sum(V, na.rm = TRUE)) |>
      dplyr::pull(V),
    pi * 2 * 1 * 0.5)
})

test_that("check oval cylinder with uneven segments", {
  # oval cross section cylinder of length 2
  cyl <- tibble::tribble(
    ~arclen, ~width, ~height,
    0, 1, 0.5,
    0.2, 1, 0.5,
    2, 1, 0.5
  )

  expect_equal(
    get_volume(cyl, arclen, width, height) |>
      dplyr::summarize(V = sum(V, na.rm = TRUE)) |>
      dplyr::pull(V),
    pi * 2 * 1 * 0.5)
})

test_that("check circular cone", {
  # formula for volume of a circular cone is
  # 1/3 pi r^2 h
  r <- 2
  h <- 4

  cyl <- tibble::tribble(
    ~arclen, ~width, ~height,
    0, r, r,
    h, 0, 0,
  )

  expect_equal(
    get_volume(cyl, arclen, width, height) |>
      dplyr::summarize(V = sum(V, na.rm = TRUE)) |>
      dplyr::pull(V),
    1/3 * pi * r^2 * h,
    tolerance = 1e-6)
})

test_that("check circular cone pointed the other way", {
  # formula for volume of a circular cone is
  # 1/3 pi r^2 h
  r <- 2
  h <- 4

  cyl <- tibble::tribble(
    ~arclen, ~width, ~height,
    0, 0, 0,
    h, r, r,
  )

  expect_equal(
    get_volume(cyl, arclen, width, height) |>
      dplyr::summarize(V = sum(V, na.rm = TRUE)) |>
      dplyr::pull(V),
    1/3 * pi * r^2 * h,
    tolerance = 1e-6)
})

test_that("check elliptical cone pointed the other way", {
  # formula for volume of an elliptical cross-section cone is
  # 1/3 pi a b h
  # https://mathworld.wolfram.com/EllipticCone.html
  a <- 2
  b <- 1
  h <- 4

  cyl <- tibble::tribble(
    ~arclen, ~width, ~height,
    0, 0, 0,
    h, a, b,
  )

  expect_equal(
    get_volume(cyl, arclen, width, height) |>
      dplyr::summarize(V = sum(V, na.rm = TRUE)) |>
      dplyr::pull(V),
    1/3 * pi * a * b * h,
    tolerance = 1e-6)
})




