test_that("simulate_data returns correct structure", {
  set.seed(42)
  df <- simulate_data(scale = 10, n = 100, label = 1)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("observations", "status", "label"))
  expect_equal(nrow(df), 100)
})

test_that("observations are non-negative and status is binary", {
  set.seed(42)
  df <- simulate_data(scale = 5, n = 200, accrual_time = 6, follow_up_time = 12)

  expect_true(all(df$observations >= 0))
  expect_true(all(df$status %in% c(0L, 1L)))
})

test_that("label column matches argument", {
  set.seed(42)
  df0 <- simulate_data(scale = 8, n = 50, label = 0)
  df1 <- simulate_data(scale = 8, n = 50, label = 1)

  expect_true(all(df0$label == 0))
  expect_true(all(df1$label == 1))
})

test_that("no censoring when scale_loss = NULL and follow_up_time = Inf", {
  set.seed(42)
  df <- simulate_data(scale = 5, n = 200)

  expect_true(all(df$status == 1))
})

test_that("administrative censoring bounds observation times", {
  set.seed(42)
  accrual <- 6
  follow_up <- 12
  df <- simulate_data(
    scale = 5, n = 500,
    accrual_time = accrual, follow_up_time = follow_up
  )

  expect_true(all(df$observations <= accrual + follow_up))
})

test_that("censor_beyond_tau = TRUE caps all observations at tau", {
  set.seed(42)
  tau <- 8
  df <- simulate_data(
    scale = 10, n = 300,
    tau = tau, censor_beyond_tau = TRUE
  )

  expect_true(all(df$observations <= tau))
})

test_that("censor_beyond_tau = FALSE does not cap observations at tau", {
  set.seed(42)
  # short tau with large scale ensures many obs exceed tau
  df <- simulate_data(scale = 20, n = 500, tau = 1, censor_beyond_tau = FALSE)

  expect_true(any(df$observations > 1))
})

test_that("heavy loss to follow-up produces mostly censored observations", {
  set.seed(42)
  df <- simulate_data(scale = 100, scale_loss = 0.01, n = 300)

  expect_gt(mean(df$status == 0), 0.95)
})

test_that("parameterisation = 2 produces valid output", {
  set.seed(42)
  df <- simulate_data(
    scale = 0.01, shape = 2, n = 100,
    parameterisation = 2,
    accrual_time = 6, follow_up_time = 12
  )

  expect_equal(nrow(df), 100)
  expect_true(all(df$observations >= 0))
  expect_true(all(df$status %in% c(0L, 1L)))
})

test_that("parameterisation = 3 produces valid output", {
  set.seed(42)
  df <- simulate_data(
    scale = 0.1, shape = 2, n = 100,
    parameterisation = 3,
    accrual_time = 6, follow_up_time = 12
  )

  expect_equal(nrow(df), 100)
  expect_true(all(df$observations >= 0))
  expect_true(all(df$status %in% c(0L, 1L)))
})

test_that("piecewise Weibull with breakpoints produces valid output", {
  set.seed(42)
  df <- simulate_data(
    scale = c(10, 5), shape = c(1, 2),
    breakpoints = 6,
    n = 200, accrual_time = 6, follow_up_time = 12
  )

  expect_equal(nrow(df), 200)
  expect_true(all(df$observations >= 0))
  expect_true(all(df$status %in% c(0L, 1L)))
})

test_that("scalar shape is expanded for piecewise scale (exponential path)", {
  set.seed(42)
  # shape = 1 triggers int_rpexp; two-piece scale exercises the expansion
  df <- simulate_data(
    scale = c(8, 4), shape = 1,
    breakpoints = 5,
    n = 100
  )

  expect_equal(nrow(df), 100)
  expect_true(all(df$observations >= 0))
  expect_true(all(df$status == 1))  # no censoring applied
})

test_that("piecewise loss to follow-up produces valid output", {
  set.seed(42)
  # scalar shape_loss = 1 expanded to match two-piece scale_loss
  df <- simulate_data(
    scale = 10, n = 200,
    scale_loss = c(20, 10), shape_loss = 1,
    breakpoints_loss = 5
  )

  expect_equal(nrow(df), 200)
  expect_true(all(df$observations >= 0))
  expect_true(all(df$status %in% c(0L, 1L)))
  # with finite loss, some censoring should occur
  expect_true(any(df$status == 0))
})

test_that("n = 1 edge case returns a single-row data frame", {
  set.seed(42)
  df <- simulate_data(scale = 5, n = 1)

  expect_equal(nrow(df), 1)
  expect_true(df$status %in% c(0L, 1L))
})
