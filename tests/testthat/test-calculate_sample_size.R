# Shared base args used across tests (plots suppressed throughout)
base_args <- list(
  scale_ctrl = 6,
  scale_trmt = 10,
  accrual_time = 6,
  follow_up_time = 3,
  tau = 4,
  scale_loss = 10,
  plot_example_data = FALSE,
  plot_design_curves = FALSE
)

# output structure ------------------------------------------------------------

test_that("calculate_sample_size returns a list containing only populated elements", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8)))

  expect_type(result, "list")
  # SS elements present (RMSTD and LRT closed form on by default)
  expect_true("Sample size for RMST difference determined by closed-form solution" %in% names(result))
  expect_true("Sample size for LRT determined by closed-form solution" %in% names(result))
  # power elements absent (power is given, so it is not computed)
  expect_false("Power for RMST difference determined by closed-form solution" %in% names(result))
  expect_false("Power for LRT determined by closed-form solution" %in% names(result))
  # simulation and Satterthwaite elements absent by default
  expect_false("RMSTD power determined by simulation" %in% names(result))
  expect_false("Satterthwaite-corrected sample size for RMST difference" %in% names(result))
})

# closed-form flags -----------------------------------------------------------

test_that("disabled closed-form flags leave SS elements absent from result", {
  result <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8,
    RMSTD_closed_form = FALSE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = FALSE
  )))

  expect_false("Sample size for RMST difference determined by closed-form solution" %in% names(result))
  expect_false("Sample size for RMST ratio determined by closed-form solution" %in% names(result))
  expect_false("Sample size for LRT determined by closed-form solution" %in% names(result))
  # RMST values also absent when no closed-form method requires them
  expect_false("RMST treatment group" %in% names(result))
  expect_false("RMST control group" %in% names(result))
})

test_that("enabled closed-form flags produce positive finite SS", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8)))

  ss_rmstd <- result[["Sample size for RMST difference determined by closed-form solution"]]
  ss_lrt   <- result[["Sample size for LRT determined by closed-form solution"]]

  expect_true(is.finite(ss_rmstd) && ss_rmstd > 0)
  expect_true(is.finite(ss_lrt)   && ss_lrt   > 0)
})

test_that("RMSTR closed-form produces a positive finite SS when enabled", {
  result <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8, RMSTR_closed_form = TRUE
  )))

  expect_true("Sample size for RMST ratio determined by closed-form solution" %in% names(result))
  ss_rmstr <- result[["Sample size for RMST ratio determined by closed-form solution"]]
  expect_true(is.finite(ss_rmstr) && ss_rmstr > 0)
})

# Satterthwaite correction ----------------------------------------------------

test_that("Satterthwaite SS is absent without correction, present and positive with it", {
  result_no_sat <- do.call(calculate_sample_size, c(base_args, list(power = 0.8, satterthwaite_corr = FALSE)))
  result_sat    <- do.call(calculate_sample_size, c(base_args, list(power = 0.8, satterthwaite_corr = TRUE)))

  expect_false("Satterthwaite-corrected sample size for RMST difference" %in% names(result_no_sat))

  expect_true("Satterthwaite-corrected sample size for RMST difference" %in% names(result_sat))
  ss_sat <- result_sat[["Satterthwaite-corrected sample size for RMST difference"]]
  expect_true(is.finite(ss_sat) && ss_sat > 0)
})

test_that("Satterthwaite SS differs from standard SS", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8, satterthwaite_corr = TRUE)))

  ss_standard <- result[["Sample size for RMST difference determined by closed-form solution"]]
  ss_sat      <- result[["Satterthwaite-corrected sample size for RMST difference"]]

  expect_false(isTRUE(all.equal(ss_standard, ss_sat)))
})

# RMST values -----------------------------------------------------------------

test_that("RMST outputs are present and ordered correctly when RMSTD_closed_form = TRUE", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8)))

  expect_true("RMST treatment group" %in% names(result))
  expect_true("RMST control group" %in% names(result))
  expect_true("RMST difference" %in% names(result))
  expect_true("RMST ratio" %in% names(result))
  expect_gt(result[["RMST treatment group"]], result[["RMST control group"]])
  expect_gt(result[["RMST difference"]], 0)
  expect_gt(result[["RMST ratio"]], 1)
})

# SS monotonicity -------------------------------------------------------------

test_that("SS increases with higher target power", {
  ss_low  <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.6, LRT_closed_form = FALSE
  )))[["Sample size for RMST difference determined by closed-form solution"]]

  ss_high <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.9, LRT_closed_form = FALSE
  )))[["Sample size for RMST difference determined by closed-form solution"]]

  expect_lt(ss_low, ss_high)
})

test_that("SS decreases with larger treatment effect", {
  ss_small_effect <- do.call(calculate_sample_size, c(
    modifyList(base_args, list(scale_trmt = 7)),
    list(power = 0.8, LRT_closed_form = FALSE)
  ))[["Sample size for RMST difference determined by closed-form solution"]]

  ss_large_effect <- do.call(calculate_sample_size, c(
    modifyList(base_args, list(scale_trmt = 15)),
    list(power = 0.8, LRT_closed_form = FALSE)
  ))[["Sample size for RMST difference determined by closed-form solution"]]

  expect_gt(ss_small_effect, ss_large_effect)
})

test_that("two-sided test requires larger SS than one-sided", {
  ss_1 <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8, sides = 1, LRT_closed_form = FALSE
  )))[["Sample size for RMST difference determined by closed-form solution"]]

  ss_2 <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8, sides = 2, LRT_closed_form = FALSE
  )))[["Sample size for RMST difference determined by closed-form solution"]]

  expect_gt(ss_2, ss_1)
})

# consistency with calculate_power --------------------------------------------

test_that("calculate_power at SS from calculate_sample_size recovers target power of 0.8", {
  n_target <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8,
    RMSTD_closed_form = TRUE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = FALSE
  )))[["Sample size for RMST difference determined by closed-form solution"]]

  pwr_result <- do.call(calculate_power, c(base_args, list(
    n = n_target,
    RMSTD_closed_form = TRUE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = FALSE
  )))

  expect_equal(
    pwr_result[["Power for RMST difference determined by closed-form solution"]],
    0.8,
    tolerance = 1e-6
  )
})

# regression values -----------------------------------------------------------

test_that("SS matches reference values", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8)))

  expect_equal(result[["Sample size for RMST difference determined by closed-form solution"]], 404.5066,  tolerance = 1e-3)
  expect_equal(result[["Sample size for LRT determined by closed-form solution"]],             294.4213,  tolerance = 1e-3)
  expect_equal(result[["RMST control group"]],   2.919497,  tolerance = 1e-5)
  expect_equal(result[["RMST treatment group"]], 3.2968,    tolerance = 1e-5)
  expect_equal(result[["RMST difference"]],      0.3773023, tolerance = 1e-5)
  expect_equal(result[["RMST ratio"]],           1.129235,  tolerance = 1e-5)
})

test_that("Satterthwaite-corrected SS matches reference value", {
  result <- do.call(calculate_sample_size, c(base_args, list(power = 0.8, satterthwaite_corr = TRUE)))

  expect_equal(
    result[["Satterthwaite-corrected sample size for RMST difference"]],
    407.652,
    tolerance = 1e-3
  )
})
