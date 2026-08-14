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

test_that("calculate_power returns a list containing only populated elements", {
  result <- do.call(calculate_power, c(base_args, list(n = 405)))

  expect_type(result, "list")
  # power elements present (RMSTD and LRT closed form on by default)
  expect_true("Power for RMST difference determined by closed-form solution" %in% names(result))
  expect_true("Power for LRT determined by closed-form solution" %in% names(result))
  # sample size elements absent (n is given, so SS not computed)
  expect_false("Sample size for RMST difference determined by closed-form solution" %in% names(result))
  expect_false("Sample size for LRT determined by closed-form solution" %in% names(result))
  # simulation elements absent (simulations off by default)
  expect_false("RMSTD power determined by simulation" %in% names(result))
  expect_false("RMSTR power determined by simulation" %in% names(result))
  expect_false("LRT power determined by simulation" %in% names(result))
})

# closed-form flags -----------------------------------------------------------

test_that("disabled closed-form flags leave power elements absent from result", {
  result <- do.call(calculate_power, c(base_args, list(
    n = 405,
    RMSTD_closed_form = FALSE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = FALSE
  )))

  expect_false("Power for RMST difference determined by closed-form solution" %in% names(result))
  expect_false("Power for RMST ratio determined by closed-form solution" %in% names(result))
  expect_false("Power for LRT determined by closed-form solution" %in% names(result))
  # RMST values also absent when no closed-form method requires them
  expect_false("RMST treatment group" %in% names(result))
  expect_false("RMST control group" %in% names(result))
})

test_that("RMSTR closed-form power is present and in [0, 1] when enabled", {
  result <- do.call(calculate_power, c(base_args, list(n = 405, RMSTR_closed_form = TRUE)))

  expect_true("Power for RMST ratio determined by closed-form solution" %in% names(result))
  pwr <- result[["Power for RMST ratio determined by closed-form solution"]]
  expect_gte(pwr, 0)
  expect_lte(pwr, 1)
})

# Satterthwaite correction ----------------------------------------------------

test_that("Satterthwaite element absent without correction, present with it", {
  result_no_sat <- do.call(calculate_power, c(base_args, list(n = 405, satterthwaite_corr = FALSE)))
  result_sat    <- do.call(calculate_power, c(base_args, list(n = 405, satterthwaite_corr = TRUE)))

  expect_false("Satterthwaite-corrected power for RMST difference" %in% names(result_no_sat))

  expect_true("Satterthwaite-corrected power for RMST difference" %in% names(result_sat))
  pwr_sat <- result_sat[["Satterthwaite-corrected power for RMST difference"]]
  expect_gte(pwr_sat, 0)
  expect_lte(pwr_sat, 1)
})

# RMST values -----------------------------------------------------------------

test_that("RMST outputs are present and ordered correctly when RMSTD_closed_form = TRUE", {
  result <- do.call(calculate_power, c(base_args, list(n = 405)))

  expect_true("RMST treatment group" %in% names(result))
  expect_true("RMST control group" %in% names(result))
  expect_true("RMST difference" %in% names(result))
  expect_true("RMST ratio" %in% names(result))
  # larger scale -> longer survival -> higher RMST
  expect_gt(result[["RMST treatment group"]], result[["RMST control group"]])
  expect_gt(result[["RMST difference"]], 0)
  expect_gt(result[["RMST ratio"]], 1)
})

# power range and monotonicity ------------------------------------------------

test_that("closed-form powers are in [0, 1]", {
  result <- do.call(calculate_power, c(base_args, list(n = 405)))

  expect_gte(result[["Power for RMST difference determined by closed-form solution"]], 0)
  expect_lte(result[["Power for RMST difference determined by closed-form solution"]], 1)
  expect_gte(result[["Power for LRT determined by closed-form solution"]], 0)
  expect_lte(result[["Power for LRT determined by closed-form solution"]], 1)
})

test_that("power increases monotonically with sample size", {
  pwr_small <- do.call(calculate_power, c(base_args, list(
    n = 100, LRT_closed_form = FALSE
  )))[["Power for RMST difference determined by closed-form solution"]]

  pwr_large <- do.call(calculate_power, c(base_args, list(
    n = 1000, LRT_closed_form = FALSE
  )))[["Power for RMST difference determined by closed-form solution"]]

  expect_lt(pwr_small, pwr_large)
})

test_that("power approaches 1 for very large sample size", {
  result <- do.call(calculate_power, c(base_args, list(
    n = 1e6, LRT_closed_form = FALSE
  )))

  expect_gt(result[["Power for RMST difference determined by closed-form solution"]], 0.999)
})

test_that("power equals alpha under null effect (scale_ctrl == scale_trmt)", {
  null_args <- modifyList(base_args, list(scale_trmt = 6))
  result <- do.call(calculate_power, c(null_args, list(n = 500, LRT_closed_form = FALSE)))

  expect_equal(
    result[["Power for RMST difference determined by closed-form solution"]],
    0.025,
    tolerance = 1e-6
  )
})

test_that("two-sided test yields lower power than one-sided for same n", {
  pwr_1 <- do.call(calculate_power, c(base_args, list(
    n = 405, sides = 1, LRT_closed_form = FALSE
  )))[["Power for RMST difference determined by closed-form solution"]]

  pwr_2 <- do.call(calculate_power, c(base_args, list(
    n = 405, sides = 2, LRT_closed_form = FALSE
  )))[["Power for RMST difference determined by closed-form solution"]]

  expect_gt(pwr_1, pwr_2)
})

# consistency with calculate_sample_size -------------------------------------

test_that("power at n from calculate_sample_size recovers target power of 0.8", {
  ss_result <- do.call(calculate_sample_size, c(base_args, list(
    power = 0.8,
    RMSTD_closed_form = TRUE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = FALSE
  )))
  n_for_power <- ss_result[["Sample size for RMST difference determined by closed-form solution"]]

  pwr_result <- do.call(calculate_power, c(base_args, list(
    n = n_for_power,
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

test_that("closed-form power matches reference values", {
  result <- do.call(calculate_power, c(base_args, list(n = 405)))

  expect_equal(result[["Power for RMST difference determined by closed-form solution"]], 0.8004779, tolerance = 1e-5)
  expect_equal(result[["Power for LRT determined by closed-form solution"]],            0.90756,   tolerance = 1e-5)
  expect_equal(result[["RMST control group"]],   2.919497,  tolerance = 1e-5)
  expect_equal(result[["RMST treatment group"]], 3.2968,    tolerance = 1e-5)
  expect_equal(result[["RMST difference"]],      0.3773023, tolerance = 1e-5)
  expect_equal(result[["RMST ratio"]],           1.129235,  tolerance = 1e-5)
})

# simulation tests ------------------------------------------------------------

test_that("RMSTD simulation produces power in [0, 1] and is present in result", {
  skip_on_cran()
  set.seed(42)
  result <- do.call(calculate_power, c(base_args, list(
    n = 405,
    RMSTD_simulation = TRUE,
    M = 10
  )))

  expect_true("RMSTD power determined by simulation" %in% names(result))
  pwr <- result[["RMSTD power determined by simulation"]]
  expect_gte(pwr, 0)
  expect_lte(pwr, 1)
})

test_that("RMSTR simulation produces power in [0, 1] and is present in result", {
  skip_on_cran()
  set.seed(42)
  result <- do.call(calculate_power, c(base_args, list(
    n = 405,
    RMSTR_simulation = TRUE,
    M = 10
  )))

  expect_true("RMSTR power determined by simulation" %in% names(result))
  pwr <- result[["RMSTR power determined by simulation"]]
  expect_gte(pwr, 0)
  expect_lte(pwr, 1)
})

test_that("LRT simulation produces power in [0, 1] and is present in result", {
  skip_on_cran()
  set.seed(42)
  result <- do.call(calculate_power, c(base_args, list(
    n = 405,
    LRT_simulation = TRUE,
    M = 10
  )))

  expect_true("LRT power determined by simulation" %in% names(result))
  pwr <- result[["LRT power determined by simulation"]]
  expect_gte(pwr, 0)
  expect_lte(pwr, 1)
})
