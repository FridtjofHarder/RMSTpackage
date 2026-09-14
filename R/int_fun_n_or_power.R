#' Internal function calculating n given power, or power given n
#'
#' @noRd
int_fun_n_or_power <- function(
    scale_ctrl,
    scale_trmt,
    scale_loss = NULL,
    shape_ctrl = 1,
    shape_trmt = 1,
    shape_loss = 1,
    breakpoints_ctrl = NULL,
    breakpoints_trmt = NULL,
    breakpoints_loss = NULL,
    accrual_time = 0,
    follow_up_time = Inf,
    tau = NULL,
    sides = 1,
    power = NULL,
    one_sided_alpha = 0.025,
    margin_RMSTD = 0,
    margin_RMSTR = 1,
    margin_LRT = 1,
    RMSTD_closed_form = TRUE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = TRUE,
    satterthwaite_corr = FALSE,
    RMSTD_simulation = FALSE, # RMSTD = RMST_trmt - RMST_ctrl = RMST_arm1 - RMST_arm0
    RMSTR_simulation = FALSE, # RMSTR = RMST_trmt / RMST_ctrl = RMST_arm1 / RMST_arm0
    LRT_simulation = FALSE,   # HR = h(trmt) / h(ctrl = h_arm1 / h_arm0)
    censor_beyond_tau = FALSE,
    M = 1,
    n = NULL,
    plot_example_data = FALSE,
    plot_design_curves = FALSE,
    parameterisation = 1
) {
  # basic definitions -----------------------------------------------------------
  if (length(shape_ctrl) == 1 & length(scale_ctrl) > 1) shape_ctrl <- rep(1, length(scale_ctrl))
  if (length(shape_trmt) == 1 & length(scale_trmt) > 1) shape_trmt <- rep(1, length(scale_trmt))
  if (length(shape_loss) == 1 & length(scale_loss) > 1) shape_loss <- rep(1, length(scale_loss))
  check_inputs(scale_ctrl = scale_ctrl,
               scale_trmt = scale_trmt,
               scale_loss = scale_loss,
               shape_ctrl = shape_ctrl,
               shape_trmt = shape_trmt,
               shape_loss = shape_loss,
               breakpoints_ctrl = breakpoints_ctrl,
               breakpoints_trmt = breakpoints_trmt,
               breakpoints_loss = breakpoints_loss,
               accrual_time = accrual_time,
               follow_up_time = follow_up_time,
               tau = tau,
               sides = sides,
               power = power,
               one_sided_alpha = one_sided_alpha,
               RMSTD_closed_form = RMSTD_closed_form,
               RMSTR_closed_form = RMSTR_closed_form,
               parameterisation = parameterisation
  )
  ss_RMSTD_closed_form <- ss_RMSTR_closed_form <- ss_LRT_closed_form <-
  pwr_RMSTD_closed_form <- pwr_RMSTR_closed_form <- pwr_LRT_closed_form <-
  pwr_RMSTD_simulated <- pwr_RMSTR_simulated <- pwr_LRT_simulated <-
  RMST_ctrl <- RMST_trmt <- True_RMSTD <- True_RMSTR <-
  ss_RMSTD_closed_form_sat <- ss_RMSTR_closed_form_sat <-
  pwr_RMSTD_closed_form_sat <- pwr_RMSTR_closed_form_sat <- NA

  # capture scales to use later for plot function
  original_scales <- list(scale_ctrl = scale_ctrl,
                          scale_trmt = scale_trmt,
                          scale_loss = scale_loss)
  # capture breakpoints to use later for plot function
  original_breakpoints <- list(breakpoints_ctrl = breakpoints_ctrl,
                               breakpoints_trmt = breakpoints_trmt,
                               breakpoints_loss = breakpoints_loss)

  breakpoints_ctrl <- normalize_breakpoints(breakpoints_ctrl)
  breakpoints_trmt <- normalize_breakpoints(breakpoints_trmt)
  breakpoints_loss <- normalize_breakpoints(breakpoints_loss)



  if (parameterisation != 1) {
    scale_ctrl <- reparameterise(parameterisation, scale_ctrl, shape_ctrl)
    scale_trmt <- reparameterise(parameterisation, scale_trmt, shape_trmt)
    scale_loss <- reparameterise(parameterisation, scale_loss, shape_loss)
  }

  if (RMSTD_closed_form || RMSTR_closed_form) { # get RMSTD and RMSTR
    RMST_ctrl <- get_theoretical_rmst(scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl, tau = tau)
    RMST_trmt <- get_theoretical_rmst(scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt, tau = tau)
    True_RMSTD <- RMST_trmt - RMST_ctrl
    True_RMSTR <- RMST_trmt / RMST_ctrl
  }
  if (margin_RMSTD != 0) {
    stopifnot(
      "Noninferiority margin of RMST difference must be below assumed RMST difference." =
        margin_RMSTD < True_RMSTD
    )
  }
  if (margin_RMSTR != 1) {
    stopifnot(
      "Noninferiority margin of RMST ratio must be below assumed RMST ratio." =
        margin_RMSTR < True_RMSTR
    )
  }
  if (margin_LRT != 1) {
    stopifnot(
      "Hazard ratio is not constant since shape parameters differ between groups.
              Noninferiority margin for hazard ratio is appropriate only when hazard ratio is constant" =
        shape_trmt == shape_ctrl
    )
    h_ctrl <- get_h(x = 1, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl)
    h_trmt <- get_h(x = 1, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt)
    true_HR <- h_trmt / h_ctrl
    stopifnot(
      "Noninferiority margin of hazard ratio must be above assumed hazard ratio." =
        margin_LRT > true_HR
    )
  }
  # closed form ----------------------------------------------------------------
  shared_args <- list(
    scale_ctrl = scale_ctrl, scale_trmt = scale_trmt, scale_loss = scale_loss,
    shape_ctrl = shape_ctrl, shape_trmt = shape_trmt, shape_loss = shape_loss,
    breakpoints_ctrl = breakpoints_ctrl, breakpoints_trmt = breakpoints_trmt,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time, follow_up_time = follow_up_time, tau = tau,
    sides = sides, power = power, alpha = one_sided_alpha
  )
  if (RMSTD_closed_form) {
    rmstd_args <- c(shared_args, list(margin = margin_RMSTD, RMST_ctrl = RMST_ctrl, RMST_trmt = RMST_trmt, contrast = "difference"))
    if(is.null(n)){ # calculate n if unspecified
      ss_RMSTD_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstd_args)
      if(satterthwaite_corr){ # use n from RMST closed form for df calculation
        ss_RMSTD_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(satterthwaite_n = ss_RMSTD_closed_form)))
      }
    }
    if(is.null(power)){ # calculate power if unspecified
      pwr_RMSTD_closed_form <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(n = n)))
      if(satterthwaite_corr){ # use specified n for df calculation
        pwr_RMSTD_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(n = n, satterthwaite_n = n)))
      }
    }
  }
  if (RMSTR_closed_form) {
    rmstr_args <- c(shared_args, list(margin = margin_RMSTR, RMST_ctrl = RMST_ctrl, RMST_trmt = RMST_trmt, contrast = "ratio"))
    if(is.null(n)){ # calculate n if unspecified
      ss_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstr_args)
      if(satterthwaite_corr){ # use n from RMST closed form for df calculation
        ss_RMSTR_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstr_args, list(satterthwaite_n = ss_RMSTR_closed_form)))
      }
    }
    if(is.null(power)){ # calculate power if unspecified
      pwr_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, c(rmstr_args, list(n = n)))
      if(satterthwaite_corr){ # use specified n for df calculation
        pwr_RMSTR_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstr_args, list(n = n, satterthwaite_n = n)))
      }
    }
  }
  if (LRT_closed_form){
    LRT_args <- c(shared_args, list(censor_beyond_tau = censor_beyond_tau, margin_LRT = margin_LRT))
    if(is.null(n)){ # calculate n if unspecified
      ss_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, LRT_args)
    }
    if(is.null(power)){ # calculate power if unspecified
      pwr_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, c(LRT_args, list(n = n)))
    }
  }
  # simulations  ---------------------------------------------------------------
  if (RMSTD_simulation || RMSTR_simulation || LRT_simulation) {
    n_per_arm <- round(n / 2)
    sim_shared <- list(
      scale_loss = scale_loss, shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time, follow_up_time = follow_up_time,
      n = n_per_arm
    )
    # bundle all per-iteration inputs so they can be sent to workers cleanly
    worker_args <- list(
      scale_trmt = scale_trmt, shape_trmt = shape_trmt, breakpoints_trmt = breakpoints_trmt,
      scale_ctrl = scale_ctrl, shape_ctrl = shape_ctrl, breakpoints_ctrl = breakpoints_ctrl,
      sim_shared = sim_shared,
      tau = tau, one_sided_alpha = one_sided_alpha,
      margin_RMSTD = margin_RMSTD, margin_RMSTR = margin_RMSTR, margin_LRT = margin_LRT,
      censor_beyond_tau = censor_beyond_tau,
      RMSTD_simulation = RMSTD_simulation,
      RMSTR_simulation = RMSTR_simulation,
      LRT_simulation = LRT_simulation
    )
    one_sim <- function(i, args) {
      simulated_data <- rbind(
        do.call(simulate_data, c(list(scale = args$scale_trmt, shape = args$shape_trmt,
                                      breakpoints = args$breakpoints_trmt, label = 1), args$sim_shared)),
        do.call(simulate_data, c(list(scale = args$scale_ctrl, shape = args$shape_ctrl,
                                      breakpoints = args$breakpoints_ctrl, label = 0), args$sim_shared))
      )
      result_i <- list(RMSTD = 0, RMSTR = 0, LRT = 0, tau_changed = FALSE)
      if (args$RMSTD_simulation || args$RMSTR_simulation) {
        tau_temp <- args$tau
        min_max <- min(
          max(simulated_data$observations[simulated_data$label == 0]),
          max(simulated_data$observations[simulated_data$label == 1])
        )
        if (min_max < args$tau) {
          tau_temp <- min_max
          result_i$tau_changed <- TRUE
        }
        result <- survRM2::rmst2(
          simulated_data$observations,
          simulated_data$status,
          simulated_data$label,
          tau = tau_temp,
          alpha = args$one_sided_alpha * 2
        )$unadjusted.result
        if (args$RMSTD_simulation) result_i$RMSTD <- as.numeric(result[1, 2] > args$margin_RMSTD)
        if (args$RMSTR_simulation) result_i$RMSTR <- as.numeric(result[2, 2] > args$margin_RMSTR)
      }
      if (args$LRT_simulation) {
        if (args$censor_beyond_tau)
          simulated_data$status[simulated_data$observations > args$tau] <- 0
        fit <- survival::coxph(survival::Surv(observations, status) ~ label, data = simulated_data)
        result_i$LRT <- as.numeric(summary(fit)$conf.int[, "upper .95"] < args$margin_LRT)
      }
      return(result_i)
    }
    sim_results <- lapply(seq_len(M), one_sim, args = worker_args)
    if (RMSTD_simulation) pwr_RMSTD_simulated <- mean(sapply(sim_results, `[[`, "RMSTD"))
    if (RMSTR_simulation) pwr_RMSTR_simulated <- mean(sapply(sim_results, `[[`, "RMSTR"))
    if (LRT_simulation)   pwr_LRT_simulated   <- mean(sapply(sim_results, `[[`, "LRT"))
    if (any(sapply(sim_results, `[[`, "tau_changed")))
      warning("tau was reduced to the minimum largest observation across groups in at least one iteration.")
  }
  # plot example data if requested ---------------------------------------------

  if (plot_design_curves | plot_example_data) {
    plot_surv(
    scale_ctrl = original_scales$scale_ctrl, # since scale parameters may have been altered by reparameterise() above
    scale_trmt = original_scales$scale_trmt,
    scale_loss = original_scales$scale_loss,
    shape_ctrl = shape_ctrl,
    shape_trmt = shape_trmt,
    shape_loss = shape_loss,
    breakpoints_ctrl = original_breakpoints$breakpoints_ctrl,
    breakpoints_trmt = original_breakpoints$breakpoints_trmt,
    breakpoints_loss = original_breakpoints$breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    censor_beyond_tau = censor_beyond_tau,
    n = n,
    parameterisation = parameterisation,
    plot_hazards = TRUE,
    plot_HR = TRUE,
    plot_reverse_KM = TRUE,
    plot_log_log = TRUE,
    plot_proportions = TRUE
    )
  }
  # returns -----------------------------------------------------------------
  result <- list(
    "Sample size for RMST difference determined by closed-form solution" = ss_RMSTD_closed_form,
    "Satterthwaite-corrected sample size for RMST difference" = ss_RMSTD_closed_form_sat,
    "Sample size for RMST ratio determined by closed-form solution" = ss_RMSTR_closed_form,
    "Satterthwaite-corrected sample size for RMST ratio" = ss_RMSTR_closed_form_sat,
    "Sample size for LRT determined by closed-form solution" = ss_LRT_closed_form,
    "Power for RMST difference determined by closed-form solution" = pwr_RMSTD_closed_form,
    "Satterthwaite-corrected power for RMST difference" = pwr_RMSTD_closed_form_sat,
    "Power for RMST ratio determined by closed-form solution" = pwr_RMSTR_closed_form,
    "Satterthwaite-corrected power for RMST ratio" = pwr_RMSTR_closed_form_sat,
    "Power for LRT determined by closed-form solution" = pwr_LRT_closed_form,
    "RMSTD power determined by simulation" = pwr_RMSTD_simulated,
    "RMSTR power determined by simulation" = pwr_RMSTR_simulated,
    "LRT power determined by simulation" = pwr_LRT_simulated,
    "RMST treatment group" = RMST_trmt,
    "RMST control group" = RMST_ctrl,
    "RMST difference" = True_RMSTD,
    "RMST ratio" = True_RMSTR
  )
  result <- Filter(Negate(is.na), result)
  return(result)
}
