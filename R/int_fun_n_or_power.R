#' Internal function calculating n given power, or power given n
#'
#' @noRd
int_fun_n_or_power <- function(scale_ctrl,
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
                               alpha = 0.025 * sides,
                               margin_RMSTD = 0,
                               margin_RMSTR = 1,
                               margin_LRT = 1,
                               RMSTD_closed_form = TRUE,
                               RMSTR_closed_form = FALSE,
                               LRT_closed_form = TRUE,
                               satter_corr = FALSE,
                               RMSTD_simulation = FALSE,
                               # RMSTD = RMST_trmt - RMST_ctrl = RMST_arm1 - RMST_arm0
                               RMSTR_simulation = FALSE,
                               # RMSTR = RMST_trmt / RMST_ctrl = RMST_arm1 / RMST_arm0
                               LRT_simulation = FALSE,
                               # HR = h(trmt) / h(ctrl = h_arm1 / h_arm0)
                               censor_beyond_tau = FALSE,
                               M = 1,
                               n = NULL,
                               plot_design_curves = FALSE,
                               parameterisation = 1,
                               output = "NULL") {
  # basic definitions -----------------------------------------------------------
  if (is.null(tau) &
      (
        RMSTD_closed_form |
        RMSTR_closed_form |
        RMSTD_simulation | RMSTR_simulation | censor_beyond_tau
      )) {
    stopifnot(
      "Trial length as sum of accrual_time plus follow_up_time must be defined if tau is undefined" =
        is.finite(accrual_time + follow_up_time)
    )
    tau <- accrual_time + follow_up_time
    warning(
      "tau was not defined and has been set to trial length as sum of accrual_time plus follow_up_time."
    )
  }
  if (length(shape_ctrl) == 1 &
      length(scale_ctrl) > 1)
    shape_ctrl <- rep(1, length(scale_ctrl))
  if (length(shape_trmt) == 1 &
      length(scale_trmt) > 1)
    shape_trmt <- rep(1, length(scale_trmt))
  if (length(shape_loss) == 1 &
      length(scale_loss) > 1)
    shape_loss <- rep(1, length(scale_loss))
  check_inputs(
    scale_ctrl = scale_ctrl,
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
    alpha = alpha,
    margin_RMSTD = margin_RMSTD,
    margin_RMSTR = margin_RMSTR,
    margin_LRT = margin_LRT,
    RMSTD_closed_form = RMSTD_closed_form,
    RMSTR_closed_form = RMSTR_closed_form,
    parameterisation = parameterisation
  )
  ss_RMSTD_closed_form <- ss_RMSTR_closed_form <- ss_LRT_closed_form <-
    pwr_RMSTD_closed_form <- pwr_RMSTR_closed_form <- pwr_LRT_closed_form <-
    pwr_RMSTD_simul <- pwr_RMSTD_simul_satter <- pwr_RMSTR_simul <- pwr_LRT_simul <-
    RMST_ctrl <- RMST_trmt <- True_RMSTD <- True_RMSTR <- pwr_RMSTR_simul_satter <-
    ss_RMSTD_closed_form_satter <- ss_RMSTR_closed_form_satter <-
    pwr_RMSTD_closed_form_satter <- pwr_RMSTR_closed_form_satter <- hazard_ratio <- pwr_RMSTD_closed_form_satter_df_from_n <-
    pwr_RMSTD_simul_satter_n <- NA

  # capture scales to use later for plot function
  original_scales <- list(scale_ctrl = scale_ctrl,
                          scale_trmt = scale_trmt,
                          scale_loss = scale_loss)
  # capture breakpoints to use later for plot function
  original_breakpoints <- list(
    breakpoints_ctrl = breakpoints_ctrl,
    breakpoints_trmt = breakpoints_trmt,
    breakpoints_loss = breakpoints_loss
  )

  breakpoints_ctrl <- normalize_breakpoints(breakpoints_ctrl)
  breakpoints_trmt <- normalize_breakpoints(breakpoints_trmt)
  breakpoints_loss <- normalize_breakpoints(breakpoints_loss)

  if (parameterisation != 1) {
    scale_ctrl <- reparameterise(parameterisation, scale_ctrl, shape_ctrl)
    scale_trmt <- reparameterise(parameterisation, scale_trmt, shape_trmt)
    scale_loss <- reparameterise(parameterisation, scale_loss, shape_loss)
  }

  if (RMSTD_closed_form ||
      RMSTR_closed_form) {
    # get RMSTD and RMSTR
    RMST_ctrl <- get_theoretical_rmst(
      scale = scale_ctrl,
      shape = shape_ctrl,
      breakpoints = breakpoints_ctrl,
      tau = tau
    )
    RMST_trmt <- get_theoretical_rmst(
      scale = scale_trmt,
      shape = shape_trmt,
      breakpoints = breakpoints_trmt,
      tau = tau
    )
    True_RMSTD <- RMST_trmt - RMST_ctrl
    True_RMSTR <- RMST_trmt / RMST_ctrl
  }

  if(shape_ctrl == shape_trmt & breakpoints_ctrl == breakpoints_trmt){
    hazard_ratio <- scale_trmt / scale_ctrl
    if(!length(unique(hazard_ratio)) == 1){
      hazard_ratio <- "Hazard ratio is not constant"
    }
  } else {
    hazard_ratio <- "Hazard ratio is not constant"
  }
  if ((RMSTD_closed_form | RMSTD_simulation) & margin_RMSTD != 0) {
    stopifnot(
      "Noninferiority margin of RMST difference must not be higher than assumed RMST difference." =
        margin_RMSTD <= True_RMSTD
    )
  }
  if ((RMSTR_closed_form | RMSTR_simulation) & margin_RMSTR != 1) {
    stopifnot(
      "Noninferiority margin of RMST ratio must not be higher than assumed RMST ratio." =
        margin_RMSTR <= True_RMSTR
    )
  }
  if (margin_LRT != 1) {
    if (shape_trmt != shape_ctrl) {
      warning(
        "Hazard ratio is not constant since shape parameters differ between groups.
              Noninferiority margin for hazard ratio is appropriate only when hazard ratio is constant"
      )
    }
    h_ctrl <- get_h(
      x = 1,
      scale = scale_ctrl,
      shape = shape_ctrl,
      breakpoints = breakpoints_ctrl
    )
    h_trmt <- get_h(
      x = 1,
      scale = scale_trmt,
      shape = shape_trmt,
      breakpoints = breakpoints_trmt
    )
    true_HR <- h_trmt / h_ctrl
    stopifnot("Noninferiority margin of hazard ratio must be above assumed hazard ratio." =
                margin_LRT > true_HR)
  }
  # closed form ----------------------------------------------------------------
  shared_args <- list(
    scale_ctrl = scale_ctrl,
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
    alpha = alpha,
    output = output
  )
  if (RMSTD_closed_form) {
    rmstd_args <- c(
      shared_args,
      list(
        margin = margin_RMSTD,
        RMST_ctrl = RMST_ctrl,
        RMST_trmt = RMST_trmt,
        contrast = "difference"
      )
    )
    if (output == "n") {
      # calculate n if unspecified
      ss_RMSTD_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstd_args)
      if (satter_corr) {
        # use n from RMST closed form for df calculation
        ss_RMSTD_closed_form_satter <- do.call(get_ss_pwr_cf_RMST, c(
          rmstd_args,
          list(
            n = ss_RMSTD_closed_form,
            satter_corr = TRUE,
            satter_df_from = "events"
          )
        ))
      }
    }
    if (output == "power") {
      # calculate power if unspecified
      pwr_RMSTD_closed_form <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(n = n)))
      if (satter_corr) {
        # calculate df from #events
        pwr_RMSTD_closed_form_satter <- do.call(get_ss_pwr_cf_RMST, c(
          rmstd_args,
          list(
            n = n,
            satter_corr = TRUE,
            satter_df_from = "events"
          )
        ))
        # use specified n for df calculation
        pwr_RMSTD_closed_form_satter_df_from_n <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(
          n = n,
          satter_corr = TRUE,
          satter_df_from = "n"
        )))
      }
    }
  }
  if (RMSTR_closed_form) {
    rmstr_args <- c(
      shared_args,
      list(
        margin = margin_RMSTR,
        RMST_ctrl = RMST_ctrl,
        RMST_trmt = RMST_trmt,
        contrast = "ratio"
      )
    )
    if (output == "n") {
      ss_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstr_args)
      if (satter_corr) {
        # use n from RMST closed form for df calculation
        ss_RMSTR_closed_form_satter <- do.call(get_ss_pwr_cf_RMST, c(
          rmstr_args,
          list(n = ss_RMSTR_closed_form, satter_corr = TRUE, satter_df_from = "events")
        ))
      }
    }
    if (output == "power") {
      pwr_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, c(rmstr_args, list(n = n)))
      if (satter_corr) {
        # use specified n for df calculation
        pwr_RMSTR_closed_form_satter <- do.call(get_ss_pwr_cf_RMST, c(rmstr_args, list(
          n = n, satter_corr = TRUE, satter_df_from = "events"
        )))
      }
    }
  }
  if (LRT_closed_form) {
    LRT_args <- c(
      shared_args,
      list(censor_beyond_tau = censor_beyond_tau, margin_LRT = margin_LRT)
    )
    if (output == "n") {
      # calculate n if unspecified
      ss_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, LRT_args)
    }
    if (output == "power") {
      # calculate power if unspecified
      pwr_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, c(LRT_args, list(n = n)))
    }
  }
  # simulations  ---------------------------------------------------------------
  if (RMSTD_simulation || RMSTR_simulation || LRT_simulation) {
    n_per_arm <- round(n / 2)
    sim_shared <- list(
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      n = n_per_arm
    )
    # bundle all per-iteration inputs so they can be sent to workers cleanly
    worker_args <- list(
      scale_trmt = scale_trmt,
      shape_trmt = shape_trmt,
      breakpoints_trmt = breakpoints_trmt,
      scale_ctrl = scale_ctrl,
      shape_ctrl = shape_ctrl,
      breakpoints_ctrl = breakpoints_ctrl,
      sim_shared = sim_shared,
      tau = tau,
      alpha = alpha,
      sides = sides,
      n = n,
      margin_RMSTD = margin_RMSTD,
      margin_RMSTR = margin_RMSTR,
      margin_LRT = margin_LRT,
      satter_corr = satter_corr,
      RMSTD_simulation = RMSTD_simulation,
      RMSTR_simulation = RMSTR_simulation,
      LRT_simulation = LRT_simulation,
      censor_beyond_tau = censor_beyond_tau
    )
      sim_results <- lapply(seq_len(M), one_sim, args = worker_args)
      if (RMSTD_simulation) {
        pwr_RMSTD_simul <- mean(sapply(sim_results, `[[`, "RMSTD"))
        if (satter_corr) {
          pwr_RMSTD_simul_satter <- mean(sapply(sim_results, `[[`, "RMSTD_satter"))
          pwr_RMSTD_simul_satter_n <- mean(sapply(sim_results, `[[`, "RMSTD_satter_n"))
        }
      }
      if (RMSTR_simulation) {
        pwr_RMSTR_simul <- mean(sapply(sim_results, `[[`, "RMSTR"))
        if (satter_corr) {
          pwr_RMSTR_simul_satter <- mean(sapply(sim_results, `[[`, "RMSTR_satter"))
        }
      }
      if (LRT_simulation)
        pwr_LRT_simul   <- mean(sapply(sim_results, `[[`, "LRT"))
      if (any(sapply(sim_results, `[[`, "tau_changed")))
        warning(
          "tau has been reduced to the minimum largest observation across groups in at least one iteration."
        )
    }
    # plot example data if requested ---------------------------------------------

    if (plot_design_curves) {
      plot_surv(
        scale_ctrl = original_scales$scale_ctrl,
        # since scale parameters may have been altered by reparameterise() above
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
      "Satterthwaite-corrected sample size for RMST difference" = ss_RMSTD_closed_form_satter,
      "Sample size for RMST ratio determined by closed-form solution" = ss_RMSTR_closed_form,
      "Satterthwaite-corrected sample size for RMST ratio" = ss_RMSTR_closed_form_satter,
      "Sample size for LRT determined by closed-form solution" = ss_LRT_closed_form,
      "Power for RMST difference determined by closed-form solution" = pwr_RMSTD_closed_form,
      "Satterthwaite-corrected power for RMST difference" = pwr_RMSTD_closed_form_satter,
      "Satterthwaite-corrected power for RMST difference with n for df" = pwr_RMSTD_closed_form_satter_df_from_n,
      "Power for RMST ratio determined by closed-form solution" = pwr_RMSTR_closed_form,
      "Satterthwaite-corrected power for RMST ratio" = pwr_RMSTR_closed_form_satter,
      "Power for LRT determined by closed-form solution" = pwr_LRT_closed_form,
      "RMSTD power determined by simulation" = pwr_RMSTD_simul,
      "Satterthwaite-corrected RMSTD power determined by simulation" = pwr_RMSTD_simul_satter,
      "Satterthwaite-corrected RMSTD power determined by simulation with n for df" = pwr_RMSTD_simul_satter_n,
      "RMSTR power determined by simulation" = pwr_RMSTR_simul,
      "Satterthwaite-corrected RMSTR power determined by simulation" = pwr_RMSTR_simul_satter,
      "LRT power determined by simulation" = pwr_LRT_simul,
      "RMST treatment group" = RMST_trmt,
      "RMST control group" = RMST_ctrl,
      "RMST difference" = True_RMSTD,
      "RMST ratio" = True_RMSTR,
      "Hazard ratio" = hazard_ratio
    )
    result <- Filter(Negate(is.na), result)
    return(result)
}

#' Internal helper functions
#'
#' used only for int_fun_n_or_power()
#'
#' @noRd
one_sim <- function(i, args) {
  simul_data <- rbind(do.call(simulate_data, c(
    list(
      scale = args$scale_trmt,
      shape = args$shape_trmt,
      breakpoints = args$breakpoints_trmt,
      label = 1
    ),
    args$sim_shared
  )), do.call(simulate_data, c(
    list(
      scale = args$scale_ctrl,
      shape = args$shape_ctrl,
      breakpoints = args$breakpoints_ctrl,
      label = 0
    ),
    args$sim_shared
  )))
  result_i <- list(
    RMSTD = 0,
    RMSTR = 0,
    RMSTD_satter = 0,
    RMSTD_satter_n = 0,
    RMSTR_satter = 0,
    LRT = 0,
    tau_changed = FALSE
  )
  if (args$RMSTD_simulation || args$RMSTR_simulation) {
    tau_temp <- args$tau
    min_max <- min(max(simul_data$observations[simul_data$label == 0]),
                   max(simul_data$observations[simul_data$label == 1]))
    if (min_max < args$tau) {
      tau_temp <- min_max
      result_i$tau_changed <- TRUE
    }
    fit <- survRM2::rmst2(
      simul_data$observations,
      simul_data$status,
      simul_data$label,
      tau = tau_temp,
      alpha = args$alpha * 2 / args$sides # keep alpha unchanged for sides = 2, double for sides = 1
    )
    result <- fit$unadjusted.result
    if (args$RMSTD_simulation) {
      if (args$sides == 1) {
        result_i$RMSTD <- as.numeric(result[1, 2] > args$margin_RMSTD) # result[1, 2] is lower limit of CI for RMSTD
        if (args$satter_corr) {
          se_diff <- sqrt(fit$RMST.arm0$rmst.var + fit$RMST.arm1$rmst.var)
          # df from #events
          satter_df_sim_events <- get_satter_df_simulation(
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "diff",
            satter_df_from = "events",
            n = args$n
          )
          print("simul df from events")
          print(satter_df_sim_events)
          half_CI_events <- se_diff * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim_events)
          lower_CI_satter <- result[1, 1] - half_CI_events # result[1, 1] is estimate of RMSTD
          result_i$RMSTD_satter <- as.numeric(lower_CI_satter > args$margin_RMSTD)
          # df from n
          satter_df_sim_n <- get_satter_df_simulation( # get df from n
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "diff",
            satter_df_from = "n",
            n = args$n
          )
          print("simul df from n")
          print(satter_df_sim_n)
          half_CI_n <- se_diff * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim_n)
          lower_CI_satter <- result[1, 1] - half_CI_n # result[1, 1] is estimate of RMSTD
          result_i$RMSTD_satter_n <- as.numeric(lower_CI_satter > args$margin_RMSTD)
        }
      }
      if (args$sides == 2) {
        result_i$RMSTD <- as.numeric(result[1, 4] < args$alpha) # result[1, 4] is p value for RMSTD
        if (args$satter_corr) {
          se_diff <- sqrt(fit$RMST.arm0$rmst.var + fit$RMST.arm1$rmst.var)
          # df from #events
          satter_df_sim_events <- get_satter_df_simulation(
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "diff",
            satter_df_from = "events",
            n = args$n)
          half_CI_events <- se_diff * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim_events)
          # get satterthwaite lower bound
          lower_CI_satter <- result[1, 1] - half_CI_events # result[1, 1] is estimate of RMSTD
          upper_CI_satter <- result[1, 1] + half_CI_events # result[1, 1] is estimate of RMSTD
          result_i$RMSTD_satter <- as.numeric(lower_CI_satter > 0 |
                                                upper_CI_satter < 0)
          # df from n
          satter_df_sim_n <- get_satter_df_simulation(
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "diff",
            satter_df_from = "n",
            n = args$n)
          half_CI_events <- se_diff * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim_n)
          lower_CI_satter <- result[1, 1] - half_CI_events # result[1, 1] is estimate of RMSTD
          upper_CI_satter <- result[1, 1] + half_CI_events # result[1, 1] is estimate of RMSTD
          result_i$RMSTD_satter_n <- as.numeric(lower_CI_satter > 0 |
                                                  upper_CI_satter < 0)

        }
      }
    }
    if (args$RMSTR_simulation) {
      if (args$sides == 1) {
        result_i$RMSTR <- as.numeric(result[2, 2] > args$margin_RMSTR) # result[2, 2] is lower limit of CI for RMSTR
        if (args$satter_corr) {
          se_log_ratio <- sqrt(
            fit$RMST.arm0$rmst.var / fit$RMST.arm0$rmst[["Est."]]^2 +
              fit$RMST.arm1$rmst.var / fit$RMST.arm1$rmst[["Est."]]^2
          )
          satter_df_sim <- get_satter_df_simulation(
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "ratio",
            n = args$n
          )
          half_CI_log <- se_log_ratio * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim)
          log_lower_CI_satter <- log(result[2, 1]) - half_CI_log # result[1, 1] is estimate of RMSTR
          result_i$RMSTR_satter <- as.numeric(log_lower_CI_satter > log(args$margin_RMSTR))
        }
      }
      if (args$sides == 2) {
        result_i$RMSTR <- as.numeric(result[2, 4] < args$alpha) # result[2, 4] is p value for RMSTR
        if (args$satter_corr) {
          se_log_ratio <- sqrt(
            fit$RMST.arm0$rmst.var / fit$RMST.arm0$rmst[["Est."]]^2 +
              fit$RMST.arm1$rmst.var / fit$RMST.arm1$rmst[["Est."]]^2
          )
          satter_df_sim <- get_satter_df_simulation(
            tau = args$tau,
            fit_object = fit,
            simul_data = simul_data,
            contrast = "ratio",
            n = args$n
          )
          half_CI_log <- se_log_ratio * get_crit_satter(alpha = args$alpha, sides = args$sides, satter_df = satter_df_sim)
          log_lower_CI_satter <- log(result[2, 1]) - half_CI_log # result[1, 1] is estimate of RMSTR
          log_upper_CI_satter <- log(result[2, 1]) + half_CI_log # result[1, 1] is estimate of RMSTR
          result_i$RMSTR_satter <- as.numeric(log_lower_CI_satter > 0 |
                                                log_upper_CI_satter < 0)
        }
      }
    }
  }
  if (args$LRT_simulation) {
    if (args$censor_beyond_tau) {
      simul_data$status[simul_data$observations > args$tau] <- 0
    }
    fit <- survival::coxph(survival::Surv(observations, status) ~ label, data = simul_data)
    log_CI_LRT <- stats::confint(fit, level = 1 - args$alpha * 2 / args$sides)
    if (args$sides == 1) {
      result_i$LRT <- log_CI_LRT[2] < log(args$margin_LRT)
    }
    if (args$sides == 2) {
      result_i$LRT <- (log_CI_LRT[2] < 0) | (log_CI_LRT[1] > 0)
    }
  }
  return(result_i)
}

get_crit_satter <- function(alpha, sides, satter_df){
  return(stats::qt(1 - alpha / sides, df = satter_df))
}

