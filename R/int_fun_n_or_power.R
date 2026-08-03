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
    RMSTR_closed_form = TRUE,
    LRT_closed_form = TRUE,
    satterthwaite_corr = FALSE,
    RMSTD_simulation = FALSE, # RMSTD = RMST_trmt - RMST_ctrl = RMST_arm1 - RMST_arm0
    RMSTR_simulation = FALSE, # RMSTR = RMST_trmt / RMST_ctrl = RMST_arm1 / RMST_arm0
    LRT_simulation = FALSE,   # HR = h(trmt) / h(ctrl = h_arm1 / h_arm0)
    censor_beyond_tau = FALSE,
    M = 1,
    n = NULL,
    plot_example_data = TRUE,
    plot_design_curves = TRUE,
    parameterisation = 1
) {
  # basic definitions -----------------------------------------------------------
  if (length(shape_ctrl) == 1 & length(scale_ctrl) > 1) shape_ctrl <- rep(shape_ctrl, length(scale_ctrl))
  if (length(shape_trmt) == 1 & length(scale_trmt) > 1) shape_trmt <- rep(shape_trmt, length(scale_trmt))
  if (length(shape_loss) == 1 & length(scale_loss) > 1) shape_loss <- rep(shape_loss, length(scale_loss))
  check_inputs(scale_ctrl = scale_ctrl,
               scale_trmt = scale_trmt,
               scale_loss = scale_loss,
               shape_ctrl = shape_ctrl,
               shape_trmt = shape_trmt,
               shape_loss = shape_loss,
               breakpoints_ctrl = breakpoints_ctrl,
               breakpoints_trmt = breakpoints_trmt,
               breakpoints_loss = breakpoints_loss,
               follow_up_time = follow_up_time,
               tau = tau,
               sides = sides,
               power = power,
               one_sided_alpha = one_sided_alpha,
               RMSTD_closed_form = RMSTD_closed_form,
               RMSTR_closed_form = RMSTR_closed_form,
               n,
               parameterisation = parameterisation
  )
  ss_RMSTD_closed_form <- ss_RMSTR_closed_form <- ss_LRT_closed_form <-
  pwr_RMSTD_closed_form <- pwr_RMSTR_closed_form <- pwr_LRT_closed_form <-
  pwr_RMSTD_simulated <- pwr_RMSTR_simulated <- pwr_LRT_simulated <-
  RMST_ctrl <- RMST_trmt <- True_RMSTD <- True_RMSTR <-
  ss_RMSTD_closed_form_sat <- ss_RMSTR_closed_form_sat <-
  pwr_RMSTD_closed_form_sat <- pwr_RMSTR_closed_form_sat <- NA

  breakpoints_ctrl <- normalize_breakpoints(breakpoints_ctrl)
  breakpoints_trmt <- normalize_breakpoints(breakpoints_trmt)
  breakpoints_loss <- normalize_breakpoints(breakpoints_loss)

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
  # reparameterise  --------------------------------------------------------------
  if (parameterisation != 1) {
    scale_ctrl <- reparameterize(parameterisation, scale_ctrl, shape_ctrl)
    scale_trmt <- reparameterize(parameterisation, scale_trmt, shape_trmt)
    scale_loss <- reparameterize(parameterisation, scale_loss, shape_loss)
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
        ss_RMSTD_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(satterthwaite_n = n)))
      }
    }
    if(is.null(power)){ # calculate power if unspecified
      pwr_RMSTD_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstd_args)
      if(satterthwaite_corr){ # use specified n for df calculation
        pwr_RMSTD_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(satterthwaite_n = n)))
      }
    }
  }
  if (RMSTR_closed_form) {
    rmstr_args <- c(shared_args, list(margin = margin_RMSTR, RMST_ctrl = RMST_ctrl, RMST_trmt = RMST_trmt, contrast = "ratio"))
    if(is.null(n)){ # calculate n if unspecified
      ss_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstd_args)
      if(satterthwaite_corr){ # use n from RMST closed form for df calculation
        ss_RMSTR_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(satterthwaite_n = n)))
      }
    }
    if(is.null(power)){ # calculate power if unspecified
      pwr_RMSTR_closed_form <- do.call(get_ss_pwr_cf_RMST, rmstd_args)
      if(satterthwaite_corr){ # use specified n for df calculation
        pwr_RMSTR_closed_form_sat <- do.call(get_ss_pwr_cf_RMST, c(rmstd_args, list(satterthwaite_n = n)))
      }
    }
  }
  if (LRT_closed_form){
    LRT_args <- c(shared_args, list(censor_beyond_tau = censor_beyond_tau, margin = margin_LRT))
    if(is.null(n)){ # calculate n if unspecified
      ss_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, LRT_args)
    }
    if(is.null(power)){ # calculate n if unspecified
      pwr_LRT_closed_form <- do.call(get_ss_pwr_cf_LRT, LRT_args)
    }
  }
  # simulations  ---------------------------------------------------------------
  if (RMSTD_simulation || RMSTR_simulation || LRT_simulation) {
    tau_changed <- FALSE
    if (RMSTD_simulation) RMSTD_simul_results <- rep(0, M)
    if (RMSTR_simulation) RMSTR_simul_results <- rep(0, M)
    if (LRT_simulation)   LRT_simul_results   <- rep(0, M)
    n_per_arm <- round(n / 2)
    sim_shared <- list(
      scale_loss = scale_loss, shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time, follow_up_time = follow_up_time,
      n = n_per_arm
    )
    for (i in 1:M) {
      simulated_data <- rbind(
        simulate_data(scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt, label = 1, !!!sim_shared),
        simulate_data(scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl, label = 0, !!!sim_shared)
      )
      if (RMSTD_simulation || RMSTR_simulation) {
        tau_temp <- tau
        min_max <- min(
          max(simulated_data$observations[simulated_data$label == 0]),
          max(simulated_data$observations[simulated_data$label == 1])
        )
        if (min_max < tau) {
          tau_temp <- min_max
          tau_changed <- TRUE
        }
        result <- survRM2::rmst2(
          simulated_data$observations,
          simulated_data$status,
          simulated_data$label,
          tau = tau_temp,
          alpha = one_sided_alpha * 2
        )$unadjusted.result
        if (RMSTD_simulation) RMSTD_simul_results[i] <- as.numeric(result[1, 2] > margin_RMSTD)
        if (RMSTR_simulation) RMSTR_simul_results[i] <- as.numeric(result[2, 2] > margin_RMSTR)
      }
      if (LRT_simulation) {
        if (censor_beyond_tau)
          simulated_data$status[simulated_data$observations > tau] <- 0
        fit <- survival::coxph(survival::Surv(observations, status) ~ label, data = simulated_data)
        LRT_simul_results[i] <- as.numeric(summary(fit)$conf.int[, "upper .95"] < margin_LRT)
      }
    }
    if (RMSTD_simulation) pwr_RMSTD_simulated <- mean(RMSTD_simul_results)
    if (RMSTR_simulation) pwr_RMSTR_simulated <- mean(RMSTR_simul_results)
    if (LRT_simulation)   pwr_LRT_simulated   <- mean(LRT_simul_results)
    if (tau_changed)
      warning("tau was reduced to the minimum largest observation across groups in at least one iteration.")
  }
  # plot example data if requested ---------------------------------------------

  if (plot_design_curves) {
    x <- NULL
    graphics::curve(
      ppweibull::ppweibull(
        x,
        rate = 1 / scale_trmt^shape_trmt,
        alpha = shape_trmt,
        t = breakpoints_trmt,
        lower.tail = FALSE
      ),
      col = "darkblue",
      xlab = "t",
      ylab = "S(t)",
      ylim = c(0, 1),
      xlim = c(0, 1.5 * tau),
      lwd = 2,
      main = "Design survival curves",
      yaxt = "n"
    )
    graphics::axis(
      2,
      at = seq(1, 0, by = -0.2),
      labels = paste0(seq(100, 0, by = -20), "%"),
      las = 1
    )
    graphics::curve(
      ppweibull::ppweibull(
        x,
        rate = 1 / scale_ctrl^shape_ctrl,
        alpha = shape_ctrl,
        t = breakpoints_ctrl,
        lower.tail = FALSE
      ),
      col = "red",
      lwd = 2,
      add = TRUE
    )
    graphics::abline(v = tau, col = "black", lwd = 2)
    graphics::text(
      x = tau,
      y = 0.1,
      pos = 4,
      labels = bquote("Time horizon " * tau * " = " * .(tau)),
      cex = .8
    )
    graphics::legend(
      "bottomleft",
      legend = c(
        paste0(
          "Treatment group with \n",
          "scale = ",
          paste(round(scale_trmt, 2), collapse = ", "),
          " and shape = ",
          paste(round(shape_trmt, 2), collapse = ", ")
        ),
        paste0(
          "Control group with \n",
          "scale = ",
          paste(round(scale_ctrl, 2), collapse = ", "),
          " and shape = ",
          paste(round(shape_ctrl, 2), collapse = ", ")
        )
      ),
      col = c("darkblue", "red"),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = .8
    )
  }
  if (plot_example_data) {
    if (is.na(n)) n <- 200
    plot_surv(
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
      censor_beyond_tau = censor_beyond_tau,
      n = round(n / 2)
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
    "Satterthwaite-corrected pwoer for RMST ratio" = pwr_RMSTR_closed_form_sat,
    "Power for LRT determined by closed-form solution" = pwr_LRT_closed_form,
    "RMSTD power determined by simulation" = pwr_RMSTD_simulated,
    "RMSTR power determined by simulation" = pwr_RMSTR_simulated,
    "LRT power determined by simulation" = pwr_LRT_simulated,
    "RMST treatment group" = RMST_trmt,
    "RMST control group" = RMST_ctrl,
    "RMST difference" = True_RMSTD,
    "RMST ratio" = True_RMSTR
  )
  return(result)
}
