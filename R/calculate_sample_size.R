#' Determines sample size or test power
#'
#' Calculates the sample size given a desired test power, or simulates the test power given a sample size. Supports tests on difference and ratio in restricted mean survival time (RMST), and log rank test (LRT). Supports superiority and non-inferiority tests.
#'
#' Sample size and power determination for both superiority and non-inferiority analysis are supported.
#' Survival curves need to be defined by \dfn{scale} and \dfn{shape} parameter, in the standard parameterization
#' defined by \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}. Sample size and power
#' can be determined for log rank test, RMST difference, and RMST ratio.
#'
#' @param scale_ctrl Specifies the \dfn{scale parameter} in the control group. Can be a scalar (weibull or exponential survival), or a vector (piecewise weibull).
#' @param scale_trmt Specifies the \dfn{scale parameter} in the treatment group. Can be a scalar (weibull or exponential survival), or a vector (piecewise weibull).
#' @param shape_ctrl Specifies the \dfn{shape parameter} in the control group. Defaults to \code{shape_ctrl = 1}, simplifying to exponential survival. If \code{length(shape_ctrl) = 1} and \code{length(shape_ctrl) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_trmt Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt = 1}, simplifying to exponential survival. If \code{length(shape_trmt) = 1} and \code{length(shape_trmt) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param breakpoints_ctrl Vector of breakpoints of the piecewise weibull distribution in the control group. Must have length of \code{scale_ctrl}and \code{shape_ctrl}. First element must be \code{0}.
#' @param breakpoints_trmt Vector of breakpoints of the piecewise weibull distribution in the treatment group. Must have length of \code{scale_trmt}and \code{shape_trmt}. First element must be \code{0}.
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param scale_loss Specifies the \dfn{scale parameter} of loss to follow-up. No loss to follow-up is assumed if undefined.
#' @param shape_loss Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential loss.
#' @param sides Sidedness of inference test, either \code{1} or \code{2}. \code{sides = 1} assumes alternative hypothesis of: \itemize{
#' \item \eqn{H_1\text{: } \text{RMST}_\text{difference} = \text{RMST}_\text{trmt} - \text{RMST}_\text{ctrl} > 0},
#' \item \eqn{H_1\text{: } \text{RMST}_\text{ratio} = \text{RMST}_\text{trmt} / \text{RMST}_\text{ctrl} > 1}, or
#' \item \eqn{H_1\text{: } \text{HR} = h(t)_\text{trmt} / h(t)_\text{ctrl} < 1}.}
#' @param power Test power with \code{power} \eqn{=1-\beta}.
#' @param one_sided_alpha \eqn{\alpha} level for one-sided inference test.
#' @param margin_RMSTD Non-inferiority margin for RMST difference. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{RMST}_\text{difference} > } \code{margin_RMSTD}, with  default \code{margin_RMSTD} \eqn{=0} simplifying to superiority test.
#' @param margin_RMSTR Non-inferiority margin for RMST ratio. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{RMST}_\text{ratio} > } \code{margin_RMSTR}, with  default \code{margin_RMSTR} \eqn{=1} simplifying to superiority test.
#' @param margin_LRT Non-inferiority margin for log rank test in terms of hazard ratio \eqn{\text{HR}}. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{HR} < } \code{margin_LRT}, with  default \code{margin_LRT} \eqn{=1} simplifying to superiority test.
#' @param RMSTD_closed_form Logical. Specifies whether to calculate sample size for RMST difference test.
#' @param RMSTR_closed_form Logical. Specifies whether to calculate sample size for RMST ratio test.
#' @param LRT_closed_form Logical. Specifies whether to calculate sample size for log rank test.
#' @param satterthwaite_corr Logical. Adds sample size calculation based on t-distributed test statistic, with degrees of freedom found by the Satterthwaite approximation using the number of events in each group. Number of events is calculated based on the sample size determined based on standard normal distribution of test statistic.
#' @param RMSTD_simulation Logical. Specifies whether to determine RMST difference test power via simulation.
#' @param RMSTR_simulation Logical. Specifies whether to determine RMST ratio test power via simulation.
#' @param LRT_simulation Logical. Specifies whether to determine log rank test power via simulation.
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored for simulations and log rank test if \code{TRUE}.
#' @param M Number of iterations when running simulation.
#' @param simulation_n Specifies sample size for simulations and for example plot.
#' @param plot_example_data Logical. Specifies whether to create a plot with example data. Plots with total sample size of \eqn{n = 100} if \code{simulation_n} is undefined.
#' @param plot_design_curves Logical. Specifies whether to plot survival curves.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#'
#' @return Returns a list with total sample sizes for each test and a test power.
#'
#' @export
#'
#' @examples
#'
#' # Sample size for superiority test with Satterthwaite approximation
#'   args_sup <- list(
#'   scale_ctrl = 6,
#'   scale_trmt = 10,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   tau = 4,
#'   scale_loss = 10,
#'   satterthwaite_corr = TRUE
#'   )
#'   result_sup <- do.call(calculate_sample_size, args = args_sup)
#'   print(result_sup)
#'
#' # Validate by running simulation with sample size previously obtained for RMST difference
#'   args_sup_sim <- args_sup
#'   args_sup_sim$RMSTD_simulation <- TRUE
#'   args_sup_sim$simulation_n <-
#'     result_sup$
#'     `Sample size for RMST difference determined by closed-form solution`
#'   result_sup_sim <- do.call(calculate_sample_size, args = args_sup_sim)
#'   print(result_sup_sim)
#'
#' # Sample size for non-inferiority test with margin
#'   args_noninf <- args_sup
#'   args_noninf$margin_LRT <- 1.3 # define noninferiority margin in terms of HR
#' # find RMST difference and RMST ratio margins equivalent to HR margin
#'   contrasts <- convert_contrast_ph(scale_ctrl = args_noninf$scale_ctrl,
#'                                    tau = args_noninf$tau,
#'                                    HR = args_noninf$margin_LRT)
#'   print(contrasts$RMSTD) # display RMST difference margin
#'   print(contrasts$RMSTR) # display RMST ratio margin
#'   args_noninf$margin_RMSTD <- contrasts$RMSTD
#'   args_noninf$margin_RMSTR <- contrasts$RMSTR
#'   result_noninf <- do.call(calculate_sample_size, args = args_noninf)
#'   print(result_noninf)
#'
#' # Assume heavy loss to follow up
#'   args_sup_loss <- args_sup
#'   args_sup_loss$scale_loss <- 2
#'   result_sup_loss <- do.call(calculate_sample_size, args = args_sup_loss)
#'   print(result_sup_loss)
#'
#' # Censure all observations past tau for LRT:
#' # eliminates information advantage of LRT over RMST based methods.
#'   args_sup_tau_cen <- args_sup
#'   args_sup_tau_cen$censor_beyond_tau <- TRUE
#'   result_sup_tau_cen <- do.call(calculate_sample_size,
#'                                 args = args_sup_tau_cen)
#'   print(result_sup_tau_cen)
#'
calculate_sample_size <- function(
  scale_ctrl,
  scale_trmt,
  shape_ctrl = 1,
  shape_trmt = 1,
  breakpoints_ctrl = 0,
  breakpoints_trmt = 0,
  accrual_time = 0,
  follow_up_time = NULL,
  tau = NULL,
  scale_loss = NULL,
  shape_loss = 1,
  sides = 1,
  power = 0.8,
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
  M = 5000,
  simulation_n = NA,
  plot_example_data = TRUE,
  plot_design_curves = TRUE,
  parameterization = 1
) {
  # basic definitions -----------------------------------------------------------
  # fill output elements in case they are not filled later on
  power_RMSTD_simulated <- power_RMSTR_simulated <- power_LRT_simulated <-
    ss_closed_form_RMSTD <- ss_closed_form_RMSTR <- ss_closed_form_LRT <-
    RMST_ctrl <- RMST_trmt <- True_RMSTD <- True_RMSTR <-
    ss_closed_form_RMSTD_sat <- ss_closed_form_RMSTR_sat <- NA

  if(length(shape_ctrl) == 1 & length(scale_ctrl) > 1){
    shape_ctrl <- rep(shape_ctrl, length(scale_ctrl))
  }
  if(length(shape_trmt) == 1 & length(scale_trmt) > 1){
    shape_trmt <- rep(shape_trmt, length(scale_trmt))
  }

  # error management --------------------------------------------------------
  stopifnot(
    # throw error when parameterization misspecified
    "Parameterization must be defined as either 1, 2, or 3." = parameterization ==
      1 ||
      parameterization == 2 ||
      parameterization == 3
  )
  # throw error when functions misspecified
  if (is.null(scale_ctrl) || is.null(scale_trmt)) {
    stop(
      "Please specify scale parameters for both treatment and survival group."
    )
  }
  # throw error when tau should have been specified
  if (!is.numeric(tau)) {
    stop("Please specify valid time horizon tau.")
  }
  if (is.null(follow_up_time)) {
    warning("Follow_up_time not specified, has been set to Inf.")
    follow_up_time <- Inf
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
        margin_RMSTD < True_RMSTR
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
      "Noninferiority margin of hazrad ratio be above assumed hazard ratio." =
        margin_LRT > true_HR
    )
  }
  if(length(scale_loss) > 1 | length(shape_loss) > 1){
    stop("Scale_loss and shape_loss must have length 1")
  }
  if(length(scale_ctrl) != length(shape_ctrl) | length(scale_trmt) != length(shape_trmt)){
    stop("Scale and shape parameter must have same length in each group")
  }
  # reparameterize if necessary --------------------------------------------------------------
  if (parameterization != 1) {
    scale_ctrl <- reparameterize(
      parameterization = parameterization,
      scale = scale_ctrl,
      shape = shape_ctrl)
    scale_trmt <- reparameterize(
      parameterization = parameterization,
      scale = scale_trmt,
      shape = shape_trmt)
    scale_loss <- reparameterize(
      parameterization = parameterization,
      scale = scale_loss,
      shape = shape_loss)
  }
  # closed form ----------------------------------------------------------------
  if (RMSTD_closed_form) {
    ss_closed_form_RMSTD <- get_ss_cf_RMSTD(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      sides = sides,
      power = power,
      alpha = one_sided_alpha,
      margin = margin_RMSTD,
      RMST_ctrl = RMST_ctrl,
      RMST_trmt = RMST_trmt
    )
    if (satterthwaite_corr) {
      ss_closed_form_RMSTD_sat <- get_ss_cf_RMSTD(
        scale_ctrl = scale_ctrl,
        scale_trmt = scale_trmt,
        shape_ctrl = shape_ctrl,
        shape_trmt = shape_trmt,
        breakpoints_ctrl = breakpoints_ctrl,
        breakpoints_trmt = breakpoints_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        tau = tau,
        scale_loss = scale_loss,
        shape_loss = shape_loss,
        sides = sides,
        power = power,
        alpha = one_sided_alpha,
        margin = margin_RMSTD,
        satterthwaite_n = ss_closed_form_RMSTD,
        RMST_ctrl = RMST_ctrl,
        RMST_trmt = RMST_trmt
      )
    }
  }
  if (RMSTR_closed_form) {
    ss_closed_form_RMSTR <- get_ss_cf_RMSTR(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      sides = sides,
      power = power,
      alpha = one_sided_alpha,
      margin = margin_RMSTR,
      RMST_ctrl = RMST_ctrl,
      RMST_trmt = RMST_trmt
    )
    if (satterthwaite_corr) {
      ss_closed_form_RMSTR_sat <- get_ss_cf_RMSTR(
        scale_ctrl = scale_ctrl,
        scale_trmt = scale_trmt,
        shape_ctrl = shape_ctrl,
        shape_trmt = shape_trmt,
        breakpoints_ctrl = breakpoints_ctrl,
        breakpoints_trmt = breakpoints_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        tau = tau,
        scale_loss = scale_loss,
        shape_loss = shape_loss,
        sides = sides,
        power = power,
        alpha = one_sided_alpha,
        margin = margin_RMSTR,
        satterthwaite_n = ss_closed_form_RMSTR,
        RMST_ctrl = RMST_ctrl,
        RMST_trmt = RMST_trmt
      )
    }
  }
  if (LRT_closed_form) {
    ss_closed_form_LRT <- get_ss_cf_LRT(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      censor_beyond_tau = censor_beyond_tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      sides = sides,
      power = power,
      alpha = one_sided_alpha,
      margin_LRT = margin_LRT
    )
  }
  # simulations  ---------------------------------------------------------------
  if (RMSTD_simulation || RMSTR_simulation || LRT_simulation) {
    tau_changed <- FALSE # indicator whether tau has been reduced in at least one iteration
    if (RMSTD_simulation) {
      RMSTD_simul_results <- rep(0, M)
    }
    if (RMSTR_simulation) {
      RMSTR_simul_results <- rep(0, M)
    }
    if (LRT_simulation) {
      LRT_simul_results <- rep(0, M)
    }
    for (i in 1:M) {
      simulated_data <- simulate_data(
        scale = scale_trmt,
        shape = shape_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss,
        n = round(
          simulation_n /
            2
        ),
        label = 1 # arm 1 = trmt
      )
      simulated_data <- rbind(
        simulated_data,
        simulate_data(
          scale = scale_ctrl,
          shape = shape_ctrl,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time,
          scale_loss = scale_loss,
          shape_loss = shape_loss,
          n = round(
            simulation_n /
              2
          ),
          label = 0 # arm 0 = ctrl
        )
      )
      if (RMSTD_simulation || RMSTR_simulation) {
        tau_temp <- tau # define tau_temp in case tau is too large
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
        lower_RMSTD <- result[1, 2]
        lower_RMSTR <- result[2, 2]
        if (RMSTD_simulation) {
          RMSTD_simul_results[i] <- as.numeric(lower_RMSTD > margin_RMSTD)
        }
        if (RMSTR_simulation) {
          RMSTR_simul_results[i] <- as.numeric(lower_RMSTR > margin_RMSTR)
        }
      }
      if (LRT_simulation) {
        if (censor_beyond_tau) {
          # censor all observations beyond tau if requested
          simulated_data$status[simulated_data$observations > tau] <- 0
        }
        fit <- survival::coxph(
          survival::Surv(observations, status) ~ label,
          data = simulated_data
        )
        LRT_simul_results[i] <-
          as.numeric(summary(fit)$conf.int[, "upper .95"] < margin_LRT)
      }
    }
    if (RMSTD_simulation) {
      power_RMSTD_simulated <- sum(RMSTD_simul_results) / M
    }
    if (RMSTR_simulation) {
      power_RMSTR_simulated <- sum(RMSTR_simul_results) / M
    }
    if (LRT_simulation) {
      power_LRT_simulated <- sum(LRT_simul_results) / M
    }
    if (tau_changed) {
      warning(
        "Warning: tau has been reduced to the minimum of the two
      largest observations in both groups in at least one simulated iteration."
      )
    }
  }
  # plot example data if requested ---------------------------------------------

  if (plot_design_curves) {
    x <- NULL
    graphics::curve(
      stats::pweibull(
        x,
        scale = scale_trmt,
        shape = shape_trmt,
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
      stats::pweibull(
        x,
        scale = scale_ctrl,
        shape = shape_ctrl,
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
          round(scale_trmt, 2),
          " and shape = ",
          round(shape_trmt, 2)
        ),
        paste0(
          "Control group with \n",
          "scale = ",
          round(scale_ctrl, 2),
          " and shape = ",
          round(shape_ctrl, 2)
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
    if (is.na(simulation_n)) simulation_n <- 200
    plot_surv(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      censor_beyond_tau = censor_beyond_tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      n = round(simulation_n / 2)
    )
  }
  # returns -----------------------------------------------------------------
  result <- list(
    "Sample size for RMST difference determined by closed-form solution" = ss_closed_form_RMSTD,
    "Satterthwaite-corrected sample size for RMST difference" = ss_closed_form_RMSTD_sat,
    "Sample size for RMST ratio determined by closed-form solution" = ss_closed_form_RMSTR,
    "Satterthwaite-corrected sample size for RMST ratio" = ss_closed_form_RMSTR_sat,
    "Sample size for LRT determined by closed-form solution" = ss_closed_form_LRT,
    "RMSTD power determined by simulation" = power_RMSTD_simulated,
    "RMSTR power determined by simulation" = power_RMSTR_simulated,
    "LRT power determined by simulation" = power_LRT_simulated,
    "RMST treatment group" = RMST_trmt,
    "RMST control group" = RMST_ctrl,
    "RMST difference" = True_RMSTD,
    "RMST ratio" = True_RMSTR
  )
  return(result)
}
