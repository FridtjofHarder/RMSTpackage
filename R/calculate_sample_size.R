#' Calculates sample size or test power
#'
#' Currently only determination of test power via simulation! Function for calculating the sample size given a desired test power, or the test power given a sample size.
#' Sample size and test power can be determined either analytically or by simulation.
#'
#'
#' Sample size and power determination for both superiory and non-inferiority analysis are supported.
#' Survival curves need to be defined by \dfn{scale} and \dfn{shape} parameter, in the standard parameterization
#' defined by \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}. Sample size and power
#' can be determined for either Cox regression or RMST-based methods (test for RMST difference and RMST ratio).
#' For comparing Cox regression with RMST based methods, all observations past time horizon \eqn{\tau}
#' may be censored if desired.
#'
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau} at which to evaluate RMST with \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param sides Sidedness of inference test, either \code{1} or \code{2}.
#' @param power Test power with \code{power} \eqn{=1-\beta}.
#' @param margin_cox Non-inferiority margin for Cox regression. \code{margin_cox} \eqn{=1} simplifies to superiority test.
#' @param margin_rmst Non-inferiority margin for RMST difference. \code{margin_rmst} \eqn{=0} simplifies to superiority test.
#' @param M Positive Integer. Number of iterations when running simulation.
#' @param plot_design_curves Boolean. Specifies whether to plot survival curves.
#' @param plot_example_data Boolean. Specifies whether to create a plot with example data .
#' @param one_sided_alpha MISSING
#' @param RMSTD_closed_form MISSING
#' @param RMSTR_closed_form MISSING
#' @param RMSTD_simulation MISSING
#' @param RMSTR_simulation MISSING
#' @param cox_ph_simulation MISSING
#' @param censor_beyond_tau MISSING
#' @param simulation_sample_size MISSING
#' @param loss_scale MISSING
#' @param loss_shape MISSING
#'
#' @return Returns a list with sample size and a test power.
#'
#' @export
#' @references
#' Yung, G., & Liu, Y. (2020). Sample size and power for the weighted log-rank test and Kaplan-Meier based tests with allowance for
#' nonproportional hazards. \emph{Biometrics}, \strong{76(3)}, 939-950.
#'
#' @examples
#' calculate_sample_size(
#' scale_trmt = 1.4,
#' scale_ctrl = 1,
#' accrual_time = 1,
#' follow_up_time = 10,
#' tau = 1,
#' RMSTD_simulation = TRUE,
#' RMSTR_simulation = TRUE,
#' cox_ph_simulation = TRUE,
#' loss_scale = 1,
#' M = 100,
#' simulation_sample_size = 100)
#'
calculate_sample_size <- function(scale_trmt,
                                  shape_trmt = 1,
                                  scale_ctrl,
                                  shape_ctrl = 1,
                                  parameterization = 1,
                                  accrual_time = 0,
                                  follow_up_time = NULL,
                                  tau = NULL,
                                  sides = 1,
                                  one_sided_alpha = 0.025,
                                  power = 0.8,
                                  margin_cox = 1,
                                  margin_rmst = 0,
                                  RMSTD_closed_form = FALSE,
                                  RMSTR_closed_form = FALSE,
                                  RMSTD_simulation = FALSE,
                                  RMSTR_simulation = FALSE,
                                  cox_ph_simulation = FALSE,
                                  censor_beyond_tau = FALSE,
                                  M = 100,
                                  simulation_sample_size,
                                  plot_design_curves = TRUE,
                                  plot_example_data = TRUE,
                                  loss_scale = NULL,
                                  loss_shape = 1){

  # error management --------------------------------------------------------

  stopifnot( # throw error when parameterization misspecified
    "parameterization must be defined as either 1, 2, or 3" = parameterization == 1 ||
      parameterization == 2 || parameterization == 3
  )

  # throw error when functions misspecified
  if(is.null(scale_trmt)||is.null(scale_ctrl)){
    stop("please specify scale parameters for both treatment and survival group")
  }

  # throw error when tau should have been specified
  if(!is.numeric(tau)){
    stop("please specify valid time horizon tau")
  }

  if(is.null(follow_up_time)) {
    warning("follow_up_time not specified, has been set to Inf")
    follow_up_time <- Inf
  }
  total_time <- accrual_time + follow_up_time

  # main function --------------------------------------------------------
  # browser()
  # convert to standard parameterization if needed
  if(parameterization == 2){
    scale_trmt <- 1/(scale_trmt^(1/shape_trmt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
  }

  if(parameterization == 3){
    scale_trmt <- 1/scale_trmt
    scale_ctrl <- 1/scale_ctrl
  }

  # simulate trial if simulations requested
  if(RMSTD_simulation || RMSTR_simulation || cox_ph_simulation){
    if(RMSTD_simulation) {RMSTD_simul_results <- rep(0, M)}
    if(RMSTR_simulation) {RMSTR_simul_results <- rep(0, M)}
    if(cox_ph_simulation){cox_ph_simul_results <- rep(0, M)}
    for (i in 1:M){
      simulated_data <- simulate_data(scale = scale_trmt,
                                      shape = shape_trmt,
                                      accrual_time = accrual_time,
                                      follow_up_time = follow_up_time,
                                      loss_scale = loss_scale, # loss is assumed to follow Weibull
                                      loss_shape = loss_shape,
                                      sample_size = round(simulation_sample_size/2),
                                      label = 1) # arm 1 = trmt
      simulated_data <- rbind(simulated_data,
                               simulate_data(scale = scale_ctrl,
                               shape = scale_ctrl,
                               accrual_time = accrual_time,
                               follow_up_time = follow_up_time,
                               loss_scale = loss_scale, # loss is assumed to follow Weibull
                               loss_shape = loss_shape,
                               sample_size = round(simulation_sample_size/2),
                               label = 0)) # arm 0 = ctrl
      # determine whether testing for RMSTD turns out positive
      if(RMSTD_simulation){
      # handle large tau by limiting tau to minmax observation
        if(min(max(simulated_data$observations[simulated_data$label == 0]),
               max(simulated_data$observations[simulated_data$label == 1])) < tau){
          tau <- min(max(simulated_data$observations[simulated_data$label == 0]),
                     max(simulated_data$observations[simulated_data$label == 1]))
        }
        result <-  survRM2::rmst2(simulated_data$observations, simulated_data$status, simulated_data$label, tau = tau,
                         alpha = one_sided_alpha * 2)$unadjusted.result
        lower <-  result[1, 2]
        RMSTD_simul_results[i] <- as.numeric(lower > -margin_rmst)
      }
      # determine whether testing for cox_ph turns out positive
      if(cox_ph_simulation){
        if(censor_beyond_tau){ # censor all observations beyond tau if requested
          simulated_data$status[simulated_data$observations > tau] <- 0
        }
        fit <- survival::coxph(survival::Surv(observations, status) ~ label, data = simulated_data)
        cox_ph_simul_results[i] <-
          as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
      }

    }
    power_RMSTD_simulated <- sum(RMSTD_simul_results)/M
    power_cox_ph_simulated <- sum(cox_ph_simul_results)/M
  }

  # sample size RMSTD by closed form

  # plot example data if requested
  if(plot_example_data){
    data_frame_ctrl <- simulate_data(scale = scale_ctrl, shape = shape_ctrl,
                                     accrual_time = accrual_time,
                                     follow_up_time = follow_up_time,
                                     loss_scale = loss_scale,
                                     loss_shape = loss_shape,
                                     sample_size = round(simulation_sample_size/2),
                                     label = 0)
    simulated_data <- rbind(data_frame_ctrl, simulate_data(scale = scale_trmt, shape = shape_trmt,
                                     accrual_time = accrual_time,
                                     follow_up_time = follow_up_time,
                                     loss_scale = loss_scale,
                                     loss_shape = loss_shape,
                                     sample_size = round(simulation_sample_size/2),
                                     label = 1))
    surv_obj <- survival::Surv(time = simulated_data$observations,
                               event = simulated_data$status)

    plot(survival::survfit(surv_obj~simulated_data$label), mark.time=T, conf.int = F, xlab = "t", ylab = "S(t)",
         col = c("red", "green"))
    graphics::legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trmt, 2), " and shape =",
                                         round(shape_trmt, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape_ctrl, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)

  }

  # plot design curves if requested
  if(plot_design_curves){
    x <- NULL
    graphics::curve(stats::pweibull(x, scale = scale_trmt, shape = shape_trmt, lower.tail = FALSE),
          col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1), xlim = c(0, 1.5*tau))
    graphics::curve(stats::pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE),
          col = "red", add = TRUE)
    graphics::abline(v = tau, col = "blue")
    graphics::text(x = tau, y = 0.1, pos = 4, labels = bquote("time horizon " * tau * " = " * .(tau)))
    graphics::legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trmt, 2), " and shape =",
                                         round(shape_trmt, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape_ctrl, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)

  }
  result <- list("RMSTD power determined by simulation" = power_RMSTD_simulated,
                 "Cox PH power determined by simulation" = power_cox_ph_simulated)
  return(result)

  }


