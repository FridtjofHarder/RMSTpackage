#' Calculate sample size
#'
#' Two approaches have been proposed for sample size calculation in the context of survival study design when choosing the difference or ratio between restricted mean survival times (RMSTs) as a contrast.
#' One approach is based on simulating trials according to prespecified parameters, the second approach is based on a closed form approximation as proposed by Yung and Liu (2020).
#' This function allows for comparing both approaches. The functions from the package \pkg{npsurvSS} have been extended to allow for sample size calculations when designing not only superiority, but also noninferioity trials.
#' Simulations for sample size calculations are based in the package \pkg{SSRMST} which natively allows for sample size calculation in superiority and noninferiority study design (Horiguchi and Uno, 2017).
#'
#' ToDo:  cure fraction
#'        admin censoring
#'        loss to FU
#'        plot example data as KM curve and as recruitment plot
#'        Cox regression univariate
#'        Replace SSRMST with own package
#'        allow for censoring of all at patient time for logRank (vergleich, ob dadurch das Resultat von RMST-Test sich verändert)
#'        zwei Fragen: 1. Ist RMST so gut wie Cox?
#'                    2. ist clsoed form so gut wie simul?
#'        Anmeldung!!!
#'        # Unterschied Zensur nach Zeithorizont, Zensur nach Studienzeit, oder Zensur zum Median
#'        RFrage: Unterschied Makuch-Simon vs LRank test sample size
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time length of accrual period
#' @param follow_up_time length of follow-up period
#' @param tau A scalar speifying the time horizon \eqn{\tau} at which to evaluate RMST with \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param sides sidedness of inference test
#' @param alpha level of \eqn{\alpha}-error
#' @param power test power with \code{power} \eqn{=1-\beta}
#' @param margin noninferiority margin. \code{margin} \eqn{>0} specifies a margin below the RMST of the control group.
#' @param simulation indicates whether to run a simulation based in the package \pkg{SSRMST}. Boolean.
#' @param M number of iterations when running simulation.
#' @param make_plot specifies whether to plot survival curves. Boolean.
#'
#' @return Returns a dataframe with a sample size and a test power. The sample size is calculated using the closed form function provided by the package
#' \pkg{npsurvSS}, the test power is obtained by running a simulation based on the package \pkg{SSRMST}. Both methods (closed form and simulation) come
#' to simiular results when the test power calculated is similar to the test power provided in the function call (default \code{power} \eqn{=0.8}).
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
                                  alpha = 0.025,
                                  power = 0.8,
                                  margin = 0,
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
                                  loss_scale, loss_shape = 1){

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
  browser()
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
                                      sample_size = round(100/2),
                                      label = 1) # arm 1 = trmt
      simulated_data <- rbind(simulated_data,
                               simulate_data(scale = scale_trmt,
                               shape = scale_ctrl,
                               accrual_time = accrual_time,
                               follow_up_time = follow_up_time,
                               loss_scale = loss_scale, # loss is assumed to follow Weibull
                               loss_shape = loss_shape,
                               sample_size = round(100/2),
                               label = 0)) # arm 0 = ctrl
      # determine whether testing for RMSTD turns out positive
      if(RMSTD_simulation){
      # handle large tau by limiting tau to minmax observation
        if(min(max(simulated_data$observations[simulated_data$arm == 0]),
               max(simulated_data$observations[simulated_data$arm == 1])) < tau){
          tau <- min(max(simulated_data$observations[simulated_data$arm == 0]),
                     max(simulated_data$observations[simulated_data$arm == 1]))
        }
        result <-  rmst2(data$time, data$status, data$arm, tau = tau,
                         alpha = one_sided_alpha * 2)$unadjusted.result
        lower <-  result[1, 2]
        RMSTD_simul_results[i] <- as.numeric(lower > -margin_rmst)
      }
      # determine whether testing for cox_ph turns out positive
      if(cox_ph_simulation){
        if(censor_beyond_tau){ # censor all observations beyond tau if requested
          simulated_data$status[simulated_data$observations > tau] <- 0
        }
        result <-  coxph(Surv(time, status) ~ arm, data = data)


      }





    }
  }

  # sample size RMSTD by simulation

  # sample size RMSTR by simulation

  # sample size Cox by simulation

  # plot example data if requested
  if(plot_example_data){
    data_frame_ctrl <- simulate_data(scale = scale_ctrl, shape = shape_ctrl,
                                     accrual_time = accrual_time,
                                     follow_up_time = follow_up_time,
                                     loss_scale = loss_scale,
                                     loss_shape = loss_shape,
                                     sample_size = simulation_sample_size,
                                     label = "ctrl")
    data_frame_trmt <- simulate_data(scale = scale_ctrl, shape = shape_ctrl,
                                     accrual_time = accrual_time,
                                     follow_up_time = follow_up_time,
                                     loss_scale = loss_scale,
                                     loss_shape = loss_shape,
                                     sample_size = simulation_sample_size,
                                     label = "trmt")
  }

  # plot design curves if requested
  if(plot_design_curves){
    curve(pweibull(x, scale = scale_trmt, shape = shape_trmt, lower.tail = FALSE),
          col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1), xlim = c(0, 1.5*tau))
    curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE),
          col = "red", add = TRUE)
    abline(v = tau, col = "blue")
    text(x = tau, y = 0.1, pos = 4, labels = paste("time horizon τ  =", tau))
    legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trmt, 2), " and shape =",
                                         round(shape_trmt, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape_ctrl, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)

  }
  result <- list("n calculated via closed form" = total_sample_size_npsurvSS, "simulated test power" = power_SSRMST)
  return(result)

  }


