#' Compare sample size calculations
#'
#' Two approaches have been proposed for sample size calculation in the context of survival study design when choosing the difference or ratio between restricted mean survival times (RMSTs) as a contrast.
#' One approach is based on simulating trials according to prespecified parameters, the second approach is based on a closed form approximation as proposed by Yung and Liu (2020).
#' This function allows for comparing both approaches. The functions from the package \pkg{npsurvSS} have been extended to allow for sample size calculations when designing not only superiority, but also noninferioity trials.
#' Simulations for sample size calculations are based in the package \pkg{SSRMST} which natively allows for sample size calculation in superiority and noninferiority study design (Horiguchi and Uno, 2017).
#'
#' @param scale_trt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trt} \eqn{=1}, simplifying to exponential survival.
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
#' result <- compare_sample_size(scale_trt = 1.4, shape_trt = 1, scale_ctrl = 1, shape_ctrl = 1, parameterization = 1, accrual_time = 1, follow_up_time = 10, tau = 1, sides = 1, alpha = 0.025, power = 0.8, margin = 0.1, simulation = TRUE, M = 100)
#'
compare_sample_size <- function(scale_trt, shape_trt = 1, scale_ctrl, shape_ctrl = 1, parameterization = 1,
                                accrual_time = 0,
                                follow_up_time = NULL,
                                tau = NULL, sides = 1,
                                alpha = 0.025, power = 0.8, margin = 0,
                                simulation = TRUE, M = 1000){

  #browser()
  if(is.null(tau)) stop("please specify horizon")
  if(is.null(follow_up_time)) {
    warning("follow_up_time not specified, has been set to horizon")
    follow_up_time <- taun
  }
  total_time <- accrual_time + follow_up_time
  # initialize df for storing sample sizes

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
  }

  RMST_trt <- integrate(pweibull, shape = shape_trt, scale = scale_trt, lower = 0,
                          upper = tau, lower.tail = FALSE)$value
  RMST_ctrl <- integrate(pweibull, shape = shape_ctrl, scale = scale_ctrl, lower = 0,
                           upper = tau, lower.tail = FALSE)$value
  Delta <- RMST_trt-RMST_ctrl
  arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time, surv_shape = shape_trt, surv_scale = 1/scale_trt, loss_scale = 0.00001, follow_time = follow_up_time)
  arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time, surv_shape = shape_ctrl,surv_scale = 1/scale_ctrl,loss_scale = 0.00001, follow_time = follow_up_time)
  design  <- npsurvSS:::calc_design(arm_trt, arm_ctrl, test = list(test = "rmst difference", milestone = tau))
  total_sample_size_npsurvSS <- sum((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(design$tsigma2) * stats::qnorm(power) )^2 /
      (design$delta + margin)^2 *c(0.5, 0.5))

  #browser()
  if(simulation){
    power_SSRMST <- SSRMST::ssrmst(ac_number = total_sample_size_npsurvSS, ac_period = accrual_time, tot_time = accrual_time + follow_up_time, tau = tau, scale0 = scale_ctrl, shape0 = shape_ctrl,
                           scale1 = scale_trt, shape1 = shape_trt, margin = margin, ntest = M)
  }
  result <- list("n calculated via closed form" = total_sample_size_npsurvSS, "simulated test power" = power_SSRMST)
  return(result)

}


