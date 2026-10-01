#' Determines tau when test power and sample size are given
#'
#' Determines the smallest time horizon \eqn{\tau} given a desired test power and assumed sample size, using the closed form solution. Supports tests on difference and ratio in restricted mean survival time (RMST), and log-rank test (LRT).
#' Supports superiority and non-inferiority tests.
#'
#' Survival can be defined by either a single Weibull function, or a piecewise exponential function.
#' Weibull survival curves need to be defined by \dfn{scale} and, optionally, \dfn{shape} parameter, in the standard parameterisation
#' given by \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}. If breakpoints in time are provided, survival and loss can be defined by a piecewise exponential function,
#' where the scale parameters correspond to piecewise hazard rates.
#'
#' @param scale_ctrl Required. Specifies the \dfn{scale parameter} in the control group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_trmt Required. Specifies the \dfn{scale parameter} in the treatment group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_loss Required. Specifies the \dfn{scale parameter} for loss to follow-up. Can be a scalar (Weibull or exponential loss), or a vector (piecewise Weibull).
#' No loss to follow-up is assumed if undefined.
#' @param shape_ctrl Specifies the \dfn{shape parameter} in the control group. Defaults to \code{shape_ctrl = 1}, simplifying to exponential survival.
#' If \code{length(shape_ctrl) = 1} and \code{length(scale_ctrl) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_trmt Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt = 1}, simplifying to exponential survival.
#' If \code{length(shape_trmt) = 1} and \code{length(scale_trmt) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_loss Specifies the \dfn{shape parameter} for loss to follow-up. Defaults to \code{shape_loss = 1}, simplifying to exponential loss.
#' If \code{length(shape_loss) = 1} and \code{length(scale_loss) > 1}, the same shape parameter will be assumed for each section of the loss distribution.
#' @param breakpoints_ctrl Vector of breakpoints of the piecewise Weibull distribution in the control group.
#' Must have length of \code{scale_ctrl} \eqn{-1} and \code{shape_ctrl} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_trmt Vector of breakpoints of the piecewise Weibull distribution in the treatment group.
#' Must have length of \code{scale_trmt} \eqn{-1} and \code{shape_trmt} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_loss Vector of breakpoints of the piecewise Weibull distribution for loss to follow-up.
#' Must have length of \code{scale_loss} \eqn{-1} and \code{shape_loss} \eqn{-1}. First element must be \code{> 0}.
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param sides Sidedness of inference test, either \code{1} or \code{2}. \code{sides = 1} assumes alternative hypothesis of: \itemize{
#' \item \eqn{H_1\text{: } \text{RMST}_\text{difference} = \text{RMST}_\text{trmt} - \text{RMST}_\text{ctrl} > 0},
#' \item \eqn{H_1\text{: } \text{RMST}_\text{ratio} = \text{RMST}_\text{trmt} / \text{RMST}_\text{ctrl} > 1}, or
#' \item \eqn{H_1\text{: } \text{HR} = h(t)_\text{trmt} / h(t)_\text{ctrl} < 1}.}
#' @param power Test power with \code{power} \eqn{=1-\beta}.
#' @param one_sided_alpha \eqn{\alpha} level for one-sided inference test.
#' @param margin_RMSTD Non-inferiority margin for RMST difference. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{RMST}_\text{difference} > } \code{margin_RMSTD}, with  default \code{margin_RMSTD} \eqn{=0} simplifying to superiority test.
#' @param margin_RMSTR Non-inferiority margin for RMST ratio. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{RMST}_\text{ratio} > } \code{margin_RMSTR}, with  default \code{margin_RMSTR} \eqn{=1} simplifying to superiority test.
#' @param margin_LRT Non-inferiority margin for log-rank test in terms of hazard ratio \eqn{\text{HR}}. Assumes alternative hypothesis of \eqn{H_1\text{: } \text{HR} < } \code{margin_LRT}, with  default \code{margin_LRT} \eqn{=1} simplifying to superiority test.
#' @param RMSTD_closed_form Logical. Specifies whether to calculate sample size for RMST difference test.
#' @param RMSTR_closed_form Logical. Specifies whether to calculate sample size for RMST ratio test.
#' @param LRT_closed_form Logical. Specifies whether to calculate sample size for log-rank test.
#' @param satterthwaite_corr Logical. Adds sample size calculation based on t-distributed test statistic, with degrees of freedom found by the Satterthwaite approximation using the number of events in each group. Number of events is calculated based on the sample size determined based on standard normal distribution of test statistic.
#' @param n Integer specifying sample size for calculating power.
#' @param plot_design_curves Logical. Specifies whether to plot survival curves.
#' @param parameterisation Define only if Weibull function is specified, not for piecewise exponential survival. One of: \itemize{
#' \item \code{parameterisation = 1}: Default. Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}},
#' \item \code{parameterisation = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterisation = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}}. This is the parameterisation used for the base \R{} function \code{stats::pweibull()}.}
#'
#' For \code{shape = 1}, the Weibull function simplifies to exponential survival with
#' \eqn{\mathrm{scale} = \mathrm{hazard}} for \code{parameterisation = 1} or \code{2}, and
#' \eqn{\mathrm{scale} = 1 / \mathrm{hazard}} for \code{parameterisation = 3}.
#'
#' @param interval Interval over which the function searches for the smallest
#'   \eqn{\tau} to achieve the desired test power. Defaults to
#'   \eqn{[1,\ \mathtt{accrual\_time} + \mathtt{follow\_up\_time}]}{[1, accrual_time + follow_up_time]}.
#'
#' @return Returns a list with \eqn{\tau} required for each test to achieve the
#' desired test power with the specified sample size given the assumed design parameters,
#' and returns RMST in each group, RMST difference, and RMST ratio. Assumes all participants
#' to be censored beyond tau for log-rank test. Returns \code{NA}
#' when no \eqn{\tau} can be found between \eqn{0} and total trial time.
#'
#' @export
#'
#' @examples
#'
#' # tau for range of different tests
#'   args_sup <- list(
#'   scale_ctrl = 0.17,
#'   scale_trmt = 0.1,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   scale_loss = 0.1,
#'   satterthwaite_corr = TRUE,
#'   n= 405)
#'   result_tau <- do.call(calculate_tau, args = args_sup)
#'   print(result_tau)

calculate_tau <- function(
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
    sides = 1,
    power = 0.8,
    one_sided_alpha = 0.025,
    margin_RMSTD = 0,
    margin_RMSTR = 1,
    margin_LRT = 1,
    RMSTD_closed_form = TRUE,
    RMSTR_closed_form = FALSE,
    LRT_closed_form = TRUE,
    satterthwaite_corr = FALSE,
    n = NULL,
    plot_design_curves = FALSE,
    parameterisation = 1,
    interval = c(1, accrual_time + follow_up_time)
) {
  tau_RMSTD_closed_form = NA
  tau_RMSTD_closed_form_sat = NA
  tau_RMSTR_closed_form = NA
  tau_RMSTR_closed_form_sat = NA
  tau_LRT_closed_form  = NA
  common_args <- list(
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
    sides = sides,
    one_sided_alpha = one_sided_alpha,
    margin_RMSTD = margin_RMSTD,
    margin_RMSTR = margin_RMSTR,
    margin_LRT = margin_LRT,
    RMSTD_closed_form = RMSTD_closed_form,
    RMSTR_closed_form = RMSTR_closed_form,
    LRT_closed_form = LRT_closed_form,
    satterthwaite_corr = satterthwaite_corr,
    n = n,
    power = 0.8, # delete later!!!
    parameterisation = parameterisation
  )

  browser()

  tau_RMSTD <- uniroot(
    function(latest_tau) do.call(get_power_diff, c(list(tau = latest_tau, which_test = 1), common_args)),
    interval = c(1, 3), check.conv = TURE
  )$root
}

get_power_diff <- function(
    scale_ctrl,
    scale_trmt,
    scale_loss,
    shape_ctrl,
    shape_trmt,
    shape_loss,
    breakpoints_ctrl,
    breakpoints_trmt,
    breakpoints_loss,
    accrual_time,
    follow_up_time,
    tau,
    sides,
    power,
    one_sided_alpha,
    margin_RMSTD,
    margin_RMSTR,
    margin_LRT,
    RMSTD_closed_form,
    RMSTR_closed_form,
    LRT_closed_form,
    satterthwaite_corr,
    n,
    plot_design_curve,
    parameterisation,
    which_test
){
  power_res <- calculate_power(
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
    one_sided_alpha = one_sided_alpha,
    margin_RMSTD = margin_RMSTD,
    margin_RMSTR = margin_RMSTR,
    margin_LRT = margin_LRT,
    RMSTD_closed_form = RMSTD_closed_form, # RMSTD = RMST_trmt - RMST_ctrl = RMST_arm1 - RMST_arm0
    RMSTR_closed_form = TRUE, # RMSTR = RMST_trmt / RMST_ctrl = RMST_arm1 / RMST_arm0
    LRT_closed_form = LRT_closed_form, # HR = h(trmt) / h(ctrl = h_arm1 / h_arm0)
    satterthwaite_corr = satterthwaite_corr,
    censor_beyond_tau = TRUE,
    n = n,
    plot_design_curves = FALSE,
    parameterisation = parameterisation)
    power_diff <- as.data.frame(power_res) - power
    return(as.numeric(power_diff[which_test]))
}

