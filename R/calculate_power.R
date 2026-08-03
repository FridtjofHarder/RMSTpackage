#' Determines test power when sample size is given
#'
#' Calculates and simulates the test power given a sample size. Supports tests on difference and ratio in restricted mean survival time (RMST), and log rank test (LRT). Supports superiority and non-inferiority tests.
#'
#' Survival curves need to be defined by \dfn{scale} and \dfn{shape} parameter, in the standard parameterisation
#' defined by \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}. If breakpoints in time are provided, survival and loss can be defined over piecewise exponential, or piecewise Weibull functions. Power
#' can be determined for log rank test, RMST difference, and RMST ratio.
#'
#' @param scale_ctrl Required. Specifies the \dfn{scale parameter} in the control group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_trmt Required. Specifies the \dfn{scale parameter} in the treatment group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_loss Required. Specifies the \dfn{scale parameter} for loss to follow-up. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull). No loss to follow-up is assumed if undefined.
#' @param shape_ctrl Specifies the \dfn{shape parameter} in the control group. Defaults to \code{shape_ctrl = 1}, simplifying to exponential survival. If \code{length(shape_ctrl) = 1} and \code{length(scale_ctrl) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_trmt Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt = 1}, simplifying to exponential survival. If \code{length(shape_trmt) = 1} and \code{length(scale_trmt) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_loss Specifies the \dfn{shape parameter} for loss to follow-up. Defaults to \code{shape_loss = 1}, simplifying to exponential loss. If \code{length(shape_loss) = 1} and \code{length(scale_loss) > 1}, the same shape parameter will be assumed for each section of the loss distribution.
#' @param breakpoints_ctrl Vector of breakpoints of the piecewise Weibull distribution in the control group. Must have length of \code{scale_ctrl} \eqn{-1} and \code{shape_ctrl} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_trmt Vector of breakpoints of the piecewise Weibull distribution in the treatment group. Must have length of \code{scale_trmt} \eqn{-1} and \code{shape_trmt} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_loss Vector of breakpoints of the piecewise Weibull distribution for loss to follow-up. Must have length of \code{scale_loss} \eqn{-1} and \code{shape_loss} \eqn{-1}. First element must be \code{> 0}.
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param sides Sidedness of inference test, either \code{1} or \code{2}. \code{sides = 1} assumes alternative hypothesis of: \itemize{
#' \item \eqn{H_1\text{: } \text{RMST}_\text{difference} = \text{RMST}_\text{trmt} - \text{RMST}_\text{ctrl} > 0},
#' \item \eqn{H_1\text{: } \text{RMST}_\text{ratio} = \text{RMST}_\text{trmt} / \text{RMST}_\text{ctrl} > 1}, or
#' \item \eqn{H_1\text{: } \text{HR} = h(t)_\text{trmt} / h(t)_\text{ctrl} < 1}.}
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
#' @param n Integer specifying sample size for calculating power.
#' @param plot_example_data Logical. Specifies whether to create a plot with example data. Plots with total sample size of \eqn{n = 100} if \code{n} is undefined.
#' @param plot_design_curves Logical. Specifies whether to plot survival curves.
#' @param parameterisation One of: \itemize{
#' \item \code{parameterisation = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterisation = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterisation = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#'
#' @return Returns a list with total sample sizes for each test and a test power.
#'
#' @export
#'
calculate_power <- function(
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
    M = 1000,
    n = NULL,
    plot_example_data = TRUE,
    plot_design_curves = TRUE,
    parameterisation = 1){
  int_fun_n_or_power(
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
    RMSTR_closed_form = RMSTR_closed_form, # RMSTR = RMST_trmt / RMST_ctrl = RMST_arm1 / RMST_arm0
    LRT_closed_form = LRT_closed_form, # HR = h(trmt) / h(ctrl = h_arm1 / h_arm0)
    satterthwaite_corr = satterthwaite_corr,
    censor_beyond_tau = censor_beyond_tau,
    plot_example_data = plot_example_data,
    plot_design_curves = plot_design_curves,
    parameterisation = parameterisation
  )

}
