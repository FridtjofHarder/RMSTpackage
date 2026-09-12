#' Determines sample size when test power is given
#'
#' Calculates the sample size given a desired test power. Supports tests on difference and ratio in restricted mean survival time (RMST), and log rank test (LRT).
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
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
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
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored for simulations and log rank test if \code{TRUE}.
#' @param plot_example_data Logical. Specifies whether to create a plot with example data. Plots with total sample size of \eqn{n = 100} if \code{n} is undefined.
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
#' @return Returns a list with total sample sizes for each test and a test power.
#'
#' @export
#'
#' @examples
#'
#' # Sample size for superiority test with Satterthwaite approximation
#'   args_sup <- list(
#'   scale_ctrl = 0.17,
#'   scale_trmt = 0.1,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   tau = 4,
#'   scale_loss = 0.1,
#'   satterthwaite_corr = TRUE
#'   )
#'   result_sup <- do.call(calculate_sample_size, args = args_sup)
#'   print(result_sup)
#'
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
#' # Assume heavy loss to follow-up
#'   args_sup_loss <- args_sup
#'   args_sup_loss$scale_loss <- 0.2
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
  power = 0.8,
  one_sided_alpha = 0.025,
  margin_RMSTD = 0,
  margin_RMSTR = 1,
  margin_LRT = 1,
  RMSTD_closed_form = TRUE,
  RMSTR_closed_form = FALSE,
  LRT_closed_form = TRUE,
  satterthwaite_corr = FALSE,
  censor_beyond_tau = FALSE,
  plot_example_data = FALSE,
  plot_design_curves = FALSE,
  parameterisation = 1
) {
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
    power = power,
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

