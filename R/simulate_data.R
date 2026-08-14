#' Simulates a trial arm
#'
#' Simulates survival of one arm, including administrative censoring and loss to follow-up. Survival function needs to be specified
#' as Weibull or exponential.
#'
#' @param scale Specifies the \dfn{scale parameter}. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_loss Specifies the \dfn{scale parameter} for loss to follow-up. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull). No loss to follow-up is assumed if undefined.
#' @param shape Specifies the \dfn{shape parameter}. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param shape_loss Specifies the \dfn{shape parameter} for loss to follow-up. Defaults to \code{shape_loss = 1}, simplifying to exponential loss. If \code{length(shape_loss) = 1} and \code{length(scale_loss) > 1}, the same shape parameter will be assumed for each section of the loss distribution.
#' @param breakpoints Vector of breakpoints of the piecewise Weibull distribution. Must have length of \code{scale} \eqn{-1} and \code{shape} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_loss Vector of breakpoints of the piecewise Weibull distribution for loss to follow-up. Must have length of \code{scale_loss} \eqn{-1} and \code{shape_loss} \eqn{-1}. First element must be \code{> 0}.
#' @param parameterisation One of: \itemize{
#' \item \code{parameterisation = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterisation = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterisation = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored if \code{TRUE}.
#' @param n Sample size.
#' @param label Group label.
#'
#' @return Data frame containing observations times, status (event = 1, censored
#'  = 0), and group label.
#'
#' @export
#'
#' @examples
#' ctrl_df <- simulate_data(
#'   scale = 6,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   scale_loss = 10,
#'   n = 50,
#'   label = 0,
#' )
#' trmt_df <- simulate_data(
#'   scale = 10,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   scale_loss = 10,
#'   n = 50,
#'   label = 1,
#' )
#' surv_df <- rbind(ctrl_df, trmt_df) # build data frame for both control and treatment group
#' head(surv_df)
#' tail(surv_df)
#'
simulate_data <- function(
  scale,
  scale_loss = NULL, # loss is assumed to follow Weibull
  shape = 1,
  shape_loss = 1,
  breakpoints = NULL,
  breakpoints_loss = NULL,
  parameterisation = 1,
  accrual_time = 0,
  follow_up_time = Inf,
  tau = NULL,
  censor_beyond_tau = FALSE,
  n,
  label = 0
) {
  # convert to standard parameterisation if needed
  if (parameterisation != 1) {
    scale <- reparameterize(
      parameterisation = parameterisation,
      scale = scale,
      shape = shape
    )
    scale_loss <- reparameterize(
      parameterisation = parameterisation,
      scale = scale_loss,
      shape = shape_loss
    )
  }
  breakpoints <- normalize_breakpoints(breakpoints)
  breakpoints_loss <- normalize_breakpoints(breakpoints_loss)
  if (length(shape) == 1 && length(scale) > 1) {
    shape <- rep(shape, length(scale))
  }
  if (length(shape_loss) == 1 && length(scale_loss) > 1) {
    shape_loss <- rep(shape_loss, length(scale_loss))
  }

  total_time <- accrual_time + follow_up_time

  # draw event times; use fast piecewise exponential sampler when shape = 1 everywhere
  if (all(shape == 1)) {
    observations <- int_rpexp(n = n, rates = 1 / scale, breakpoints = breakpoints)
  } else {
    observations <- ppweibull::rpweibull(n = n, rate = 1 / scale^shape, alpha = shape, t = breakpoints)
  }
  status <- rep(1, n)

  # censor observations if loss to follow-up is defined
  if (!is.null(scale_loss)) {
    if (all(shape_loss == 1)) {
      loss_to_follow_up <- int_rpexp(n = n, rates = 1 / scale_loss, breakpoints = breakpoints_loss)
    } else {
      loss_to_follow_up <- ppweibull::rpweibull(n = n, rate = 1 / scale_loss^shape_loss, alpha = shape_loss, t = breakpoints_loss)
    }
    status[loss_to_follow_up < observations] <- 0
    observations <- pmin(observations, loss_to_follow_up)
  }

  # censor observations if total time is not Inf
  if (follow_up_time != Inf) {
    admin_loss <- total_time - stats::runif(n = n, max = accrual_time)
    status[admin_loss < observations] <- 0
    observations <- pmin(observations, admin_loss)
  }

  # censor observations beyond tau if requested
  if (!is.null(tau) && censor_beyond_tau) {
    status[tau < observations] <- 0
    observations <- pmin(observations, tau)
  }

  # construct data frame with observation time, status, and group label
  data_df <- data.frame(observations, status, label)

  return(data_df)
}
