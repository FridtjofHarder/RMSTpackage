#' Simulates a trial arm
#'
#' Simulates survival of one arm, including administrative censoring and loss to follow-up. Survival function needs to be fully defined
#' as a Weibull or exponential distribution.
#'
#' Details
#'
#' @param scale A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time Specifies the accrual time during which participants are recruited.
#' @param follow_up_time Specifies the follow-up time after accrual.
#' @param loss_scale A scalar \eqn{>0} specifying the \dfn{scale parameter} of loss to follow-up. No loss to follow-up is assumed if undefined.
#' @param loss_shape A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential loss.
#' @param sample_size Sample size.
#' @param label Group label.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau}.
#' @param censor_beyond_tau Boolean. All observations past tau \eqn{\tau} censored if \code{c(TRUE)}.
#'
#' @return data frame containing observations times, status (event = 1, censored
#'  = 1), and group label
#'
#' @export
#'
#' @examples
#' simulate_data(
#'   scale = 1,
#'   accrual_time = 1,
#'   follow_up_time = 1,
#'   loss_scale = 1,
#'   sample_size = 100,
#'   label = 0,
#' )
simulate_data <- function(
    scale,
    shape = 1,
    parameterization = 1,
    accrual_time = 0,
    follow_up_time = Inf,
    loss_scale = NULL, # loss is assumed to follow Weibull
    loss_shape = 1,
    sample_size,
    label = 0,
    tau = NULL,
    censor_beyond_tau = FALSE) {
  # convert to standard parameterization if needed
  if (parameterization != 1){
    scale <- reparameterize(parameterization = parameterization,
                            scale = scale,
                            shape = shape)
    loss_scale <- reparameterize(parameterization = parameterization,
                            scale = loss_scale,
                            shape = loss_shape)
  }

  total_time <- accrual_time + follow_up_time

  # draw event times vom weibull distribution
  observations <- stats::rweibull(n = sample_size, shape = shape, scale = scale)
  status <- rep(1, sample_size)

  # censor observations if loss to follow up is defined
  if (!is.null(loss_scale)) {
    loss_to_follow_up <- stats::rweibull(
      n = sample_size,
      shape = loss_shape,
      scale = loss_scale
    )
    status[loss_to_follow_up < observations] <- 0
    observations <- pmin(observations, loss_to_follow_up)
  }

  # censor observations if total time is not Inf
  if (follow_up_time != Inf) {
    admin_loss <- total_time - stats::runif(n = sample_size, max = accrual_time)
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
