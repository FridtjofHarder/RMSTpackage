#' Simulates a trial arm
#'
#' Simulates survival of one arm, including administrative censoring and loss to follow-up. Survival function needs to be specified
#' as Weibull or exponential.
#'
#' Details
#'
#' @param scale Specifies the \dfn{scale parameter} in the treatment group.
#' @param shape Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored if \code{TRUE}.
#' @param n Sample size.
#' @param loss_scale Specifies the \dfn{scale parameter} of loss to follow-up. No loss to follow-up is assumed if undefined.
#' @param loss_shape Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential loss.
#' @param label Group label.
#'
#' @return Data frame containing observations times, status (event = 1, censored
#'  = 0), and group label
#'
#' @export
#'
#' @examples
#' ctrl_df <- simulate_data(
#'   scale = 6,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   loss_scale = 10,
#'   n = 50,
#'   label = 0,
#' )
#' trmt_df <- simulate_data(
#'   scale = 10,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   loss_scale = 10,
#'   n = 50,
#'   label = 1,
#' )
#' surv_df <- rbind(ctrl_df, trmt_df) # build data frame for both control and treatment group
#' head(surv_df)
#' tail(surv_df)
#'
#'
simulate_data <- function(
    scale,
    shape = 1,
    parameterization = 1,
    accrual_time = 0,
    follow_up_time = Inf,
    tau = NULL,
    censor_beyond_tau = FALSE,
    n,
    loss_scale = NULL, # loss is assumed to follow Weibull
    loss_shape = 1,
    label = 0
  ) {
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

  # draw event times from weibull distribution
  observations <- stats::rweibull(n = n, shape = shape, scale = scale)
  status <- rep(1, n)

  # censor observations if loss to follow up is defined
  if (!is.null(loss_scale)) {
    loss_to_follow_up <- stats::rweibull(
      n = n,
      shape = loss_shape,
      scale = loss_scale
    )
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
