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
#' @param plot_data Boolean. Plots Kaplan Meier curve of simulated data.
#' @param plot_recruitment Boolean. Plots recruitment plot.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau}.
#' @param plot_curves Boolean. Plots design curves if \code{TRUE}.
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
#'   plot_curves = TRUE,
#'   plot_data = TRUE,
#'   plot_recruitment = TRUE
#' )
simulate_data <- function(scale,
                          shape = 1,
                          parameterization = 1,
                          accrual_time = 0,
                          follow_up_time = Inf,
                          loss_scale = NULL, # loss is assumed to follow Weibull
                          loss_shape = 1,
                          sample_size,
                          label = 0,
                          plot_curves = FALSE,
                          plot_data = FALSE,
                          plot_recruitment = FALSE,
                          tau = NULL,
                          censor_beyond_tau = FALSE) {
  # convert to standard parameterization if needed
  if (parameterization == 2) {
    scale <- 1 / (scale^(1 / shape))
  }
  if (parameterization == 3) {
    scale <- 1 / scale
  }
  total_time <- accrual_time + follow_up_time

  # draw event times vom weibull distribution
  observations <- stats::rweibull(n = sample_size, shape = shape, scale = scale)
  status <- rep(1, sample_size)

  # censor observations if loss to follow up is defined
  if (!is.null(loss_scale)) {
    loss_to_follow_up <- stats::rweibull(n = sample_size, shape = loss_shape, scale = loss_scale)
    status[loss_to_follow_up < observations] <- 0
    observations <- pmin(observations, loss_to_follow_up)
  }

  # censor observations if total time is not Inf
  if (follow_up_time != Inf) {
    admin_loss <- total_time - stats::runif(n = sample_size, max = accrual_time)
    status[admin_loss < observations] <- 0
    observations <- pmin(observations, admin_loss)
  }

  # construct data frame with observation time, status, and group label
  data_df <- data.frame(observations, status, label)

  # plot design curves if requested
  if (plot_curves) {
    x <- NULL
    graphics::curve(stats::pweibull(x, scale = scale, shape = shape, lower.tail = FALSE),
      col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1)
    )
    graphics::legend("bottomleft",
      legend = paste0(
        "design curve with \n", "scale =",
        round(scale, 2), " and shape =",
        round(shape, 2)
      ),
      col = c("green"), lty = 1, y.intersp = 1.5, bty = "n", cex = 0.8
    )
  }

  # plot KM curve if requested
  if (plot_data) {
    surv_obj <- survival::Surv(time = observations, event = status)
    plot(survival::survfit(surv_obj ~ 1), mark.time = T, conf.int = F, xlab = "t", ylab = "S(t)")
  }

  # produce recruitment plot if requested
  if (plot_recruitment) {
    if (100 < sample_size) {
      sample_size <- round(seq(from = 1, to = sample_size, length.out = 100))
    }
    recruitment <- total_time - admin_loss # recruitment time in study time
    stop <- recruitment + observations # last observation in study time
    df <- data.frame(recruitment, stop)
    df_sorted <- df[order(recruitment), ]
    plot(
      x = df_sorted$recruitment, y = 1:sample_size,
      xlim = c(0, total_time), ylim = c(0, sample_size),
      xlab = "study time", ylab = "participants ordered by entry into study"
    )
    graphics::title("Individual observation trails in study time")
    graphics::segments(
      x0 = df_sorted$recruitment, y0 = 1:length(recruitment),
      x1 = df_sorted$stop
    )
  }

  return(data_df)
}

#### just for testing
