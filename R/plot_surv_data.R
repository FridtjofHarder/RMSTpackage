#' Plots example survival data
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param xlim Defaults to \code{c(0, 1.5*tau)}.
#' @param ylim Range of y-axis as survival percentages. Defaults to \code{c(0, 100)}.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau}. The time horizon will be marked in the figure by a vertical line.
#' @param n An integer \eqn{>0} specifying the total sample size. Will be increased to the next even number if uneven. Group sample sizes are assumed to be equal.
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period.
#' @param censor_beyond_tau Boolean. All observations past tau \eqn{\tau} censored if \code{c(TRUE)}.
#' @param loss_scale Scale of Weibull distributed loss to follow-up.
#' @param loss_shape Shape of Weibull distributed loss to follow-up.
#' @param plot_reverse_KM Boolean. Will plot a reverse KM curve if  \code{c(TRUE)}.
#'
#' @export
#'
#' @examples
#' plot_surv_data(
#' scale_trmt = 1.4,
#' scale_ctrl = 1,
#' accrual_time = 1,
#' follow_up_time = 10,
#' tau = 1,
#' loss_scale = 1,
#' n = 100)
plot_surv_data <- function(scale_trmt, scale_ctrl, shape_trmt = 1,
                           shape_ctrl = 1, parameterization = 1, tau = NULL,
                           xlim = NULL, ylim = c(0, 100), n = NULL,
                           accrual_time = 0, follow_up_time = NULL,
                           censor_beyond_tau = FALSE, loss_scale = NULL,
                           loss_shape = 1, plot_reverse_KM = FALSE){
  data_frame_ctrl <- simulate_data(scale = scale_ctrl, shape = shape_ctrl,
                                   accrual_time = accrual_time,
                                   follow_up_time = follow_up_time,
                                   loss_scale = loss_scale,
                                   loss_shape = loss_shape,
                                   sample_size = round(n / 2),
                                   label = 0)
  simulated_data <- rbind(data_frame_ctrl, simulate_data(scale = scale_trmt, shape = shape_trmt,
                                                         accrual_time = accrual_time,
                                                         follow_up_time = follow_up_time,
                                                         loss_scale = loss_scale,
                                                         loss_shape = loss_shape,
                                                         sample_size = round(n / 2),
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

  # reverse KM
  if (plot_reverse_KM){
    print("hello world")
  }

}
