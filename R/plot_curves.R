#' Draws survival curves and marks time horizon tau
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param xlim Defaults to \code{c(0, 1.5*tau)}.
#' @param ylim Range of y-axis as survival percentages. Defaults to \code{c(0, 100)}.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau}. The time horizon will be marked in the figure by a vertical line.
#'
#' @export
#'
#' @examples
#' plot_surv_curves(scale_trmt = 1, scale_ctrl = 0.9, tau = 1) # xlim set to c(0, 1.5*tau)
#' plot_surv_curves(scale_trmt = 1, scale_ctrl = 0.9, shape_ctrl = 1.1,
#' shape_trmt = 0.9) # xlim set automatically
#' plot_surv_curves(scale_trmt = 1, scale_ctrl = 0.9, shape_ctrl = 1.1,
#' shape_trmt = 0.9, xlim = c(0.5, 0.6), ylim = (0.5, 0.6) # xlim set automatically
plot_survival <- function(scale_trmt, scale_ctrl, shape_trmt = 1,
                             shape_ctrl = 1, parameterization = 1, tau = NULL,
                             xlim = NULL, ylim = c(0, 100), n = 100,
                          plot_curves = TRUE, plot_data = TRUE){

  x <- NULL
  if (!is.null(tau) && is.null(xlim)) {xlim <- c(0, 1.5*tau)}
graphics::curve(100*stats::pweibull(x, scale = scale_trmt, shape = shape_trmt,
                                lower.tail = FALSE),
                col = "darkblue", xlab = "t", ylab = "S(t)/%", ylim = ylim,
                xlim = xlim, main = "Design survival curves", yaxt = "n",
                lwd = 2)
axis(2, at = seq(ylim[1], ylim[2], round(diff(ylim)/10)),
     labels = paste(seq(ylim[1], ylim[2], round(diff(ylim)/10)), "%", sep = ""),
     las = 1)
graphics::curve(100*stats::pweibull(x, scale = scale_ctrl, shape = shape_ctrl,
                                lower.tail = FALSE),
                col = "red", add = TRUE, lwd = 2)
if (!is.null(tau)){ # mark tau if defined
graphics::abline(v = tau, col = "black", lwd = 2)
graphics::text(x = tau, y = 0.1, pos = 4,
               labels = bquote("Time horizon " * tau * " = " * .(tau)))
}
graphics::legend("bottomleft", legend=c(paste0("Treatment group with \n", "scale = ",
                                               round(scale_trmt, 2), " and shape = ",
                                               round(shape_trmt, 2)),
                                        paste0("Control group with \n", "scale = ",
                                               round(scale_ctrl, 2), " and shape = ",
                                               round(shape_ctrl, 2))),
                 col=c("darkblue", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 1)
}

#' Plots example survival data
#'
#' @param scale_trmt
#' @param scale_ctrl
#' @param shape_trmt
#' @param shape_ctrl
#' @param parameterization
#' @param tau
#' @param xlim
#' @param ylim
#' @param n
#' @param accrual_time
#'
#' @returns
#' @export
#'
#' @examples
plot_surv_data <- function(scale_trmt, scale_ctrl, shape_trmt = 1,
                           shape_ctrl = 1, parameterization = 1, tau = NULL,
                           xlim = NULL, ylim = c(0, 100), n = NULL,
                           accrual_time = 0, follow_up_time = NULL,
                           censor_beyond_tau = FALSE, loss_scale = NULL,
                           loss_shape = 1){
  data_frame_ctrl <- simulate_data(scale = scale_ctrl, shape = shape_ctrl,
                                   accrual_time = accrual_time,
                                   follow_up_time = follow_up_time,
                                   loss_scale = loss_scale,
                                   loss_shape = loss_shape,
                                   sample_size = round(simulation_sample_size/2),
                                   label = 0)
  simulated_data <- rbind(data_frame_ctrl, simulate_data(scale = scale_trmt, shape = shape_trmt,
                                                         accrual_time = accrual_time,
                                                         follow_up_time = follow_up_time,
                                                         loss_scale = loss_scale,
                                                         loss_shape = loss_shape,
                                                         sample_size = round(simulation_sample_size/2),
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

}
