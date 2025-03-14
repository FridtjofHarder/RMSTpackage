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
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period.
#' @param tau A scalar \eqn{>0} specifying the time horizon \eqn{\tau}. The time horizon will be marked in the figure by a vertical line.
#' @param plot_data Boolean. Will plot example data if \code{c(TRUE)}.
#' @param n An integer \eqn{>0} specifying the total sample size when plotting example sata is requested.
#' @param plot_hazard Bloolean. Will plot hazard if \code{c(TRUE)}.
#' @param plot_RMSTD_over_tau Boolean. Will plot RMST\eqn{\tau}) if \code{c(TRUE)}.
#'
#' @export
#'
#' @examples
#' plot_curves(
#' scale_trmt = 1.4,
#' scale_ctrl = 1,
#' accrual_time = 1,
#' follow_up_time = 10,
#' tau = 1)
plot_curves <- function(scale_trmt, scale_ctrl,
                        shape_trmt = 1,
                        shape_ctrl = 1,
                        parameterization = 1,
                        xlim = NULL, ylim = c(0, 100),
                        accrual_time = 0,
                        follow_up_time = NULL,
                        tau = NULL,
                        plot_data = TRUE,
                        n = 100,
                        plot_hazard = TRUE,
                        plot_RMSTD_over_tau = FALSE){

  x <- NULL
  if (!is.null(tau) && is.null(xlim)) {xlim <- c(0, 1.5*tau)}
  graphics::curve(100*stats::pweibull(x, scale = scale_trmt, shape = shape_trmt,
                                      lower.tail = FALSE),
                  col = "darkblue", xlab = "t", ylab = "S(t) / %", ylim = ylim,
                  xlim = xlim, main = "Design survival curves", yaxt = "n",
                  lwd = 2)
  graphics::axis(2, at = seq(ylim[1], ylim[2], round(diff(ylim)/10)),
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

  # plot hazard curves if desired

  if(plot_hazard){

    graphics::curve(hazard(x, shape = shape_trmt, scale = scale_trmt),
                    col = "darkblue", xlab = "t", ylab = "h(t)",
                    xlim = xlim, main = "Hazard rates", lwd = 2, ylim = c(0, 2)) ## muss dynamisch angepasst werden!
    graphics::curve(hazard(x, shape = shape_ctrl, scale = scale_ctrl),
                    col = "red", lwd = 2, add = TRUE)
    graphics::legend("topleft", legend=c("Hazard in treatment group",
                                         "Hazard in control group"),
                     col=c("darkblue", "red"), lty=1:1,
                     y.intersp = 1.5, bty = "n", cex = 1)

  # plot RMST and RMSTD over tau
  }
  # browser()
  if (plot_RMSTD_over_tau) {
    # draw RMST(tau) for trmt
    graphics::curve(RMST_over_tau(x, shape = shape_trmt, scale = scale_trmt),
                    col = "darkblue", xlab = expression(tau), ylab = "RMST; RMSTD",
                    main=expression(RMST(tau) ~ "and" ~ RMSTD(tau)),
                    lwd = 2, ylim = c(-1, 2), xlim = c(0, 5))
    # draw RMST(tau) for ctrl
    graphics::curve(RMST_over_tau(x, shape = shape_ctrl, scale = scale_ctrl),
                    col = "red", lwd = 2, add = TRUE)
    # draw RMSTD(tau)
    graphics::curve(RMSTD_over_tau(x, shape_ctrl = shape_ctrl, scale_ctrl = scale_ctrl,
                                   shape_trmt = shape_trmt, scale_trmt = scale_trmt),
                    col = "black", lwd = 2, add = TRUE)
    # legend
    graphics::legend("topleft", legend=c(expression(RMST(tau) ~ "in treatment group"),
                                         expression(RMST(tau) ~ "in control group"),
                                         expression(RMSTD(tau))),
                     col=c("darkblue", "red", "black"), lty=1:1,
                     y.intersp = 1.5, bty = "n", cex = 1)

  }



}

# Weibull hazard function
hazard <- function(t, shape, scale){ # Funktion noch auslagern?
  return( (shape / scale) * (t / scale) ^ (shape -1) )
}

# RMST over tau
RMST_over_tau <- function(tau, shape = 1, scale) {
  sapply(tau, function(tau) {
    stats::integrate(stats::pweibull,
                     shape = shape,
                     scale = scale, lower = 0,
                     upper = tau,
                     lower.tail = F)$value
  })
}

RMSTD_over_tau <-  function(tau, shape_trmt = 1, scale_trmt, shape_ctrl = 1, scale_ctrl){
  sapply(tau, function(tau){
    RMST_over_tau(tau = tau, shape = shape_trmt, scale = scale_trmt) -
      RMST_over_tau(tau = tau, shape = shape_ctrl, scale = scale_ctrl)
  })

}

