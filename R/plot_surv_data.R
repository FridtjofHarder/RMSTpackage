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
#' @param plot_log_log Boolean. Will plot a log-log plot for assessing proportionality of hazards if \code{c(TRUE)}.
#' @param plot_extended Boolean. Will produce extended plots (..)
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
#' loss_scale = 0.2,
#' n = 100,
#' plot_reverse_KM = TRUE,
#' plot_log_log = TRUE,
#' plot_extended = TRUE)
plot_surv_data <- function(scale_trmt, scale_ctrl, shape_trmt = 1,
                           shape_ctrl = 1, parameterization = 1, tau = NULL,
                           xlim = NULL, ylim = c(0, 100), n = NULL,
                           accrual_time = 0, follow_up_time = Inf,
                           censor_beyond_tau = FALSE, loss_scale = NULL,
                           loss_shape = 1, plot_reverse_KM = FALSE,
                           plot_log_log = FALSE,
                           plot_extended){

  # create data_frame_ctrl
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

  if (!is.null(tau) && is.null(xlim)) {xlim <- c(0, 1.5*tau)} # set xlim if undefined

  # plot KM estimator
  par(mar=c(5,6,4,1)+.1)
  plot(survival::survfit(surv_obj~simulated_data$label), mark.time=T,
       conf.int = F, xlab = "t", ylab = expression(hat(S)(t) ~ "in %"),
       col = c("red", "darkblue"),
       xlim = xlim, ylim = c(0, 1), lwd = 2,
       main = "Kaplan Meier estimators for treatment and control group",
       yaxt = "n")

  axis(2, at = seq(1, 0, by = -0.2), labels = paste0(seq(100, 0, by = -20), "%"), las = 1)

  # mark tau if defined
  if (!is.null(tau)){ # mark tau if defined
    graphics::abline(v = tau, col = "black", lwd = 2)
    graphics::text(x = tau, y = 0.1, pos = 4,
                   labels = bquote("Time horizon " * tau * " = " * .(tau)))
  }

  # draw design curves
  curve(pweibull(x, shape = shape_trmt, scale = scale_trmt, lower.tail = FALSE),
        from = xlim[1], to = xlim[2],
        add = TRUE, col = "darkblue", lwd = 2, lty = 2)
  curve(pweibull(x, shape = shape_ctrl, scale = scale_ctrl, lower.tail = FALSE),
        from = xlim[1], to = xlim[2],
        add = TRUE, col = "red", lwd = 2, lty = 2)

  # create legend
  graphics::legend("topright", legend=c(paste0("Treatment group with \n", "scale = ",
                                                 round(scale_trmt, 2), " and shape = ",
                                                 round(shape_trmt, 2)),
                                          paste0("Control group with \n", "scale = ",
                                                 round(scale_ctrl, 2), " and shape = ",
                                                 round(shape_ctrl, 2))),
                   col=c("darkblue", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 1)


  # reverse KM
  if (plot_reverse_KM){
    # reverse indicator for event and censure
    simulated_data_reverse <- simulated_data
    simulated_data_reverse$status <- as.numeric(simulated_data$status == 0)

    surv_obj_reverse <- survival::Surv(time = simulated_data_reverse$observations,
                               event = simulated_data_reverse$status)

    # create plot reverse KM
    plot(survival::survfit(surv_obj_reverse~simulated_data_reverse$label), mark.time=T,
         conf.int = F, xlab = "t", ylab = "Censure-free observations in %", col = c("red", "darkblue"),
         xlim = xlim, ylim = c(0, 1), lwd = 2, yaxt = "n",
         main = "Reverse Kaplan Meier estimators for treatment and control group")
    axis(2, at = seq(0, 1, by = 0.2), labels = paste0(seq(0, 100, by = 20), "%"))

    # draw reverse design curves
    # p(uncensored) = p(not lost to FU) X p(not lost to admin)
    if(!is.null(loss_scale) && !is.null(loss_shape)){
      p_uncensored <- NULL
    } else{
      p_uncensored <- NULL
    }

    # mark tau if defined
    if (!is.null(tau)){ # mark tau if defined
      graphics::abline(v = tau, col = "black", lwd = 2)
      graphics::text(x = tau, y = 0.1, pos = 4,
                     labels = bquote("Time horizon " * tau * " = " * .(tau)))
    }

    # create legend
    graphics::legend("bottomleft", legend=c("Treatment group",
                                            "Control group"),
                     col=c("darkblue", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 1)
  }
  #browser()
  if (plot_log_log){
    plot(survival::survfit(surv_obj~simulated_data$label), fun = "cloglog",
         xlab = "t",
         ylab = "log-log survival", main = "Log-log curves for assessing proportionality of hazards",
         col = c("red", "darkblue"), lwd = 2)
    graphics::legend("topleft", legend=c(paste0("Treatment group with \n", "scale = ",
                                                 round(scale_trmt, 2), " and shape = ",
                                                 round(shape_trmt, 2)),
                                          paste0("Control group with \n", "scale = ",
                                                 round(scale_ctrl, 2), " and shape = ",
                                                 round(shape_ctrl, 2))),
                     col=c("darkblue", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 1)
  }
  browser()
  # define arms in npsurvSS
  arm_npsurvSS_0 <- npsurvSS::create_arm(size = 1,
                                         accr_time = accrual_time,
                                         follow_time = follow_up_time,
                                         surv_scale = 1/scale_ctrl,
                                         surv_shape = shape_ctrl,
                                         loss_scale = loss_scale,
                                         loss_shape = loss_shape)
  arm_npsurvSS_1 <- npsurvSS::create_arm(size = 1,
                                         accr_time = accrual_time,
                                         follow_time = follow_up_time,
                                         surv_scale = 1/scale_trmt,
                                         surv_shape = shape_trmt,
                                         loss_scale = loss_scale,
                                         loss_shape = loss_shape)

  if(plot_extended){
    # admin censoring
    curve(100 * p_not_censored_admin_npsurvSS(x, arm_npsurvSS = arm_npsurvSS_0),
         xlim = c(0, follow_up_time + accrual_time),
         ylim = c(0, 100),
         xlab = "t", ylab = "Proportion remaining in %",
         main = "Extended plot differentiating causes for censoring",
         lwd = 2, col = "lightgrey"
         )
    # pts not lost to lost to FU
    curve(100 * npsurvSS::ploss(q = x, arm = arm_npsurvSS_0, lower.tail = FALSE),
          lwd = 2, col = "darkgrey",
          add = TRUE
          )
    # pts not lost to censoring
    curve(100 * p_not_being_censored_npsurvSS(x, arm_npsurvSS = arm_npsurvSS_0),
          lwd = 2, col = "black",
          add = TRUE)

    # pts not lost to neither censoring nor event ctrl
    curve(100 * p_neither_censored_nor_event_npsurvSS(x, arm_npsurvSS = arm_npsurvSS_0),
          lwd = 2, col = "#EA95BA",
          add = TRUE)

    # pts not lost to neither censoring nor event trmt
    curve(100 * p_neither_censored_nor_event_npsurvSS(x, arm_npsurvSS = arm_npsurvSS_1),
          lwd = 2, col = "steelblue1",
          add = TRUE)

    curve(100 * pweibull(x, shape = shape_ctrl, scale = scale_ctrl, lower.tail = FALSE),
          add = TRUE, col = "red", lwd = 2)

    curve(100 * pweibull(x, shape = shape_trmt, scale = scale_trmt, lower.tail = FALSE),
          add = TRUE, col = "darkblue", lwd = 2)

    graphics::legend("topright",
                     legend=c("Share not lost to admnistrative censoring",
                              "Share not lost to FU",
                              "Share lost neither to FU nor to administrative censoring",
                              "Share in control group lost to neither events nor to censoring",
                              "Share in treatment group lost to neither events nor to censoring",
                              "Survival in control group",
                              "Survival in treatment group"),
                     col=c("lightgrey", "darkgrey", "black", "#EA95BA",
                           "steelblue1", "red", "darkblue"),
                     lty=1:1, y.intersp = 1.5, bty = "n", cex = 1)
  }

}

# aux function to calculate probabilty of not being censored by admin censoring in case of linear accrual
p_not_censored_admin_npsurvSS <- function(x, arm_npsurvSS){
  return(npsurvSS::paccr(arm = arm_npsurvSS, q = arm_npsurvSS$total_time - x))
}

# aux function to calculate probabilty of not being censored
p_not_being_censored_npsurvSS <- function(x, arm_npsurvSS){
  return(npsurvSS::paccr(arm = arm_npsurvSS, q = arm_npsurvSS$total_time - x) *
    npsurvSS::ploss(q = x, arm = arm_npsurvSS, lower.tail = FALSE))
}
# aux function to calculate probabilty of not being censored by either admin or loss to FU
p_neither_censored_nor_event_npsurvSS <- function(x, arm_npsurvSS){
  return(npsurvSS::paccr(arm = arm_npsurvSS, q = arm_npsurvSS$total_time - x) *
           npsurvSS::ploss(q = x, arm = arm_npsurvSS, lower.tail = FALSE) *
           npsurvSS::psurv(q = x, arm = arm_npsurvSS, lower.tail = FALSE))
}
