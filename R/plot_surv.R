#' Plots example survival data, survival functions, and other characteristics
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
#' \item \code{parameterization = 3}: Specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
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
#' @param plot_recruitment Boolean. Plots recruitment plot.
#' @param plot_extended Boolean. Will produce extended plots (..)

#'
#' @export
#'
#' @examples
#' plot_surv(
#'   scale_trmt = 1.4,
#'   scale_ctrl = 1,
#'   accrual_time = 1,
#'   follow_up_time = 10,
#'   tau = 1,
#'   loss_scale = 0.2,
#'   n = 100,
#'   plot_reverse_KM = TRUE,
#'   plot_log_log = TRUE,
#'   plot_extended = TRUE
#' )
plot_surv <- function(
    scale_trmt,
    scale_ctrl,
    shape_trmt = 1,
    shape_ctrl = 1,
    parameterization = 1,
    tau = NULL,
    xlim = NULL,
    ylim = c(0, 100),
    n = NULL,
    accrual_time = 0,
    follow_up_time = Inf,
    censor_beyond_tau = FALSE,
    loss_scale = NULL,
    loss_shape = 1,
    plot_reverse_KM = FALSE,
    plot_log_log = FALSE,
    plot_recruitment = FALSE,
    plot_extended = FALSE
    ) {
  browser()

  if (parameterization != 1){
    scale_trmt <- reparameterize(parameterization = parameterization,
                                 scale = scale_trmt,
                                 shape = shape_trmt)
    scale_ctrl <- reparameterize(parameterization = parameterization,
                                 scale = scale_ctrl,
                                 shape = shape_ctrl)
    loss_scale <- reparameterize(parameterization = parameterization,
                                 scale = loss_scale,
                                 shape = loss_shape)
  }

  # create data_frame_ctrl

  data_frame_ctrl <- simulate_data(
    scale = scale_ctrl,
    shape = shape_ctrl,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    loss_scale = loss_scale,
    loss_shape = loss_shape,
    sample_size = round(n / 2),
    label = 0
  )
  simulated_data <- rbind(
    data_frame_ctrl,
    simulate_data(
      scale = scale_trmt,
      shape = shape_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      loss_scale = loss_scale,
      loss_shape = loss_shape,
      sample_size = round(n / 2),
      label = 1
    )
  )
  surv_obj <- survival::Surv(
    time = simulated_data$observations,
    event = simulated_data$status
  )

  if (!is.null(tau) && is.null(xlim)) {
    xlim <- c(0, 1.5 * tau)
  } # set xlim if undefined
  # plot KM estimator
  graphics::par(mar = c(5, 6, 4, 1) + .1)
  plot(
    survival::survfit(surv_obj ~ simulated_data$label),
    mark.time = T,
    conf.int = F,
    xlab = "t",
    ylab = expression(hat(S)(t) ~ "in %"),
    col = c("red", "darkblue"),
    xlim = xlim,
    ylim = c(0, 1),
    lwd = 2,
    main = "Kaplan Meier estimators for treatment and control group",
    yaxt = "n"
  )

  graphics::axis(
    2,
    at = seq(1, 0, by = -0.2),
    labels = paste0(seq(100, 0, by = -20), "%"),
    las = 1
  )

  # mark tau if defined
  if (!is.null(tau)) {
    # mark tau if defined
    graphics::abline(v = tau, col = "black", lwd = 2)
    graphics::text(
      x = tau,
      y = 0.1,
      pos = 4,
      labels = bquote("Time horizon " * tau * " = " * .(tau))
    )
  }

  # draw design curves
  graphics::curve(
    stats::pweibull(
      x,
      shape = shape_trmt,
      scale = scale_trmt,
      lower.tail = FALSE
    ),
    from = xlim[1],
    to = xlim[2],
    add = TRUE,
    col = "darkblue",
    lwd = 2,
    lty = 2
  )
  graphics::curve(
    stats::pweibull(
      x,
      shape = shape_ctrl,
      scale = scale_ctrl,
      lower.tail = FALSE
    ),
    from = xlim[1],
    to = xlim[2],
    add = TRUE,
    col = "red",
    lwd = 2,
    lty = 2
  )

  # create legend
  graphics::legend(
    "topright",
    legend = c(
      paste0(
        "Treatment group with \n",
        "scale = ",
        round(scale_trmt, 2),
        " and shape = ",
        round(shape_trmt, 2)
      ),
      paste0(
        "Control group with \n",
        "scale = ",
        round(scale_ctrl, 2),
        " and shape = ",
        round(shape_ctrl, 2)
      )
    ),
    col = c("darkblue", "red"),
    lty = 1:1,
    y.intersp = 1.5,
    bty = "n",
    cex = 1
  )

  # reverse KM
  if (plot_reverse_KM) {
    # reverse indicator for event and censure
    simulated_data_reverse <- simulated_data
    simulated_data_reverse$status <- as.numeric(simulated_data$status == 0)

    surv_obj_reverse <- survival::Surv(
      time = simulated_data_reverse$observations,
      event = simulated_data_reverse$status
    )

    # create plot reverse KM
    plot(
      survival::survfit(surv_obj_reverse ~ simulated_data_reverse$label),
      mark.time = T,
      conf.int = F,
      xlab = "t",
      ylab = "Censure-free observations in %",
      col = c("red", "darkblue"),
      xlim = xlim,
      ylim = c(0, 1),
      lwd = 2,
      yaxt = "n",
      main = "Reverse Kaplan Meier estimators for treatment and control group"
    )
    graphics::axis(
      2,
      at = seq(0, 1, by = 0.2),
      labels = paste0(seq(0, 100, by = 20), "%")
    )

    graphics::text(
      x = graphics::par("usr")[2] - 0.01,
      y = graphics::par("usr")[4] - 0.01,
      labels = "+ indicates event",
      adj = c(1, 1),
      cex = 0.9
    )

    # draw reverse design curves

    # mark tau if defined
    if (!is.null(tau)) {
      # mark tau if defined
      graphics::abline(v = tau, col = "black", lwd = 2)
      graphics::text(
        x = tau,
        y = 0.1,
        pos = 4,
        labels = bquote("Time horizon " * tau * " = " * .(tau))
      )
    }

    # create legend
    graphics::legend(
      "bottomleft",
      legend = c(
        "Treatment group",
        "Control group"
      ),
      col = c("darkblue", "red"),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = 1
    )
  }
  if (plot_log_log) {

    plot(
      survival::survfit(surv_obj ~ simulated_data$label),
      xlim = c(min(simulated_data$observations), ),
      fun = "cloglog",
      xlab = "t",
      ylab = "-ln[-ln(survival)]",
      main = "Log-log curves for assessing proportionality of hazards",
      col = c("red", "darkblue"),
      lwd = 2
    )
    graphics::legend(
      "topleft",
      legend = c(
        paste0(
          "Treatment group with \n",
          "scale = ",
          round(scale_trmt, 2),
          " and shape = ",
          round(shape_trmt, 2)
        ),
        paste0(
          "Control group with \n",
          "scale = ",
          round(scale_ctrl, 2),
          " and shape = ",
          round(shape_ctrl, 2)
        )
      ),
      col = c("darkblue", "red"),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = 1
    )
  }
  if (plot_recruitment) {
    indices <- 1:sample_size # indices for y-axis
    sample_size_plot <- sample_size
    # if sample size is above 100, only display sample of 100
    if (100 < sample_size) {
      indices <- round(seq(from = 1, to = sample_size, length.out = 100))
      sample_size_plot <- 100
    }
    if (total_time == Inf) {
      xlim <- c(0, max(observations))
    } else {
      xlim <- c(0, total_time)
    }
    if (follow_up_time == Inf) {
      recruitment <- rep(0, sample_size)
    } else {
      recruitment <- (total_time - admin_loss)[indices]
    } # recruitment time in study time
    stop <- recruitment + observations[indices] # last observation in study time
    df <- data.frame(recruitment, stop)
    df_sorted <- df[order(recruitment), ]
    plot(
      x = df_sorted$recruitment,
      y = 1:sample_size_plot,
      xlim = xlim,
      ylim = c(0, sample_size_plot),
      xlab = "Study time",
      ylab = "Participants ordered by entry into study"
    )
    graphics::title("Individual observation trails in study time")
    graphics::segments(
      x0 = df_sorted$recruitment,
      y0 = 1:length(recruitment),
      x1 = df_sorted$stop
    )
  }
  if (plot_extended) {
    # admin censoring
    if (follow_up_time != Inf && is.null(xlim)) {
      xlim <- c(0, follow_up_time + accrual_time)
    }
    graphics::curve(
      100 *
        sapply(
          x,
          get_p_not_lost_admin,
          follow_up_time = follow_up_time,
          accrual_time = accrual_time
        ),
      xlim = xlim,
      ylim = c(0, 100),
      xlab = "t",
      ylab = "Proportion remaining in %",
      main = "Extended plot differentiating causes for censoring",
      lwd = 2,
      col = "lightgrey"
    )
    # pts not lost to FU
    if (is.null(loss_scale)) {
      graphics::abline(h = 100, col = "darkgrey", lwd = 2)
    } else {
      graphics::curve(
        100 *
          stats::pweibull(
            q = x,
            shape = loss_shape,
            scale = loss_scale,
            lower.tail = FALSE
          ),
        lwd = 2,
        col = "darkgrey",
        add = TRUE
      )
    }
    # pts not lost to censoring
    graphics::curve(
      100 *
        sapply(
          x,
          get_p_not_censored,
          follow_up_time = follow_up_time,
          accrual_time = accrual_time,
          loss_scale = loss_scale,
          loss_shape = loss_shape
        ),
      lwd = 2,
      col = "black",
      add = TRUE
    )

    # pts not lost to neither censoring nor event ctrl
    graphics::curve(
      100 *
        sapply(
          x,
          get_p_at_risk,
          scale = scale_ctrl,
          shape = shape_ctrl,
          loss_scale = loss_scale,
          loss_shape = loss_shape,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        ),
      lwd = 2,
      col = "#EA95BA",
      add = TRUE
    )

    # pts not lost to neither censoring nor event trmt
    graphics::curve(
      100 *
        sapply(
          x,
          get_p_at_risk,
          scale = scale_trmt,
          shape = shape_trmt,
          loss_scale = loss_scale,
          loss_shape = loss_shape,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        ),
      lwd = 2,
      col = "steelblue1",
      add = TRUE
    )

    graphics::curve(
      100 *
        stats::pweibull(
          x,
          shape = shape_ctrl,
          scale = scale_ctrl,
          lower.tail = FALSE
        ),
      add = TRUE,
      col = "red",
      lwd = 2
    )

    graphics::curve(
      100 *
        stats::pweibull(
          x,
          shape = shape_trmt,
          scale = scale_trmt,
          lower.tail = FALSE
        ),
      add = TRUE,
      col = "darkblue",
      lwd = 2
    )

    graphics::legend(
      "topright",
      legend = c(
        "Share not lost to admnistrative censoring",
        "Share not lost to FU",
        "Share lost neither to FU nor to administrative censoring",
        "Share in control group lost to neither events nor to censoring",
        "Share in treatment group lost to neither events nor to censoring",
        "Survival in control group",
        "Survival in treatment group"
      ),
      col = c(
        "lightgrey",
        "darkgrey",
        "black",
        "#EA95BA",
        "steelblue1",
        "red",
        "darkblue"
      ),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = 1
    )
    # plot RMST(tau), RMSTD(tau), and RMSTR(tau)
    browser()
    if (follow_up_time == Inf){
      upper_limit <- 1 * tau
      graphics::curve(
        100 *
          sapply(
            x,
            get_theoretical_rmst,
            scale = scale_ctrl,
            shape = shape_ctrl),
        lwd = 2,
        col = "#EA95BA"
      )
    }
  }
}

utils::globalVariables(c("x")) # prevents warnings on undefined variables when running devtools::check()
