#' Produces extensive surival plots
#'
#' Plots example survival data, survival functions, reverse KM plot, recruitment plot, and censoring functions differentiating causes for censoring.
#'
#' @param scale_trmt Specifies the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl Specifies the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored if \code{TRUE}.
#' @param n Specifies the total sample size. Increases to next even number if uneven. Group sample sizes are assumed to be equal.
#' @param loss_scale Specifies the \dfn{scale parameter} of loss to follow-up. No loss to follow-up is assumed if undefined.
#' @param loss_shape Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential loss.
#' @param plot_reverse_KM Logical. Will plot a reverse KM curve if \code{c(TRUE)}, indicating censure-free follow-up.
#' @param plot_log_log Logical. Will plot a log-log plot for assessing proportionality of hazards if \code{TRUE}.
#' @param plot_recruitment Logical. Plots recruitment plot indicating time between recruitment and last observation in study time. Will
#' plot a representative sample of size \eqn{100} if \code{n} \eqn{> 100}.
#' @param plot_extended Logical. Will produce extended plots differentiating causes of loss to follow up.
#' @param xlim Range of plot x-axis. Defaults to \code{c(0, 1.5*tau)}.
#' @param ylim Range of plot y-axis as survival percentages. Defaults to \code{c(0, 100)}.
#'
#' @export
#'
#' @examples
#'
#' # plot full range of plots with sample size n = 1000
#' args_plot <- list(scale_trmt = 10,
#' scale_ctrl = 6,
#' accrual_time = 6,
#' follow_up_time = 3,
#' tau = 4,
#' loss_scale = 10,
#' n = 1000,
#' plot_reverse_KM = TRUE,
#' plot_log_log = TRUE,
#' plot_recruitment = TRUE,
#' plot_extended = TRUE)
#' do.call(plot_surv, args = args_plot)
#'
#' # crossing survival curves
#' args_cross <- args_plot
#' args_cross$shape_ctrl <- .7
#' args_cross$shape_trmt <- 1.3
#' do.call(plot_surv, args = args_cross)
#'
plot_surv <- function(
    scale_trmt,
    scale_ctrl,
    shape_trmt = 1,
    shape_ctrl = 1,
    parameterization = 1,
    accrual_time = 0,
    follow_up_time = Inf,
    tau = NULL,
    censor_beyond_tau = FALSE,
    n = NULL,
    loss_scale = NULL,
    loss_shape = 1,
    plot_reverse_KM = FALSE,
    plot_log_log = FALSE,
    plot_recruitment = FALSE,
    plot_extended = FALSE,
    xlim = NULL,
    ylim = c(0, 100)
    ) {

# reparameterize ----------------------------------------------------------
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

# simulate data and create surv object ------------------------------------
  data_frame_ctrl <- simulate_data(
    scale = scale_ctrl,
    shape = shape_ctrl,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    loss_scale = loss_scale,
    loss_shape = loss_shape,
    n = round(n / 2),
    label = 0
  )
  data_frame_trmt <- simulate_data(
    scale = scale_trmt,
    shape = shape_trmt,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    loss_scale = loss_scale,
    loss_shape = loss_shape,
    n = round(n / 2),
    label = 1
  )
  simulated_data <- rbind(
    data_frame_ctrl, data_frame_trmt)

  if(censor_beyond_tau){ # censor all observations beyond tau if requested
    simulated_data$status[simulated_data$observations > tau] <- 0
    simulated_data$observations[simulated_data$observations > tau] <- tau
  }

  surv_obj <- survival::Surv(
    time = simulated_data$observations,
    event = simulated_data$status
  )

# plot --------------------------------------------------------------------

  if (!is.null(tau) && is.null(xlim)) { # define xlim in relation to tau if not specified
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

  if (!is.null(tau)) {
    # mark tau if defined
    graphics::abline(v = tau, col = "black", lwd = 2)
    graphics::text(
      x = tau,
      y = 0.1,
      pos = 4,
      labels = bquote("Time horizon " * tau * " = " * .(tau)),
      cex = 0.8
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
    "bottomleft",
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
    cex = .8
  )

  # reverse KM
  if (plot_reverse_KM) {
    # reverse indicator for event and censure
    simulated_data_reverse <- simulated_data
    simulated_data_reverse$status <- as.numeric(simulated_data$status == 0) # reverse status

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
      cex = 0.8
    )

    # mark tau if defined
    if (!is.null(tau)) {
      # mark tau if defined
      graphics::abline(v = tau, col = "black", lwd = 2)
      graphics::text(
        x = tau,
        y = 0.1,
        pos = 4,
        labels = bquote("Time horizon " * tau * " = " * .(tau)),
        cex = .8
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
      cex = .8
    )
  }
  if (plot_log_log) {

    plot(
      survival::survfit(surv_obj ~ simulated_data$label),
      xlim = c(min(simulated_data$observations), max(simulated_data$observations)),
      fun = "cloglog",
      xlab = "t",
      ylab = expression(-log(-log(hat(S)(t)))),
      main = "Parallel log-log curves indicate proportional hazards",
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
      cex = 0.8
    )
  }
  if (plot_recruitment) {
    # if sample size is above 100, only display sample of 100. Else display all.
    if (100 < n) {
      indices <- round(seq(from = 1, to = n, length.out = 100))
      sample_size_recruitment <- 100
    } else {
      indices <- 1:n
      sample_size_recruitment <- n
    }
    df_recruitment <- simulated_data[indices, ]
    df_recruitment$accrual_timepoint <- stats::runif(sample_size_recruitment, min = 0, max = accrual_time)
    df_recruitment$last_observation <- df_recruitment$accrual_timepoint + df_recruitment$observations
    df_recruitment <- df_recruitment[order(df_recruitment$accrual_timepoint), ]

    make_color <- function(label, status){ # fun generates colors according to group label and status
      base_col = ifelse(label == 0, "red", "darkblue")
      opacity  <- ifelse(status == 1, 1, 0.3)
      grDevices::adjustcolor(base_col, alpha.f = opacity)
    }
    cols <- mapply(make_color, df_recruitment$label, df_recruitment$status)

    plot( # prepare empty plot
      NA, NA,
      xlim = c(min(df_recruitment$accrual_timepoint), max(df_recruitment$last_observation)),
      ylim = c(0, sample_size_recruitment),
      xlab = "Time",
      main = "Time from recruitment to last observation in study time",
      ylab = NA,
      yaxt = "n",          # no default y ticks
    )
    graphics::segments(
      x0 = c(df_recruitment$accrual_timepoint),
      y0 = 1:100,
      x1 = c(df_recruitment$last_observation),
      y1 = 1:100,
      lwd = 2,
      col = cols
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
      main = "Extended plot on survival and censoring",
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
      cex = 0.8
    )
    # plot RMST(tau), RMSTD(tau), and RMSTR(tau)
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
