#' Produces extensive survival plots
#'
#' Plots example survival data, survival functions, reverse KM plot, recruitment plot, and censoring functions differentiating causes for censoring.
#'
#' @param scale_ctrl Specifies the \dfn{scale parameter} in the control group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_trmt Specifies the \dfn{scale parameter} in the treatment group. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull).
#' @param scale_loss Specifies the \dfn{scale parameter} for loss to follow-up. Can be a scalar (Weibull or exponential survival), or a vector (piecewise Weibull). No loss to follow-up is assumed if undefined.
#' @param shape_ctrl Specifies the \dfn{shape parameter} in the control group. Defaults to \code{shape_ctrl = 1}, simplifying to exponential survival. If \code{length(shape_ctrl) = 1} and \code{length(scale_ctrl) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_trmt Specifies the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt = 1}, simplifying to exponential survival. If \code{length(shape_trmt) = 1} and \code{length(scale_trmt) > 1}, the same shape parameter will be assumed for each section of the survival function.
#' @param shape_loss Specifies the \dfn{shape parameter} for loss to follow-up. Defaults to \code{shape_loss = 1}, simplifying to exponential loss. If \code{length(shape_loss) = 1} and \code{length(scale_loss) > 1}, the same shape parameter will be assumed for each section of the loss distribution.
#' @param breakpoints_ctrl Vector of breakpoints of the piecewise Weibull distribution in the control group. Must have length of \code{scale_ctrl} \eqn{-1} and \code{shape_ctrl} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_trmt Vector of breakpoints of the piecewise Weibull distribution in the treatment group. Must have length of \code{scale_trmt} \eqn{-1} and \code{shape_trmt} \eqn{-1}. First element must be \code{> 0}.
#' @param breakpoints_loss Vector of breakpoints of the piecewise Weibull distribution for loss to follow-up. Must have length of \code{scale_loss} \eqn{-1} and \code{shape_loss} \eqn{-1}. First element must be \code{> 0}.
#' @param accrual_time Length of accrual period.
#' @param follow_up_time Length of follow-up period. Set to \code{Inf} if unspecified.
#' @param tau Specifies the time horizon \eqn{\tau} at which to evaluate \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param censor_beyond_tau Logical. All observations past \eqn{\tau} are censored if \code{TRUE}.
#' @param plot_hazards Logical. Will plot hazard rates.
#' @param plot_HR Logical. Will plot the hazard ratio.
#' @param plot_reverse_KM Logical. Will plot a reverse KM curve, indicating censoring-free follow-up.
#' @param plot_log_log Logical. Will plot a log-log plot for assessing proportionality of hazards if \code{TRUE}.
#' @param plot_proportions Logical. Will plot the proportions of participants being lost to events, censporing, or are still among the risk set.
#' @param xlim Range of plot x-axis. Defaults to \code{c(0, 1.5*tau)}.
#' @param ylim Range of plot y-axis as survival percentages. Defaults to \code{c(0, 100)}.
#' @param parameterisation Define only if Weibull function is specified, not for piecewise exponential survival. One of: \itemize{
#' \item \code{parameterisation = 1}: Default. Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}},
#' \item \code{parameterisation = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterisation = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}}. This is the parameterisation used for the base \R{} function \code{stats::pweibull()}.}
#'
#' @export
#'
#' @examples
#'
#' args_plot <- list(
#'   scale_ctrl = 0.17,
#'   scale_trmt = 0.1,
#'   accrual_time = 6,
#'   follow_up_time = 3,
#'   tau = 4,
#'   scale_loss = 0.1,
#'   plot_hazards = TRUE,
#'   plot_HR = TRUE,
#'   plot_reverse_KM = TRUE,
#'   plot_log_log = TRUE,
#'   plot_proportions = TRUE
#' )
#' do.call(plot_surv, args = args_plot)
#'
#' # crossing survival curves
#' args_cross <- args_plot
#' args_cross$shape_ctrl <- .7
#' args_cross$shape_trmt <- 1.3
#' do.call(plot_surv, args = args_cross)
#'
#' # piecewise exponential
#' args_pex <- args_plot
#' args_pex$scale_ctrl <- c(0.17, 0.4, 0.5)
#' args_pex$breakpoints_ctrl <- c(1, 2)
#' do.call(plot_surv, args = args_pex)
#'
plot_surv <- function(
  scale_ctrl,
  scale_trmt,
  scale_loss = NULL,
  shape_ctrl = 1,
  shape_trmt = 1,
  shape_loss = 1,
  breakpoints_ctrl = NULL,
  breakpoints_trmt = NULL,
  breakpoints_loss = NULL,
  accrual_time = 0,
  follow_up_time = Inf,
  tau = NULL,
  censor_beyond_tau = FALSE,
  plot_hazards = TRUE,
  plot_HR = TRUE,
  plot_reverse_KM = TRUE,
  plot_log_log = TRUE,
  plot_proportions = TRUE,
  xlim = NULL,
  ylim = c(0, 100),
  parameterisation = 1
) {

# error management --------------------------------------------------------
  browser()
  if (length(shape_ctrl) == 1 & length(scale_ctrl) > 1) shape_ctrl <- rep(1, length(scale_ctrl))
  if (length(shape_trmt) == 1 & length(scale_trmt) > 1) shape_trmt <- rep(1, length(scale_trmt))
  if (length(shape_loss) == 1 & length(scale_loss) > 1) shape_loss <- rep(1, length(scale_loss))
  check_inputs(scale_ctrl = scale_ctrl,
               scale_trmt = scale_trmt,
               scale_loss = scale_loss,
               shape_ctrl = shape_ctrl,
               shape_trmt = shape_trmt,
               shape_loss = shape_loss,
               breakpoints_ctrl = breakpoints_ctrl,
               breakpoints_trmt = breakpoints_trmt,
               breakpoints_loss = breakpoints_loss,
               accrual_time = accrual_time,
               follow_up_time = follow_up_time,
               parameterisation = parameterisation,
               sides = 1
  )
  breakpoints_ctrl <- normalize_breakpoints(breakpoints_ctrl)
  breakpoints_trmt <- normalize_breakpoints(breakpoints_trmt)
  breakpoints_loss <- normalize_breakpoints(breakpoints_loss)



  # capture scales to use later for plot function
  original_scales <- list(scale_ctrl = scale_ctrl,
                          scale_trmt = scale_trmt,
                          scale_loss = scale_loss)
  # reparameterise ----------------------------------------------------------
  if (parameterisation != 1) {
    scale_ctrl <- reparameterise(
      parameterisation = parameterisation,
      scale = scale_ctrl,
      shape = shape_ctrl
    )
    scale_trmt <- reparameterise(
      parameterisation = parameterisation,
      scale = scale_trmt,
      shape = shape_trmt
    )
    scale_loss <- reparameterise(
      parameterisation = parameterisation,
      scale = scale_loss,
      shape = shape_loss
    )
  }

    if (!is.null(tau) && is.null(xlim)) { # define xlim in relation to tau if not specified
    xlim <- c(0, 1.5 * tau)
  }
  # shared helpers ------------------------------------------------------------
  x_grid <- seq(xlim[1], xlim[2], length.out = 1000)

  percent_y_scale <- ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = paste0(seq(0, 100, by = 20), "%")
  )

  # Adds a vertical line and label marking tau, if tau is defined.
  tau_layers <- function(tau) {
    if (is.null(tau)) return(NULL)
    list(
      ggplot2::geom_vline(xintercept = tau, color = "black", linewidth = 1),
      ggplot2::annotate(
        "text",
        x = tau,
        y = 0.1,
        hjust = 0,
        label = paste0("'Time horizon ' * tau * ' = ' * ", tau),
        parse = TRUE,
        size = 3
      )
    )
  }

  # plot survival--------------------------------------------------------------------

  ctrl_label <- paste0(
    "Control group with \n",
    "scale = ", paste(round(original_scales$scale_ctrl, 2), collapse = ", "),
    " and shape = ", paste(round(shape_ctrl, 2), collapse = ", ")
  )
  trmt_label <- paste0(
    "Treatment group with \n",
    "scale = ", paste(round(original_scales$scale_trmt, 2), collapse = ", "),
    " and shape = ", paste(round(shape_trmt, 2), collapse = ", ")
  )

  surv_data <- data.frame(t = x_grid)
  surv_data$Control <- ppweibull::ppweibull(q = x_grid, alpha = shape_ctrl, rate = scale_ctrl^shape_ctrl, t = breakpoints_ctrl, lower.tail = FALSE)
  surv_data$Treatment <- ppweibull::ppweibull(q = x_grid, alpha = shape_trmt, rate = scale_trmt^shape_trmt, t = breakpoints_trmt, lower.tail = FALSE)
  surv_data_long <- stats::reshape(
    surv_data,
    varying = c("Control", "Treatment"),
    v.names = "S",
    timevar = "group",
    times = c(ctrl_label, trmt_label),
    direction = "long"
  )

  p_surv <- ggplot2::ggplot(surv_data_long, ggplot2::aes(x = t, y = S, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    tau_layers(tau) +
    ggplot2::scale_color_manual(values = c(stats::setNames("red", ctrl_label), stats::setNames("darkblue", trmt_label))) +
    percent_y_scale +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
    ggplot2::labs(x = "t", y = "S(t) in %", color = NULL, title = "Survival functions for treatment and control group") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  print(p_surv)

  # plot HR --------------------------------------------------------------------

  if (plot_HR) {
    hr_data <- data.frame(
      t = x_grid,
      HR = get_h(x = x_grid, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt) /
        get_h(x = x_grid, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl)
    )

    p_hr <- ggplot2::ggplot(hr_data, ggplot2::aes(x = t, y = HR)) +
      ggplot2::geom_line(color = "black", linewidth = 1.2) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(x = "t", y = "Hazard ratio") +
      ggplot2::theme_bw()

    print(p_hr)
  }

  # plot hazards--------------------------------------------------------------------

  if (plot_hazards) {
    hazard_ctrl <- get_h(x = x_grid, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl)
    hazard_trmt <- get_h(x = x_grid, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt)
    ylim_hazards <- range(c(hazard_ctrl, hazard_trmt), finite = TRUE)

    hazard_data <- data.frame(t = x_grid)
    hazard_data$Control <- hazard_ctrl
    hazard_data$Treatment <- hazard_trmt
    hazard_data_long <- stats::reshape(
      hazard_data,
      varying = c("Control", "Treatment"),
      v.names = "hazard",
      timevar = "group",
      times = c("Hazard rate in control group", "Hazard rate in treatment group"),
      direction = "long"
    )

    p_hazard <- ggplot2::ggplot(hazard_data_long, ggplot2::aes(x = t, y = hazard, color = group)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_color_manual(values = c(
        "Hazard rate in control group" = "red",
        "Hazard rate in treatment group" = "darkblue"
      )) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim_hazards) +
      ggplot2::labs(x = "t", y = "Hazard rates", color = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    print(p_hazard)
  }

  # plot loglog--------------------------------------------------------------------
  if (plot_log_log) {
    loglog_ctrl <- -log(-log(ppweibull::ppweibull(q = x_grid, alpha = shape_ctrl, rate = scale_ctrl^shape_ctrl, t = breakpoints_ctrl, lower.tail = FALSE)))
    loglog_trmt <- -log(-log(ppweibull::ppweibull(q = x_grid, alpha = shape_trmt, rate = scale_trmt^shape_trmt, t = breakpoints_trmt, lower.tail = FALSE)))
    ylim_loglog <- range(c(loglog_ctrl, loglog_trmt), finite = TRUE)

    loglog_data <- data.frame(t = x_grid)
    loglog_data$Control <- loglog_ctrl
    loglog_data$Treatment <- loglog_trmt
    loglog_data_long <- stats::reshape(
      loglog_data,
      varying = c("Control", "Treatment"),
      v.names = "loglog",
      timevar = "group",
      times = c("Control group", "Treatment group"),
      direction = "long"
    )

    p_loglog <- ggplot2::ggplot(loglog_data_long, ggplot2::aes(x = t, y = loglog, color = group)) +
      ggplot2::geom_line(linewidth = 1) +
      tau_layers(tau) +
      ggplot2::scale_color_manual(values = c("Control group" = "red", "Treatment group" = "darkblue")) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim_loglog) +
      ggplot2::labs(x = "t", y = "-log-log(S(t))", color = NULL, title = "-log-log plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    print(p_loglog)
  }

  # potential follow-up -----------------------------------------------------

  follow_up_data <- data.frame(t = x_grid)
  follow_up_data$p_not_censored <- sapply(
    x_grid,
    function(xi) get_p_not_censored(
      x = xi,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss
    )
  )

  p_follow_up <- ggplot2::ggplot(follow_up_data, ggplot2::aes(x = t, y = p_not_censored)) +
    ggplot2::geom_line(linewidth = 1) +
    tau_layers(tau) +
    percent_y_scale +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
    ggplot2::labs(x = "t", y = "Proportion under observation in %", title = "Potential follow-up") +
    ggplot2::theme_bw()

  print(p_follow_up)

  # plot stacked area charts of proportions --------------------------------------------

  if (plot_proportions) {
    band_levels <- c(
      "Under observation",
      "Lost to event",
      "Lost to follow-up",
      "Lost to administrative censoring"
    )
    band_colors <- c(
      "Under observation" = "#BA1650",
      "Lost to event" = "#F5C700",
      "Lost to follow-up" = "#FFD3F0",
      "Lost to administrative censoring" = "#B7E7FC"
    )

    # Number of points to keep per band for plotting. The competing-risk
    # probabilities themselves are still computed on the full, fine x_area
    # grid for numerical accuracy; only the plotted curve is thinned, since
    # rendering geom_ribbon() with hundreds of thousands of points per band
    # is dramatically slower than the underlying computation (unlike base
    # graphics::polygon(), which draws large point counts cheaply).
    max_plot_points <- 2000

    # Converts cumulative competing-risk probabilities into stacked bands,
    # thinned to at most max_plot_points points for fast rendering.
    build_proportion_bands <- function(x, probs) {
      cum0 <- rep(0, length(x))
      cum1 <- probs$p_obs
      cum2 <- cum1 + probs$p_event
      cum3 <- cum2 + probs$p_loss
      cum4 <- cum3 + probs$p_admin

      if (length(x) > max_plot_points) {
        idx <- unique(c(1, round(seq(1, length(x), length.out = max_plot_points)), length(x)))
        x <- x[idx]
        cum0 <- cum0[idx]; cum1 <- cum1[idx]; cum2 <- cum2[idx]; cum3 <- cum3[idx]; cum4 <- cum4[idx]
      }

      rbind(
        data.frame(t = x, ymin = cum0, ymax = cum1, band = band_levels[1]),
        data.frame(t = x, ymin = cum1, ymax = cum2, band = band_levels[2]),
        data.frame(t = x, ymin = cum2, ymax = cum3, band = band_levels[3]),
        data.frame(t = x, ymin = cum3, ymax = cum4, band = band_levels[4])
      )
    }

    plot_proportion_bands <- function(band_data, title, subtitle) {
      band_data$band <- factor(band_data$band, levels = band_levels)
      ggplot2::ggplot(band_data, ggplot2::aes(x = t, ymin = ymin, ymax = ymax, fill = band)) +
        ggplot2::geom_ribbon() +
        ggplot2::scale_fill_manual(values = band_colors, breaks = band_levels) +
        percent_y_scale +
        ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
        ggplot2::labs(x = "t", y = "Proportion in %", fill = NULL, title = title, subtitle = subtitle) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
    }

    x_area <- seq(xlim[1], xlim[2], by = 0.001)

    # ctrl group: plot stacked area chart -------------------------------------------------
    probs_ctrl <- get_competing_risk_probs(
      x = x_area,
      scale_ctrl = scale_ctrl,
      shape_ctrl = shape_ctrl,
      breakpoints_ctrl = breakpoints_ctrl,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time
    )
    p_ctrl_proportions <- plot_proportion_bands(
      build_proportion_bands(x_area, probs_ctrl),
      title = "Control group:",
      subtitle = "Proportion of Subjects by Event and Censoring Type"
    )
    print(p_ctrl_proportions)

    # trmt group: plot stacked area chart -------------------------------------------------
    probs_trmt <- get_competing_risk_probs(
      x = x_area,
      scale_ctrl = scale_trmt,
      shape_ctrl = shape_trmt,
      breakpoints_ctrl = breakpoints_trmt,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time
    )
    p_trmt_proportions <- plot_proportion_bands(
      build_proportion_bands(x_area, probs_trmt),
      title = "Treatment group:",
      subtitle = "Proportion of Subjects by Event and Censoring Type"
    )
    print(p_trmt_proportions)
  }

  invisible(NULL)
}


# helper function for stacked area chart ----------------------------------

#' Compute competing-risks probabilities for event and loss to follow-up
#'
#' @noRd
get_competing_risk_probs <- function(
    x,
    scale_ctrl,
    shape_ctrl,
    breakpoints_ctrl,
    scale_loss,
    shape_loss,
    breakpoints_loss,
    accrual_time,
    follow_up_time
) {

  # Create a fine grid for numerical integration
  t_grid <- x
  dt <- t_grid[2] - t_grid[1]
  t_grid_admin_loss <- seq(0, accrual_time + follow_up_time, by = dt)

  # Hazards on the grid
  hE <- get_h(t_grid, scale = scale_ctrl, shape = shape_ctrl,
              breakpoints = breakpoints_ctrl)
  if(is.null(scale_loss)){
    hL <- rep(0, length(t_grid))} else{
  hL <- get_h(t_grid, scale = scale_loss, shape = shape_loss,
              breakpoints = breakpoints_loss)
    }
  hA_temp <- rep(0, length(t_grid_admin_loss))
  hA_temp[(follow_up_time / dt) : ((follow_up_time + accrual_time) / dt)] <-
    1 / (accrual_time - seq(0, accrual_time, by = dt))
  # t_grid may extend beyond accrual_time + follow_up_time (e.g. when the
  # plotting xlim exceeds the trial length); pad with 0 there instead of
  # letting out-of-range indexing silently introduce NAs, since there is no
  # further administrative censoring hazard once the trial has ended.
  hA <- rep(0, length(t_grid))
  n_copy <- min(length(hA_temp), length(t_grid))
  hA[seq_len(n_copy)] <- hA_temp[seq_len(n_copy)]

  max_finite <- max(hE[is.finite(hE)]) # in case some element in h_all is Inf
  hE[is.infinite(hE)] <- max_finite
  max_finite <- max(hL[is.finite(hL)]) # in case some element in h_all is Inf
  hE[is.infinite(hL)] <- max_finite


  # Cumulative overall hazard H_all(t) = int_0^t [hE(u) + hL(u)] du
  h_all <- hE + hL + hA


  H_all <- cumsum(h_all) * dt

  # Overall survival (no event, no loss) S_all(t) = exp(-H_all(t))
  S_all <- exp(-H_all)

  # Integrand for CIFs: h_k(t) * S_all(t)
  integrand_E <- hE * S_all
  integrand_L <- hL * S_all
  integrand_A <- hA * S_all

  # Cumulative incidence functions (numerical integration)
  F_event_grid <- cumsum(integrand_E) * dt
  F_loss_grid  <- cumsum(integrand_L) * dt
  F_admin_grid  <- cumsum(integrand_A) * dt

  # Interpolate back to the requested x values
  p_obs    <- stats::approx(t_grid, S_all,        xout = x, rule = 2)$y
  p_event  <- stats::approx(t_grid, F_event_grid, xout = x, rule = 2)$y
  p_loss   <- stats::approx(t_grid, F_loss_grid,  xout = x, rule = 2)$y
  p_admin  <- stats::approx(t_grid, F_admin_grid,  xout = x, rule = 2)$y

  # Small numerical correction: ensure they sum to 1
  total <- p_obs + p_event + p_loss + p_admin
  p_obs   <- p_obs   / total
  p_event <- p_event / total
  p_loss  <- p_loss  / total
  p_admin <- p_admin / total

  list(
    p_obs   = p_obs,
    p_event = p_event,
    p_loss  = p_loss,
    p_admin = p_admin
  )
}

utils::globalVariables(c("HR", "hazard", "loglog", "p_not_censored", "ymin", "ymax", "band"))
