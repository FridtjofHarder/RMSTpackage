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
#' @param plot_HR Logical. Will plot hazard ratio.
#' @param plot_reverse_KM Logical. Will plot a reverse KM curve, indicating censoring-free follow-up.
#' @param plot_log_log Logical. Will plot a log-log plot for assessing proportionality of hazards if \code{TRUE}.
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
#' # plot full range of plots with sample size n = 1000
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
  n = NULL,
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
  # plot survival--------------------------------------------------------------------

  graphics::par(mar = c(5, 6, 4, 1) + .1)
    plot(
      NA,
      xlab = "t",
      ylab = expression(S(t) ~ "in %"),
      xlim = xlim,
      ylim = c(0, 1),
      main = "Survival functions for treatment and control group",
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
    ppweibull::ppweibull(q = x, alpha = shape_ctrl, rate = scale_ctrl^shape_ctrl, t = breakpoints_ctrl, lower.tail = FALSE),
    from = xlim[1],
    to = xlim[2],
    add = TRUE,
    col = "red",
    lwd = 2,
    lty = 1
  )
  graphics::curve(
    ppweibull::ppweibull(q = x, alpha = shape_trmt, rate = scale_trmt^shape_trmt, t = breakpoints_trmt, lower.tail = FALSE),
    from = xlim[1],
    to = xlim[2],
    add = TRUE,
    col = "darkblue",
    lwd = 2,
    lty = 1
  )
  graphics::legend(
    "topright",
    legend = c(
      paste0(
        "Control group with \n",
        "scale = ",
        paste(round(original_scales$scale_ctrl, 2), collapse = ", "),
        " and shape = ",
        paste(round(shape_ctrl, 2), collapse = ", ")
      ),
      paste0(
        "Treatment group with \n",
        "scale = ",
        paste(round(original_scales$scale_trmt, 2), collapse = ", "),
        " and shape = ",
        paste(round(shape_trmt, 2), collapse = ", ")
      )
    ),
    col = c("red", "darkblue"),
    lty = 1:1,
    y.intersp = 1.5,
    bty = "n",
    cex = .8
  )
  # plot HR --------------------------------------------------------------------

  if(plot_HR){
    graphics::curve(
      get_h(x = x, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt) /
        get_h(x = x, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl),
      from = xlim[1],
      to = xlim[2],
      col = "black",
      lwd = 4,
      lty = 1,
      ylab = "Hazard ratio",
      xlab = "t"
    )
  }
  # plot hazards--------------------------------------------------------------------

  if(plot_hazards){
    x_grid <- seq(xlim[1], xlim[2], length.out = 1000)

    hazard_ctrl <- get_h(
      x = x_grid,
      scale = scale_ctrl,
      shape = shape_ctrl,
      breakpoints = breakpoints_ctrl
    )

    hazard_trmt <- get_h(
      x = x_grid,
      scale = scale_trmt,
      shape = shape_trmt,
      breakpoints = breakpoints_trmt
    )

    ylim_hazards <- range(
      c(hazard_ctrl, hazard_trmt),
      finite = TRUE
    )

  graphics::curve(
    get_h(x = x, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl),
    from = xlim[1],
    to = xlim[2],
    ylim = ylim_hazards,
    col = "red",
    lwd = 4,
    lty = 1,
    ylab = "Hazard rates"
  )
  graphics::curve(
    get_h(x = x, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt),
    from = xlim[1],
    to = xlim[2],
    add = TRUE,
    col = "darkblue",
    lwd = 4,
    lty = 1
  )

  graphics::legend(
    "topright",
    legend = c("Hazard rate in control group",
               "Hazard rate in treatment group"),
    col = c("red", "darkblue"),
    lty = 1:1,
    y.intersp = 1.5,
    bty = "n",
    cex = .8
  )
  }

  # plot loglog--------------------------------------------------------------------
  if (plot_log_log) {

    x_grid <- seq(xlim[1], xlim[2], length.out = 1000)
    loglog_ctrl <- -log(-log(ppweibull::ppweibull(q = x_grid, alpha = shape_ctrl, rate = scale_ctrl^shape_ctrl, t = breakpoints_ctrl, lower.tail = FALSE)))
    loglog_trmt <- -log(-log(ppweibull::ppweibull(q = x_grid, alpha = shape_trmt, rate = scale_trmt^shape_trmt, t = breakpoints_trmt, lower.tail = FALSE)))

    ylim_loglog <- range(
      c(loglog_ctrl, loglog_trmt),
      finite = TRUE
    )
    graphics::par(mar = c(5, 6, 4, 1) + .1)
    plot(
      NA,
      xlab = "t",
      ylab = "-log-log(S(t))",
      xlim = xlim,
      ylim = ylim_loglog,
      main = "-log-log plot",
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
      -log(-log(ppweibull::ppweibull(q = x, alpha = shape_ctrl, rate = scale_ctrl^shape_ctrl, t = breakpoints_ctrl, lower.tail = FALSE))),
      from = xlim[1],
      to = xlim[2],
      add = TRUE,
      col = "red",
      lwd = 2,
      lty = 1
    )
    graphics::curve(
      -log(-log(ppweibull::ppweibull(q = x, alpha = shape_trmt, rate = scale_trmt^shape_trmt, t = breakpoints_trmt, lower.tail = FALSE))),
      from = xlim[1],
      to = xlim[2],
      add = TRUE,
      col = "darkblue",
      lwd = 2,
      lty = 1
    )
    graphics::legend(
      "topright",
      legend = c("Control group", "Treatment group"),
      col = c("red", "darkblue"),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = .8
    )
  }

  # potential follow-up -----------------------------------------------------
  graphics::par(mar = c(5, 6, 4, 1) + .1)
  plot(
    NA,
    xlab = "t",
    ylab = "Proportion under observation in %",
    xlim = xlim,
    ylim = c(0, 1),
    main = "Potential follow-up",
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
  my_x <- seq(xlim[1], xlim[2], length.out = 1000)
  my_y <- sapply(
    my_x,
    function(xi) get_p_not_censored(
      x = xi,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss
    )
  )

  graphics::lines(my_x, my_y,
    lwd = 2,
    lty = 1
  )


# ctrl group: plot stacked area chart ctrl -------------------------------------------------
  if(plot_proportions){
  graphics::par(mar = c(5, 6, 4, 1) + .1)
  plot(
    NA,
    xlab = "t",
    ylab = "Proportion in %",
    xlim = xlim,
    ylim = c(0, 1),
    main = "Control group:\nProportion of Subjects by Event and Censoring Type",
    yaxt = "n"
  )

  graphics::axis(
    2,
    at = seq(1, 0, by = -0.2),
    labels = paste0(seq(100, 0, by = -20), "%"),
    las = 1
  )
  x <- seq(xlim[1], xlim[2], by = 0.001)
  probs <- get_competing_risk_probs(
    x = x,
    scale_ctrl = scale_ctrl,
    shape_ctrl = shape_ctrl,
    breakpoints_ctrl = breakpoints_ctrl,
    scale_loss = scale_loss,
    shape_loss = shape_loss,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time
  )
  cum0 <- rep(0, length(x))
  cum1 <- probs$p_obs
  cum2 <- cum1 + probs$p_event
  cum3 <- cum2 + probs$p_loss
  cum4 <- cum3 + probs$p_admin

  # Bottom: under observation
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum1, rev(cum0)),
    col = "#BA1650",
    border = NA
  )

  # Event
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum2, rev(cum1)),
    col = "#F5C700",
    border = NA
  )

  # Loss to follow-up
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum3, rev(cum2)),
    col = "#FFD3F0",
    border = NA
  )

  # Administrative censoring (will be invisible here since p_admin = 0)
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum4, rev(cum3)),
    col = "#B7E7FC",
    border = NA
  )

  legend(
    "topright",
    legend = c(
      "Under observation",
      "Lost to event",
      "Lost to follow-up",
      "Lost to administrative censoring"
    ),
    fill = c("#BA1650", "#F5C700", "#FFD3F0", "#B7E7FC"),
    bty = "n"
  )

  # trmt group: plot stacked area chart ctrl -------------------------------------------------
  graphics::par(mar = c(5, 6, 4, 1) + .1)
  plot(
    NA,
    xlab = "t",
    ylab = "Proportion in %",
    xlim = xlim,
    ylim = c(0, 1),
    main = "Treatment group:\nProportion of Subjects by Event and Censoring Type",
    yaxt = "n"
  )

  graphics::axis(
    2,
    at = seq(1, 0, by = -0.2),
    labels = paste0(seq(100, 0, by = -20), "%"),
    las = 1
  )
  x <- seq(xlim[1], xlim[2], by = 0.001)
  probs <- get_competing_risk_probs(
    x = x,
    scale_ctrl = scale_trmt,
    shape_ctrl = shape_trmt,
    breakpoints_ctrl = breakpoints_trmt,
    scale_loss = scale_loss,
    shape_loss = shape_loss,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time
  )
  cum0 <- rep(0, length(x))
  cum1 <- probs$p_obs
  cum2 <- cum1 + probs$p_event
  cum3 <- cum2 + probs$p_loss
  cum4 <- cum3 + probs$p_admin

  # Bottom: under observation
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum1, rev(cum0)),
    col = "#BA1650",
    border = NA
  )

  # Event
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum2, rev(cum1)),
    col = "#F5C700",
    border = NA
  )

  # Loss to follow-up
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum3, rev(cum2)),
    col = "#FFD3F0",
    border = NA
  )

  # Administrative censoring (will be invisible here since p_admin = 0)
  graphics::polygon(
    x = c(x, rev(x)),
    y = c(cum4, rev(cum3)),
    col = "#B7E7FC",
    border = NA
  )

  legend(
    "topright",
    legend = c(
      "Under observation",
      "Lost to event",
      "Lost to follow-up",
      "Lost to administrative censoring"
    ),
    fill = c("#BA1650", "#F5C700", "#FFD3F0", "#B7E7FC"),
    bty = "n"
  )
  }
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
  hA <- hA_temp[1:length(t_grid)]

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
  p_obs    <- approx(t_grid, S_all,        xout = x, rule = 2)$y
  p_event  <- approx(t_grid, F_event_grid, xout = x, rule = 2)$y
  p_loss   <- approx(t_grid, F_loss_grid,  xout = x, rule = 2)$y
  p_admin  <- approx(t_grid, F_admin_grid,  xout = x, rule = 2)$y

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

utils::globalVariables(c("x")) # prevents warnings on undefined variables when running devtools::check()
