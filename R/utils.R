#' Internal helper functions
#'
#' Small utilities used only inside the package.
#'
#' @noRd

# calculate true RMST for weibull function
get_theoretical_rmst <- function(tau, scale, shape) {
  stats::integrate(function(y) {
    stats::pweibull(y, shape, scale, lower.tail = FALSE)
  }, lower = 0, upper = tau)$value
}

# calculate true hazard
get_h <- function(x, scale, shape) {
  h <- stats::dweibull(shape = shape, scale = scale, x = x) /
    stats::pweibull(
      shape = shape,
      scale = scale,
      q = x,
      lower.tail = FALSE
    )
  return(h)
}

# calculate probability of not lost to administrative censoring
get_p_not_lost_admin <- function(x, follow_up_time, accrual_time) {
  if (x <= follow_up_time || follow_up_time == Inf) {
    return(1)
  } else {
    return(max(((follow_up_time + accrual_time - x) / accrual_time
    ), 0))
  }
}

get_p_not_lost_admin_v <- Vectorize(get_p_not_lost_admin)

# calculate p(not being censored) as product of p(not lost to admin. censoring) *
# p (not lost to FU). Return only p(not being censored) if no loss to FU.
get_p_not_censored <- function(x,
                               follow_up_time,
                               accrual_time,
                               loss_scale,
                               loss_shape) {
  if (is.null(loss_scale) || is.null(loss_shape)) {
    return(
      get_p_not_lost_admin(
        x = x,
        follow_up_time = follow_up_time,
        accrual_time = accrual_time
      )
    )
  } else
    return (
      get_p_not_lost_admin(
        x = x,
        follow_up_time = follow_up_time,
        accrual_time = accrual_time
      ) *    stats::pweibull(
        q = x,
        shape = loss_shape,
        scale = loss_scale,
        lower.tail = FALSE
      )
    )
}

# calculate probability at risk and capture loss_scale = NULL
get_p_at_risk <- function(x,
                          scale,
                          shape,
                          loss_scale,
                          loss_shape,
                          accrual_time,
                          follow_up_time) {
  return(
    stats::pweibull(
      q = x,
      # 1st term: survival
      shape = shape,
      scale = scale,
      lower.tail = FALSE
    ) *
      get_p_not_censored(
        x = x,
        follow_up_time = follow_up_time,
        # 2nd term: loss to censoring by any cause
        accrual_time = accrual_time,
        loss_scale = loss_scale,
        loss_shape = loss_shape
      )
  )
}

# calculate true sigma2 for RMST
get_sigma2_rmst <- function(tau,
                            scale,
                            shape,
                            loss_scale,
                            loss_shape,
                            accrual_time,
                            follow_up_time) {
  inner <- function(x) {
    sapply(x, function(x1) {
      stats::integrate(function(x2) {
        stats::pweibull(x2,
                        scale = scale,
                        shape = shape,
                        lower.tail = FALSE)
      }, lower = x1, upper = tau)$value
    })
  }
  stats::integrate(function(x) {
    inner(x)^2 * get_h(x, scale = scale, shape = shape) /
      sapply(
        X = x,
        get_p_at_risk,
        scale = scale,
        shape = shape,
        accrual_time = accrual_time,
        loss_shape = loss_shape,
        loss_scale = loss_scale,
        follow_up_time = follow_up_time
      )
  }, lower = 0, upper = tau)$value
}

# get sample size by closed-form margin still missing
get_ss_cf <- function(sides = 2,
                      power = 0.8,
                      alpha = 0.05,
                      scale_ctrl,
                      shape_ctrl,
                      scale_trmt,
                      shape_trmt,
                      tau,
                      loss_scale,
                      loss_shape,
                      follow_up_time,
                      accrual_time,
                      margin = 0) {
  delta <- get_theoretical_rmst(tau = tau,
                                scale = scale_trmt,
                                shape = shape_trmt) -
    get_theoretical_rmst(tau = tau,
                         scale = scale_ctrl,
                         shape = shape_ctrl)
  sigma2_ctrl <- get_sigma2_rmst(
    tau = tau,
    scale = scale_ctrl,
    shape = shape_ctrl,
    # my method
    loss_scale = loss_scale,
    loss_shape = loss_shape,
    follow_up_time = follow_up_time,
    accrual_time = accrual_time
  )
  sigma2_trmt <- get_sigma2_rmst(
    tau = tau,
    scale = scale_trmt,
    shape = shape_trmt,
    # my method
    loss_scale = loss_scale,
    loss_shape = loss_shape,
    follow_up_time = follow_up_time,
    accrual_time = accrual_time
  )
  sigma2 <- sigma2_ctrl / 0.5 + sigma2_trmt / 0.5
  return((
    sqrt(sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(sigma2) * stats::qnorm(power)
  )^2 /
    (delta - margin)^2)
}

# RMST over tau
RMST_over_tau <- function(tau, shape = 1, scale) {
  sapply(tau, function(tau) {
    stats::integrate(stats::pweibull,
                     shape = shape,
                     scale = scale, lower = 0,
                     upper = tau,
                     lower.tail = FALSE
    )$value
  })
}

# RMSTD over tau
RMSTD_over_tau <- function(tau, shape_trmt = 1, scale_trmt, shape_ctrl = 1, scale_ctrl) {
  sapply(tau, function(tau) {
    RMST_over_tau(tau = tau, shape = shape_trmt, scale = scale_trmt) -
      RMST_over_tau(tau = tau, shape = shape_ctrl, scale = scale_ctrl)
  })
}
