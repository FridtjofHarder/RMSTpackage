#' Internal helper functions
#'
#' Small utilities used only inside the package.
#'
#' @noRd

# get design parameters ---------------------------------------------------

# calculate true RMST for Weibull function
get_theoretical_rmst <- function(scale, shape, breakpoints, tau) {
  stats::integrate(
    function(y) ppweibull::ppweibull(q = y, rate = 1 / scale^shape, alpha = shape, t = breakpoints, lower.tail = FALSE),
    lower = 0,
    upper = tau
  )$value
}

# calculate true hazard
get_h <- function(x, scale, shape, breakpoints) {
  h <- ppweibull::dpweibull(x = x, alpha = shape, rate = 1 / scale^shape, t = breakpoints) /
    ppweibull::ppweibull(q = x, rate = 1 / scale^shape, alpha = shape, t = breakpoints, lower.tail = FALSE)
  return(h)
}

# get probability distributions -------------------------------------------

# calculate probability of not lost to administrative censoring
get_p_not_lost_admin <- function(x, accrual_time, follow_up_time) {
  if (x <= follow_up_time || follow_up_time == Inf) {
    return(1)
  } else {
    return(max(((follow_up_time + accrual_time - x) / accrual_time), 0))
  }
}

# calculate p(not being censored) as product of p(not lost to admin. censoring) *
# p(not lost to follow-up). Return only p(not being censored) if no loss to follow-up.
get_p_not_censored <- function(
  x,
  accrual_time,
  follow_up_time,
  scale_loss,
  shape_loss,
  breakpoints_loss
) {
  if (is.null(scale_loss) || is.null(shape_loss)) {
    return(
      get_p_not_lost_admin(
        x = x,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      )
    )
  } else {
    return(
      get_p_not_lost_admin(
        x = x,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) *
        ppweibull::ppweibull(q = x, rate = 1 / scale_loss^shape_loss, alpha = shape_loss, t = breakpoints_loss, lower.tail = FALSE)
    )
  }
}

# calculate probability at risk and capture scale_loss = NULL
get_p_at_risk <- function(
  x,
  scale,
  scale_loss,
  shape,
  shape_loss,
  breakpoints,
  breakpoints_loss,
  accrual_time,
  follow_up_time
) {
  return(
    ppweibull::ppweibull(q = x, rate = 1 / scale^shape, alpha = shape, t = breakpoints, lower.tail = FALSE) *
      get_p_not_censored(
        x = x,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss,
        breakpoints_loss = breakpoints_loss
      )
  )
}

# get event density function f(t) * p(not censored) = -d/dt S(t) * p(not censored)
get_density <- function(
  x,
  scale,
  scale_loss,
  shape,
  shape_loss,
  breakpoints,
  breakpoints_loss,
  accrual_time,
  follow_up_time
) {
  ppweibull::dpweibull(x = x, alpha = shape, rate = 1 / scale^shape, t = breakpoints) *
    get_p_not_censored(
      x = x,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints_loss = breakpoints_loss
    )
}

# get sigma2 and delta--------------------------------------------------------------

# calculate true sigma2 for RMST
get_sigma2_rmst <- function(
  scale,
  scale_loss,
  shape,
  shape_loss,
  breakpoints,
  breakpoints_loss,
  accrual_time,
  follow_up_time,
  tau
) {
  inner <- function(x) {
    sapply(x, function(x1) {
      stats::integrate(
        function(x2) ppweibull::ppweibull(q = x2, rate = 1 / scale^shape, alpha = shape, t = breakpoints, lower.tail = FALSE),
        lower = x1,
        upper = tau
      )$value
    })
  }
  stats::integrate(
    function(x) {
      inner(x)^2 *
        get_h(x, scale = scale, shape = shape, breakpoints = breakpoints) /
        sapply(
          X = x,
          get_p_at_risk,
          scale = scale,
          scale_loss = scale_loss,
          shape = shape,
          shape_loss = shape_loss,
          breakpoints = breakpoints,
          breakpoints_loss = breakpoints_loss,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        )
    },
    lower = 0,
    upper = tau
  )$value
}

# calculate true sigma2 for LRT
get_sigma2_LRT <- function(scale_ctrl,
                           scale_trmt,
                           scale_loss = NULL,
                           shape_ctrl = 1,
                           shape_trmt = 1,
                           shape_loss = 1,
                           breakpoints_ctrl,
                           breakpoints_trmt,
                           breakpoints_loss,
                           accrual_time = 0,
                           follow_up_time = NULL,
                           tau = NULL,
                           censor_beyond_tau = FALSE) {
  if (censor_beyond_tau) {
    total_time <- tau
  } else {
    total_time <- accrual_time + follow_up_time
  }
  sigma2 <- stats::integrate(
    Vectorize(function(x) {
      get_p_at_risk(
        x,
        scale = scale_ctrl,
        scale_loss = scale_loss,
        shape = shape_ctrl,
        shape_loss = shape_loss,
        breakpoints = breakpoints_ctrl,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) * get_p_at_risk(
        x,
        scale = scale_trmt,
        scale_loss = scale_loss,
        shape = shape_trmt,
        shape_loss = shape_loss,
        breakpoints = breakpoints_trmt,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) /
        (
          get_p_at_risk(
            x,
            scale = scale_ctrl,
            scale_loss = scale_loss,
            shape = shape_ctrl,
            shape_loss = shape_loss,
            breakpoints = breakpoints_ctrl,
            breakpoints_loss = breakpoints_loss,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time
          ) + get_p_at_risk(
            x,
            scale = scale_trmt,
            scale_loss = scale_loss,
            shape = shape_trmt,
            shape_loss = shape_loss,
            breakpoints = breakpoints_trmt,
            breakpoints_loss = breakpoints_loss,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time
          )
        )^2 *
        (get_density(
          x,
          scale = scale_ctrl,
          scale_loss = scale_loss,
          shape = shape_ctrl,
          shape_loss = shape_loss,
          breakpoints = breakpoints_ctrl,
          breakpoints_loss = breakpoints_loss,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        ) +
          get_density(
            x,
            scale = scale_trmt,
            scale_loss = scale_loss,
            shape = shape_trmt,
            shape_loss = shape_loss,
            breakpoints = breakpoints_trmt,
            breakpoints_loss = breakpoints_loss,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time
          )) / 2
    }),
    lower = 0,
    upper = total_time
  )$value
  return(sigma2)
}

# new delta calculation: Schoenfeld instead of asymptotic npsurvSS
get_delta_LRT <- function(scale_ctrl,
                          scale_trmt,
                          scale_loss = NULL,
                          shape_ctrl = 1,
                          shape_trmt = 1,
                          shape_loss = 1,
                          breakpoints_ctrl,
                          breakpoints_trmt,
                          breakpoints_loss = breakpoints_loss,
                          accrual_time = 0,
                          follow_up_time = NULL,
                          tau = NULL,
                          censor_beyond_tau = FALSE,
                          margin_LRT = 1) {
  if (censor_beyond_tau) {
    total_time <- tau
  } else {
    total_time <- accrual_time + follow_up_time
  }
  delta_LRT <- stats::integrate(Vectorize(function(x) {
    (log(get_h(x, scale = scale_trmt, shape = shape_trmt, breakpoints = breakpoints_trmt)) -
      log(get_h(x, scale = scale_ctrl, shape = shape_ctrl, breakpoints = breakpoints_ctrl)) -
      log(margin_LRT)) *
      get_p_at_risk(
        x,
        scale = scale_trmt,
        scale_loss = scale_loss,
        shape = shape_trmt,
        shape_loss = shape_loss,
        breakpoints = breakpoints_trmt,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) *
      get_p_at_risk(
        x,
        scale = scale_ctrl,
        scale_loss = scale_loss,
        shape = shape_ctrl,
        shape_loss = shape_loss,
        breakpoints = breakpoints_ctrl,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) /
      (get_p_at_risk(
        x,
        scale = scale_trmt,
        scale_loss = scale_loss,
        shape = shape_trmt,
        shape_loss = shape_loss,
        breakpoints = breakpoints_trmt,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) +
        get_p_at_risk(
          x,
          scale = scale_ctrl,
          scale_loss = scale_loss,
          shape = shape_ctrl,
          shape_loss = shape_loss,
          breakpoints = breakpoints_ctrl,
          breakpoints_loss = breakpoints_loss,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        ))^2 *
      (get_density(
        x,
        scale = scale_ctrl,
        scale_loss = scale_loss,
        shape = shape_ctrl,
        shape_loss = shape_loss,
        breakpoints = breakpoints_ctrl,
        breakpoints_loss = breakpoints_loss,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time
      ) +
        get_density(
          x,
          scale = scale_trmt,
          scale_loss = scale_loss,
          shape = shape_trmt,
          shape_loss = shape_loss,
          breakpoints = breakpoints_trmt,
          breakpoints_loss = breakpoints_loss,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time
        )) / 2
  }), lower = 0, upper = total_time)$value
  return(delta_LRT)
}

# get sample size closed form ---------------------------------------------

# get sample size by closed-form for RMSTD
get_ss_cf_RMSTD <- function(
  scale_ctrl,
  scale_trmt,
  scale_loss,
  shape_ctrl,
  shape_trmt,
  shape_loss,
  breakpoints_ctrl,
  breakpoints_trmt,
  breakpoints_loss,
  accrual_time,
  follow_up_time,
  tau,
  sides = 2,
  power = 0.8,
  alpha = 0.05,
  margin = 0,
  satterthwaite_n = NA,
  RMST_ctrl,
  RMST_trmt
) {
  sigma2_ctrl <- get_sigma2_rmst(
    scale = scale_ctrl,
    scale_loss = scale_loss,
    shape = shape_ctrl,
    shape_loss = shape_loss,
    breakpoints = breakpoints_ctrl,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau
  )
  sigma2_trmt <- get_sigma2_rmst(
    scale = scale_trmt,
    scale_loss = scale_loss,
    shape = shape_trmt,
    shape_loss = shape_loss,
    breakpoints = breakpoints_trmt,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau
  )
  sigma2 <- sigma2_ctrl / 0.5 + sigma2_trmt / 0.5
  if (!is.na(satterthwaite_n)) {
    df <- get_satterthwaite_df(
      scale_ctrl = scale_ctrl,
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
      tau = tau,
      satterthwaite_n = satterthwaite_n,
      sigma2_ctrl = sigma2_ctrl,
      sigma2_trmt = sigma2_trmt
    )
    return(
      sigma2 * (stats::qt(1 - alpha / sides, df) + stats::qt(power, df))^2 /
        (RMST_trmt - RMST_ctrl - margin)^2
    )
  }
  return(
    sigma2 * (stats::qnorm(1 - alpha / sides) + stats::qnorm(power))^2 /
      (RMST_trmt - RMST_ctrl - margin)^2
  )
}

# get sample size by closed-form for RMSTR
get_ss_cf_RMSTR <- function(
  scale_ctrl,
  scale_trmt,
  scale_loss,
  shape_ctrl,
  shape_trmt,
  shape_loss,
  breakpoints_ctrl,
  breakpoints_trmt,
  breakpoints_loss,
  accrual_time,
  follow_up_time,
  tau,
  sides = 2,
  power = 0.8,
  alpha = 0.05,
  margin = 1,
  satterthwaite_n = NA,
  RMST_ctrl,
  RMST_trmt
) {
  sigma2_ctrl <- get_sigma2_rmst(
    scale = scale_ctrl,
    scale_loss = scale_loss,
    shape = shape_ctrl,
    shape_loss = shape_loss,
    breakpoints = breakpoints_ctrl,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau
  )
  sigma2_trmt <- get_sigma2_rmst(
    scale = scale_trmt,
    scale_loss = scale_loss,
    shape = shape_trmt,
    shape_loss = shape_loss,
    breakpoints = breakpoints_trmt,
    breakpoints_loss = breakpoints_loss,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau
  )
  sigma2 <- sigma2_ctrl / .5 / RMST_ctrl^2 +
    sigma2_trmt / .5 / RMST_trmt^2
  if (!is.na(satterthwaite_n)) {
    df <- get_satterthwaite_df(
      scale_ctrl = scale_ctrl,
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
      tau = tau,
      satterthwaite_n = satterthwaite_n,
      sigma2_ctrl = sigma2_ctrl,
      sigma2_trmt = sigma2_trmt
    )
    return(
      sigma2 * (stats::qt(1 - alpha / sides, df) + stats::qt(power, df))^2 /
        (log(RMST_trmt / RMST_ctrl) - log(margin))^2
    )
  }
  return(
    sigma2 * (stats::qnorm(1 - alpha / sides) + stats::qnorm(power))^2 /
      (log(RMST_trmt / RMST_ctrl) - log(margin))^2
  )
}

# get sample size by closed-form for LRT
get_ss_cf_LRT <- function(
  scale_ctrl,
  scale_trmt,
  scale_loss,
  shape_ctrl,
  shape_trmt,
  shape_loss,
  breakpoints_ctrl,
  breakpoints_trmt,
  breakpoints_loss,
  accrual_time,
  follow_up_time,
  tau,
  censor_beyond_tau,
  sides = 1,
  power = 0.8,
  alpha = 0.025,
  margin_LRT = 1
) {
  sigma2 <- get_sigma2_LRT(
    scale_ctrl = scale_ctrl,
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
    tau = tau,
    censor_beyond_tau = censor_beyond_tau
  )
  delta <- get_delta_LRT(
    scale_ctrl = scale_ctrl,
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
    tau = tau,
    censor_beyond_tau = censor_beyond_tau,
    margin_LRT = margin_LRT
  )
  return(
    sigma2 * (stats::qnorm(1 - alpha / sides) + stats::qnorm(power))^2 /
      delta^2
  )
}

# misc ------------------------------------------------------------------

# prepend 0 to a breakpoints vector if not already present; return 0 for NULL
normalize_breakpoints <- function(x) {
  if (is.null(x)) return(0)
  if (x[1] != 0) c(0, x) else x
}

reparameterize <- function(parameterisation, scale, shape) {
  if (is.null(scale)) {
    return(NULL)
  }
  if (parameterisation == 2) {
    return(scale^shape)
  }
  if (parameterisation == 3) {
    return(1 / scale)
  }
}

get_satterthwaite_df <- function(scale_ctrl, scale_trmt,
                                 scale_loss, shape_ctrl, shape_trmt,
                                 shape_loss, breakpoints_ctrl, breakpoints_trmt,
                                 breakpoints_loss,
                                 accrual_time, follow_up_time, tau,
                                 satterthwaite_n, sigma2_ctrl, sigma2_trmt) {
  events_ctrl <- stats::integrate(
    Vectorize(function(x) get_density(
      x,
      scale = scale_ctrl,
      scale_loss = scale_loss,
      shape = shape_ctrl,
      shape_loss = shape_loss,
      breakpoints = breakpoints_ctrl,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time
    )),
    lower = 0, upper = tau
  )$value * satterthwaite_n
  events_trmt <- stats::integrate(
    Vectorize(function(x) get_density(
      x,
      scale = scale_trmt,
      scale_loss = scale_loss,
      shape = shape_trmt,
      shape_loss = shape_loss,
      breakpoints = breakpoints_trmt,
      breakpoints_loss = breakpoints_loss,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time
    )),
    lower = 0, upper = tau
  )$value * satterthwaite_n
  df <- (sigma2_ctrl / events_ctrl + sigma2_trmt / events_trmt)^2 /
    ((sigma2_ctrl / events_ctrl)^2 / (events_ctrl - 1) + (sigma2_trmt / events_trmt)^2 / (events_trmt - 1))
  return(df)
}

# check inputs
check_inputs <- function(scale_ctrl = NULL,
                         scale_trmt = NULL,
                         scale_loss = NULL,
                         shape_ctrl = NULL,
                         shape_trmt = NULL,
                         shape_loss = NULL,
                         breakpoints_ctrl = NULL,
                         breakpoints_trmt = NULL,
                         breakpoints_loss = NULL,
                         follow_up_time = NULL,
                         tau = NULL,
                         sides = NULL,
                         power = NULL,
                         one_sided_alpha = NULL,
                         RMSTD_closed_form = FALSE,
                         RMSTR_closed_form = FALSE,
                         parameterisation = NULL){
  if (is.null(scale_ctrl) || is.null(scale_trmt)) {
    stop(
      "Please specify scale parameters for both treatment and survival group."
    )
  }
  if(length(scale_ctrl) != length(shape_ctrl) || length(scale_trmt) != length(shape_trmt) || length(scale_loss) != length(shape_loss)){
    stop("Scale and shape parameter must have same length in each group")
  }
  stopifnot("first element in breakpoint vectors must be larger than 0" =
              (is.null(breakpoints_ctrl[1]) ||  breakpoints_ctrl[1] > 0) &&
              (is.null(breakpoints_trmt[1]) ||  breakpoints_trmt[1] > 0) &&
              (is.null(breakpoints_loss[1]) ||  breakpoints_loss[1] > 0))
  stopifnot("breakpoints must be in increasing order" =
              is.null(breakpoints_ctrl[1]) || all(diff(breakpoints_ctrl) > 0) &&
              is.null(breakpoints_trmt[1]) || all(diff(breakpoints_trmt) > 0) &&
              (is.null(breakpoints_loss[1]) ||  breakpoints_loss[1] > 0))
  if (follow_up_time == Inf) {
    warning("follow_up_time not specified, no administrative censoring will be applied.")
  }
  if (RMSTD_closed_form || RMSTR_closed_form) {
    stopifnot("Please specify valid time horizon tau > 0." = tau > 0 && !is.null(tau))
  }
  stopifnot(
    "Parameterization must be defined as either 1, 2, or 3." = parameterisation == 1 || parameterisation == 2 || parameterisation == 3
  )
  stopifnot("sides must be set to either 1 or 2" = any(sides == c(1, 2)))
  stopifnot("one_sides_alpha must be larger than 0 and smaller than 1" = (0 < one_sided_alpha & one_sided_alpha < 1))
  stopifnot("power must be larger than 0 and smaller than 1" = (0 < power & power < 1))
}
