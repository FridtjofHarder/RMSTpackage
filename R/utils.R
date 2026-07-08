#' Internal helper functions
#'
#' Small utilities used only inside the package.
#'
#' @noRd

# get design parameters ---------------------------------------------------

# calculate true RMST for weibull function
get_theoretical_rmst <- function(scale, shape, breakpoints, tau) {
  stats::integrate(
    function(y) my_pew_surv(q = y, scale = scale, shape = shape, breakpoints = breakpoints),
    lower = 0,
    upper = tau
  )$value
}

# calculate true hazard
get_h <- function(x, scale, shape, breakpoints) {
  h <- my_pew_dens(x = x, scale = scale, shape = shape, breakpoints = breakpoints) /
    my_pew_surv(q = x, scale = scale, shape = shape, breakpoints = breakpoints)
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
# p (not lost to FU). Return only p(not being censored) if no loss to FU.
get_p_not_censored <- function(
  x,
  accrual_time,
  follow_up_time,
  scale_loss,
  shape_loss,
  breakpoints
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
        my_pew_surv(q = x, scale = scale_loss, shape = shape_loss, breakpoints = 0)
    )
  }
}

# calculate probability at risk and capture scale_loss = NULL
get_p_at_risk <- function(
  x,
  scale,
  shape,
  breakpoints,
  accrual_time,
  follow_up_time,
  scale_loss,
  shape_loss
) {
  return(
    my_pew_surv(q = x, scale = scale, shape = shape, breakpoints = breakpoints) *
      get_p_not_censored(
        x = x,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss,
        breakpoints = breakpoints
      )
  )
}

# get event density function f(t) * p(not censored) = -d/dt S(t) * p(not censored)
get_density <- function(
  x,
  scale,
  shape,
  breakpoints,
  accrual_time,
  follow_up_time,
  scale_loss,
  shape_loss
) {
  my_pew_dens(x = x, scale = scale, shape = shape, breakpoints = breakpoints) *
    get_p_not_censored(
      x = x,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      breakpoints = breakpoints
    )
}

# get sigma2 and delta--------------------------------------------------------------

# calculate true sigma2 for RMST
get_sigma2_rmst <- function(
  scale,
  shape,
  breakpoints,
  accrual_time,
  follow_up_time,
  tau,
  scale_loss,
  shape_loss
) {
  inner <- function(x) {
    sapply(x, function(x1) {
      stats::integrate(
        function(x2) my_pew_surv(q = x2, scale = scale, shape = shape, breakpoints = breakpoints),
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
          shape = shape,
          breakpoints = breakpoints,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time,
          scale_loss = scale_loss,
          shape_loss = shape_loss
        )
    },
    lower = 0,
    upper = tau
  )$value
}

# calculate true sigma2 for LRT
get_sigma2_LRT <- function(scale_ctrl,
                           scale_trmt,
                           shape_ctrl = 1,
                           shape_trmt = 1,
                           breakpoints_ctrl,
                           breakpoints_trmt,
                           accrual_time = 0,
                           follow_up_time = NULL,
                           tau = NULL,
                           censor_beyond_tau = FALSE,
                           scale_loss = NULL,
                           shape_loss = 1) {
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
        shape = shape_ctrl,
        breakpoints = breakpoints_ctrl,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) * get_p_at_risk(
        x,
        scale = scale_trmt,
        shape = shape_trmt,
        breakpoints = breakpoints_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) /
        (
          get_p_at_risk(
            x,
            scale = scale_ctrl,
            shape = shape_ctrl,
            breakpoints = breakpoints_ctrl,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time,
            scale_loss = scale_loss,
            shape_loss = shape_loss
          ) + get_p_at_risk(
            x,
            scale = scale_trmt,
            shape = shape_trmt,
            breakpoints = breakpoints_trmt,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time,
            scale_loss = scale_loss,
            shape_loss = shape_loss
          )
        )^2 *
        (get_density(
          x,
          scale = scale_ctrl,
          shape = shape_ctrl,
          breakpoints = breakpoints_ctrl,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time,
          scale_loss = scale_loss,
          shape_loss = shape_loss
        ) +
          get_density(
            x,
            scale = scale_trmt,
            shape = shape_trmt,
            breakpoints = breakpoints_trmt,
            accrual_time = accrual_time,
            follow_up_time = follow_up_time,
            scale_loss = scale_loss,
            shape_loss = shape_loss
          )) / 2
    }),
    lower = 0,
    upper = total_time
  )$value
  return(sigma2)
}

# new delta calculation: schoenfeld instead of asymptotic npsurvSS
get_delta_LRT <- function(scale_ctrl,
                          scale_trmt,
                          shape_ctrl = 1,
                          shape_trmt = 1,
                          breakpoints_ctrl,
                          breakpoints_trmt,
                          accrual_time = 0,
                          follow_up_time = NULL,
                          tau = NULL,
                          censor_beyond_tau = FALSE,
                          scale_loss = NULL,
                          shape_loss = 1,
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
        shape = shape_trmt,
        breakpoints = breakpoints_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) *
      get_p_at_risk(
        x,
        scale = scale_ctrl,
        shape = shape_ctrl,
        breakpoints = breakpoints_ctrl,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) /
      (get_p_at_risk(
        x,
        scale = scale_trmt,
        shape = shape_trmt,
        breakpoints = breakpoints_trmt,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) +
        get_p_at_risk(
          x,
          scale = scale_ctrl,
          shape = shape_ctrl,
          breakpoints = breakpoints_ctrl,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time,
          scale_loss = scale_loss,
          shape_loss = shape_loss
        ))^2 *
      (get_density(
        x,
        scale = scale_ctrl,
        shape = shape_ctrl,
        breakpoints = breakpoints_ctrl,
        accrual_time = accrual_time,
        follow_up_time = follow_up_time,
        scale_loss = scale_loss,
        shape_loss = shape_loss
      ) +
        get_density(
          x,
          scale = scale_trmt,
          shape = shape_trmt,
          breakpoints = breakpoints_trmt,
          accrual_time = accrual_time,
          follow_up_time = follow_up_time,
          scale_loss = scale_loss,
          shape_loss = shape_loss
        )) / 2
  }), lower = 0, upper = total_time)$value
  return(delta_LRT)
}

# get sample size closed form ---------------------------------------------

# get sample size by closed-form for RMSTD
get_ss_cf_RMSTD <- function(
  scale_ctrl,
  scale_trmt,
  shape_ctrl,
  shape_trmt,
  breakpoints_ctrl,
  breakpoints_trmt,
  accrual_time,
  follow_up_time,
  tau,
  scale_loss,
  shape_loss,
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
    shape = shape_ctrl,
    breakpoints = breakpoints_ctrl,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss
  )
  sigma2_trmt <- get_sigma2_rmst(
    scale = scale_trmt,
    shape = shape_trmt,
    breakpoints = breakpoints_trmt,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss
  )
  sigma2 <- sigma2_ctrl / 0.5 + sigma2_trmt / 0.5
  if (!is.na(satterthwaite_n)) {
    df <- get_satterthwaite_df(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      satterthwaite_n = satterthwaite_n,
      sigma2_ctrl = sigma2_ctrl,
      sigma2_trmt = sigma2_trmt
    )
    return(
      (sqrt(sigma2) *
        stats::qt(1 - alpha / sides, df) +
        sqrt(sigma2) * stats::qt(power, df))^2 /
        (RMST_trmt - RMST_ctrl - margin)^2
    )
  }
  return(
    (sqrt(sigma2) *
      stats::qnorm(1 - alpha / sides) +
      sqrt(sigma2) * stats::qnorm(power))^2 /
      (RMST_trmt - RMST_ctrl - margin)^2
  )
}

# get sample size by closed-form for RMSTR
get_ss_cf_RMSTR <- function(
  scale_ctrl,
  scale_trmt,
  shape_ctrl,
  shape_trmt,
  breakpoints_ctrl,
  breakpoints_trmt,
  accrual_time,
  follow_up_time,
  tau,
  scale_loss,
  shape_loss,
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
    shape = shape_ctrl,
    breakpoints = breakpoints_ctrl,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss
  )
  sigma2_trmt <- get_sigma2_rmst(
    scale = scale_trmt,
    shape = shape_trmt,
    breakpoints = breakpoints_trmt,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss
  )
  sigma2 <- sigma2_ctrl / .5 / RMST_ctrl^2 +
    sigma2_trmt / .5 / RMST_trmt^2
  if (!is.na(satterthwaite_n)) {
    df <- get_satterthwaite_df(
      scale_ctrl = scale_ctrl,
      scale_trmt = scale_trmt,
      shape_ctrl = shape_ctrl,
      shape_trmt = shape_trmt,
      breakpoints_ctrl = breakpoints_ctrl,
      breakpoints_trmt = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      tau = tau,
      scale_loss = scale_loss,
      shape_loss = shape_loss,
      satterthwaite_n = satterthwaite_n,
      sigma2_ctrl = sigma2_ctrl,
      sigma2_trmt = sigma2_trmt
    )
    return(
      (sqrt(sigma2) *
        stats::qt(1 - alpha / sides, df) +
        sqrt(sigma2) * stats::qt(power, df))^2 /
        (log(RMST_trmt / RMST_ctrl) - log(margin))^2
    )
  }
  return(
    (sqrt(sigma2) *
      stats::qnorm(1 - alpha / sides) +
      sqrt(sigma2) * stats::qnorm(power))^2 /
      (log(RMST_trmt / RMST_ctrl) - log(margin))^2
  )
}

# get sample size by closed-form for LRT
get_ss_cf_LRT <- function(
  scale_ctrl,
  scale_trmt,
  shape_ctrl,
  shape_trmt,
  breakpoints_ctrl,
  breakpoints_trmt,
  accrual_time,
  follow_up_time,
  tau,
  censor_beyond_tau,
  scale_loss,
  shape_loss,
  sides = 1,
  power = 0.8,
  alpha = 0.025,
  margin_LRT = 1
) {
  sigma2 <- get_sigma2_LRT(
    scale_ctrl = scale_ctrl,
    scale_trmt = scale_trmt,
    shape_ctrl = shape_ctrl,
    shape_trmt = shape_trmt,
    breakpoints_ctrl = breakpoints_ctrl,
    breakpoints_trmt = breakpoints_trmt,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    censor_beyond_tau = censor_beyond_tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss
  )
  delta <- get_delta_LRT(
    scale_ctrl = scale_ctrl,
    scale_trmt = scale_trmt,
    shape_ctrl = shape_ctrl,
    shape_trmt = shape_trmt,
    breakpoints_ctrl = breakpoints_ctrl,
    breakpoints_trmt = breakpoints_trmt,
    accrual_time = accrual_time,
    follow_up_time = follow_up_time,
    tau = tau,
    censor_beyond_tau = censor_beyond_tau,
    scale_loss = scale_loss,
    shape_loss = shape_loss,
    margin_LRT = margin_LRT
  )
  return(
    (sqrt(sigma2) *
      stats::qnorm(1 - alpha / sides) +
      sqrt(sigma2) * stats::qnorm(power))^2 /
      (delta)^2
  )
}

# misc ------------------------------------------------------------------

# reparameterize non-standard parameterizations
reparameterize <- function(parameterization, scale, shape) {
  if (is.null(scale)) {
    return(NULL)
  }
  if (parameterization == 2) {
    return(scale^shape)
  }
  if (parameterization == 3) {
    return(1 / scale)
  }
}

# get satterthwaite degrees of freedom
get_satterthwaite_df <- function(scale_ctrl, scale_trmt,
                                 shape_ctrl, shape_trmt,
                                 breakpoints_ctrl, breakpoints_trmt,
                                 accrual_time, follow_up_time, tau,
                                 scale_loss, shape_loss,
                                 satterthwaite_n, sigma2_ctrl, sigma2_trmt) {
  events_ctrl <- stats::integrate(
    Vectorize(function(x) get_density(
      x,
      scale = scale_ctrl,
      shape = shape_ctrl,
      breakpoints = breakpoints_ctrl,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss
    )),
    lower = 0, upper = tau
  )$value * satterthwaite_n
  events_trmt <- stats::integrate(
    Vectorize(function(x) get_density(
      x,
      scale = scale_trmt,
      shape = shape_trmt,
      breakpoints = breakpoints_trmt,
      accrual_time = accrual_time,
      follow_up_time = follow_up_time,
      scale_loss = scale_loss,
      shape_loss = shape_loss
    )),
    lower = 0, upper = tau
  )$value * satterthwaite_n
  df <- (sigma2_ctrl / events_ctrl + sigma2_trmt / events_trmt)^2 /
    ((sigma2_ctrl / events_ctrl)^2 / (events_ctrl - 1) + (sigma2_trmt / events_trmt)^2 / (events_trmt - 1))
  return(df)
}

# my weibulls -------------------------------------------------------------
my_pew_surv <- function(q, scale, shape, breakpoints) {
  n_intervals <- length(scale)
  shape <- rep_len(shape, n_intervals)
  if (n_intervals == 1) {
    S <- exp(- (q / scale[1])^shape[1])
  } else {
    S <- vapply(q, function(t) {
      H_t <- 0
      for (i in seq_len(n_intervals)) {
        left <- breakpoints[i]
        if (i < n_intervals) {
          right <- breakpoints[i + 1]
        } else {
          right <- Inf
        }
        if (t > left) {
          time_in_interval <- min(t, right) - left
          H_t <- H_t + (time_in_interval / scale[i])^shape[i]
        }
        if (t <= right) break
      }
      exp(-H_t)
    }, numeric(1))
  }
  return(S)
}

my_pew_dens <- function(x, scale, shape, breakpoints) {
  n_intervals <- length(scale)
  shape <- rep_len(shape, n_intervals)
  if (n_intervals == 1) {
    S_t <- exp(- (x / scale[1])^shape[1])
    h_t <- (shape[1] / scale[1]) * (x / scale[1])^(shape[1] - 1)
    f_t <- h_t * S_t
  } else {
    f_t <- vapply(x, function(t) {
      H_t <- 0
      h_t <- NA_real_
      for (i in seq_len(n_intervals)) {
        left <- breakpoints[i]
        if (i < n_intervals) {
          right <- breakpoints[i + 1]
        } else {
          right <- Inf
        }
        if (t > left) {
          time_in_interval <- min(t, right) - left
          H_t <- H_t + (time_in_interval / scale[i])^shape[i]
          if (t <= right) {
            h_t <- (shape[i] / scale[i]) * (time_in_interval / scale[i])^(shape[i] - 1)
            break
          }
        }
      }
      S_t <- exp(-H_t)
      h_t * S_t
    }, numeric(1))
  }
  return(f_t)
}

my_pew_rand <- function(n, scale, shape, breakpoints) {
  n_intervals <- length(scale)
  shape <- rep_len(shape, n_intervals)
  out <- numeric(n)
  for (j in seq_len(n)) {
    u <- stats::runif(1)
    target_H <- -log(u)
    if (n_intervals == 1) {
      out[j] <- scale[1] * target_H^(1 / shape[1])
    } else {
      H_cum <- 0
      for (i in seq_len(n_intervals)) {
        left <- breakpoints[i]
        if (i < n_intervals) {
          right <- breakpoints[i + 1]
          interval_length <- right - left
          H_next <- H_cum + (interval_length / scale[i])^shape[i]
          if (target_H <= H_next) {
            out[j] <- left + scale[i] * (target_H - H_cum)^(1 / shape[i])
            break
          } else {
            H_cum <- H_next
          }
        } else {
          out[j] <- left + scale[i] * (target_H - H_cum)^(1 / shape[i])
          break
        }
      }
    }
  }
  return(out)
}
