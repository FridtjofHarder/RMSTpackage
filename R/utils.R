#' Internal helper functions
#'
#' Small utilities used only inside the package.
#'
#' @noRd

# calculate true RMST for weibull function
get_theoretical_rmst <- function(tau, scale, shape){
  stats::integrate(function(y) stats::pweibull(y, shape, scale, lower.tail = FALSE),
                   lower = 0, upper = tau)$value
}

# calculate true hazard
get_h <- function(x, scale, shape){
  h <- dweibull(shape = shape, scale = scale, x = x) /
    pweibull(shape = shape, scale = scale, q = x, lower.tail = FALSE)
  return (h)
}

# calculate probability at risk
get_p_at_risk <- function(x, scale, shape, loss_scale, loss_shape, accrual_time, follow_up_time){
  return(get_p_loss_to_FU(x, follow_up_time = follow_up_time, accrual_time = accrual_time) *
           stats::pweibull(x, shape = shape, scale = scale, lower.tail = FALSE) *
           stats::pweibull(x, shape = loss_shape, scale = loss_scale, lower.tail = FALSE))
}

# calculate probability of loss to FU
get_p_loss_to_FU <- function(x, follow_up_time, accrual_time){
  if (x <= follow_up_time) return(1) else
    return (max(((follow_up_time + accrual_time - x) / accrual_time), 0))
}

# calculate true sigma2 for RMST
get_sigma2_rmst <- function(tau, scale, shape, loss_scale, loss_shape,
                            accrual_time, follow_up_time) {
  inner <- function(x) {
    sapply(x, function(x1) stats::integrate(function(x2) pweibull(x2, scale = scale, shape = shape, lower.tail = FALSE),
                                            lower=x1,
                                            upper=tau)$value
    )
  }
  stats::integrate(function(x) inner(x)^2 * get_h(x, scale = scale, shape = shape) /
                     sapply(X = x, get_p_at_risk, scale = scale, shape = shape, accrual_time = accrual_time,
                            loss_shape = loss_shape, loss_scale = loss_scale, follow_up_time = follow_up_time),
                   lower=0,
                   upper=tau)$value
}

# get sample size by closed-form margin still missing
get_ss_cf <- function(sides = 2, power = 0.8, alpha = 0.05, scale_ctrl,
                      shape_ctrl, scale_trmt, shape_trmt, tau, loss_scale, loss_shape,
                      follow_up_time, accrual_time, margin = 0){
  delta <- get_theoretical_rmst(tau = tau, scale = scale_trmt, shape = shape_trmt) -
    get_theoretical_rmst(tau = tau, scale = scale_ctrl, shape = shape_ctrl)
  sigma2_ctrl <- get_sigma2_rmst(tau = tau, scale = scale_ctrl, shape = shape_ctrl, # my method
                                 loss_scale = loss_scale, loss_shape = loss_shape,
                                 follow_up_time = follow_up_time, accrual_time = accrual_time)
  sigma2_trmt <- get_sigma2_rmst(tau = tau, scale = scale_trmt, shape = shape_trmt, # my method
                                 loss_scale = loss_scale, loss_shape = loss_shape,
                                 follow_up_time = follow_up_time, accrual_time = accrual_time)
  sigma2 <- sigma2_ctrl / 0.5 + sigma2_trmt / 0.5
  return ((sqrt(sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(sigma2) * stats::qnorm(power))^2 /
            (delta-margin)^2)
}

# just some tests to compare my functions with npsurvSS DELETE LATER!!!!!!!!!!!

# define params ----------------------------------------------------------------
scale_trmt <- 1.5
scale_ctrl <- 1
shape_trmt <- 1.5
shape_ctrl <- 1.5
follow_up_time <- 1
accrual_time <- 1
loss_scale <- 1.5
loss_shape <- 1.5
tau <-  teval <-  1.9
sides <- 2
alpha <- 0.05
power <- 0.8

# define arms for npsurvSS -----------------------------------------------------
arm_trmt <- npsurvSS::create_arm(size=1,
                                 accr_time = accrual_time,
                                 follow_time = follow_up_time,
                                 surv_scale = 1 / scale_trmt,
                                 surv_shape = shape_trmt,
                                 loss_scale = 1 / loss_scale,
                                 loss_shape = loss_shape
)

arm_ctrl <- npsurvSS::create_arm(size=1,
                      accr_time = accrual_time,
                      follow_time = follow_up_time,
                      surv_scale = 1 / scale_ctrl,
                      surv_shape = shape_ctrl,
                      loss_scale = 1 / loss_scale,
                      loss_shape = loss_shape
                      )

# get sigma2 -------------------------------------------------------------------
# compare sigma2

(sigma2_npsurvSS_ctrl <- npsurvSS:::sigma2j_rmst(arm = arm_ctrl, teval = tau))
(sigma2_npsurvSS_trmt <- npsurvSS:::sigma2j_rmst(arm = arm_trmt, teval = tau))
(sigma2_npsurvSS <-  sigma2_npsurvSS_ctrl / 0.5 + sigma2_npsurvSS_trmt / 0.5)

(sigma2_ctrl <- get_sigma2_rmst(tau = tau, scale = scale_ctrl, shape = shape_ctrl, # my method
                                loss_scale = loss_scale, loss_shape = loss_shape,
                                follow_up_time = follow_up_time, accrual_time = accrual_time))

(sigma2_trmt <- get_sigma2_rmst(tau = tau, scale = scale_trmt, shape = shape_trmt, # my method
                                loss_scale = loss_scale, loss_shape = loss_shape,
                                follow_up_time = follow_up_time, accrual_time = accrual_time))

(sigma2  <- sigma2_ctrl / 0.5 + sigma2_trmt / 0.5)

# get n ------------------------------------------------------------------------
npsurvSS::size_two_arm(arm_ctrl, arm_trmt, list(test="rmst difference", milestone=tau),
             sides = 2, alpha = .05) # npsurvSS native

theoretical_delta <- get_theoretical_rmst(tau = tau, scale = scale_trmt, shape = shape_trmt) -
  get_theoretical_rmst(tau = tau, scale = scale_ctrl, shape = shape_ctrl) # get Delta from indivual RMSTs

get_ss_cf(sides = sides, scale_ctrl = scale_ctrl,
          shape_ctrl = shape_ctrl, scale_trmt = scale_trmt, shape_trmt = shape_trmt,
          tau = tau, loss_scale = loss_scale, loss_shape = loss_shape,
          follow_up_time = follow_up_time, accrual_time = accrual_time)

sides <- 2

(sqrt(sigma2) * stats::qnorm(1 - .05 / 2) + sqrt(sigma2) * stats::qnorm(.8))^2 /
  theoretical_delta^2


