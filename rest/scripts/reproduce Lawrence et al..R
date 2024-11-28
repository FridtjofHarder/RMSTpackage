# reproduce Lawrence, J., Qiu, J., Bai, S., & Hung, H. J. (2019). Difference in restricted mean survival time: small sample distribution and asymptotic relative efficiency. Statistics in Biopharmaceutical Research, 11(1), 61-66.
# acc. to paper: alpha ~ 0.025 for both COX and RMST-D

if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survRM2')) install.packages('survRM2'); library('survRM2')

# functions ---------------------------------------------------------------

power_simul_noninf_cox <- function(scale_ctrl = NULL, shape_ctrl = 1, scale_trt = NULL,
                                   shape_trt = 1, scale_loss = NULL, shape_loss = 1,
                                   accrual_time = NULL, follow_up_time = NULL,
                                   tau = NULL, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                   margin_hr = 1, n = NULL, M = 100,
                                   parameterization = 1, censor_beyond_tau = FALSE){
  positive <- rep(NA, M)

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
    scale_loss <- 1/(scale_loss^(1/shape_loss))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
    scale_loss <- 1/scale_loss
  }

  n_per_group <- round(n/2)
  total_time <- accrual_time + follow_up_time

  for (i in 1:M){

    events <- c(rweibull(n_per_group, scale = scale_ctrl, shape = shape_ctrl),
                rweibull(n_per_group, scale = scale_trt, shape = shape_trt))

    loss_to_fu  <- rweibull(2*n_per_group, scale = scale_loss, shape = shape_loss)
    admin_loss  <- total_time - runif(2*n_per_group, max = accrual_time)

    status <- as.numeric(events <= loss_to_fu & events <= admin_loss)
    time <- pmin(events, loss_to_fu, admin_loss)
    arm <- c(rep(0, n_per_group), rep(1, n_per_group)) #  0 = ctrl, 1 = trt

    if(censor_beyond_tau){
      status[time>tau] <- 0
    }

    data <- data.frame(time, status, arm)

    result <-  coxph(Surv(time, status) ~ arm, data = data)
    #browser()
    upper <-  summary(result)$conf.int[, "upper .95"]

    if(upper < margin_hr){positive[i] <- 1}
    else positive[i] <- 0
  }

  print(sum(is.na(positive)))
  power <-  sum(positive, na.rm = TRUE)/sum(!is.na(positive))
  return(power)

}

power_simul_noninf_rmst <- function(scale_ctrl = NULL, shape_ctrl = 1, scale_trt = NULL,
                                    shape_trt = 1, scale_loss = NULL, shape_loss = 1,
                                    accrual_time = NULL, follow_up_time = NULL,
                                    tau = NULL, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                    margin_rmst = 0, n = NULL, M = 100,
                                    parameterization = 1, handling_large_tau = "adjust tau"){
  positive <- rep(NA, M)

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
    scale_loss <- 1/(scale_loss^(1/shape_loss))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
    scale_loss <- 1/scale_loss
  }

  n_per_group <- round(n/2)
  total_time <- accrual_time + follow_up_time
  counter_large_tau <- 0

  for (i in 1:M){

    events <- c(rweibull(n_per_group, scale = scale_ctrl, shape = shape_ctrl),
                rweibull(n_per_group, scale = scale_trt, shape = shape_trt))

    loss_to_fu  <- rweibull(2*n_per_group, scale = scale_loss, shape = shape_loss)
    admin_loss  <- total_time - runif(2*n_per_group, max = accrual_time)

    status <- as.numeric(events <= loss_to_fu & events <= admin_loss)
    time <- pmin(events, loss_to_fu, admin_loss)
    arm <- c(rep(0, n_per_group), rep(1, n_per_group)) #  0 = ctrl, 1 = trt
    data <- data.frame(time, status, arm)
    #browser()
    if(max(time[arm == 0]) <= tau || max(time[arm == 1]) <= tau){
      counter_large_tau <- counter_large_tau + 1
      if(handling_large_tau == "adjust tau") {tau <-  min(max(time[arm == 0]), max(time[arm == 1]))}
      if(handling_large_tau == "discard") next
    }

    result <-  rmst2(data$time, data$status, data$arm, tau = tau,
                     alpha = one_sided_alpha * 2)$unadjusted.result
    lower <-  result[1, 2]
    positive[i] <- if(lower > -margin_rmst){positive[i] <- 1}
    else positive[i] <- 0
  }

  print(sum(is.na(positive)))
  print(counter_large_tau)
  power <-  sum(positive, na.rm = TRUE)/sum(!is.na(positive))
  return(power)
}

ss_closed_form_noninf_rmst <- function(scale_trt, shape_trt = 1, scale_ctrl,
                                       shape_ctrl = 1, parameterization = 1,
                                       accrual_time = 0,
                                       follow_up_time = NULL,
                                       tau = NULL, sides = 1,
                                       alpha = 0.025, power = 0.8, margin = 0,
                                       plot_curves = TRUE,
                                       plot_example_data = TRUE,
                                       scale_loss = 10000, shape_loss = 1)
{
  if(is.null(tau)) stop("please specify horizon")
  if(is.null(follow_up_time)) {
    warning("follow_up_time not specified, has been set to horizon")
    follow_up_time <- tau}

  total_time <- accrual_time + follow_up_time

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
    scale_loss <- 1/(scale_loss^(1/shape_loss))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
    scale_loss <- 1/(scale_loss^(1/shape_loss))
  }

  arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                   surv_shape = shape_trt, surv_scale = 1/scale_trt,
                                   loss_scale = 1/scale_loss, loss_shape = shape_loss,
                                   follow_time = follow_up_time)
  arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                   surv_shape = shape_ctrl,surv_scale = 1/scale_ctrl,
                                   loss_scale = 1/scale_loss, loss_shape = shape_loss,
                                   follow_time = follow_up_time)

  design  <- npsurvSS:::calc_design(arm_trt, arm_ctrl, test = list(test = "rmst difference", milestone = tau)) # get Delta and sigma2 via package npsurvSS
  sample_size <- sum((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(design$tsigma2) * stats::qnorm(power) )^2 /
                       (design$delta + margin)^2 *c(0.5, 0.5))
  return(sample_size)

}

# scenario 1 ---------------------------------------------------------------

total_sample_size <-  109
sample_size_ctrl <-  ceiling(109/2)
sample_size_trt <- floor(109/2)

scale_ctrl <- 1
shape_ctrl <- 1

scale_trt <- 1
shape_trt <- 1

scale_loss <- 0.00001

accrual_time <- 0.00000001
follow_up_time <- 100

tau = 1

arm_ctrl <- create_arm(
  size = 0.7*sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = 0.7*sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = 1),
  alpha = 0.025,
  sides = 1
)

curve(pweibull(x, scale = 1/scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 5, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = 1/scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 5, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc1: survival curves for treatment and control group")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                                   scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                                   tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                   margin_hr = 1, n = total_sample_size, M = 1000,
                                   parameterization = 3, censor_beyond_tau = TRUE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                                   scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                                    tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                    margin_rmst = 0, n = total_sample_size, M = 1000,
                                    parameterization = 3, handling_large_tau = "adjust tau")

# result: false error rate = 0.025 in both tests.


# scenario 2 --------------------------------------------------------------

total_sample_size <-  109
size_per_group <- round(total_sample_size/2)

scale_ctrl <- 1
shape_ctrl <- 1

scale_trt <- 1
shape_trt <- 3

scale_loss <- 0.00001

accrual_time <- 0.00000001
follow_up_time <- 1

tau <- 1

arm_ctrl <- create_arm(
  size = size_per_group,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = size_per_group,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm( # estimates power at 0.187 instead of 0.867. WHY???
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm( # comes to expected result of power ~0.874
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = tau),
  alpha = 0.025,
  sides = 1
)

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt, shape_ctrl = shape_ctrl, shape_trt = shape_trt,
                       scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = 1, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 1, n = total_sample_size, M = 10000,
                       parameterization = 3, censor_beyond_tau = T, test = T) # estimates power similar to closed form, but very different from paper, at 0.207 vs. 0.867. WHY??

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt, shape_ctrl = shape_ctrl, shape_trt = shape_trt, # comes to expected result of power ~0.874
                        scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 0, n = total_sample_size, M = 1000,
                        parameterization = 3, handling_large_tau = "adjust tau")

curve(pweibull(x, scale = 1/scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 5, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = 1/scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 5, add = TRUE, col = "green")
legend("topright", legend= c(bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )^ .(shape_ctrl)), bquote(trt: S(t) == exp( - .(scale_trt) *t )^ .(shape_trt)))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc2: survival curves for treatment and control group")

# scenario 3 --------------------------------------------------------------

total_sample_size <- 773
sample_size_ctrl <-  ceiling(total_sample_size/2)
sample_size_trt <- floor(total_sample_size/2)

scale_ctrl <- 1
shape_ctrl <- 1

scale_trt <- 1
shape_trt <- 1

scale_loss <- 0.00001

accrual_time <- 0.00000001
follow_up_time <- 100

tau = 1

arm_ctrl <- create_arm(
  size = sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = 1),
  alpha = 0.025,
  sides = 1
)

curve(pweibull(x, scale = 1/scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 5, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = 1/scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 5, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc1: survival curves for treatment and control group")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                       scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 1, n = total_sample_size, M = 10000,
                       parameterization = 3, censor_beyond_tau = TRUE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                        scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 0, n = total_sample_size, M = 10000,
                        parameterization = 3, handling_large_tau = "adjust tau")

# result: false positive rate = 0.025 in both tests in closed form, as well as simulation .
# code identical to scenario 1 except for sample size.

# scenario 4 --------------------------------------------------------------

total_sample_size <- 773
sample_size_ctrl <-  ceiling(total_sample_size/2)
sample_size_trt <- floor(total_sample_size/2)

scale_ctrl <- 1
shape_ctrl <- 1

scale_trt <- 1
shape_trt <- 1.5

scale_loss <- 0.00001

accrual_time <- 0.00000001
follow_up_time <- 1

tau = 1

arm_ctrl <- create_arm(
  size = sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = 1),
  alpha = 0.025,
  sides = 1
)

curve(pweibull(x, scale = 1/scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 5, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = 1/scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 5, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc1: survival curves for treatment and control group")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt, shape_trt = shape_trt, shape_ctrl = shape_ctrl,
                       scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 1, n = total_sample_size, M = 10000,
                       parameterization = 3, censor_beyond_tau = TRUE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt, shape_trt = shape_trt, shape_ctrl = shape_ctrl,
                        scale_loss = scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 0, n = total_sample_size, M = 1000,
                        parameterization = 3, handling_large_tau = "adjust tau")

## result: same as in paper for RMST, but much lower power for LR / Cox regression (0.2 vs 0.8.) WHY??? similar as in sc. 2


# scenario 6 --------------------------------------------------------------

total_sample_size <- 109
sample_size_ctrl <-  ceiling(total_sample_size/2)
sample_size_trt <- floor(total_sample_size/2)

scale_ctrl <- 1 # follows R standard parameterization
shape_ctrl <- 2

scale_trt <- 1 # follows R standard parameterization
shape_trt <- 2

scale_loss <- 0.00000001

accrual_time <- 0.00000001
follow_up_time <- 2

margin_rmst <- -0.1147

tau = 1

arm_ctrl <- create_arm(
  size = sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = 1/scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = 1/scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = tau),
  alpha = 0.025,
  sides = 1
)

ss_closed_form_noninf_rmst(scale_trt = scale_trt, shape_trt = shape_trt, scale_ctrl = scale_ctrl,
                                       shape_ctrl = shape_ctrl, parameterization = 3,
                                       accrual_time = accrual_time,
                                       follow_up_time = follow_up_time,
                                       tau = tau, sides = 1,
                                       alpha = 0.025, power = 0.6, margin = -0.1147,
                                       plot_curves = TRUE,
                                       plot_example_data = TRUE,
                                       scale_loss = 0.0000001, shape_loss = 1)


curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 510, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 510, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc1: survival curves for treatment and control group")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                       scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 1, n = total_sample_size, M = 1000,
                       parameterization = 1, censor_beyond_tau = TRUE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                        scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 0, n = total_sample_size, M = 1000,
                        parameterization = 1, handling_large_tau = "adjust tau")

# result as exptected, both in closed form as well as in simulations


# scenario 10 -------------------------------------------------------------

total_sample_size <- 524
sample_size_ctrl <-  ceiling(total_sample_size/2)
sample_size_trt <- floor(total_sample_size/2)

scale_ctrl <- 8500 # follows R standard parameterization
shape_ctrl <- 1

scale_trt <- 8500 # follows R standard parameterization
shape_trt <- 1

scale_loss <- 0.00000001

accrual_time <- 510-475
follow_up_time <- 475

tau = 500

arm_ctrl <- create_arm(
  size = sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = 1/scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = 1/scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = tau),
  alpha = 0.025,
  sides = 1
)

curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 510, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 510, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc10: survival curves for treatment and control group")

print("tau = 1 censored F")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                       scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = 1*tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 8500/3686.91, n = total_sample_size, M = 1000,
                       parameterization = 1, censor_beyond_tau = FALSE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                        scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 18, n = total_sample_size, M = 10000,
                        parameterization = 1, handling_large_tau = "adjust tau")

ss_closed_form_noninf_rmst(scale_trt = scale_trt, shape_trt = 1, scale_ctrl = scale_ctrl,
                                       shape_ctrl = 1, parameterization = 1,
                                       accrual_time = accrual_time,
                                       follow_up_time = follow_up_time,
                                       tau = tau, sides = 1,
                                       alpha = 0.025, power = 0.856, margin = 18,
                                       plot_curves = TRUE,
                                       plot_example_data = TRUE,
                                       scale_loss = 1/scale_loss, shape_loss = 1)

# result: power for RMST and Cox almost identical to paper (0.8547 vs 0.856 in paper, resp. ~0.6112 vs. 0.615 for Cox).


# scenario 12 -------------------------------------------------------------

total_sample_size <- 18160
sample_size_ctrl <-  ceiling(total_sample_size/2)
sample_size_trt <- floor(total_sample_size/2)

scale_ctrl <- 8500 # follows R standard parameterization
shape_ctrl <- 1

scale_trt <- 10500 # follows R standard parameterization
shape_trt <- 1

scale_loss <- 0.00000001

accrual_time <- 510-475
follow_up_time <- 475

tau = 500

arm_ctrl <- create_arm(
  size = sample_size_ctrl,
  accr_time = accrual_time,
  surv_shape = shape_ctrl,
  surv_scale = 1/scale_ctrl,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

arm_trt <- create_arm(
  size = sample_size_trt,
  accr_time = accrual_time,
  surv_shape = shape_trt,
  surv_scale = 1/scale_trt,
  loss_scale = scale_loss,
  follow_time = follow_up_time,
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "weighted logrank"),
  alpha = 0.025,
  sides = 1
)

power_two_arm(
  arm0 = arm_ctrl,
  arm1 = arm_trt,
  test = list(test = "rmst difference", milestone = tau),
  alpha = 0.025,
  sides = 1
)

curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE), from = 0, to = 510, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = scale_trt, shape = shape_trt, lower.tail = FALSE), from = 0, to = 510, add = TRUE, col = "green")
legend("topright", legend= c(bquote(trt: S(t) == exp( - .(scale_trt) *t )), bquote(ctrl: S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "green"), lty=1, cex=0.8)
title("sc1: survival curves for treatment and control group")

power_simul_noninf_cox(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                       scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                       tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                       margin_hr = 1, n = total_sample_size, M = 1000,
                       parameterization = 1, censor_beyond_tau = TRUE)

power_simul_noninf_rmst(scale_ctrl = scale_ctrl, scale_trt = scale_trt,
                        scale_loss = 1/scale_loss, accrual_time = accrual_time, follow_up_time = follow_up_time,
                        tau = tau, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                        margin_rmst = 0, n = total_sample_size, M = 1000,
                        parameterization = 1, handling_large_tau = "adjust tau")



# result as exptected, both in closed form as well as in simulations











