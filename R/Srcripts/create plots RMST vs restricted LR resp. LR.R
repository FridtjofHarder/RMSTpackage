# creates plot with sample size for RMST difference against tau and sample size for RMST ratio against tau. Vertical line with sample size for unrestricted log rank test for reference.
# creates plot for survival curves

if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')

# functions ---------------------------------------------------------------

ss_closed_form_lr <- function(scale_trt, shape_trt = 1, scale_ctrl,
                              shape_ctrl = 1, parameterization = 1,
                              accrual_time = 0,
                              follow_up_time = NULL,
                              tau = Inf, sides = 1,
                              alpha = 0.025, power = 0.8, margin_hr = 1,
                              plot_curves = TRUE,
                              plot_example_data = TRUE,
                              scale_loss = 10000, shape_loss = 1)
{
  total_time <- accrual_time + follow_up_time
  #browser()
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

  if(tau>arm_ctrl$total_time){tau <- arm_ctrl$total_time}
  #browser()
  sigma2  <- stats::integrate(function(x) 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) *
                                0.5 * npsurvSS:::prob_risk(arm_trt, x) /
                                ( 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) + 0.5 * npsurvSS:::prob_risk(arm_trt, x) )^2 *
                                ( 0.5 * npsurvSS:::dens_event(arm_ctrl, x) + 0.5 * npsurvSS:::dens_event(arm_trt, x)),
                              lower=0,
                              upper=tau)$value
  delta <- stats::integrate(function(x) (1 / 0.5 / npsurvSS:::prob_risk(arm_ctrl, x) + 1 / 0.5 / npsurvSS:::prob_risk(arm_trt, x)) ^ (-1) *
                              ( hsurv(x, arm_trt) - hsurv(x, arm_ctrl) ),
                            lower=0,
                            upper=tau)$value
  out     <- ( sqrt(sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(sigma2) * stats::qnorm(power) )^2 /
    delta^2 *
    c(0.5, 0.5)
  return(sum(out))
}


# script ------------------------------------------------------------------

scale_ctrl <-  1.5*.0533
scale_trt <-  0.0533
accrual_time <-  .0000001 # since cannot be set to 0
follow_up_time <- 1000 # error when set to Inf
scale_loss <-  .00000000000006 # since cannot be set to 0

arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_scale = scale_trt,
                                 loss_scale = scale_loss, follow_time = follow_up_time)
arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_scale = scale_ctrl,
                                 loss_scale = scale_loss, follow_time = follow_up_time)

LR_sample_size <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt)['n']

tau_vector = seq(1, 100, 0.1)
RMST_sample_size_difference <- RMST_sample_size_ratio <- sample_size_restricted_LRT<- NULL

for(i in 1:length(tau_vector)){
  RMST_sample_size_difference[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst difference", milestone = tau_vector[i]))['n']
  RMST_sample_size_ratio[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst ratio", milestone = tau_vector[i]))['n']
  sample_size_restricted_LRT[i] <- ss_closed_form_lr(scale_trt = scale_trt, scale_ctrl = scale_ctrl, parameterization = 3, accrual_time = accrual_time,
                                                  tau = tau_vector[i], follow_up_time = follow_up_time, scale_loss = scale_loss)

}

plot(x = tau_vector, y = RMST_sample_size_difference, type = 'l', log = 'y', xlab = "tau", ylab = "sample size RMST", col = "blue")
lines(x = tau_vector, y = RMST_sample_size_ratio, col = "green")
lines(x = tau_vector, y =  sample_size_restricted_LRT, col = "yellow")
legend("topright", legend=c("sample size unrestricted LR test closed form", "sample size RMST difference against tau", "sample size RMST ratio against tau", "sample size LR test censored beyond tau"),
       col=c("red", "blue", "green", "yellow"), lty=1, cex=0.8)
abline(h = LR_sample_size, col = "red")
title("n for RMST difference and ratio against tau, and for unrestricted LR test.")

# plot survival curves
curve(pweibull(x, scale = 1/scale_ctrl, shape = 1, lower.tail = FALSE), from = 0, to = 30, col = "red", xlab = "t", ylab = "S(t)")
curve(pweibull(x, scale = 1/scale_trt, shape = 1, lower.tail = FALSE), from = 0, to = 30, add = TRUE, col = "green")
legend("topright", legend= c(bquote(S(t) == exp( - .(scale_trt) *t )), bquote(S(t) == exp( - .(scale_ctrl) *t )))
       , col=c("red", "blue", "green", "yellow"), lty=1, cex=0.8)
title("survival curves for treatment and control group")

#### closed form function for sample size calculation log rank test with censoring beyond tau

ss_closed_form_lr <- function(scale_trt, shape_trt = 1, scale_ctrl,
                              shape_ctrl = 1, parameterization = 1,
                              accrual_time = 0,
                              follow_up_time = NULL,
                              tau = Inf, sides = 1,
                              alpha = 0.025, power = 0.8, margin_hr = 1,
                              plot_curves = TRUE,
                              plot_example_data = TRUE,
                              scale_loss = 10000, shape_loss = 1)
{
  total_time <- accrual_time + follow_up_time
  #browser()
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

  if(tau>arm_ctrl$total_time){tau <- arm_ctrl$total_time}
  #browser()
  sigma2  <- stats::integrate(function(x) 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) *
                                0.5 * npsurvSS:::prob_risk(arm_trt, x) /
                                ( 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) + 0.5 * npsurvSS:::prob_risk(arm_trt, x) )^2 *
                                ( 0.5 * npsurvSS:::dens_event(arm_ctrl, x) + 0.5 * npsurvSS:::dens_event(arm_trt, x)),
                              lower=0,
                              upper=tau)$value
  delta <- stats::integrate(function(x) (1 / 0.5 / npsurvSS:::prob_risk(arm_ctrl, x) + 1 / 0.5 / npsurvSS:::prob_risk(arm_trt, x)) ^ (-1) *
                              ( hsurv(x, arm_trt) - hsurv(x, arm_ctrl) ),
                            lower=0,
                            upper=tau)$value
  out     <- ( sqrt(sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(sigma2) * stats::qnorm(power) )^2 /
    delta^2 *
    c(0.5, 0.5)
  return(sum(out))
}



