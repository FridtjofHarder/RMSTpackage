scale_ctrl <-  1.5*.0533
scale_trt <-  0.0533
shape_ctrl <- 1
shape_trt <- 1
tau <-  14
accrual_time <-  .0000001#17
follow_up_time <- 1000
shape_loss <- 1
scale_loss <-  .00000000000006
margin_rmst <- 0
convert_contrast(scale_trt = NULL,
                             shape_trt = 1,
                             scale_ctrl = scale_trt,
                             shape_ctrl = 1,
                             parameterization = 1,
                             HR = NULL,
                             median_diff = NULL,
                             percentile_diff = NULL,
                             percentile = NULL,
                             survival_diff = NULL,
                             t = NULL,
                             RMSTD = margin_rmst,
                             RMSTR = NULL,
                             output = "HR",
                             tau = tau,
                             plot_curves = TRUE)
M <- 100000

ss_closed_form_cox(scale_trt = scale_trt, shape_trt = 1, scale_ctrl = scale_ctrl,
                    shape_ctrl = 1, parameterization = 1,
                    accrual_time = accrual_time,
                    follow_up_time = follow_up_time,
                    tau = tau, sides = 1,
                    alpha = 0.025, power = 0.8, margin_hr = 1,
                    plot_curves = TRUE,
                    plot_example_data = TRUE,
                    scale_loss = scale_loss, shape_loss = 1)

# already known: rmsamplesize() is equivalent to rmsamplesizeequv(rmstDiffLower = 0. rmstDiffUpper = Inf)
# already known: my closed form, simulation and rmsamplesizeequv() are equivalent for noninf, but not for admin censoring, probably because
# rmsamplesizequiv assumes poisson accrual
# problem: simulation stops when last obersvation too small

rmsamplesizeequiv(beta = 0.2, alpha = 0.025, milestone = tau,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             lambda1 = scale_trt,
             lambda2 = scale_ctrl,
             gamma1 = scale_loss,
             gamma2 = scale_loss,
             accrualDuration = accrual_time,
             followupTime = follow_up_time, fixedFollowup = FALSE,
             rmstDiffLower = -margin,
             rmstDiffUpper = 1000)

n_closed_form <- ss_closed_form_noninf_rmst(scale_trt = scale_trt, scale_ctrl = scale_ctrl, accrual_time = accrual_time,
                           follow_up_time = follow_up_time, tau = tau, margin = margin, shape_loss = shape_loss, scale_loss = scale_loss,
                           parameterization = 3)
n_closed_form

power_simul_noninf_rmst(scale_trt = scale_trt, scale_ctrl = scale_ctrl, accrual_time = accrual_time,
                           follow_up_time = follow_up_time, tau = tau, margin = margin, shape_loss = shape_loss, scale_loss = scale_loss,
                           parameterization = 3, M = M, n = n_closed_form, handling_large_tau = "discard")

arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_trt, surv_scale = scale_trt,
                                 loss_scale = scale_loss, follow_time = follow_up_time)
arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_ctrl,surv_scale = scale_ctrl,
                                 loss_scale = scale_loss, follow_time = follow_up_time)
size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, test = list(test="rmst difference", milestone = tau))

# cox/LR

lrsamplesizeequiv(beta = 0.2, alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             lambda2 = scale_trt,
             lambda1 = scale_ctrl,
             accrualDuration = accrual_time,
             followupTime = follow_up_time, fixedFollowup = FALSE,
             hazardRatioLower = 0.00001, hazardRatioUpper = 1)

lrsamplesizeequiv(alpha = 0.05, hazardRatioLower = 0.4, hazardRatioUpper = 1.8,
                  allocationRatioPlanned = 1, accrualTime = seq(0, 8),
                  accrualIntensity = 26/9*seq(1, 9),
                  lambda1 = 0.06,
                  lambda2 = 0.05,
                  accrualDuration = NA,
                  followupTime = 300, fixedFollowup = T) ### not working, gives error constantly!!!

power_simul_noninf_cox(scale_trt = scale_trt, scale_ctrl = scale_ctrl, accrual_time = accrual_time,
                        follow_up_time = follow_up_time, tau = tau, margin_hr = 1, shape_loss = shape_loss, scale_loss = scale_loss,
                        parameterization = 3, M = M, n = 369, censor_beyond_tau = FALSE)

arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_trt, surv_scale = scale_trt,
                                 loss_scale = scale_loss, follow_time = follow_up_time)
arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_ctrl,surv_scale = scale_ctrl,
                                 loss_scale = scale_loss, follow_time = follow_up_time)

size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt)

#### test function npsurvSS hacked

scale_test <- 1
loss_scale_test <- 0.0000001
accrual_time_test <- 1
follow_up_time_test <- 1
arm0 <- arm_ctrl
arm1 <- arm_trt
p0 <- p1 <- 0.5

arm_test  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time_test,
                                 surv_scale = scale_test,
                                 loss_scale = loss_scale_test, follow_time = follow_up_time_test)

sigma2_1  <- stats::integrate(function(x)
                              p0 * npsurvSS:::prob_risk(arm0, x) * p1 * npsurvSS:::prob_risk(arm1, x) /
                              ( p0 * npsurvSS:::prob_risk(arm0, x) + p1 * npsurvSS:::prob_risk(arm1, x) )^2 *
                              ( p0 * npsurvSS:::dens_event(arm0, x) + p1 * npsurvSS:::dens_event(arm1, x)),
                            lower=0,
                            upper=arm0$total_time)$value

delta <- stats::integrate(function(x) (1 / 0.5 / npsurvSS:::prob_risk(arm0, x) + 1 / 0.5 / npsurvSS:::prob_risk(arm1, x)) ^ (-1) *
                            ( hsurv(x, arm1) - hsurv(x, arm0) ),
                          lower=0,
                          upper=arm0$total_time)$value

ss_closed_form_cox(scale_trt, shape_trt = 1, scale_ctrl,
                   shape_ctrl = 1, parameterization = 3,
                   accrual_time = accrual_time,
                   follow_up_time = follow_up_time,
                   tau = tau, sides = 1,
                   alpha = 0.025, power = 0.8, margin_hr = 1,
                   plot_curves = TRUE,
                   plot_example_data = TRUE,
                   scale_loss = scale_loss, shape_loss = 1)

size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt)
arm_ctrl$total_time

sapply(seq(1, 20), ss_closed_form_cox, scale_trt = scale_trt, shape_trt = 1, scale_ctrl = scale_ctrl,
       shape_ctrl = 1, parameterization = 3, accrual_time = accrual_time,
       follow_up_time = follow_up_time,
       sides = 1,
       alpha = 0.025, power = 0.8, margin_hr = 1,
       plot_curves = TRUE,
       plot_example_data = TRUE,
       scale_loss = scale_loss, shape_loss = 1)

power_simul_noninf_cox(scale_ctrl = scale_ctrl, shape_ctrl = 1, scale_trt = scale_trt,
                                   shape_trt = 1, scale_loss = scale_loss, shape_loss = 1,
                                   accrual_time = accrual_time, follow_up_time = follow_up_time,
                                   total_time = NULL,
                                   tau = 2, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                   margin_hr = 1, n = 2*784.5, M = 100000,
                                   parameterization = 3, censor_beyond_tau = TRUE)

size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt)
size_rmst <- function(arm_ctrl, arm_trt, tau) size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst difference", milestone = tau))
sapply(seq(1, 50), size_rmst, arm_ctrl = arm_ctrl, arm_trt = arm_trt)

size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst difference", milestone = 100))

follow_up_time_vector <- seq(1, 1000, 10)
sample_size_lr_vector <- rep(NA, length(follow_up_time_vector))

for(i in 1:length(follow_up_time)){
  arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                   surv_shape = shape_trt, surv_scale = scale_trt,
                                   loss_scale = scale_loss, follow_time = follow_up_time_vector[i])
  arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                   surv_shape = shape_ctrl,surv_scale = scale_ctrl,
                                   loss_scale = scale_loss, follow_time = follow_up_time_vector[i])
  sample_size_lr_vector[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt)['n']
  follow_up_time_vector[i]

}

arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_trt, surv_scale = scale_trt,
                                 loss_scale = scale_loss, follow_time = 100)
arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                 surv_shape = shape_ctrl,surv_scale = scale_ctrl,
                                 loss_scale = scale_loss, follow_time = 100)

tau_vector = seq(1, 100, 0.1)
RMST_sample_size_difference <- RMST_sample_size_ratio <- NULL

for(i in 1:length(tau_vector)){
 RMST_sample_size_difference[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst difference", milestone = tau_vector[i]))['n']
}

for(i in 1:length(tau_vector)){
  RMST_sample_size_ratio[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst ratio", milestone = tau_vector[i]))['n']
}

LR_sample_size <- min(sample_size_lr_vector, na.rm = TRUE) # total sample size for LR test with follow up time=901, shapes = 1, no loss, scale_ctrl <-  1.5*.0533, scale_trt <-  0.0533. = 195.975.

plot(x = tau_vector, y = RMST_sample_size_difference, type = 'l', log = 'y', xlab = "tau", ylab = "sample size RMST", col = "blue")
lines(x = tau_vector, y = RMST_sample_size_ratio, col = "green")
legend("topright", legend=c("sample size unrestricted LR test closed form", "sample size RMST difference by tau", "sample size RMST ratio by tau"),
       col=c("red", "blue", "green"), lty=1, cex=0.8)
abline(h = LR_sample_size, col = "red")

### scenario 4 Lawrence

events_ctrl <- ctrl <- rweibull(ceiling(109/2), shape = 1, scale = 1)
events_trt <- ctrl <- rweibull(floor(109/2), shape = 3, scale = 1)
events <- c(events_ctrl, events_trt)
arm <- c(rep(0, length(events_ctrl)), rep(1, length(events_trt)))
status <- rep(1, 109)
status[events>=1] <- 0
events[status == 0] <- 1
dataframe <- data.frame(events = events, status = status, arm = arm)
new_dataframe <- dataframe[status == 1,]
result <-  coxph(Surv(events, status) ~ arm, data = new_dataframe)
summary(result)
summary(first_result)
second_result <- result
