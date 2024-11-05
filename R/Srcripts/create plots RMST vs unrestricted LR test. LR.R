# creates plot with sample size for RMST difference against tau and sample size for RMST ratio against tau. Vertical line with sample size for unrestricted log rank test for reference.

if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')

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
RMST_sample_size_difference <- RMST_sample_size_ratio <- NULL

for(i in 1:length(tau_vector)){
  RMST_sample_size_difference[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst difference", milestone = tau_vector[i]))['n']
  RMST_sample_size_ratio[i] <- size_two_arm(arm0 = arm_ctrl, arm1 = arm_trt, list(test = "rmst ratio", milestone = tau_vector[i]))['n']
}

plot(x = tau_vector, y = RMST_sample_size_difference, type = 'l', log = 'y', xlab = "tau", ylab = "sample size RMST", col = "blue")
lines(x = tau_vector, y = RMST_sample_size_ratio, col = "green")
legend("topright", legend=c("sample size unrestricted LR test closed form", "sample size RMST difference against tau", "sample size RMST ratio against tau"),
       col=c("red", "blue", "green"), lty=1, cex=0.8)
abline(h = LR_sample_size, col = "red")
title("n for RMST difference and ratio against tau, and for unrestricted LR test.")
