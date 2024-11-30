# compare RMST with Cox-PH regarding different patterns of administrative censorins

# load required packages
if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survRM2')) install.packages('survRM2'); library('survRM2')

#set parameters
parameterization = 1 # S(t)=exp[-(t/scale)^shape] = Parameterization 1
size_ctrl <- size_trmt <- 100
iterations <- 1000 # iterations for simulation
scale_ctrl <- 1
scale_trmt <- 1.5
follow_up <- 1
accrual <- 1

# function for simulating Cox PH
cox_power_superiority <- function(iterations = 1000,
                                  arm_ctrl,
                                  arm_trmt,
                                  tau){
  counter_cox_ph_positive <- rep(NA, iterations) # initialize counter of positive tests

  # Cox simulated with all observations censored beyond tau
  for(i in 1:iterations){
    arm_ctrl_sample <- simulate_arm(arm_ctrl, label = 0) #draw samples without restriction
    arm_trmt_sample <- simulate_arm(arm_trmt, label = 1)
    data <- data.frame(time = c(arm_ctrl_sample$time.surv,
                                arm_trmt_sample$time.surv),
                       arm = c(rep(0, size_ctrl), rep(1, size_trmt)),
                       status = as.numeric(
                         c(arm_ctrl_sample$time.surv,
                           arm_trmt_sample$time.surv) < tau
                       )) # create data frame of observation times, arm label, and status. All observations beyond tau are censored.
    fit <- coxph(Surv(time, status) ~ arm, data = data) # arm_trmt has sign. better survival, then upper CI limit of CI will be < 1
    counter_cox_ph_positive[i] <-
      as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
    # without restriction at tau
  }
  return(sum(counter_cox_ph_positive)/iterations)
}

# define arms in npsurvSS
arm_ctrl <- create_arm(size=size_ctrl, accr_time=accrual, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                       loss_scale=0.00000001, follow_time=follow_up)
arm_trmt <- create_arm(size=size_trmt, accr_time=accrual, surv_scale=1/scale_trmt,
                       loss_scale=0.00000001, follow_time=follow_up)

# 1: tau is set directly before admin censoring sets in ---------------------------
tau <- 1

# power RMST
power_two_arm(arm_ctrl, arm_trmt, list(test="rmst difference", milestone=tau))

# power Cox
cox_power_superiority(arm_ctrl <- arm_ctrl,
                      arm_trmt <- arm_trmt,
                      tau <- tau)

# 2: admin censoring sets in before tau. tau is set at follow-up + 1/2 recruitment period. ---------

tau <- 1.5

# power RMST
power_two_arm(arm_ctrl, arm_trmt, list(test="rmst difference", milestone=tau))

# power Cox
cox_power_superiority(arm_ctrl <- arm_ctrl,
                      arm_trmt <- arm_trmt,
                      tau <- tau)

# 3: tau set at end follow-up + recruitment period -------------------------------------------------

tau <- 2
# define curves in npsurvSS

arm_ctrl <- create_arm(size=size_ctrl, accr_time=accrual, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                       loss_scale=0.00000001, follow_time=follow_up)
arm_trmt <- create_arm(size=size_trmt, accr_time=accrual, surv_scale=1/scale_trmt,
                       loss_scale=0.00000001, follow_time=follow_up)
# power RMST
power_two_arm(arm_ctrl, arm_trmt, list(test="rmst difference", milestone=tau))

# power Cox
cox_power_superiority(arm_ctrl <- arm_ctrl,
                      arm_trmt <- arm_trmt,
                      tau <- tau)
