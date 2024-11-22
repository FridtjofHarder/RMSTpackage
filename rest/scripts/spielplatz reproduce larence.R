# lawrence ------------------------------

# install and initialize all needed packages
if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('ggplot2')) install.packages('survRM2'); library('survRM2')

# simulate RMST after method by lawrence and Coxph for scenario2, esp. for Cox with:
# 1. follow_time = 1
# 2. Follow_time = inf but manually censored at 1
# 3. Follow_time = inf w/o censoring

# define parameters
size_ctrl <- size_trmt <- 55
iterations <- 1000 # iterations for simulation
scale_ctrl <- 1
scale_trmt <- 1
shape_ctrl <- 1
shape_trmt <- 3 # should be 3, all else 1
tau <- 1 # acc. to lawrence, all obervations beyond tau are censored

# define arms without restriction to tau in npsurvSS
arm_ctrl_unrestricted <- create_arm(size=size_ctrl, accr_time=.00000001, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                                         surv_shape = shape_ctrl, loss_scale=0.00000001, follow_time=Inf)
arm_trmt_unrestricted <- create_arm(size=size_trmt, accr_time=.00000001, surv_scale=1/scale_trmt,
                                         surv_shape = shape_trmt, loss_scale=0.00000001, follow_time=Inf)
# dismiss warning "either follow_time nor total_time were defined"

# plot curves
ggplot() + xlim(0, 2) + ylim(0, 1) +
  geom_function(fun = npsurvSS::psurv,
                args = list(arm = arm_ctrl_unrestricted,
                            lower.tail = FALSE),
                aes(color = "control")) +
  geom_function(fun = npsurvSS::psurv,
                args = list(arm = arm_trmt_unrestricted,
                            lower.tail = FALSE),
                aes(color = "treatment"))

counter_cox_ph_positive <- rep(NA, iterations)
# Cox simulated with all observations censored beyond tau
for(i in 1:iterations){
  arm_ctrl_sample <- simulate_arm(arm_ctrl_unrestricted, label = 0) #draw samples without restriction
  arm_trmt_sample <- simulate_arm(arm_trmt_unrestricted, label = 1)
  data <- data.frame(time = c(arm_ctrl_sample$time.surv,
                              arm_trmt_sample$time.surv),
                     arm = c(rep(0, 55), rep(1, 55)),
                     status = as.numeric(
                       c(arm_ctrl_sample$time.surv,
                         arm_trmt_sample$time.surv)<1
                     )) # create data frame of observation times, arm label, and status. All observations beyond tau are censored.
  fit <- coxph(Surv(time, status) ~ arm, data = data) # if arm1 has better survival, then upper bound of CI will be below 1
  counter_cox_ph_positive[i] <-
    as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
}

sum(counter_cox_ph_positive)/iterations

  # simulate 2 arms with restriction by limiting follow-up to t = 1
  arm0_restricted_to_tau_simul <- simulate_arm(arm0_restricted_to_tau, label = 0)
  arm1_restricted_to_tau_simul <- simulate_arm(arm1_restricted_to_tau, label = 1)
  data <- data.frame(time = c(arm0_restricted_to_tau_simul$time.surv,
                              arm1_restricted_to_tau_simul$time.surv),
                     arm = c(rep(0, 55), rep(1, 55)),
                     status = 1)
  fit <- coxph(Surv(time, status) ~ arm, data = data) # if arm1 has better survival, then upper bound of CI will be below 1
  counter_cox_ph_restricted_to_tau[i] <-
    as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
}

# Cox cimulated with FU restricted to tau
# Cox simulated with unrestricted FU

counter_cox_ph_unrestricted <- rep(NA, iterations)
counter_cox_ph_restricted_to_tau <- rep(NA, iterations)
counter_cox_ph_restricted_manually <- rep(NA, iterations)

for(i in 1:iterations){
  # simulate 2 arms practically w/o restriction (admin censoring at t=10)
  arm0_unrestricted_simul <- simulate_arm(arm0_unrestricted, label = 0)
  arm1_unrestricted_simul <- simulate_arm(arm1_unrestricted, label = 1)
  data <- data.frame(time = c(arm0_unrestricted_simul$time.surv,
                              arm1_unrestricted_simul$time.surv),
                     arm = c(rep(0, 55), rep(1, 55)),
                     status = 1)
  fit <- coxph(Surv(time, status) ~ arm, data = data) # if arm1 has better survival, then upper bound of CI will be below 1
  counter_cox_ph_unrestricted[i] <-
    as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)

  # simulate 2 arms with restriction by limiting follow-up to t = 1
  arm0_restricted_to_tau_simul <- simulate_arm(arm0_restricted_to_tau, label = 0)
  arm1_restricted_to_tau_simul <- simulate_arm(arm1_restricted_to_tau, label = 1)
  data <- data.frame(time = c(arm0_restricted_to_tau_simul$time.surv,
                              arm1_restricted_to_tau_simul$time.surv),
                     arm = c(rep(0, 55), rep(1, 55)),
                     status = 1)
  fit <- coxph(Surv(time, status) ~ arm, data = data) # if arm1 has better survival, then upper bound of CI will be below 1
  counter_cox_ph_restricted_to_tau[i] <-
    as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
}
paste("test power of LRT without restricting survival is ",
      sum(counter_cox_ph_unrestricted)/iterations)
paste("test power of LRT with follow-up = 1 is ",
      sum(counter_cox_ph_restricted_to_tau)/iterations)

# closed-form power calculated with npsurvSS, restricted to tau = 1 by stopping FU at tau = 1.
arm_ctrl_restricted_to_tau <- create_arm(size=size_ctrl, accr_time=.00001, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                                     surv_shape = shape_ctrl, loss_scale=0.00001, follow_time=tau)
arm_trmt_restricted_to_tau <- create_arm(size=size_trmt, accr_time=.00001, surv_scale=1/scale_trmt,
                                     surv_shape = shape_trmt, loss_scale=0.00001, follow_time=tau)
power_two_arm(arm_ctrl_restricted_to_tau, arm_trmt_restricted_to_tau) # cox power acc. to package npsurvSS, when FU stops at tau = 1.
