### test spieltplatz wiederköschen

# package example, but accr_time lowere dfrom 6 to 0.00001 -----------------------
arm0 <- create_arm(size=120, accr_time=.00001, surv_scale=0.05, loss_scale=0.005, follow_time=12)
arm1 <- create_arm(size=120, accr_time=.00001, surv_scale=0.03, loss_scale=0.005, follow_time=12)
power_two_arm(arm0, arm1)
power_two_arm(arm0, arm1, list(test="rmst ratio", milestone=12))

# power increases for both when surv_scale in arm1 decreases

# varying shape instead of scale

arm0 <- create_arm(size=120, accr_time=.00001, surv_scale=0.05, loss_scale=0.005, follow_time=12)
arm1 <- create_arm(size=120, accr_time=.00001, surv_scale=0.03, loss_scale=0.005, follow_time=12)
power_two_arm(arm0, arm1)
power_two_arm(arm0, arm1, list(test="rmst ratio", milestone=12))

# lawrence ------------------------------

# simulate RMST after method by lawrence and Coxph for scenario2, esp. for Cox with:
# 1. follow_time = 1
# 2. Follow_time = inf but manually censored at 1
# 3. Follow_time = inf w/o censoring

arm0_restricted_to_tau <- create_arm(size=55, accr_time=.00001, surv_scale=1, surv_shape = 1, loss_scale=0.00001, follow_time=1)
arm1_restricted_to_tau <- create_arm(size=54, accr_time=.00001, surv_scale=1, surv_shape = 3, loss_scale=0.00001, follow_time=1)
power_two_arm(arm0_restricted_to_tau, arm1_restricted_to_tau)
power_two_arm(arm0, arm1, list(test="rmst difference", milestone=1))
arm0_unrestricted <- create_arm(size=55, accr_time=.00001, surv_scale=1, surv_shape = 1, loss_scale=0.00001, follow_time=10)
arm1_unrestricted <- create_arm(size=54, accr_time=.00001, surv_scale=1, surv_shape = 3, loss_scale=0.00001, follow_time=10)


for(1 in 1:10000){
  # simulate 2 arms practically w/o restriction (admin censoring at t=10)
  arm0_unrestricted <- simulate_arm(arm0_unrestricted, label = 0)
  arm1_unrestricted <- simulate_arm(arm1_unrestricted, label = 1)
  data <- data.frame(time = c(arm0_unrestricted$time.surv, arm1_unrestricted$time.surv),
                     arm = c(rep(0, 55), rep(1, 54)),
                     status = 1)
  fit <- coxph(Surv(time, status) ~ arm, data = data)
  z_value <- summary(fit)$coefficients[1, "z"] # if arm1 has better survival, then upper bound of CI will be below 1


  # manually restrict unrestricted arms to tau = 1

  # simulate 2 arms restricted to tau = 1 by limiting follow up to t = 1
  str(summary(fit))
  str(summary(fit)$conf.int)
  int <- summary(fit)$conf.int




}
