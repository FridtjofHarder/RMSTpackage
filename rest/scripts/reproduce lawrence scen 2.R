# lawrence ------------------------------

# install and initialize all needed packages
if (!require('npsurvSS')) install.packages('npsurvSS'); library('npsurvSS')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')

# define parameters
size_ctrl <- size_trmt <- 55
iterations <- 10000 # iterations for simulation
scale_ctrl <- 3
scale_trmt <- 1
shape_ctrl <- 1
shape_trmt <- 1
tau <- 1 # acc. to lawrence, all obervations beyond tau are censored

# define arms without restriction to tau in npsurvSS
# dismiss warning "either follow_time nor total_time were defined"
arm_ctrl_unrestricted <- create_arm(size=size_ctrl, accr_time=.00000001, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                                         surv_shape = shape_ctrl, loss_scale=0.00000001, follow_time=Inf)
arm_trmt_unrestricted <- create_arm(size=size_trmt, accr_time=.00000001, surv_scale=1/scale_trmt,
                                         surv_shape = shape_trmt, loss_scale=0.00000001, follow_time=Inf)

# plot design survival curves
ggplot() + xlim(0, 2) + ylim(0, 1) +
  geom_function(fun = npsurvSS::psurv,
                args = list(arm = arm_ctrl_unrestricted,
                            lower.tail = FALSE),
                aes(color = "control")) +
  geom_function(fun = npsurvSS::psurv,
                args = list(arm = arm_trmt_unrestricted,
                            lower.tail = FALSE),
                aes(color = "treatment")) +
  labs(title = "Design survival curves", color = "condition") +
  xlab("t") + ylab("S(t)") +
  geom_vline(xintercept = tau) +
  scale_color_manual(values = c("control" = "red", "treatment" = "blue")) +
  geom_text(aes(label = paste("\u03c4", "=", tau)), x = tau+.1, y = 0.9)

# plot hazards

ggplot() + xlim(0, 2) + ylim(0, 10) +
  geom_function(fun = npsurvSS::hsurv,
                args = list(arm = arm_ctrl_unrestricted),
                aes(color = "control")) +
  geom_function(fun = npsurvSS::hsurv,
                args = list(arm = arm_trmt_unrestricted),
                aes(color = "treatment")) +
  labs(title = "Hazards of design survival curves", color = "condition") +
  xlab("t") + ylab("h(t)") +
  geom_vline(xintercept = tau) +
  scale_color_manual(values = c("control" = "red", "treatment" = "blue")) +
  geom_text(aes(label = paste("\u03c4", "=", tau)), x = tau+.1, y = 0.9)

counter_cox_ph_positive <- rep(NA, iterations) # initialize counter of positive tests

# Cox simulated with all observations censored beyond tau
for(i in 1:iterations){
  arm_ctrl_sample <- simulate_arm(arm_ctrl_unrestricted, label = 0) #draw samples without restriction
  arm_trmt_sample <- simulate_arm(arm_trmt_unrestricted, label = 1)
  data <- data.frame(time = c(arm_ctrl_sample$time.surv,
                              arm_trmt_sample$time.surv),
                     arm = c(rep(0, size_ctrl), rep(1, size_trmt)),
                     status = as.numeric(
                       c(arm_ctrl_sample$time.surv,
                         arm_trmt_sample$time.surv) < tau
                     )) # create data frame of observation times, arm label, and status. All observations beyond tau are censored.
  # set all times beyond tau equal to tau
  data$time[data$time > tau] <- tau
  fit <- coxph(Surv(time, status) ~ arm, data = data) # arm_trmt has sign. better survival, then upper CI limit of CI will be < 1
  counter_cox_ph_positive[i] <-
    as.numeric(summary(fit)$conf.int[, 'upper .95'] < 1)
  # without restriction at tau
}

paste("test power is: ",sum(counter_cox_ph_positive)/iterations) # print test power

# closed-form power calculated with npsurvSS, restricted to tau = 1 by stopping FU at tau = 1.
arm_ctrl_restricted_to_tau <- create_arm(size=size_ctrl, accr_time=.00001, surv_scale=1/scale_ctrl, # parameterization in npsurvSS follows "second alternative" in wikipedia
                                     surv_shape = shape_ctrl, loss_scale=0.00001, follow_time=tau)
arm_trmt_restricted_to_tau <- create_arm(size=size_trmt, accr_time=.00001, surv_scale=1/scale_trmt,
                                     surv_shape = shape_trmt, loss_scale=0.00001, follow_time=tau)
power_two_arm(arm_ctrl_restricted_to_tau, arm_trmt_restricted_to_tau)

# plot survival curves from last iteration as example. All survivors are censored at tau
surv_obj <- Surv(time = data$time, event = data$status == 1)
kmfit = survfit(surv_obj ~ data$arm)
plot(kmfit, col = c("red", "blue"), xlab = "t", ylab = "S(t)", mark.time = TRUE)

