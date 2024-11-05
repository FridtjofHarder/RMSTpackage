require(biostat3)  
require(rstpm2)  
require(dplyr)   
require(ggplot2)
require(npsurvSS)
require(geepack)
require(pseudo)

scale_0 <- c(.4) # set scale for exponential dis of survival in control
interval_0 <- c(0, Inf) # interval for respective scales
interval_1 <- c(0, Inf)
evaluation_point <- 8
follow_up <- 2 
accrual_time <- 7
total_length <- accrual_time + follow_up
use_KM = TRUE
use_flexsurv = TRUE # takes some time
use_splines = TRUE
use_pseudo = TRUE
pseudo_size <- 5000
m <- 1000 # sample size in MC-simulation
M <- 50 # no. of iterations
loss_to_FU_scale_0 <- .1 # define Weibull scale loss to follow up. Use 0 if no loss assumed.
loss_to_FU_scale_1 <- .1
loss_to_FU_shape_0 <- 1 # define Weibull shape loss to follow up. Use 0 if no loss assumed.
loss_to_FU_shape_1 <- 1
set_seed <- TRUE # set seed for reproducible results

if (set_seed) {set.seed(42)}


HR <- seq(from = 1, to = 5, by = .1)
RMSTR <- rep(NA, 41)
for (i in 1:41){
scale_1 <- scale_0 / HR[i]
integration_limit <- min(total_length, evaluation_point)
arm_0 <- create_arm(size = m, accr_time = accrual_time, surv_scale = scale_0, 
                    surv_interval = interval_0, loss_scale = loss_to_FU_scale_0, 
                    loss_shape = loss_to_FU_shape_0, follow_time = follow_up)
arm_1 <- create_arm(size = m, accr_time = accrual_time, surv_scale = scale_1, 
                    surv_interval = interval_1, loss_scale = loss_to_FU_scale_1, 
                    loss_shape = loss_to_FU_shape_1, follow_time = follow_up)
n_npsurvSS_0 <- size_two_arm(arm_0, arm_1, list(test="rmst difference", 
                                                milestone = integration_limit))

sim_arm_0 <- simulate_arm(arm_0, label = 0)
sim_arm_1 <- simulate_arm(arm_1, label = 1)
sim_trial <- simulate_trial(arm0 = arm_0, arm1 = arm_1)

status_0 <- status_1 <- rep(1, m)

event_times_0 <- rsurv(m, arm_0, include_cured = F)
observation_stops_0 <- runif(m, follow_up, total_length)
status_0[event_times_0 > observation_stops_0] <- 0
surv_times_0 <- pmin(event_times_0, observation_stops_0)

event_times_1 <- rsurv(m, arm_1, include_cured = F)
observation_stops_1 <- runif(m, follow_up, total_length)
status_1[event_times_1 > observation_stops_1] <- 0
surv_times_1 <- pmin(event_times_1, observation_stops_1)

if (loss_to_FU_scale_0 > 0) {
  loss_arm_0 <- create_arm(size = 1, accr_time=accrual_time, 
                           surv_scale = loss_to_FU_scale_0, 
                           surv_shape = loss_to_FU_shape_0, loss_scale = 0, 
                           follow_time = follow_up)
  loss_to_FU_times_0 <- rsurv(m, loss_arm_0, include_cured = F)
  status_0[surv_times_0 > loss_to_FU_times_0] <- 0
  surv_times_0 <- pmin(surv_times_0, loss_to_FU_times_0)
}

if (loss_to_FU_scale_1 > 0) {
  loss_arm_1 <- create_arm(size = 1, accr_time=accrual_time, 
                           surv_scale = loss_to_FU_scale_1, 
                           surv_shape = loss_to_FU_shape_1, loss_scale = 0, 
                           follow_time = follow_up)
  loss_to_FU_times_1 <- rsurv(m, loss_arm_1, include_cured = F)
  status_1[surv_times_1 > loss_to_FU_times_1] <- 0
  surv_times_1 <- pmin(surv_times_1, loss_to_FU_times_1)
} 

integration_limit <- min(evaluation_point, max(surv_times_0), 
                         max(surv_times_0))

surv_obj_0 <- Surv(surv_times_0, status_0)
surv_obj_1 <- Surv(surv_times_1, status_1)

times_column <- c(surv_times_0, surv_times_1)
group_column <- c(rep(0, m), rep(1, m))
status_column <- c(status_0, status_1)
survival_frame <- data.frame(times_column, group_column, status_column)
fit_obj <- stpm2(Surv(times_column, status_column) ~ group_column, df = 3, 
                 data=survival_frame, tvc = list(group_column = 1))

summary_0 <- predict(fit_obj, newdata=data.frame(group_column = 0, 
                                                 times_column = integration_limit), type = "rmst", 
                     se.fit = TRUE)
RMST_splines_0 <- summary_0[[1]]


summary_1 <- predict(fit_obj, newdata=data.frame(group_column = 1, 
                                                 times_column = integration_limit), type = "rmst", 
                     se.fit = TRUE)
RMST_splines_1 <- summary_1[[1]]

RMSTR[i] <- RMST_splines_1/RMST_splines_0
}
## RMSTR, RMSTD, HR are defined.
#############################################################
par(mfrow=c(2,1))
x <- 0:(100 * (follow_up + accrual_time)) / 100
plot(HR, RMSTR)
curve(x*(1 - exp(-scale_0*integration_limit/x))/(1 - exp(-scale_0*integration_limit)), from=1, to=5, , xlab="HR", ylab="RMSTR", add = T) ## x = HR
title("RMSTR against HR for hazard_0 = .4, integration limit = 8, \n exponential model")

plot(survfit(surv_obj_0~1), mark.time=T, conf.int = F, xlab = "t", ylab = "S(t)", 
     xlim = c(0, total_length))
lines(survfit(surv_obj_1~1), mark.time=T, col="red", conf.int = F)
p_0 <- psurv(x, arm_0, lower.tail=F)
lines(x, p_0, col="black", lty = "dotted")
p_1 <- psurv(x, arm_1, lower.tail=F)
lines(x, p_1, col="red", lty = "dotted")
title("hazard ratio = 5, lambda_0 = .4")




