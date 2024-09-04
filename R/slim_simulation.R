# non-exported function. Slim simulation
# allows for simulating power with loss to FU and admiin censoring. Benchmark with ssrmst
# compare drawing for one group each vs drawing for both

scale_ctrl <- 1
scale_trt <- 1.6
shape_ctrl <- 1
shape_trt <- 1
tau <- 0.75
accrual_time <- 0.5
follow_up_time <- 100
total_time <- accrual_time + follow_up_time
margin <-  0
scale_loss <- 10
shape_loss <- 1
seed <- 42
ntest <- 1000
n <- 250
one_sided_alpha <- 0.025
M <- 100

sample_ctrl <- create_arm(size = 1, accr_time = accrual_time, surv_shape = shape_ctrl, surv_scale = 1/scale_ctrl, follow_time = follow_up_time,
                          loss_scale = 0.0000001)

sample_trt <- create_arm(size = 1, accr_time = accrual_time, surv_shape = shape_trt, surv_scale = 1/scale_trt, follow_time = follow_up_time,
                          loss_scale = 0.0000001)
n_npsurvSS_RMST <- size_two_arm(sample_ctrl, sample_trt, test = list(test = "rmst difference", milestone = tau))["n"]
n_npsurvSS_LR <- size_two_arm(sample_ctrl, sample_trt, )["n"]

ssrmst(ac_period = accrual_time, ac_number = n_npsurvSS, tot_time = accrual_time + follow_up_time, scale0 = scale_ctrl, scale1 = scale_trt,
       shape0 = shape_ctrl, shape1 = shape_trt, seed = seed, ntest = 386, margin = margin, tau = tau)

slim_simulation <- function(scale_ctrl = NULL, shape_ctrl = 1, scale_trt = NULL,
                            shape_trt = 1, scale_loss = NULL, shape_loss = NULL,
                            accrual_time = NULL, follow_up_time = NULL,
                            total_time = NULL,
                            tau = NULL, alpha = 0.025, sides = 1, power = 0.8,
                            margin = 0, n_RMST = NULL, n_LR = NULL, M = M){
  positive_RMST <- NULL
  positive_LR <- NULL
  for (i in 1:M){
    print(M)
    n_per_group <- round(n_npsurvSS/2)
    events <- c(rweibull(n_per_group, scale = scale_ctrl, shape = shape_ctrl),
                rweibull(n_per_group, scale = scale_trt, shape = shape_trt))
    loss_to_fu <- rep(Inf, 2*n_per_group)
    #loss_to_fu  <- rweibull(2*round(n_npsurvSS/2), scale = scale_loss, shape = shape_loss)
    admin_loss  <- total_time - runif(2*n_per_group, max = accrual_time)

    status <- as.numeric(events <= loss_to_fu & events <= admin_loss)
    time <- pmin(events, loss_to_fu, admin_loss)
    arm <- c(rep(0, n_per_group), rep(1, n_per_group)) #  0 = ctrl, 1 = trt
    data <- data.frame(time, status, arm)
    result_RMST <-  rmst2(data$time, data$status, data$arm, tau = tau,
                 alpha = one_sided_alpha * 2)$unadjusted.result
    lower <-  result[1, 2]
    positive_RMST[i] <- if(lower > -margin){positive[i] <- 1}
    else positive_RMST[i] <- 0



  }
  return(sum(positive)/M)

  positive_RMST <-   slim_simulation(scale_ctrl = scale_ctrl, shape_ctrl = shape_ctrl, scale_trt = scale_trt,
           shape_trt = shape_trt, scale_loss = scale_loss, shape_loss = shape_loss,
           accrual_time = accrual_time, follow_up_time = follow_up_time,
           total_time = total_time,
           tau = tau, alpha = 0.025, sides = 1, power = 0.8,
           margin = 0, n_RMST = n_npsurvSS_RMST, n_LR = n_npsurvSS_LR, M = M)
  power_RMST <- sum(positive)

}

events[events<=losstofu]

rmsamplesize(beta = 0.2, alpha = 0.025, milestone = 18,
             accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
             lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = 22,
             followupTime = 18, fixedFollowup = FALSE)

# compare npsurvSS vs lrstat.
size_two_arm(sample_ctrl, sample_trt, test = list(test = "rmst difference", milestone = tau))["n"]

rmsamplesize(beta = 0.2,
            alpha = 0.025,
            milestone = 18,
            allocationRatioPlanned = 1, accrualTime = seq(0, 5),
            accrualIntensity = 1/9*seq(1, 6),
            lambda1 = 0.0533,
            lambda2 = 1.5*0.0533,
            accrualDuration = 5,
            followupTime = 18, fixedFollowup = TRUE, round = 0)

sample_ctrl <- create_arm(size = 1, accr_time = 5, surv_scale = 1.5*0.0533, follow_time =18,
                          loss_scale = 0.0000001)

sample_trt <- create_arm(size = 1, accr_time = 5, surv_scale = 0.0533, follow_time =18,
                          loss_scale = 0.0000001)

size_two_arm(sample_ctrl, sample_trt, test = list(test = "rmst difference", milestone = 18))["n"]

