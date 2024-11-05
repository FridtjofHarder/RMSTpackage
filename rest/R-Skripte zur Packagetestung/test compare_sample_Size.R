scale_trt <- 1.4
shape_trt <- 1
scale_ctrl <- 1
shape_ctrl <- 1
accrual_time <- 1
follow_up_time <- 10
tau <- 1
sides <- 1
alpha <- 0.025
power <- 0.8
margin <- 0.1
closed_form <- T
simulation <- T
M <- 100

compare_sample_size(scale_trt = scale_trt, shape_trt = shape_trt, scale_ctrl = scale_ctrl, shape_ctrl = shape_ctrl,
                                accrual_time = accrual_time,
                                follow_up_time = follow_up_time,
                                tau = tau, sides = sides,
                                alpha = alpha, power = power, margin = -margin, closed_form = closed_form,
                                simulation = simulation, M = M)

test_arm <- create_arm(size = 1, accr_time = 1, surv_shape = 1.1, surv_scale = 1.2, loss_scale = 0.000001, follow_time = 100)
psurv(1, test_arm, lower.tail = F)
pweibull(1, shape = 1.1, scale = 1/1.2, lower.tail = F)

ssrmst(ac_rate = 117, ac_period = accrual_time, tot_time = accrual_time + follow_up_time, tau = tau, scale0 = scale_ctrl, shape0 = shape_ctrl,
       scale1 = scale_trt, shape1 = shape_trt, margin = margin, ntest = 1000)


# if(npsurvSS_RMSTD | npsurvSS_LRT_unrestricted | npsurvSS_LRT_restricted){
#   npsurvSS_arm_t <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_t, surv_shape = shape_t,
#                                loss_scale = 0.0000000000000000000000001, loss_shape = 1, total_time = total_time)
#   npsurvSS_arm_c <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_c, surv_shape = shape_c,
#                                loss_scale = 0.0000000000000000000000001, loss_shape = 1, follow_time = follow_up_time)
#   if(npsurvSS_RMSTD){
#     sample_sizes_df["npsurvSS RMSTD"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "rmst difference", milestone = horizon), power = 0.8, alpha = alpha, sides = sides)['n'])
#   }
#   if(npsurvSS_LRT_unrestricted){
#     sample_sizes_df["npsurvSS LRT unrestricted"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
#   }
#   if(npsurvSS_LRT_restricted){
#     npsurvSS_arm_t$follow_time <- npsurvSS_arm_c$follow_time <- horizon
#     npsurvSS_arm_t$total_time <- npsurvSS_arm_c$total_time <- horizon + accrual_time
#     sample_sizes_df["npsurvSS LRT restricted"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
#   }
#   print(sample_sizes_df)
# }
# if(RMSTdesign_closed_form & margin == 0){
#   sample_sizes_df["RMSTdesign closed form"] <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = F, method = 'tau_star',
#                                                              two.sided = sides == 2, alpha = alpha)$n)
#   print(sample_sizes_df)
# }
# if(RMSTdesign_simulation & margin == 0){
#   sample_sizes_df["RMSTdesign simulation"]  <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = T, method = 'tau_star',
#                                                              two.sided = sides == 2, alpha = alpha, M = 10)$n)
#   print(sample_sizes_df)
# }
# if(powerRMST & npsurvSS_RMSTD){
#   sample_sizes_df["powerRMST power"] <-  powerRMST(n = sample_sizes_df[["npsurvSS RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=margin, one_sided_alpha=alpha/sides, seed=NULL, M = 10)$power
#   print(sample_sizes_df)
# }
# if(ssrmst & npsurvSS_RMSTD){
#   sample_sizes_df["ssrmst power"] <- ssrmst(ac_number = sample_sizes_df[["npsurvSS RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=0,
#                                             one_sided_alpha=alpha/sides, ntest=10)$power1
#   print(sample_sizes_df)
# }
# if(survmixer & margin == 0){
#   sample_sizes_df["survmixer"] <-  round(survm_samplesize(
#     ascale0_r = 1/scale_c,
#     ascale0_nr = 1/scale_c,
#     ascale1_r = 1/scale_t,
#     ascale1_nr = 1/scale_t,
#     delta_p = 0,
#     p0 = 1,
#     m0_r = 1,
#     m0_nr = 1,
#     diffm_r = 1,
#     diffm_nr = 1,
#     S0_r = 1,
#     S0_nr = 1,
#     diffS_r = 1,
#     diffS_nr = 1,
#     Delta_r = 1,
#     Delta_nr = 1,
#     ascale_cens = 1000,
#     tau = horizon,
#     bshape0 = 1,
#     bshape1 = 1,
#     alpha = alpha/sides,
#     beta = 0.2,
#     set_param = 0
#   )[1, 2])
# }
# print(sample_sizes_df)
# return(sample_sizes_df)
#                               }
#


