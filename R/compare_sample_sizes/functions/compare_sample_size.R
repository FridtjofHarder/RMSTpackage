#compares sample sizes for superiority trials and noninferiority trials, RMST-difference vs. LRT. For different scales but same shapes.



compare_sample_size <- function(scale_t = 0.01,
                                scale_c = 0.02,
                                shape_t = 1,
                                shape_c = 1,
                                accrual_time = 0.00000001,
                                follow_up_time = 12,
                                horizon,
                                sides = 1,
                                alpha = 0.025,
                                power = 0.8,
                                margin = 0,
                                closed_form_noninf_RMSTD = T,
                                closed_form_noninf_LRT_unrestricted = T,
                                closed_form_noninf_LRT_restricted = T,
                                RMSTdesign_closed_form = T,
                                RMSTdesign_simulation = T,
                                powerRMST = T,
                                ssrmst = T,
                                survmixer = T
                                ){
  total_time <- accrual_time + follow_up_time
  # initialize df for storing sample sizes
  sample_sizes_df <- data.frame(matrix(ncol = 8))
  rownames(sample_sizes_df) <- "n or power"
  colnames(sample_sizes_df) <- c("closed_form_noninf RMSTD", "closed_form_noninf LRT unrestricted", "closed_form_noninf LRT restricted", "RMSTdesign closed form", "RMSTdesign simulation", "powerRMST power", "ssrmst power", "survmixer")
  
  if(closed_form_noninf_RMSTD | closed_form_noninf_LRT_unrestricted | closed_form_noninf_LRT_restricted){
    closed_form_noninf_arm_t <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_t, surv_shape = shape_t,
                                 loss_scale = 0.0000000000000000000000001, loss_shape = 1, total_time = total_time)
    closed_form_noninf_arm_c <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_c, surv_shape = shape_c,
                                 loss_scale = 0.0000000000000000000000001, loss_shape = 1, follow_time = follow_up_time)
    if(closed_form_noninf_RMSTD){
      sample_sizes_df["closed_form_noninf RMSTD"] <- round(closed_form_noninf(closed_form_noninf_arm_t, closed_form_noninf_arm_c, test = list(test = "rmst difference", milestone = horizon), power = 0.8, alpha = alpha, sides = sides, margin = margin)['n'])
    }
    if(closed_form_noninf_LRT_unrestricted & margin == 0){
      sample_sizes_df["closed_form_noninf LRT unrestricted"] <- round(size_two_arm(closed_form_noninf_arm_t, closed_form_noninf_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
    }
    if(closed_form_noninf_LRT_restricted & margin == 0){
      closed_form_noninf_arm_t$follow_time <- closed_form_noninf_arm_c$follow_time <- horizon
      closed_form_noninf_arm_t$total_time <- closed_form_noninf_arm_c$total_time <- horizon + accrual_time
      sample_sizes_df["closed_form_noninf LRT restricted"] <- round(size_two_arm(closed_form_noninf_arm_t, closed_form_noninf_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
    }
    
  }
  if(RMSTdesign_closed_form & margin == 0){
    sample_sizes_df["RMSTdesign closed form"] <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = F, method = 'tau_star', 
                                                                  two.sided = sides == 2, alpha = alpha)$n)
    
  }
  if(RMSTdesign_simulation & margin == 0){
    sample_sizes_df["RMSTdesign simulation"]  <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = T, method = 'tau_star', 
                                                                  two.sided = sides == 2, alpha = alpha, M = 10)$n)
    
  }
  if(powerRMST & closed_form_noninf_RMSTD){
    sample_sizes_df["powerRMST power"] <-  powerRMST(n = sample_sizes_df[["closed_form_noninf RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=margin, one_sided_alpha=alpha/sides, seed=NULL, M = 100)$power
    
  }
  if(ssrmst & closed_form_noninf_RMSTD){
    sample_sizes_df["ssrmst power"] <- ssrmst(ac_number = sample_sizes_df[["closed_form_noninf RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=0, 
                       one_sided_alpha=alpha/sides, ntest=100)$power1
    
  }
  if(survmixer & margin == 0){
    sample_sizes_df["survmixer"] <-  round(survm_samplesize(
      ascale0_r = 1/scale_c,
      ascale0_nr = 1/scale_c,
      ascale1_r = 1/scale_t,
      ascale1_nr = 1/scale_t,
      delta_p = 0,
      p0 = 1,
      m0_r = 1,
      m0_nr = 1,
      diffm_r = 1,
      diffm_nr = 1,
      S0_r = 1,
      S0_nr = 1,
      diffS_r = 1,
      diffS_nr = 1,
      Delta_r = 1,
      Delta_nr = 1,
      ascale_cens = 1000,
      tau = horizon,
      bshape0 = 1,
      bshape1 = 1,
      alpha = alpha/sides,
      beta = 0.2,
      set_param = 0
    )[1, 2])
  }
  
  return(sample_sizes_df)
}