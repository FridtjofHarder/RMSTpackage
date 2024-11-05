## install the following packages + ancillaries: npsurvSS; rmstdesign, powerrmst, ssrmst
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, numDeriv, npsurvSS, SSRMST, survmixer, ggplot2, hesim, survival)

install_github("anneae/RMSTdesign")
library(RMSTdesign)
library(numDeriv)
library(devtools)
library(npsurvSS)
library(survival)


install.packages("SSRMST")
install_github("pauknemj/survWMST")
library(survWMST)

library(SSRMST)

scale_t = 0.01
scale_c = 0.02
horizon = 5
accrual_time = 0.00000001
follow_up <- 12
total_time <- accrual_time + follow_up
margin = 0

# npsurvSS
npsurvSS_arm_t <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_t, surv_shape = 1,
                             loss_scale = 0.0000000000000000000000001, loss_shape = 1, total_time = total_time)
npsurvSS_arm_c <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_c, surv_shape = 1,
                             loss_scale = 0.0000000000000000000000001, loss_shape = 1, total_time = total_time)
n_npsurvSS <- size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "rmst difference", milestone = horizon), power = 0.8, alpha = 0.025, sides = 1)
n_npsurvSS['n']


#RMSTdesign
RMSTdesign_arm_t <- survdef(haz = scale_t)
RMSTdesign_arm_c <- survdef(haz = scale_c)
plotsurvdef(survdefC = RMSTdesign_arm_c, survdefT = RMSTdesign_arm_t, xupper = 3)
n_RMSTdesign_analytic <- RMSTpow(survdefT = RMSTdesign_arm_t, survdefC = RMSTdesign_arm_c, k1 = accrual_time, k2 = follow_up, power = 0.8, tau = horizon, sim = F, M = 1000, method = 'tau_star', 
                                 two.sided = T, alpha = 0.05)
n_RMSTdesign_analytic
n_RMSTdesign_numeric<- RMSTpow(RMSTdesign_arm_t, RMSTdesign_arm_c, k1 = accrual_time, k2 = follow_up, power = 0.8, tau = horizon, sim = T, M = 1000, method = 'tau_star', 
                               two.sided = T, alpha = 0.05)
n_RMSTdesign_numeric

#powerRMST
n_powerRMST <-  powerRMST(n = 1200, ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=0, one_sided_alpha=0.025, seed=NULL)
n_powerRMST

#ssrmst
n_ssrmst <- ssrmst(ac_number = 1278, ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=0, 
                   one_sided_alpha=0.025, ntest=100) ## not working...????
n_ssrmst

#survmixer
n_survmixer <- survm_samplesize(
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
  all_ratio = 0.5,
  alpha = 0.025,
  beta = 0.2,
  set_param = 0
)
survmixture_f(100, ascale_r = 1/scale_t, ascale_nr = 1, bshape = 1, p = 1)
psurv(q = 100, npsurvSS_arm_t, F, F)

#my simulation
my_weibull <- function(scale, shape = 1, x){
  exp(-scale*x^shape)
}