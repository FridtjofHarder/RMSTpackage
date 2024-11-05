# required functions. Place folder "compare_sample_sizes" in working directory for sourcing.
source("compare_sample_sizes/functions/testfunctions2.R")
source("compare_sample_sizes/functions/powerRMST.R")
source("compare_sample_sizes/functions/closed_form_noninf.R")
source("compare_sample_sizes/functions/compare_sample_size.R")

# install and load all required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, numDeriv, npsurvSS, SSRMST, survmixer, ggplot2, hesim, survival, RMSTdesign)

# creates dataframe, in which the sample size calculations from 8 different methods/packages are compared
sample_sizes_df <- compare_sample_size(scale_t = 0.01,
                                       scale_c = 0.02,
                                       accrual_time = 1,
                                       follow_up_time = 12,
                                       horizon = 5,
                                       sides = 1,
                                       alpha = 0.025,
                                       power = 0.8,
                                       margin = 0, # non-inferiority not yet working since not supported by all method
                                       closed_form_noninf_RMSTD = T,
                                       closed_form_noninf_LRT_unrestricted = T, # log rank test utilizing information from all events
                                       closed_form_noninf_LRT_restricted = T, # log rank test with observations censored after tau
                                       RMSTdesign_closed_form = T,
                                       RMSTdesign_simulation = T,
                                       powerRMST = T,
                                       ssrmst = T,
                                       survmixer = T
)