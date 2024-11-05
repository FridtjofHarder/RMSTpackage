########## get required packages
packages <- c("survRM2", "SSRMST", "bpcp", "npsurvSS")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

########## functions

# get n for one sample according to Makuch Simon
MakuchSimon<-function(p0,p1,delta,alpha=.9,power=.8){ # get n for one sample accodring to Makuch Simon
  Zalpha<-qnorm(alpha)
  Zbeta<-qnorm(power)
  
  n<-(p0*(1-p0)+p1*(1-p1))*((Zalpha+Zbeta)/(delta-(p0-p1)))^2
  print("For each treatment group")
  return(ceiling(n))
}

# finds scale parameter lambda for exponential survival function from S(t)
lambda_from_surv_prob <- function(p, t){ # finds scale parameter lambda for exponential survival function from S(t)
  lambda <- -log(p)/t
  return(lambda)
}

# define exponential survival function
survival_expo <- function(x, lambda){exp(-lambda*x)}

# translates pointwise survival margin to RMST difference margin between two exponential functions
RMST_margin_from_point_margin <- function(point_margin, p, t){
  lambda <- lambda_from_surv_prob(p, t) # get scale parameter lambda for exponential survival function
  lambda_marginal_function <- lambda_from_surv_prob(p+point_margin, t) # get lambda for exponential survival function passing through S(t) + pointwise margin
  RMST <- integrate(survival_expo, lambda = lambda, lower = 0, upper = t) # get RMST
  RMST_marginal_function <- integrate(survival_expo, lambda = lambda_marginal_function, lower = 0, upper = t) # get RMST for marginal survival function
  RMST_margin <- RMST_marginal_function$value - RMST$value # get RMST margin
  return(RMST_margin)
}

# estimates power for non-inferioriy trial based on RMST margin and simple exponential survival, so far without censoring, estimates also power for logrank and for Cox
simple_power_simulation_noninf <- function(tau, n0, iterations = 10000, alpha_level = .05, beta = .8, 
                                           lambda = NA, survival = NA,
                                           RMST = T, survival_diff = T, hr = T,
                                           RMST_margin = NA, survival_margin = NA, hr_margin = NA){
  count_RMST <- count_survdiff <- count_hr <- rep(NA, iterations) # prepare counter
  status <- rep(1, 2*n0) # prepare status indicator with 1 = event
  arm <- c(rep(1, n0), rep(0, n0)) # prepare arm as group variable
  
  for(i in 1:iterations){ # begin Monte Carlo
    sample_non_inferior <- rexp(n0, rate = lambda) # sample from exponential survival for presumed non-inferior arm
    sample_inferior <- rexp(n0, rate = lambda) # sample from exponential survival for comparison arm
    time <- c(sample_non_inferior, sample_inferior) # concatenate both vectors
    
    if(RMST == T){
    result <- rmst2(time, status, arm, tau, alpha = alpha_level*2) # estimate RMST and standard errors
    lower_RMST <- result$unadjusted.result["RMST (arm=1)-(arm=0)", "lower .90"]
    }
    
    if(lower_RMST > RMST_margin){ # count instances in which lower bound is above margin
      count_RMST[i] <- 1
    } else{
      count_RMST[i] <- 0
    }
    
    if(survival_diff == TRUE){
      df1 <- data.frame(time = sample_non_inferior, status = status[1:n0])
      df0 <- data.frame(time = sample_inferior, status = status[1:n0])
      survival_se1 <- summary(survfit(Surv(time, status) ~ 1, data = df1))$std.err[Position(function(x) x >= tau, sort(sample_non_inferior))]
      survival_se0 <- summary(survfit(Surv(time, status) ~ 1, data = df0))$std.err[Position(function(x) x >= tau, sort(sample_inferior))]
      survival1 <- summary(survfit(Surv(time, status) ~ 1, data = df1))$surv[Position(function(x) x >= tau, sort(sample_non_inferior))]
      survival0 <- summary(survfit(Surv(time, status) ~ 1, data = df0))$surv[Position(function(x) x >= tau, sort(sample_inferior))]
      survivaldiff <- survival1 - survival0
      lower_surv_diff <- survivaldiff-qnorm(1-alpha_level)*sqrt(survival_se1^2+survival_se0^2)
      }
    if(lower_surv_diff > survival_margin){ # count instances in which lower bound is above margin
      count_survdiff[i] <- 1
    } else{
      count_survdiff[i] <- 0
    }
    
    if(hr == T){
      status1 <- as.numeric(sample_non_inferior <= 36)
      status0 <- as.numeric(sample_inferior <= 36)
      time <- c(sample_non_inferior, sample_inferior)
      
      df <- data.frame(time = time, status = c(status1, status0), arm = arm)
      fit <- coxph(Surv(time, status) ~ arm, data = df)
      lower_hr <- exp(confint(fit, level = .9))[1]
    }
    if(lower_hr > hr_margin){ # count instances in which lower bound is above margin
      count_hr[i] <- 1
    } else{
      count_hr[i] <- 0
    }
    
  }    
  power_RMST <- sum(count_RMST)/iterations # power is share of instances with lower bound above margin
  power_survdiff <- sum(count_survdiff)/iterations
  power_hr <- sum(count_hr)/iterations
  
  return(list(power_RMST = power_RMST, power_survdiff = power_survdiff, power_hr = power_hr))
}

# closed form functions for sample size based on RMST

prob_risk <- function(arm, teval) {
  psurv(teval, arm, lower.tail=F) *
    ploss(teval, arm, lower.tail=F) *
    paccr(pmin(arm$accr_time, arm$total_time-teval), arm)
}

sigma2j_rmst <- function(arm, teval) {
  inner <- function(x) {
    sapply(x, function(x1) stats::integrate(function(x2) psurv(x2, arm, lower.tail=F),
                                            lower=x1,
                                            upper=teval)$value
    )
  }
  stats::integrate(function(x) inner(x)^2 * hsurv(x, arm) / prob_risk(arm, x),
                   lower=0,
                   upper=teval)$value
}

########## commands

# calculate n for one sample according to Makuch Simon from S(t) and margin
n0 <- MakuchSimon(.93,.93,0.055,alpha=.95,power=.80)

# get lambda and lambda marginal for exponantial survival function from S(t)
lambda <- lambda_from_surv_prob(.93, 36)
lambda_marginal <- lambda_from_surv_prob(.875, 36)
hr_margin <- lambda/lambda_marginal

# translate pointwise marin in S(t) to RMST margin
RMST_margin <- RMST_margin_from_point_margin(-.055, .93, 36)

RMST_marginal_function <- integrate(survival_expo, lambda = lambda_marginal, lower = 0, upper = 36)

RMST <- integrate(survival_expo, lambda = lambda, lower = 0, upper = 36)

# estimate power with sample size calculated via Makuch Simon. Change no. of iterations if desired.
power <- simple_power_simulation_noninf(lambda = lambda, RMST_margin = RMST_margin, tau = 36, n0 = n0, iterations = 1000, alpha_level = .05, beta = .8, survival_margin = -.055, hr_margin = hr_margin)

# power is below .8. Repeat with larger sample size:
n_increased <- n0 + 80
power_increased_n <- simple_power_simulation_noninf(lambda = lambda, RMST_margin = RMST_margin, tau = 36, n0 = n_increased, iterations = 100, alpha_level = .05, beta = .8, survival_margin = -.055, hr_margin = hr_margin)

# repeat with function from package SSRMST:
power_ssrmst <- ssrmst(ac_rate = n0*2, ac_period = 1, tot_time = 40, tau = 36, shape0 = 1, shape1 = 1, scale0 = 1/lambda, scale1 = 1/lambda, margin = -RMST_margin, one_sided_alpha = .1, ntest = 100)
# with larger sample size:
power_ssrmst_increased_n <- ssrmst(ac_rate = n_increased*2, ac_period = 1, tot_time = 40, tau = 36, shape0 = 1, shape1 = 1, scale0 = 1/lambda, scale1 = 1/lambda, margin = -RMST_margin, one_sided_alpha = .05, ntest = 100)

# closed-form approximation of sample size for RMST

arm1 <- create_arm(size=1, accr_time=.0001, surv_scale=lambda, loss_scale=0.0000000001, follow_time=37) # define arm via npsurvss package

sigma2 <- sigma2j_rmst(arm1, teval = 36) #get sigma^2
alpha <- .05 # set alpha level
sample_size_closed_RMST <- 2*(qnorm(.8)+qnorm(1-alpha))^2*sigma2/RMST_margin^2 # calculate sample size from sigma^2 and margin

# closed-form approximation of sample size for Cox regression

events_required <- function(lambda, lambda_marginal){
  margin <- log(lambda/lambda_marginal)
  return((qnorm(1-.05) + qnorm(.8))^2 / (.25 * margin^2))
}



########## conclusion:
# in the situation presented, larger sample sizes are needed when non-inferiority analysis is performed based on RMST compared to pointwise survival.

########## ToDo: LogRank test simulation und analytisch, dazu Formel anschauen, S(t) difference simulation

########## Frage: arm 1 und 2 richtig?

# analysic_var_simple <- function(lambda, tau){
#   expo_surv <- function(x, lambda){exp(-lambda * x)}
#   inner_integral <- function(lambda, t, tau){
#     integrate(expo_surv, lambda, lower = t, upper = tau, subdivisions = 20000)
#   }
#   outer_function <- function(lambda, t, tau){
#     inner_integral(lambda, t, tau)$value^2 * lambda / expo_surv(lambda, t)
#   }
#   return(
#     integrate(outer_function, tau, lambda, lower = 0, upper = tau, subdivisions = 20000)
#   )
# }

#### nüpsurvss ##########

### meins ##########
sigma2 <- sigma2j_rmst(arm1, teval = 36)
alpha <- .05
2*(qnorm(.8)+qnorm(1-alpha))^2*sigma2/RMST_margin^2
