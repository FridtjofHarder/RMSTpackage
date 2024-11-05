# Sample size determination code: 
## includes 1. function (powerRMST) to compute power under a given sample size with other specified parameters.
## includes 2. iterative process for sample size determination for a target power


# 1. powerRMST():  function to determine power under a specified sample size   
## Takes as inputs a sample size (n), accrual period (ac_period), length of trial (tot_time), time horizon (tau)
## Weibull scale parameter for Control arm(scale0),Weibull scale parameter for Experimental arm(scale1),
## Weibull shape parameter (shape), non-inferiority margin (margin), allocation ratio (allocation1),
## alpha level (one_sided_alpha), seed (seed). 


powerRMST <-
  function(n, ac_period, tot_time, tau, scale0, scale1, shape=1, margin=0, allocation1=0.5, one_sided_alpha=0.025, seed=NULL, M=20000){
    ac_rate = n / ac_period
    
    if (ac_period > tot_time){
      n0 = round(ac_rate*tot_time*(1-allocation1))
      n1 = round(ac_rate*tot_time*allocation1)
    }
    if (ac_period <= tot_time){
      n0 = round(ac_rate*ac_period*(1-allocation1))
      n1 = round(ac_rate*ac_period*allocation1)
    }
    
    ###--- test (main part) ---
    answer     = NULL
    check      = NULL
    event_arm0 = NULL
    event_arm1 = NULL
    
    if (!is.null(seed)){
      set.seed(seed)
    }
    
    for (w in 1:M){
      
      ##-- data frame --
      E             = rweibull(n0, shape, scale0)
      C             = tot_time - runif(n0, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(0,n0)
      data0         = data.frame(time, status, arm)
      ind0          = data0$status==1
      
      
      E             = rweibull(n1, shape, scale1)
      C             = tot_time - runif(n1, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(1,n1)
      data1         = data.frame(time, status, arm)
      ind1          = data1$status==1
      
      data   = rbind(data0, data1)
      data   = data[data$time>0, ]
      
      data$status[which(data$time>tau)] <- 0
      data$time[which(data$time>tau)] <- as.numeric(tau)
      
      res.km <- summary(survfit(Surv(time, status)~arm, data=data), rmean=tau, print.rmean=T)
      
      rmstC <- res.km$table[1,5]
      rmstC_SE <- res.km$table[1,6]
      rmstE <- res.km$table[2,5]
      rmstE_SE <- res.km$table[2,6]
      
      z <- qnorm(1-one_sided_alpha)
      lower <- rmstE-rmstC - z*sqrt(rmstC_SE^2+rmstE_SE^2)
      
      if (lower >= margin){
        answer[w] = 1
      } else {
        answer[w] = 0
      }
    }
    
    ###--- power ---
    power = sum(answer)/M
    
    
    ###--- output ---
    out = matrix(0,1,3)
    
    out[1,1] = n0+n1
    out[1,2] = n0
    out[1,3] = n1
    
    
    rownames(out) = c("Sample size")	
    colnames(out) = c("Total", "arm0", "arm1")
    
    Z = list()
    
    Z$result      = out
    Z$power       = power
    Z$ac_rate     = ac_rate
    Z$ac_period   = ac_period
    Z$tot_time    = tot_time
    Z$margin      = margin
    Z$tau         = tau
    
    Z
  }


# 2. Iterative process for sample size determination for target power 
## initialize sample size at 0 then add increments of 1000 until target power exceeded. Take previous iteration value then 
## continue by adding increments of 100 then increments of 10 to hone in on target power.
## This takes as input a .txt file (inputs) with 1 row per clinical trial containing 
## input variables as defined in the ssrmst_IW function above.


ninit<-0

repeat{
  if (powerRMST(n=ninit+1000, ac_period=inputs$ac_period[i], tot_time=inputs$tot_time[i], tau=inputs$tau[i], scale0=inputs$scale[i], scale1=inputs$scale[i], shape=inputs$shape[i], margin=inputs$margin_presp_conv[i], allocation1=inputs$allocation1[i], one_sided_alpha=inputs$one_sided_alpha[i])$power > inputs$design_power[i]) {break}
  ninit <- ninit+1000
}

repeat{
  if (powerRMST(n=ninit+100, ac_period=inputs$ac_period[i], tot_time=inputs$tot_time[i], tau=inputs$tau[i], scale0=inputs$scale[i], scale1=inputs$scale[i], shape=inputs$shape[i], margin=inputs$margin_presp_conv[i], allocation1=inputs$allocation1[i], one_sided_alpha=inputs$one_sided_alpha[i])$power > inputs$design_power[i]) {break}
  ninit <- ninit+100
}

repeat{
  if (powerRMST(n=ninit+10, ac_period=inputs$ac_period[i], tot_time=inputs$tot_time[i], tau=inputs$tau[i], scale0=inputs$scale[i], scale1=inputs$scale[i], shape=inputs$shape[i], margin=inputs$margin_presp_conv[i], allocation1=inputs$allocation1[i], one_sided_alpha=inputs$one_sided_alpha[i])$power > inputs$design_power[i]) {break}
  ninit <- ninit+10
}

answer$ss_margin_presp_conv[i] <- ninit+10
answer$pw_margin_presp_conv[i] <- powerRMST(n=ninit+10, ac_period=inputs$ac_period[i], tot_time=inputs$tot_time[i], tau=inputs$tau[i], scale0=inputs$scale[i], scale1=inputs$scale[i], shape=inputs$shape[i], margin=inputs$margin_presp_conv[i], allocation1=inputs$allocation1[i], one_sided_alpha=inputs$one_sided_alpha[i])$power

