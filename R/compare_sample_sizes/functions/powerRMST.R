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
