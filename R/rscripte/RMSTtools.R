# RMST Tools
#
# collection of function to deal with RMST
######################################################################
require(OurTools)
require(npsurvSS)
require(survival)
require(geepack)
require(pseudo)


# Defining a trial design and simulating survival data is based on the
# package "npsurvSS"

#------------------------------------------------------------------------
# IllustrateTrial illustrated tral arms defined with "npsurvSS"
#   - true survival curve(s)
#   - true distribution(s) of observation time
#------------------------------------------------------------------------

IllustrateTrial<-function(armC, armE = NULL, OBS = FALSE){
  
  tt<-seq(0,armC$total_time,length.out = 200)
  # Ticks on Time axis
  Tmax<-pretty(tt, n=10)
  
  if (!OBS){ # Survival function
    plot(tt,1-psurv(tt,armC),type = "l", 
         lwd=2, col = 4, las=1, cex.axis = 1.3, cex.lab = 1.3,
         xaxp=c(0,Tmax[length(Tmax)],length(Tmax)-1),
         ylim=c(0,1), yaxp=c(0,1,10),
         xlab = "Time",ylab = "S(t)")
    
    if ( !is.null(armE) ) {
      lines( tt, 1-psurv(tt,armE),lwd=2, col = 2) 
      legend("bottomleft",c("armE", "armC"), lwd=3, col = c(4,2))
    }
  } else { # Distribution of observation times
    
    plot(armC$total_time-tt,paccr(tt,armC),type = "l", main = "Distribution of observation times",
         lwd=2, col = 1, las=1, cex.axis = 1.5, cex.lab = 1.5, 
         xaxp=c(0,Tmax[length(Tmax)],length(Tmax)-1),
         ylim=c(0,1), yaxp=c(0,1,10),
         xlab = "Time",ylab = "Proportion of patients under observation")
    lines(armC$total_time-tt,
          paccr(tt,armC)*(1-ploss(armC$total_time-tt,armC)), 
          col=4, lwd=2)
    
    if ( !is.null(armE) ) {
      lines( armE$total_time-tt, 
             paccr(tt,armE)*(1-ploss(armE$total_time-tt,armE)),
             lwd=2, col = 2)
      
      legend("bottomleft",c("Ideal without loss","armC with loss", "armE with loss"), lwd=3, col = c(1,4,2))
    }
  }
  
  
}



# IllustrateTrial(arm_0, OBS = F)
# IllustrateTrial(arm_0, OBS = T)
# IllustrateTrial(arm_0, arm_1, OBS = F)
# IllustrateTrial(arm_0, arm_1, OBS = T)


#------------------------------------------------------------------------
# Calculate the true RMST (check if there is better function in package?)
TrueRMST<-function(HORIZON,ARM){
  trueRMST<-function(HORIZON,ARM){ # for a single integral
    hh<-function(tt){1-psurv(tt, ARM)}
    res<-integrate(hh, lower = 0, upper = HORIZON)$value
    return(res)
  }
  # Vectorize
  res<-numeric(0)
  for ( i in 1:length(HORIZON) ) res<-c(res,trueRMST(HORIZON[i],ARM))
  return(res)
}

# TrueRMST(arm_0$total_time,arm_0)
# TrueRMST(arm_0$total_time,arm_1)
#------------------------------------------------------------------------
# Calculate the true RSDST (check if there is better function in package?)
TrueRSDST<-function(HORIZON,ARM){
  trueRMST<-function(HORIZON,ARM){ # for a single integral
    hh<-function(tt){1-psurv(tt, ARM)}
    res<-integrate(hh, lower = 0, upper = HORIZON)$value
    return(res)
  }
  trueRSDST<-function(HORIZON,ARM){ # for a single integral

    VS<-function(tt){tt*(1-psurv(tt, ARM))}
    VAR<-function(H){integrate(VS,0,H, subdivisions = 1000)$value}
    
    res<-sqrt(2*VAR(HORIZON)-trueRMST(HORIZON,ARM)^2) 

    return(res)
  }
  # Vectorize
  res<-numeric(0)
  for ( i in 1:length(HORIZON) ) res<-c(res,trueRSDST(HORIZON[i],ARM))
  return(res)
}
#TrueRSDST(1:72, arm_0)




#------------------------------------------------------------------------
# illustrate RMSTs as a function of the time horizon
#  or the difference of the ratio
#------------------------------------------------------------------------

IllustrateRMST<-function(armC, armE, TYPE = "RMSTdiff"){
  
  if (TYPE == "RMST"){
    Ymax<-max(TrueRMST(armC$total_time,armC),TrueRMST(armE$total_time,armE))
    Ymax<-pretty(c(0,Ymax), n=10)
    tt<-seq(.2,armC$total_time,by = .2)
    Tmax<-pretty(c(0,tt), n=10)
    plot( c(0,tt), c(0,TrueRMST(tt,armC)),
          type = "l",lwd=3,col=4,las=1, cex.axis = 1.3, cex.lab = 1.3,
          xaxp=c(0,Tmax[length(Tmax)],length(Tmax)-1),
          ylim=c(0,Ymax[length(Ymax)]), yaxp=c(0,Ymax[length(Ymax)],length(Ymax)-1),
          ylab=" Restricted Mean Survival Time", xlab = "Time horizon")
    lines( c(0,tt), c(0,TrueRMST(tt,armE)),lwd=3,col=2)
    abline(v=armC$total_time, lty = 2, lwd=2)
  }
  
  if (TYPE == "RMSTdiff") {   
    Ymax<-max(abs(TrueRMST(armC$total_time,armC)-TrueRMST(armE$total_time,armE)))
    Ymax<-pretty(c(0,Ymax), n=6)
    tt<-seq(.2,armC$total_time,by = .2)
    Tmax<-pretty(c(0,tt), n=10)
    plot( c(0,tt), c(0, TrueRMST(tt,armE)) - c(0,TrueRMST(tt,armC)),
          type = "l",lwd=3,col=4,las=1, cex.axis = 1.3, cex.lab = 1.3,
          xaxp=c(0,Tmax[length(Tmax)],length(Tmax)-1),
          ylim=c(-Ymax[length(Ymax)],Ymax[length(Ymax)]), 
          yaxp=c(-Ymax[length(Ymax)],2*Ymax[length(Ymax)],length(Ymax)-1),
          ylab="RMST difference", xlab = "Time horizon")
    abline(h = 0, lty = 2, lwd=2)
    
  }
  if (TYPE == "RMSTratio"){
    Ymax<-max((TrueRMST(armC$total_time,armE)/TrueRMST(armE$total_time,armC)))
    Ymax<-pretty(c(2-Ymax,Ymax), n=6)
    tt<-seq(.2,armC$total_time,by = .2)
    Tmax<-pretty(c(0,tt), n=10)
    plot( c(0,tt), c(1, TrueRMST(tt,armE)) / c(1,TrueRMST(tt,armC)),
          type = "l",lwd=3, col=4, las=1, cex.axis = 1.3, cex.lab = 1.3,
          xaxp=c(0,Tmax[length(Tmax)],length(Tmax)-1),
          ylim=c(2-Ymax[length(Ymax)],Ymax[length(Ymax)]), 
          yaxp=c(2-Ymax[length(Ymax)],Ymax[length(Ymax)],length(Ymax)-1),
          ylab="RMST ratio", xlab = "Time horizon")
    abline(h = 1, lty = 2, lwd=2)
  }
}

# IllustrateRMST(arm_0, arm_1, TYPE = "RMST")
# IllustrateRMST(arm_0, arm_1, TYPE = "RMSTdiff")
# IllustrateRMST(arm_0, arm_1, TYPE = "RMSTratio")



#------------------------------------------------------------------------
# Simulate data from a npsurvSS defined arm
#------------------------------------------------------------------------

SimulateArm<-function(arm, N, PLOT = T){
  Arm0<-simulate_arm(arm)
  for ( i in 2:N){Arm0<-rbind(Arm0,simulate_arm(arm))}
  Arm0$Time<- pmin(Arm0$time.obs, Arm0$time.surv, 
                   Arm0$time.loss, arm$total_time)
  Arm0obj<-Surv(Arm0$Time, Arm0$censor)
  if (PLOT){plot(Arm0obj, mark.time = T, lwd=2)
    tt<-seq(0,arm$total_time,length.out = 200)
    lines( tt, 1-psurv(tt,arm),lwd=2, col = 2) }
  return(Arm0obj)
}

#------------------------------------------------------------------------
####  Functions to extract and compare restricted mean survival times
####  Based on pseudo-value method as recomended by Royston & Parmer
###   DH 2014-08-15

# require(geepack)
# require(pseudo)

# ####  using Survfit
RestMeanSurv<-function(Sob,t0){
  res<-c(NA , NA)
  try (m1<-summary(survfit(Sob~1), print.rmean=TRUE,rmean=t0), silent = T)
  if (exists("m1")) res<-m1$table[5:6]
  return(res)
}
#RestMeanSurv(Sob2,12)

#### using pseudo-values and geese
RestrictedMeanSurv<-function(Sob,t0){
  futime<-Sob[,1]
  fustat<-Sob[,2]
  pseudo = pseudomean(time=futime, event=fustat,tmax=t0)
  id<-1:length(futime)
  mod<-summary(fit <- geese(pseudo ~ 1,id=id,
                            jack = TRUE, family=gaussian, 
                            corstr="independence", scale.fix=FALSE))
  return(mod$mean[c(1,3)])
}  
#RestrictedMeanSurv(Sob1,12)

#### using pseudo-values and geese

CompareRestrictedMeanSurv<-function(Sob1,Sob2,t0){
  futime<-c(Sob1[,1],Sob2[,1])
  fustat<-c(Sob1[,2],Sob2[,2])
  gr<-c(rep(0,dim(Sob1)[1]),rep(1,dim(Sob2)[1]))
  pseudo = pseudomean(time=futime, event=fustat,tmax=t0)
  id<-1:length(futime)
  mod<-summary(fit <- geese(pseudo ~ gr,id=id,
                            jack = TRUE, family=gaussian, 
                            corstr="independence", scale.fix=FALSE))
  
  res<-unlist(c(    (mod$mean[2,c(1,3)]),
                    (mod$mean[2,4:5])))
  return(res)
}
#CompareRestrictedMeanSurv(Sob1,Sob2,12)

