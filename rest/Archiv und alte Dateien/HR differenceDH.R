internal_function <- function(t, HR_margin, shape, scale){
  result <- exp(-HR_margin*(t/scale)^shape)
  return(result) ## returns S(t) for given t, HR_margin, shape, and scale
}
internal_function<-Vectorize(internal_function)

HR_margin_to_RMST_diff <- function(HR_margin, shape, scale, tau, 
                                   PLOT = T, DIGITS=3){
  if(PLOT){
    tt<- (0:(100*tau))/100
    plot(tt, internal_function (tt, HR_margin, shape, scale), ylim=c(0,1),
         type = "l", col="red", lwd=3, xlab="time", ylab="proportion event free")
    lines(tt, internal_function (tt, 1, shape, scale),col="darkblue",lwd=3)
    abline(v=tau,lty=3)
  }
  
  TauArm1<-internal_function (tau, HR_margin, shape, scale) # return S(tau)
  TauArm2<-internal_function (tau, 1, shape, scale) # return S(tau)
  TauDiff<-TauArm1-TauArm2 # return S(t)_0 - S(t)_1
  Arm1<-integrate(internal_function, HR_margin, shape, scale, lower = 0,upper = tau)$value # return RMST for tau
  Arm2<-integrate(internal_function, HR_margin = 1, shape, scale, lower = 0,upper = tau)$value # return RMST for tau
  Diff<-Arm1-Arm2 # return RMSTD
  RelDiff<-Diff/max(Arm1,Arm2)
  result <- data.frame( HR_margin=HR_margin, TauArm1=TauArm1, TauArm2=TauArm2,TauDiff=TauDiff, 
                        RMST1=Arm1, RMST2=Arm2, Diff=Diff,RelDiff=RelDiff)
  result<-round(result, digits=DIGITS)
  return(result)
}


HR_margin_to_RMST_diff(log(.93-.055)/log(.93),1.2,26.7,10)

log(.93-.055)/log(.93)

