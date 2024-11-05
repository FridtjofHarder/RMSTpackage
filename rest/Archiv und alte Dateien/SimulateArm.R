# SimulateArm

arm_0 <- create_arm(size = 1, 
                    accr_time = 12, 
                    follow_time = 60,
                    surv_cure = .3,
                    surv_scale = 1/30, 
                    surv_shape = 2,
                    loss_shape = 1, 
                    loss_scale = .002)

IllustrateTrial(arm_0, OBS = F)
IllustrateTrial(arm_0, OBS = T)

hh<-12:72

#TrueRMST(hh, arm_0)

#TrueRSDST(hh, arm_0)

# SS<-smooth.spline(c(0,Fit$time), c(1,Fit$surv))
# S<-function(tt){predict(SS, tt, deriv = 0)$y}


RMST<-function(H) {integrate(S, 0, H, subdivisions = 1000)$value}
RMST<-Vectorize(RMST)

# Make t* KM function
VS<-function(t){t*S(t)}
VAR<-function(H){integrate(VS,0,H, subdivisions = 1000)$value}

RSDST<-function(H){sqrt(2*VAR(H)-RMST(H)^2) }
RSDST<-Vectorize(RSDST)


# Estimates directly from data

Empir<-function(H, N){
  ff<-RestMeanSurv(Obj,H)
  ff[2]<-ff[2]*sqrt(N)
  return(c(as.numeric(ff),as.numeric(ff[2]/ff[1])))
}
Empir<-Vectorize(Empir)



EmpirP<-function(H, N){
  ff<-RestrictedMeanSurv(Obj,H)
  ff[2]<-ff[2]*sqrt(N)
  return(c(as.numeric(ff),as.numeric(ff[2]/ff[1])))
}
EmpirP<-Vectorize(EmpirP)


Nevent<-function(h){ sum(Fit$n.event[Fit$time <= h])}
Nevent<-Vectorize(Nevent)
Nevent(hh)



par(mfrow=c(2,1))
N=500
Obj<-SimulateArm(arm = arm_0, N = N, PLOT = T) 
Fit<-survfit(Obj~1)



# MakeKMfunction

SS<-smooth.spline(c(0,Fit$time), c(1,Fit$surv))
S<-function(tt){predict(SS, tt, deriv = 0)$y}

lines(hh, S(hh), col = 4, lwd=3)
#S1<-splinefun(c(0,Fit$time), c(1,Fit$surv),method="hyman")
# lines(1:72, S1(1:72), col=3, lwd=3)
lines(Fit$time, cumsum(Fit$n.event)/N)
#plot(Fit$time, cumsum(Fit$n.event))

plot(hh,Empir(hh,N)[2,]/(RSDST(hh)), type = "l", ylim=c(.8,1.2), xlim=c(0,72))
abline(h=1)
lines(hh,Empir(hh,N)[2,]/(TrueRSDST(hh, arm_0)), col = 4 )



#--------------  compare SE estimates  
rr<- unlist(c(RestrictedMeanSurv(Obj, 12) [2], RestMeanSurv(Obj, 12)[2]))
for (h in hh){
 rr<- rbind(rr,
            unlist(c(RestrictedMeanSurv(Obj, h) [2], RestMeanSurv(Obj, h)[2]))
 )
}
rr<-rr[-1,]

NiceHist(rr[,1]-rr[,2])  # geese slightly larger
#---------------------------------------------------------





# TODO Simulation How many events needed for acceptable estimate?
hh<-(3:18)*4
Nsim=1000
Nsample=c( 10*(2:9), 100*(1:5))
SS<-smooth.spline(c(0,Fit$time), c(1,Fit$surv))
S<-function(tt){predict(SS, tt, deriv = 0)$y}
res<-data.frame( Sim= rep(0,length(hh)), N = rep(N,length(hh)), h=hh, Nevent=Nevent(hh),ERMST=E[1,],ERSDST=E[2,],
                 RMST = RMST(hh), RSDST = RSDST(hh), TrueRMST = TrueRMST(hh, arm_0))

for (i in 1:Nsim){
  
  for (N in Nsample){
    Obj<-SimulateArm(arm = arm_0, N = N, PLOT = F) 
    Fit<-survfit(Obj~1)
    SS<-smooth.spline(c(0,Fit$time), c(1,Fit$surv))
    S<-function(tt){predict(SS, tt, deriv = 0)$y}
    E<-Empir(hh,N)
    
    res<-rbind(res,
      data.frame( Sim= i, N = rep(N,length(hh)), h=hh, Nevent=Nevent(hh),
                  ERMST=E[1,],ERSDST=E[2,],
                  RMST = RMST(hh), RSDST = RSDST(hh), TrueRMST = TrueRMST(hh, arm_0))
    )
  }
}

res<-res[res$Sim > 0, ]

save(res, file = "RMSTbyNbyHsimulation.Rdata")

load(file = "RMSTbyNbyHsimulation.Rdata")

Nsample=c( 10*(2:9), 100*(1:5))
colnames(res)

InAB<-function(X,L,U){ !is.na(X) & !is.na(L) & !is.na(U) & X >= L & X <= U}

res$ESE<-res$ERSDST/sqrt(res$N)
res$Eupper<-res$ERMST + qnorm(0.975) * res$ESE
res$Elower<-res$ERMST - qnorm(0.975) * res$ESE
res$Ecover<-InAB(res$TrueRMST, res$Elower, res$Eupper)

res$EupperT<-res$ERMST + qt(0.975, df = res$Nevent-1) * res$ESE
res$ElowerT<-res$ERMST - qt(0.975, df = res$Nevent-1) * res$ESE
res$EcoverT<-InAB(res$TrueRMST, res$ElowerT, res$EupperT)

#----------------------------------------------------------------

yy<-aggregate(res$Ecover, by = list(res$N), FUN = mean)
colnames(yy)<-c("N","PCover")

yy1<-aggregate(res$EcoverT, by = list(res$N), FUN = mean)
colnames(yy1)<-c("N","PCoverT")

ee<-aggregate(res$Nevent, by = list(res$N), FUN = mean)
colnames(ee)<-c("N","MeanEvent")

yyy<-MergeWithMaster(yy,yy1, BY = "N")
yyy<-MergeWithMaster(yyy,ee, BY = "N")

plot(yyy$N, yyy$PCover, log = "x", ylim=c(.90,1), pch=19, type = "b")
abline(h=.95)
points(yyy$N, yyy$PCoverT, pch=20, col=4)
lines(yyy$N, yyy$PCoverT, pch=20, col=4)


#----------------------------------------------------------------------
zz<-aggregate(res$Ecover, by = list(res$N, res$h ), FUN = mean)
colnames(zz)<-c("N", "h" , "PCover")

zz1<-aggregate(res$EcoverT, by = list(res$N, res$h ), FUN = mean)
colnames(zz1)<-c("N", "h","PCoverT")

zee<-aggregate(res$Nevent, by = list(res$N, res$h ), FUN = mean)
colnames(zee)<-c("N", "h","MeanEvent")

zzz<-merge(zz,zz1)
colnames(zzz)
zzz<-merge(zzz,zee)

# plot(zzz$MeanEvent, zzz$PCoverT)
# points(zzz$MeanEvent, zzz$PCover, col = 2, pgh=19)
# abline(h = .95)



plot(zzz$h[ zzz$N == 20 ], zzz$PCover[ zzz$N == 20 ], 
     col = 0, type = "l", ylim=c(.7,1), main = "normal approximation",
     xlab = "time horizon h", ylab=" estimated covering probability")
abline(h = .95)
cc=1
for ( n in Nsample) {
  lines(zzz$h[zzz$N == n ],zzz$PCover[zzz$N == n ], col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)


plot(zzz$h[ zzz$N == 20 ], zzz$PCoverT[ zzz$N == 20 ], 
     col = 0, type = "l", ylim=c(.7,1), main = "t df=E-1 approximation",
     xlab = "time horizon h", ylab=" estimated covering probability")
abline(h = .95)
cc=1
for ( n in Nsample) {
  lines(zzz$h[zzz$N == n ],zzz$PCoverT[zzz$N == n ], col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)



plot( 5:300, qt(0.975, df=5:300), type = "l", col = 4)
abline(h=qnorm(0.975), col = 2, lwd = 2)

plot( 5:100, qt(0.975, df=5:100)/qnorm(0.975), 
      type = "l", col = 4, ylim=c(1, 1.3))
abline(h = 1, col=2, lwd = 2)
#-----------------------------------------------------------------
# Bias

bb<-aggregate((res$ERMST - res$TrueRMST)/res$TrueRMST, by = list(res$N, res$h ), FUN = mean)
colnames(bb)<-c("N", "h", "Ebias")
bb1<-aggregate((res$RMST - res$TrueRMST)/res$TrueRMST, by = list(res$N, res$h ), FUN = mean)
colnames(bb1)<-c("N", "h", "Sbias")
bbb<- merge(bb, bb1)


#par(mfrow=c(2,1))
plot(bbb$h[ bbb$N == 20 ], bbb$Ebias[ bbb$N == 20 ], 
     col = 0, type = "l", ylim=c(-.1,.02), main = "E bias",
     xlab = "time horizon h", ylab="e Bias")
abline(h = 0)
cc=1
for ( n in Nsample) {
  lines(bbb$h[ bbb$N == n ], bbb$Ebias[ bbb$N == n ], col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)

plot(bbb$h[ bbb$N == 20 ], bbb$Sbias[ bbb$N == 20 ], 
     col = 0, type = "l", ylim=c(-.3,.2), main = "S bias",
     xlab = "time horizon h", ylab="e Bias")
abline(h = 0)
cc=1
for ( n in Nsample) {
  lines(bbb$h[ bbb$N == n ], bbb$Sbias[ bbb$N == n ], col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)

#  Spline through KM not good!
#-----------------------------------------------------------------
# MSE

bb<-aggregate(((res$ERMST - res$TrueRMST)/res$TrueRMST)^2, by = list(res$N, res$h ), FUN = mean)
colnames(bb)<-c("N", "h", "E_MSE")
bb1<-aggregate(((res$RMST - res$TrueRMST)/res$TrueRMST)^2, by = list(res$N, res$h ), FUN = mean)
colnames(bb1)<-c("N", "h", "S_MSE")
bbb<- merge(bb, bb1)


#par(mfrow=c(2,1))
plot(bbb$h[ bbb$N == 20 ], sqrt(bbb$E_MSE[ bbb$N == 20 ]), 
     col = 0, type = "l",  main = "E_MSE", ylim=c(-.3,.2),
     xlab = "time horizon h", ylab="E_MSE")
abline(h = 0)
cc=1
for ( n in Nsample) {
  lines(bbb$h[ bbb$N == n ], sqrt(bbb$E_MSE[ bbb$N == n ]), col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)

plot(bbb$h[ bbb$N == 20 ], sqrt(bbb$S_MSE[ bbb$N == 20 ]), 
     col = 0, type = "l", ylim=c(-.3,.2), main = "S_MSE",
     xlab = "time horizon h", ylab="S_MSE")
abline(h = 0)
cc=1
for ( n in Nsample) {
  lines(bbb$h[ bbb$N == n ], sqrt(bbb$S_MSE[ bbb$N == n ]), col = cc)
  cc=cc+1
}
legend("bottomright",paste("N =" ,rev(Nsample)), col=13:1, lwd=2)

