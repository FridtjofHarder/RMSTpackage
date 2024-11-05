rm(list=ls())
source("RMSTtools.R")

scale_0 <- 1/30
scale_1 <- 1/20
arm_0 <- create_arm(size = 100, 
                    accr_time = 1, 
                    total_time = 2,
                    surv_cure = 0,
                    surv_scale = scale_0, 
                    surv_shape = 1,
                    loss_shape = 1, 
                    loss_scale = .01)

arm_1 <- create_arm(size = 1, 
                    accr_time = 36, 
                    follow_time = 36,
                    surv_cure = 0,
                    surv_scale = scale_1, 
                    surv_shape = 1,
                    loss_shape = 1, 
                    loss_scale = .1)


IllustrateTrial(arm_0, OBS = F)
IllustrateTrial(arm_0, OBS = T)

hh=1:72

# inverse cumulated distribution of ideal observation times
B<-function(tt, armC){ paccr(armC$total_time-tt, armC)}
plot(hh,B(hh, arm_0), type = "l", lwd=2,
     xlab="Time horizon h",
     ylab="Proportion under observation (ignoring events)")



# inverse cumulated distribution of time to drop-out
D<-function(tt, armC){ 1-ploss(tt,armC) }
lines(hh, D(hh, arm_0), col =2, lwd=2)

# inverse cumulated distribution of actual observatio times (ignoring delays)
Fobs<-function(tt, armC){ (B(tt,armC)) *(D(tt,armC)) }
lines(hh, Fobs(hh, arm_0), col =4, lwd=2)
# Checked: same as in IllustrateTrial(arm_0, OBS = T)
legend("bottomleft",c("Ideal Observation", "Drop-outs", "Combined"),
       lwd=2, col=c(1,2,4))

# Calculate the density of Fobs from the underlying densities
# of B and D

plot(hh, daccr(arm_0$total_time-hh, arm_0), type ="l", lwd=2)
lines(hh, dloss(hh, arm_0), lwd=2, col=2)

# Differentiate 1-Fobs
dobs<-function(tt,armC){ 
  daccr(armC$total_time-tt, armC) + 
    dloss(tt, armC) -
    daccr(armC$total_time-tt, armC) * (1-D(tt, armC)) -
    dloss(tt, armC) * (1-B(tt, armC))
  }
plot(hh, dobs(hh, arm_0))
#dobs(60,arm_0)
# Check should integrate to 1
# XX<-function(hh){dobs(hh, arm_0) }
# integrate(XX,0,72) # Nice!




# E_ideal
plot(hh,psurv(hh,arm_0), ylim=c(0,1),type="l",col=4, lwd=2,
     xlab="Time horizon h",
     ylab="Proportion of patients with event")
# E_observed in censored patients
H<-function(h){ (psurv(h,arm_0))*dobs(h, arm_0)}
E_cens<-function(h){integrate(H,0,h, subdivisions = 100)$value}
E_cens<-Vectorize(E_cens)
E_cens(c(34,72))

H_1<-function(h){ (psurv(h,arm_1))*dobs(h, arm_1)}
E_cens_1<-function(h){integrate(H_1,0,h, subdivisions = 100)$value}
E_cens_1<-Vectorize(E_cens_1)
E_cens_1(c(34,72))

# E_observed
lines(hh,  E_cens(hh) + psurv(hh,arm_0)*Fobs(hh,arm_0), col=2, lwd=2)
legend("bottomright", c("Proportion events without censoring",
                        "Proportion observed events with censoring"),
       lwd=2, col=c(4,2))


# proportion Events expected to be observed / Events with ideal observation
plot(hh,  (E_cens(hh) + psurv(hh,arm_0)*Fobs(hh,arm_0))/psurv(hh,arm_0), 
     col=3, lwd=2, type="l",ylim=c(0,1),
     xlab="Time horizon h",
     ylab="Proportion of observable events")

# ?? 1/ Inflation factor

### plot n acc to npsurvSS
n_npsurvSS <- rep(NA, length(hh))
for (i in 1:length(hh)){
n_npsurvSS[i] <- size_two_arm(arm_0, arm_1, test = list(test = "rmst ratio", milestone = hh[i]))[[1]]}

plot(hh, n_npsurvSS, ylim = c(0, 5000), type = "l")
legend("topright", legend = c("npsurvSS", "R+P no censoring", "R+P modified by observable proportion"), col = c("black", "red", "blue"), lty = 1)

ratio_0 <- (E_cens(hh) + psurv(hh,arm_0)*Fobs(hh,arm_0))/psurv(hh,arm_0)
ratio_1 <- (E_cens_1(hh) + psurv(hh,arm_1)*Fobs(hh,arm_1))/psurv(hh,arm_1)

### plot n without censoring acc to R+P
n_RP_no_cens <- rep(NA, length(hh))
n_RP_cens <- rep(NA, length(hh))
for (i in 1:length(hh)){
  true_RSTSD_0 <- TrueRSDST(hh[i], arm_0)
  true_RSTSD_1 <- TrueRSDST(hh[i], arm_1)
  ratio_0 <- (E_cens(hh) + psurv(hh,arm_0)*Fobs(hh,arm_0))/psurv(hh,arm_0)
  ratio_1 <- (E_cens_1(hh) + psurv(hh,arm_1)*Fobs(hh,arm_1))/psurv(hh,arm_1)
  n_RP_no_cens[i] <- ((qnorm(1-.05/2) + qnorm(.8))^2)*(true_RSTSD_0^2 + true_RSTSD_1^2)/(TrueRMST(hh[i], arm_0) - TrueRMST(hh[i], arm_1))^2
  n_RP_cens[i] <- ((qnorm(1-.05/2) + qnorm(.8))^2)*(true_RSTSD_0^2*1/ratio_0[i]^2 + true_RSTSD_1^2*1/ratio_1[i]^2)/(TrueRMST(hh[i], arm_0) - TrueRMST(hh[i], arm_1))^2
}
lines(hh, n_RP_no_cens, type = "l", col = "red")
lines(hh, n_RP_cens, type = "l", col = "blue")






### get RSTSD acc to R+P