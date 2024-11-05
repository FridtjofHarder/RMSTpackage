# install.packages("npsurvSS")
# library(npsurvSS)

### define parameters
accr_time <- 8
total_time <- 11
surv_cure <- 0
surv_shape1 <- 1
surv_scale1 <- 1
surv_shape2 <- 1
surv_scale2 <- .7
loss_shape <- 1
loss_scale <- 0
milestone <- 1


### create arms
arm0 <- create_arm(1, accr_time = accr_time, surv_shape = surv_shape1, surv_scale = surv_scale1, loss_shape = loss_shape, loss_scale = loss_scale, total_time = total_time)
arm1 <- create_arm(1, accr_time = accr_time, surv_shape = surv_shape2, surv_scale = surv_scale2, loss_shape = loss_shape, loss_scale = loss_scale, total_time = total_time) 

 
p1  <- arm1$size / (arm0$size + arm1$size)
p0  <- 1 - p1

theoretical_rmst <- function(x, arm) {
  stats::integrate(function(y) psurv(y, arm, lower.tail=F),
                   lower=0,
                   upper=x)$value
}

trueRMST <- function(x, arm) {
  stats::integrate(function(y) psurv(y, arm, lower.tail=F),
                   lower=0,
                   upper=x)$value
}

rmst0 <- theoretical_rmst(milestone, arm0) ### theoretical rmst
rmst1 <- theoretical_rmst(milestone, arm1) ### theoretical rmst

theoretical_delta   <- rmst0 - rmst1

### get sigma

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

### get prob_risk

prob_risk <- function(arm, teval) {
  psurv(teval, arm, lower.tail=F) *
    ploss(teval, arm, lower.tail=F) *
    paccr(pmin(arm$accr_time, arm$total_time-teval), arm)
}

### get theoretical sigma
sigma2_0 <- sigma2j_rmst(arm0, milestone)
sigma2_1 <- sigma2j_rmst(arm1, milestone)
sigma2  <- sigma2_0 / p0 + sigma2_1 / p1
out     <- (sqrt(sigma2) * stats::qnorm(1 - .05 / 2) + sqrt(sigma2) * stats::qnorm(.8) )^2 /
  theoretical_delta^2 * c(p0, p1)

### get size two arm

sides = 1
out     <- (sqrt(sigma2) * stats::qnorm(1 - .05 / sides) + sqrt(sigma2) * stats::qnorm(.8) )^2 /
  theoretical_delta^2 * c(p0, p1)

trueRSDST<-function(HORIZON,ARM){ # for a single integral
    
    VS<-function(tt){tt*(1-psurv(tt, ARM))}
    VAR<-function(H){integrate(VS,0,H, subdivisions = 1000)$value}
    
    res<-sqrt(2*VAR(HORIZON)-trueRMST(HORIZON,ARM)^2) 
    
    return(res)}
  

#TrueRSDST(1:72, arm_0)
trueRSDST <- trueRSDST(milestone, arm0)
sqrt(sigma2_0)

Hasenclevers_sigma <- rep(NaN, 90)
adjusted_sigma <- rep(NaN, 90)
sigma2_ratio <- rep(NaN, 90)
incremental_milestone <- rep(NaN, 90)
for (i in 1:90)
{
  incremental_milestone[i] <- milestone + i/10 - 0.1
  Hasenclevers_sigma[i] <- TrueRSDST(incremental_milestone[i], arm0)
  adjusted_sigma[i] <- sigma2j_rmst(arm0, incremental_milestone[i])
  sigma2_ratio[i] <- adjusted_sigma[i]/Hasenclevers_sigma[i]^2
}

plot(incremental_milestone, sigma2_ratio, main = "uncertainty dependent on administrative censoring", type = 'l', xlab = "milestone", ylab = "ratio adjusted sigma² to theorectical sigma²")
abline(v = total_time - accr_time, col = "red")
legend("bottomright", legend=c("administrative censoring starts", "increase of uncertainty"),
       col=c("red", "black"), lty=1:2, cex=0.8)

    
    