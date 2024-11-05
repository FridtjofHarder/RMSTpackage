St_margin_to_RMST_diff <- function(margin, lambda, tau, DIGITS = 3){
  if (exp(-lambda*tau) + margin <= 0){
    stop("Error: S(t) - margin is <= 0, marginal S(t) can therefore not be determined.
         Please choose an earlier time horizon tau or a smaller hazard lambda.")
  }
  lambda_margin <- -log(exp(-lambda*tau) + margin)/tau
  
  S <- function(t, lambda){return(exp(-lambda * t))}
  S_tau <- S(tau, lambda)
  S_tau_margin <- S(tau, lambda_margin)
  S_tau_diff <- S_tau_margin - S_tau
  RMST_tau <- integrate(S, lambda = lambda, lower = 0, upper = tau)$value
  RMST_tau_2 <- integrate(S, lambda = lambda_margin, lower = 0, upper = tau)$value
  RMST_margin <-  RMST_tau_2 - RMST_tau
  Relative_RMST_margin <- RMST_tau_2/RMST_tau
  result <- data.frame(survival_margin=margin, S_tau_ctrl=S_tau, S_tau_margin=S_tau_margin,S_tau_diff=S_tau_diff, 
                       RMST_tau_control=RMST_tau, RMST_tau_marginal=RMST_tau_2, Diff=RMST_margin,
                       Relative_RMST_margin=Relative_RMST_margin, lambda_margin=lambda_margin)
  result<-round(result, digits=DIGITS)
  return(result)
}

get_lambda_from_St <- function(survival_at_tau, tau){
  lambda <- -log(survival_at_tau)/tau
  return(lambda)
}

margin <- -.055
tau <- 36
survival_at_tau <- .95
lambda <- get_lambda_from_St(survival_at_tau, tau)
df_results <- St_margin_to_RMST_diff(margin, lambda, tau)

tt<- (0:(100*tau))/100
plot(tt, exp(-lambda*tt), ylim=c(0,1),
     type = "l", col="red", lwd=3, xlab="time", ylab="proportion event free")
lines(tt, exp(-df_results$lambda_margin*tt),col="darkblue",lwd=3)
abline(v=tau,lty=3)

dropout_rate <- .1
Median_FU <- 66
enrollment_duration <- 130
Max_FU <- 150





