tau = 1

S0_t <- function(t) pweibull(t, shape = 1, scale = 1, lower.tail = FALSE)
S1_t <- function(t, logHR = 0) S0_t(t)^exp(logHR)

RMST0 <- function(tau = 1) integrate(S0_t, lower = 0, upper = tau)$value
RMST1 <- function(tau = 1, HR = 1) integrate(S1_t, HR = HR, lower = 0, upper = tau)$value
RMSTD <- function(tau = 1, HR = 1) RMST1(tau = tau) - RMST2 (tau = tau, HR = HR)

F_logHR <- function(tau = 1, logHR = 0) {
  # browser()
  print(S0_t(t = tau))
  print(S1_t(t = tau, logHR = logHR))
  return( 2 / sqrt( 0.5 * (1-S0_t(t = tau)) + 0.5 * (1-S1_t(t = tau, logHR = logHR))) )
  }

tS0_t <- function(t) t*S0_t(t)
tS1_t <- function(t, logHR = 0) t*S1_t(t, logHR = logHR)

var_RMST0 <- 2 * integrate(tS0_t, lower = 0, upper = tau)$value -
  integrate(S0_t, lower = 0, upper = tau)$value^2

var_RMST1 <- function(tau = 1, logHR = 0){
  2 * integrate(tS1_t, logHR = logHR, lower = 0, upper = tau)$value -
    integrate(S1_t, logHR = logHR, lower = 0, upper = tau)$value^2
}

F_RMSTD <- function(tau = 1, logHR = 0){

  sqrt( var_RMST0 + var_RMST1(tau = tau, logHR = logHR) )
}


