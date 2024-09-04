#' Title
#'
#' @param scale_ctrl
#' @param shape_ctrl
#' @param scale_trt
#' @param shape_trt
#' @param scale_loss
#' @param shape_loss
#' @param accrual_time
#' @param follow_up_time
#' @param total_time
#' @param tau
#' @param alpha
#' @param sides
#' @param power
#' @param margin_rmst
#' @param n_RMST
#' @param n_LR
#' @param M
#'
#' @return
#' @export
#'
#' @examples
power_simul_noninf_rmst <- function(scale_ctrl = NULL, shape_ctrl = 1, scale_trt = NULL,
                            shape_trt = 1, scale_loss = NULL, shape_loss = 1,
                            accrual_time = NULL, follow_up_time = NULL,
                            tau = NULL, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                            margin_rmst = 0, n = NULL, M = 100,
                            parameterization = 1, handling_large_tau = "adjust tau"){
  positive <- rep(NA, M)

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
    scale_loss <- 1/(scale_loss^(1/shape_loss))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
    scale_loss <- 1/scale_loss
  }

  n_per_group <- round(n/2)
  total_time <- accrual_time + follow_up_time
  counter_large_tau <- 0

  for (i in 1:M){

    events <- c(rweibull(n_per_group, scale = scale_ctrl, shape = shape_ctrl),
                rweibull(n_per_group, scale = scale_trt, shape = shape_trt))

    loss_to_fu  <- rweibull(2*n_per_group, scale = scale_loss, shape = shape_loss)
    admin_loss  <- total_time - runif(2*n_per_group, max = accrual_time)

    status <- as.numeric(events <= loss_to_fu & events <= admin_loss)
    time <- pmin(events, loss_to_fu, admin_loss)
    arm <- c(rep(0, n_per_group), rep(1, n_per_group)) #  0 = ctrl, 1 = trt
    data <- data.frame(time, status, arm)
    #browser()
    if(max(time[arm == 0]) <= tau || max(time[arm == 1]) <= tau){
      counter_large_tau <- counter_large_tau + 1
      if(handling_large_tau == "adjust tau") {tau <-  min(max(time[arm == 0]), max(time[arm == 1]))}
      if(handling_large_tau == "discard") next
    }

    result <-  rmst2(data$time, data$status, data$arm, tau = tau,
                 alpha = one_sided_alpha * 2)$unadjusted.result
    lower <-  result[1, 2]
    positive[i] <- if(lower > -margin_rmst){positive[i] <- 1}
    else positive[i] <- 0
  }

  print(sum(is.na(positive)))
  print(counter_large_tau)
  power <-  sum(positive, na.rm = TRUE)/sum(!is.na(positive))
  return(power)

}


