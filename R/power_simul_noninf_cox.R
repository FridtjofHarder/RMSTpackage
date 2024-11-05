#' Title
#' Attention: HR>1 indicates superior treatment group survival
#' @return
#' @export
#'
#' @examples
power_simul_noninf_cox <- function(scale_ctrl = NULL, shape_ctrl = 1, scale_trt = NULL,
                                    shape_trt = 1, scale_loss = NULL, shape_loss = 1,
                                    accrual_time = NULL, follow_up_time = NULL,
                                    tau = NULL, one_sided_alpha = 0.025, sides = 1, power = 0.8,
                                    margin_hr = 1, n = NULL, M = 100,
                                    parameterization = 1, censor_beyond_tau = FALSE){
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

  for (i in 1:M){
    #browser()
    events <- c(rweibull(n_per_group, scale = scale_ctrl, shape = shape_ctrl),
                rweibull(n_per_group, scale = scale_trt, shape = shape_trt))

    loss_to_fu  <- rweibull(2*n_per_group, scale = scale_loss, shape = shape_loss)
    admin_loss  <- total_time - runif(2*n_per_group, max = accrual_time)

    status <- as.numeric(events <= loss_to_fu & events <= admin_loss)
    time <- pmin(events, loss_to_fu, admin_loss)
    arm <- c(rep(0, n_per_group), rep(1, n_per_group)) #  0 = ctrl, 1 = trt

    if(censor_beyond_tau){
      status[time>tau] <- 0
    }

    data <- data.frame(time, status, arm)

    result <-  coxph(Surv(time, status) ~ arm, data = data)
    #browser()
    upper <-  summary(result)$conf.int[, "upper .95"]

    if(upper < margin_hr){positive[i] <- 1}
    else positive[i] <- 0
  }

  print(sum(is.na(positive)))
  power <-  sum(positive, na.rm = TRUE)/sum(!is.na(positive))
  return(power)

}


