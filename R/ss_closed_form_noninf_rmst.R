#' Title
#'
#' @return
#' @export
#'
#' @examples
ss_closed_form_noninf_rmst <- function(scale_trt, shape_trt = 1, scale_ctrl,
                                       shape_ctrl = 1, parameterization = 1,
                                       accrual_time = 0,
                                       follow_up_time = NULL,
                                       tau = NULL, sides = 1,
                                       alpha = 0.025, power = 0.8, margin = 0,
                                       plot_curves = TRUE,
                                       plot_example_data = TRUE,
                                       scale_loss = 10000, shape_loss = 1)
  {
  if(is.null(tau)) stop("please specify horizon")
  if(is.null(follow_up_time)) {
    warning("follow_up_time not specified, has been set to horizon")
    follow_up_time <- tau}

    total_time <- accrual_time + follow_up_time

    if(parameterization == 2){
      scale_trt <- 1/(scale_trt^(1/shape_trt))
      scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
      scale_loss <- 1/(scale_loss^(1/shape_loss))
    }

    if(parameterization == 3){
      scale_trt <- 1/scale_trt
      scale_ctrl <- 1/scale_ctrl
      scale_loss <- 1/(scale_loss^(1/shape_loss))
    }

    arm_trt  <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                     surv_shape = shape_trt, surv_scale = 1/scale_trt,
                                     loss_scale = 1/scale_loss, loss_shape = shape_loss,
                                     follow_time = follow_up_time)
    arm_ctrl <- npsurvSS::create_arm(size = 1, accr_time = accrual_time,
                                     surv_shape = shape_ctrl,surv_scale = 1/scale_ctrl,
                                     loss_scale = 1/scale_loss, loss_shape = shape_loss,
                                     follow_time = follow_up_time)

    design  <- npsurvSS:::calc_design(arm_trt, arm_ctrl, test = list(test = "rmst difference", milestone = tau)) # get Delta and sigma2 via package npsurvSS
    sample_size <- sum((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(design$tsigma2) * stats::qnorm(power) )^2 /
                                        (design$delta + margin)^2 *c(0.5, 0.5))
    return(sample_size)

}


