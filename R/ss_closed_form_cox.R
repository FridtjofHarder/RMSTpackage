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
#' @param one_sided_alpha
#' @param sides
#' @param power
#' @param margin_hr
#' @param parameterization
#' @param censor_beyond_tau
#'
#' @return
#' @export
#'
#' @examples
ss_closed_form_cox <- function(scale_trt, shape_trt = 1, scale_ctrl,
                               shape_ctrl = 1, parameterization = 1,
                               accrual_time = 0,
                               follow_up_time = NULL,
                               tau = Inf, sides = 1,
                               alpha = 0.025, power = 0.8, margin_hr = 1,
                               plot_curves = TRUE,
                               plot_example_data = TRUE,
                               scale_loss = 10000, shape_loss = 1)
{
  total_time <- accrual_time + follow_up_time
  #browser()
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

  if(tau>arm_ctrl$total_time){tau <- arm_ctrl$total_time}
  #browser()
  sigma2  <- stats::integrate(function(x) 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) *
                                0.5 * npsurvSS:::prob_risk(arm_trt, x) /
                                ( 0.5 * npsurvSS:::prob_risk(arm_ctrl, x) + 0.5 * npsurvSS:::prob_risk(arm_trt, x) )^2 *
                                ( 0.5 * npsurvSS:::dens_event(arm_ctrl, x) + 0.5 * npsurvSS:::dens_event(arm_trt, x)),
                              lower=0,
                              upper=tau)$value
  delta <- stats::integrate(function(x) (1 / 0.5 / npsurvSS:::prob_risk(arm_ctrl, x) + 1 / 0.5 / npsurvSS:::prob_risk(arm_trt, x)) ^ (-1) *
                              ( hsurv(x, arm_trt) - hsurv(x, arm_ctrl) ),
                            lower=0,
                            upper=tau)$value
  out     <- ( sqrt(sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(sigma2) * stats::qnorm(power) )^2 /
    delta^2 *
    c(0.5, 0.5)
  return(sum(out))
}

