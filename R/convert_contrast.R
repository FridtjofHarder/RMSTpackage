#' Converts contrast to RMSTD or RMSTR
#'
#' This function returns the RMST difference or RMST ratio.
#'
#' This is a text providing some details of this function. This is \strong{strong}, this is \bold{bold}.
#'\code{code}, \preformatted{preformatted}, \kbd{keyboard-characters}. \samp{samp}
#'\verb{verb}, \pkg{package_name}, \env{environment_variable}, \option{option},\command{command_name}
#'\dfn{term}, \abbr{abbr}
#'
#' @param scale_trt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parametrization = 3}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param median_diff A scalar specifying the difference in time \eqn{\Delta t} at which each group has survival \eqn{S(t) = 0.5}. \code{median_diff} \eqn{> 0} suggests superior survival in the treatment group.
#' @param percentile_diff A scalar specifying the difference in time t \eqn{\Delta t} at which each group has the survival \eqn{S(t) =} \code{percentile}.
#' @param percentile A scalar specifying at which percentile to evaluate \code{percentile_diff}.
#' @param survival_diff A scalar specifying the survival difference in survival \eqn{\Delta S(t)} between treatment and control group.
#' @param output A string specifying whether to calculate the difference in RMSTs (\code{"RMSTD"}) or the ratio between RMSTs (\code{"RMSTR"}).
#' @param HR A scalar specifying a hazard ratio. \code{HR} \eqn{< 1} suggests a lower hazard in the treatment group than in the control group.
#' @param tau A scalar speifying the time horizon \eqn{\tau} at which to evaluate RMST with \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param t A scalar specifying the time \eqn{t} at which to evaluate \code{survival_diff}.
#' @param extended_output
#'
#' @return Returns either the difference (RMSTD) or ratio (RMSTR) in between RMSTs of treatment group and control group. An RMSTD \eqn{> 0} or an RMSTR \eqn{> 1} indicate a larger RMST in the treatment group.
#' @export
#'
#' @examples
#' RMSTD <- convert_contrast(2)
convert_contrast <- function(scale_trt = NULL, shape_trt = 1, scale_ctrl = NULL, shape_ctrl = 1,
                             parameterization = 1, HR = NULL, median_diff = NULL, percentile_diff = NULL, percentile = NULL,
                             survival_diff = NULL, t = NULL, output = "RMSTD", tau = NULL, extended_output = F){
  # convert from HR

  # throws errors if inputs are misspecified
  #stopifnot("error: more than one contrast, or no contrast were specified, or contrast specified as inappropriate data type. Please specify exactly one contrast as a scalar"
  #          = 1!=sum(is.numeric(HR, median_diff, percentile_diff, survival_diff)))
  stopifnot("output needs to be specified as either 'RMSTD' or 'RMSTR'" = output == "RMSTD"||output == "RMSTR")
  stopifnot("parameterization must be defined as either 1, 2, or 3" = parameterization == 1||parameterization == 2||parameterization == 3)
  stopifnot("please specify exactly one scale parameter" = xor(is.null(scale_trt), is.null(scale_ctrl)))

  if(!is.null(scale_trt)&&shape_ctrl!=1 || !is.null(scale_ctrl)&&shape_trt!=1){stop("please specify shape parameter only for the group for which scale has been specified")}

  if(shape_trt == 1){
    shape_trt <- shape_ctrl
  }
  if(shape_ctrl == 1){
    shape_trt <- shape_trt
  }

  if(parameterization == 2){
    scale_trt <- 1/(scale_trt^(1/shape_trt))
    scale_ctrl <- 1/(scale_ctrl^(1/shape_ctrl))
  }

  if(parameterization == 3){
    scale_trt <- 1/scale_trt
    scale_ctrl <- 1/scale_ctrl
  }

  if(!is.null(HR)){
    if(!is.null(scale_trt)){
      scale_ctrl <- scale_trt*HR
      } else{
      scale_trt <- scale_ctrl/HR
    }
  }

  if(!is.null(median_diff)){
    if(!is.null(scale_trt)){
      median_time_trt <- -log(0.5)^(1/shape_trt)*scale_trt
      median_time_ctrl <- median_time_trt - median_diff
      scale_ctrl <- median_time_ctrl*(-log(0.5))^(-1/shape_ctrl)
      } else{
      median_time_ctrl <- -log(0.5)^(1/shape_ctrl)*scale_ctrl
      median_time_trt <- median_time_ctrl + median_diff
      scale_trt <- median_time_trt*(-log(0.5))^(-1/shape_trt)
      }
  }

  if(!is.null(percentile_diff)){
    if(!is.null(scale_trt)){
      median_time_trt <- -log(percentile)^(1/shape_trt)*scale_trt
      median_time_ctrl <- median_time_trt - median_diff
      scale_ctrl <- median_time_ctrl*(-log(percentile))^(-1/shape_ctrl)
      } else{
      median_time_ctrl <- -log(percentile)^(1/shape_ctrl)*scale_ctrl
      median_time_trt <- median_time_ctrl + median_diff
      scale_trt <- median_time_trt*(-log(percentile))^(-1/shape_trt)
    }
  }

  if(!is.null(survival_diff)){
    if(!is.null(scale_trt)){
      survival_trt <- pweibull(t, scale_trt, shape_trt, lower.tail = F)
      survival_ctrl <- survival_trt - survival_diff
      if(survival_ctrl <= 0){stop("survival in treatment group minus survival difference is equal to or less than zero")}
      scale_ctrl <- t*(-log(survival_ctrl))^(-1/shape_ctrl)
    } else{
      survival_ctrl <- pweibull(t, scale_ctrl, shape_ctrl, lower.tail = F)
      survival_trt <- survival_ctrl + survival_diff
      if(survival_trt <= 0){stop("survival in control group minus survival difference is equal to or less than zero")}
      scale_trt <- t*(-log(survival_trt))^(-1/shape_trt)
    }
  }

  RMST_trt <- integrate(pweibull, shape = shape_trt, scale = scale_trt, lower = 0,
                        upper = tau, lower.tail = F)$value
  RMST_ctrl <- integrate(pweibull, shape = shape_ctrl, scale = scale_ctrl, lower = 0,
                         upper = tau, lower.tail = F)$value

  if(output == "RMSTD"){return(RMST_trt-RMST_ctrl)}
  if(output == "RMSTR"){return(RMST_trt/RMST_ctrl)}
}
