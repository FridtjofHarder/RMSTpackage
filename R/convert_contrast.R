#' Converts contrast to RMSTD or RMSTR
#'
#' This function converts different contrasts into either the difference or the ratio of restricted mean surival times (RMSTs).
#' Convertible contrasts are: \itemize{
#' \item hazard ratio (HR).
#' \item median difference: the difference in times \eqn{t} at which each survival curve passes \eqn{S(t) = 0.5}.
#' \item percentile difference: the difference in times \eqn{t} at which each survival curve passes a certain \eqn{S(t)}. Simplifies to the median difference when the percentile \eqn{0.5} is chosen. TODO: discuss the term "percentile"! Should multiply by 100?]
#' \item survival difference: difference between two survival curves \eqn{S(t)} at a specified time \eqn{t}.}
#' Survival is assumed to follow Weibull distributions in both control and treatment group with a constant HR, implying the same shape parameter in treatment and control group.
#'
#'ToDo:
#'Output als df - besprechen mit DH
#'#'
#' Function details (missing)
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
#' @param RMSTD A scalar specifying the RMSTD between control group and treatment group. Allows for converting RMSTD to HR.
#' @param RMSTR A scalar specifying the RMSTR between control group and treatment group. Allows for converting RMSTD to HR.
#' @param plot_curves specifying whether to create a plot. Boolean.
#'
#' @return Returns either the difference (RMSTD) or ratio (RMSTR) between RMSTs of treatment group and control group. An RMSTD \eqn{> 0} or an RMSTR \eqn{> 1} indicate a larger RMST in the treatment group.
#' @import stats
#' @export
#'
#' @examples
#' RMSTD <- convert_contrast(scale_trt = 1, shape_trt = 1, HR = 0.9,
#' output = "RMSTD", tau = 1)
#' RMSTR <- convert_contrast(scale_trt = 1, shape_trt = 1, survival_diff = 0.01,
#' t = 1, output = "RMSTR", tau = 1)
#' RMSTR <- convert_contrast(scale_trt = 1, shape_trt = 1, survival_diff = 0.01,
#' t = 1, output = "RMSTR", tau = 1)
convert_contrast <- function(scale_trt = NULL,
                             shape_trt = 1,
                             scale_ctrl = NULL,
                             shape_ctrl = 1,
                             parameterization = 1,
                             HR = NULL,
                             median_diff = NULL,
                             percentile_diff = NULL,
                             percentile = NULL,
                             survival_diff = NULL,
                             t = NULL,
                             RMSTD = NULL,
                             RMSTR = NULL,
                             output = "RMSTD",
                             tau = NULL,
                             plot_curves = TRUE) {
  # convert from HR
  # throws errors if inputs are misspecified
  #browser()
  stopifnot(
    "error: more than one contrast, or no contrast were specified, or contrast specified is of inappropriate data type. Please specify exactly one contrast as a scalar"
    = 1 == sum(
      !is.null(HR),
      !is.null(median_diff),
      !is.null(percentile_diff),
      !is.null(survival_diff),
      !is.null(RMSTR),
      !is.null(RMSTD)
    )
  )

  stopifnot(
    "output needs to be specified as either 'RMSTD' or 'RMSTR' or 'HR'" = output == "RMSTD" ||
      output == "RMSTR" || output == "HR"
  )
  stopifnot(
    "parameterization must be defined as either 1, 2, or 3" = parameterization == 1 ||
      parameterization == 2 || parameterization == 3
  )
  stopifnot("please specify exactly one scale parameter" = xor(is.null(scale_trt), is.null(scale_ctrl)))

  # throws error if shape parameter is specified for the group, for which scale is unspecified
  if (!is.null(scale_trt) &&
      shape_ctrl != 1 ||
      !is.null(scale_ctrl) &&
      shape_trt != 1) {
    stop("please specify shape parameter only for the group for which scale has been specified")
  }

  # sets shape parameter equal
  if (shape_trt == 1) {
    shape_trt <- shape_ctrl
  }
  # sets shape parameter equal
  if (shape_ctrl == 1) {
    shape_ctrl <- shape_trt
  }
  # converts scale parameter wrt parameterization
  if (parameterization == 2) {
    if (!is.null(scale_trt)) scale_trt <- 1 / (scale_trt ^ (1 / shape_trt))
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / (scale_ctrl ^ (1 / shape_ctrl))
  }
  # converts scale parameter wrt parameterization
  if (parameterization == 3) {
    if (!is.null(scale_trt)) scale_trt <- 1 / scale_trt
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / scale_ctrl
  }
  # specifies previously undefined scale parameter if HR is chosen as input contrast
  if (!is.null(HR)) {
    if (!is.null(scale_trt)) {
      scale_ctrl <- scale_trt * HR
    } else{
      scale_trt <- scale_ctrl / HR
    }
  }
  # if median difference is chosen as input parameter, median in one group is
  # calculated from scale and shape parameters, and input from other group is
  # calculated from median of first group minus median difference. Undefined
  # Scale and shape parameters are calculated from median.
  if (!is.null(median_diff)) {
    if (!is.null(scale_trt)) {
      median_time_trt <- (-log(0.5)) ^ (1 / shape_trt) * scale_trt
      median_time_ctrl <- median_time_trt - median_diff
      scale_ctrl <- median_time_ctrl * (-log(0.5)) ^ (-1 / shape_ctrl)
    } else{
      median_time_ctrl <- (-log(0.5)) ^ (1 / shape_ctrl) * scale_ctrl
      median_time_trt <- median_time_ctrl + median_diff
      scale_trt <- median_time_trt * (-log(0.5)) ^ (-1 / shape_trt)
    }
  }
  # procedure similar to median difference
  if (!is.null(percentile_diff)) {
    if (!is.null(scale_trt)) {
      median_time_trt <- (-log(percentile)) ^ (1 / shape_trt) * scale_trt
      median_time_ctrl <- median_time_trt - median_diff
      scale_ctrl <- median_time_ctrl * (-log(percentile)) ^ (-1 / shape_ctrl)
    } else{
      median_time_ctrl <- (-log(percentile)) ^ (1 / shape_ctrl) * scale_ctrl
      median_time_trt <- median_time_ctrl + median_diff
      scale_trt <- median_time_trt * (-log(percentile)) ^ (-1 / shape_trt)
    }
  }
  # procedure similar to median difference
  if (!is.null(survival_diff)) {
    if (!is.null(scale_trt)) {
      survival_trt <- pweibull(t, scale_trt, shape_trt, lower.tail = F)
      survival_ctrl <- survival_trt - survival_diff
      if (survival_ctrl <= 0) {
        stop(
          "survival in treatment group minus survival difference is equal to or less than zero"
        )
      }
      scale_ctrl <- t * (-log(survival_ctrl)) ^ (-1 / shape_ctrl)
    } else{
      survival_ctrl <- pweibull(t, scale_ctrl, shape_ctrl, lower.tail = F)
      survival_trt <- survival_ctrl + survival_diff
      if (survival_trt <= 0) {
        stop("survival in control group minus survival difference is equal to or less than zero")
      }
      scale_trt <- t * (-log(survival_trt)) ^ (-1 / shape_trt)
    }
  }
  # calculate RMSTR or RMSTD.
  if (output == "RMSTD" || output == "RMSTR"){
    RMST_trt <- integrate(
      pweibull,
      shape = shape_trt,
      scale = scale_trt,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMST_ctrl <- integrate(
      pweibull,
      shape = shape_ctrl,
      scale = scale_ctrl,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
  }

  # calculate RMST in unspecified group when RMSTR or RMSTD are given as input
  if (xor(!is.null(RMSTD), !is.null(RMSTR))) {
    # when either RMSTR or RMSTD defined
    if (!is.null(scale_trt)) { #if scale_ctrl is undefined, calculate it from RMST of trt and RMSTD/RMSTR
      RMST_trt <- integrate(
        pweibull,
        shape = shape_trt,
        scale = scale_trt,
        lower = 0,
        upper = tau,
        lower.tail = F
      )$value
      if (!is.null(RMSTD)) {
        RMST_ctrl <- RMST_trt - RMSTD
      }
      if (!is.null(RMSTR)) {
        RMST_ctrl <- RMST_trt / RMSTR
      }
      scale_ctrl_temp <- uniroot(find_root_weibull, shape = shape_ctrl, tau = tau, RMST = RMST_ctrl,
                                 lower = 0.000001,
                                 upper = 100000, tol = 0.0001)$root
    }
    if (!is.null(scale_ctrl)) {
      RMST_ctrl <- integrate(
        pweibull,
        shape = shape_ctrl,
        scale = scale_ctrl,
        lower = 0,
        upper = tau,
        lower.tail = F
      )$value
      if (!is.null(RMSTD)) {
        RMST_trt <- RMST_ctrl + RMSTD
      }
      if (!is.null(RMSTR)) {
        RMST_trt <- RMST_ctrl * RMSTR
      }

      scale_trt <- uniroot(find_root_weibull, shape = shape_trt, tau = tau,
                           RMST = RMST_trt, lower = 0.000001,
                           upper = 100000, tol = 0.0001)$root
    }
    if (is.null(scale_ctrl)) scale_ctrl <- scale_ctrl_temp
  }
  if(output == "HR") {return((scale_ctrl/scale_trt)^shape_trt)}

  if(plot_curves){
    curve(pweibull(x, scale = scale_trt, shape = shape_trt, lower.tail = FALSE),
          col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1), xlim = c(0, 1.5*tau))
    curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE),
          col = "red", add = TRUE)
    abline(v = tau, col = "blue")
    text(x = tau, y = 0.1, pos = 4, labels = paste("time horizon τ  =", tau))
    legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trt, 2), " and shape =",
                                         round(shape_trt, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape_ctrl, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)


  }
  if (output == "RMSTD") {
    return(RMST_trt - RMST_ctrl)
  }
  if (output == "RMSTR") {
    return(RMST_trt / RMST_ctrl)
  }
  # results_df <- data.fame(scale_trt, shape_trt, scale_ctrl, shape_ctrl, tau)
}



find_root_weibull <- function(unknown_scale, shape, tau, RMST){
  integrate(pweibull, shape = shape,  scale = unknown_scale,  lower = 0,  upper = tau,  lower.tail = F)$value-RMST
}
