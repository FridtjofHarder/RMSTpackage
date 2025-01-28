#' Converts contrast between RMSTD, RMSTR, and HR
#'
#' Function for converting different contrasts into either the difference or the ratio of restricted mean survival times (RMSTs), the hazard ratio (HR),
#' and other contrasts.
#'
#' Convertible contrasts are: \itemize{
#' \item hazard ratio (HR).
#' \item median difference: the difference in times \eqn{t} at which each survival curve passes \eqn{S(t) = 0.5 = 50 \%}.
#' \item percentile difference: the difference in times \eqn{t} at which each survival curve passes a certain \eqn{S(t)}. Simplifies to the median difference when the percentile \eqn{50} is chosen.
#' \item survival difference: difference between two survival curves \eqn{S(t)} at a specified time \eqn{t}.}
#' Survival is assumed to follow Weibull distributions in both control and treatment group with a constant HR,
#' implying the same shape parameter in treatment and control group.
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape A scalar \eqn{>0} specifying the \dfn{shape parameter} in both groups. Shape
#' parameters in both groups are assumed to be equal.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param median_diff A scalar specifying the difference in time \eqn{\Delta t} at which each group has survival \eqn{S(t) = 0.5}. \code{median_diff} \eqn{> 0} suggests superior survival in the treatment group.
#' @param percentile_diff A scalar specifying the difference in time \eqn{\Delta t} at which each group has the survival \eqn{S(t) =} \code{percentile}.
#' @param percentile A scalar specifying at which percentile of survival to evaluate \code{percentile_diff}.
#' \code{percentile_diff} is equal to \code{median_diff} when \code{percentile}  \eqn{=50}.
#' @param survival_diff A scalar specifying the survival difference in survival \eqn{\Delta S(t)} between treatment and control group.
#' @param HR A scalar specifying a hazard ratio. \code{HR} \eqn{< 1} suggests a lower hazard in the treatment group than in the control group.
#' @param tau A scalar specifying the time horizon \eqn{\tau} at which to evaluate RMST with \eqn{\mathrm{RMST} = \int_{0}^{\tau}S(t) \,dt}.
#' @param t A scalar specifying the time \eqn{t} at which to evaluate \code{survival_diff}.
#' @param RMSTD A scalar specifying the RMSTD between control group and treatment group. Allows for converting RMSTD to HR.
#' @param RMSTR A scalar specifying the RMSTR between control group and treatment group. Allows for converting RMSTD to HR.
#' @param plot_curves Boolean. Creates a plot if \code{TRUE}
#'
#' @return Returns a dataframe containing the inputs, the difference (RMSTD) and ratio (RMSTR) in RMST,
#' the HR between treatment group and control group, the median survival time difference
#' An RMSTD \eqn{> 0} or an RMSTR \eqn{> 1} indicate a larger RMST in the treatment group,
#' a HR \eqn{< 1} indicates lower hazard in the treatment group.
#'
#' @export
#'
#' @examples
#' RMSTD <- convert_contrast(scale_trmt = 1, shape = 1, HR = 0.9,
#' output = "RMSTD", tau = 1)
#' RMSTR <- convert_contrast(scale_trmt = 1, shape = 1, survival_diff = 0.01,
#' t = 1, output = "RMSTR", tau = 1)
#' HR <- convert_contrast(scale_trmt = 1, shape = 1, RMSTD = 0.028,
#' t = 1, output = "HR", tau = 1) # inserting RMSTD from first example
convert_contrast <- function(scale_trmt = NULL,
                             scale_ctrl = NULL,
                             parameterization = 1,
                             HR = NULL,
                             median_diff = NULL,
                             percentile_diff = NULL,
                             percentile = NULL,
                             survival_diff = NULL,
                             RMSTD = NULL,
                             RMSTR = NULL,
                             tau = NULL,
                             plot_curves = TRUE) {


# error management --------------------------------------------------------

   stopifnot( # throw error when contrast misspecified
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

  stopifnot( # throw error when output misspecified
    "output needs to be specified as either 'RMSTD' or 'RMSTR' or 'HR'" = output == "RMSTD" ||
      output == "RMSTR" || output == "HR"
  )
  stopifnot( # throw error when parameterization misspecified
    "parameterization must be defined as either 1, 2, or 3" = parameterization == 1 ||
      parameterization == 2 || parameterization == 3
  )
  # throw error when not exactly one scale
  stopifnot("please specify exactly one scale parameter. Use calculate_contrast() to determine the RMSTD, RMSTR, or HR for fully specified survival functions"
            = xor(is.null(scale_trmt), is.null(scale_ctrl)))

  # throws error if shape parameter is specified for the group, for which scale is unspecified
  if (!is.null(scale_trmt) &&
      shape != 1 ||
      !is.null(scale_ctrl) &&
      shape != 1) {
    stop("please specify shape parameter only for the group for which scale has been specified")
  }


# main function --------------------------------------------------------------------

  # converts scale parameter wrt parameterization
  if (parameterization == 2) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / (scale_trmt ^ (1 / shape))
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / (scale_ctrl ^ (1 / shape))
  }
  # converts scale parameter wrt parameterization
  if (parameterization == 3) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / scale_trmt
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / scale_ctrl
  }
  # specifies previously undefined scale parameter if HR is chosen as input contrast
  if (!is.null(HR)) {
    if (!is.null(scale_trmt)) {
      scale_ctrl <- scale_trmt * HR
    } else{
      scale_trmt <- scale_ctrl / HR
    }
  }
  # if median difference is chosen as input parameter, median in one group is
  # calculated from scale and shape parameters, and input from other group is
  # calculated from median of first group minus median difference. Undefined
  # Scale and shape parameters are calculated from median.
  if (!is.null(median_diff)) {
    if (!is.null(scale_trmt)) {
      median_time_trmt <- (-log(0.5)) ^ (1 / shape) * scale_trmt
      median_time_ctrl <- median_time_trmt - median_diff
      scale_ctrl <- median_time_ctrl * (-log(0.5)) ^ (-1 / shape)
    } else{
      median_time_ctrl <- (-log(0.5)) ^ (1 / shape) * scale_ctrl
      median_time_trmt <- median_time_ctrl + median_diff
      scale_trmt <- median_time_trmt * (-log(0.5)) ^ (-1 / shape)
    }
  }
  # procedure similar to median difference
  if (!is.null(percentile_diff)) {
    if (!is.null(scale_trmt)) {
      median_time_trmt <- (-log(percentile)) ^ (1 / shape) * scale_trmt
      median_time_ctrl <- median_time_trmt - median_diff
      scale_ctrl <- median_time_ctrl * (-log(percentile)) ^ (-1 / shape)
    } else{
      median_time_ctrl <- (-log(percentile)) ^ (1 / shape) * scale_ctrl
      median_time_trmt <- median_time_ctrl + median_diff
      scale_trmt <- median_time_trmt * (-log(percentile)) ^ (-1 / shape)
    }
  }
  # procedure similar to median difference
  if (!is.null(survival_diff)) {
    if (!is.null(scale_trmt)) {
      survival_trmt <- stats::pweibull(t, scale_trmt, shape, lower.tail = F)
      survival_ctrl <- survival_trmt - survival_diff
      if (survival_ctrl <= 0) {
        stop(
          "survival in treatment group minus survival difference is equal to or less than zero"
        )
      }
      scale_ctrl <- t * (-log(survival_ctrl)) ^ (-1 / shape)
    } else{
      survival_ctrl <- stats::pweibull(t, scale_ctrl, shape, lower.tail = F)
      survival_trmt <- survival_ctrl + survival_diff
      if (survival_trmt <= 0) {
        stop("survival in control group minus survival difference is equal to or less than zero")
      }
      scale_trmt <- t * (-log(survival_trmt)) ^ (-1 / shape)
    }
  }
  # calculate RMSTR or RMSTD.
  if (output == "RMSTD" || output == "RMSTR"){
    RMST_trmt <- stats::integrate(
      stats::pweibull,
      shape = shape,
      scale = scale_trmt,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMST_ctrl <- stats::integrate(
      stats::pweibull,
      shape = shape,
      scale = scale_ctrl,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
  }

  # calculate RMST in unspecified group when RMSTR or RMSTD are given as input
  if (xor(!is.null(RMSTD), !is.null(RMSTR))) {
    # when either RMSTR or RMSTD defined
    if (!is.null(scale_trmt)) { #if scale_ctrl is undefined, calculate it from RMST of trmt and RMSTD/RMSTR
      RMST_trmt <- stats::integrate(
        stats::pweibull,
        shape = shape,
        scale = scale_trmt,
        lower = 0,
        upper = tau,
        lower.tail = F
      )$value
      if (!is.null(RMSTD)) {
        RMST_ctrl <- RMST_trmt - RMSTD
      }
      if (!is.null(RMSTR)) {
        RMST_ctrl <- RMST_trmt / RMSTR
      }
      scale_ctrl_temp <- stats::uniroot(find_root_weibull, shape = shape, tau = tau, RMST = RMST_ctrl,
                                 lower = 0.000001,
                                 upper = 100000, tol = 0.0001)$root
    }
    if (!is.null(scale_ctrl)) {
      RMST_ctrl <- stats::integrate(
        stats::pweibull,
        shape = shape,
        scale = scale_ctrl,
        lower = 0,
        upper = tau,
        lower.tail = F
      )$value
      if (!is.null(RMSTD)) {
        RMST_trmt <- RMST_ctrl + RMSTD
      }
      if (!is.null(RMSTR)) {
        RMST_trmt <- RMST_ctrl * RMSTR
      }

      scale_trmt <- stats::uniroot(find_root_weibull, shape = shape, tau = tau,
                           RMST = RMST_trmt, lower = 0.000001,
                           upper = 100000, tol = 0.0001)$root
    }
    if (is.null(scale_ctrl)) scale_ctrl <- scale_ctrl_temp
  }
  if(output == "HR") {return((scale_ctrl/scale_trmt)^shape)}

  if(plot_curves){
    x <- NULL
    graphics::curve(stats::pweibull(x, scale = scale_trmt, shape = shape, lower.tail = FALSE),
          col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1), xlim = c(0, 1.5*tau))
    graphics::curve(stats::pweibull(x, scale = scale_ctrl, shape = shape, lower.tail = FALSE),
          col = "red", add = TRUE)
    graphics::abline(v = tau, col = "blue")
    graphics::text(x = tau, y = 0.1, pos = 4, labels = bquote("time horizon " * tau * " = " * .(tau)))
    graphics::legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trmt, 2), " and shape =",
                                         round(shape, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)


  results_df <- data.fame(scale_trmt, shape, scale_ctrl, shape, tau)
  return(results_df)
  }
}

find_root_weibull <- function(unknown_scale, shape, tau, RMST){
  stats::integrate(stats::pweibull, shape = shape,  scale = unknown_scale,  lower = 0,  upper = tau,  lower.tail = F)$value-RMST
}
