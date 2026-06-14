#' Calculates or converts contrast between RMSTD, RMSTR, and HR, and others
#'
#' Calculates or converts contrasts into either the difference or the ratio of restricted mean survival times (RMSTs), the hazard ratio (HR),
#' or other contrasts. Assumes Weibull distributed survival under proportional hazards.
#'
#' Calculates a range of different contrasts between two survival curves
#' given their scale parameters and shape parameter. Alternatively only one scale parameter may be supplied together with one contrast (e.g. hazard
#' ratio) in order to calculate the missing scale parameter and the remaining contrasts.
#' Survival curves are assumed to be Weibull functions with identical shape parameter.
#' Specify either one scale parameter and one contrast, or specify both scale parameters.
#'
#'
#' Supported contrasts are: \itemize{
#' \item Hazard ratio (HR).
#' \item Median difference: the difference in times \eqn{\Delta t} at which each survival curve passes \eqn{S(t) = 0.5 = 50 \%}.
#' \item Percentile difference: the difference in times \eqn{\Delta t} at which each survival curve passes a certain \eqn{S(t)}.
#' Simplifies to the median difference when the percentile \eqn{50} is chosen.
#' \item Survival difference: difference between two survival curves \eqn{S(t)} at a specified time \eqn{\tau}.}
#' Survival is assumed to follow Weibull distributions in both control and treatment group with a constant HR,
#' implying the same shape parameter in treatment and control group.
#'
#' @param scale_trmt Specifies the \dfn{scale parameter} in the treatment group.
#' @param scale_ctrl Specifies the \dfn{scale parameter} in the treatment group.
#' @param shape Specifies the \dfn{shape parameter} in both groups.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}},
#' \item \code{parameterization = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param RMSTD Specifies the RMSTD between control group and treatment group. Allows for converting RMSTD to HR.
#' @param RMSTR Specifies the RMSTR between control group and treatment group. Allows for converting RMSTD to HR.
#' @param tau Specifies the time horizon \eqn{\tau} where evaluate RMST with \eqn{\mathrm{RMST(\tau)} = \int_{0}^{\tau}S(t) \,dt}, and where to evaluate \eqn{\Delta S(\tau)}.
#' @param HR Specifies a hazard ratio with \eqn{\text{HR} = h(t)_{\text{trmt}} / h(t)_{\text{ctrl}}}.
#' @param median_diff Specifies the difference of median survival times, with \cr
#' \code{median_diff} \eqn{= t_{\text{median, trmt}} - t_{\text{median, ctrl}} = \{t : S_\text{trmt}(t) = 0.5\} - \{t : S_\text{ctrl}(t) = 0.5\}}.
#' @param percentile_diff Specifies the difference in time \eqn{\Delta t} at which each group has the survival \eqn{S(t) =} \code{percentile} with \cr
#' \code{percentile_diff} \eqn{= t_{\text{percentile, trmt}} - t_{\text{percentile, ctrl}} = \{t : S_\text{trmt}(t) = \code{percentile}\} - \{t : S_\text{ctrl}(t) = \code{percentile}\}}.
#' @param percentile Specifies at which percentile of survival to evaluate \code{percentile_diff}.
#' \code{percentile_diff} is equal to \code{median_diff} when \code{percentile}  \eqn{=50}.
#' @param survival_diff Specifies the survival difference \eqn{\Delta S(\tau) = S(t)_{trmt} - S(t)_{ctrl}}.
#' @param plot_curves Logical. Creates a plot if \code{TRUE}.
#'
#' @return Returns a list containing scale and shape parameters, hazard ratio, median and percentile difference, survival difference, RMST in treatment and control group,
#' and RMST difference and ratio.
#'
#' @export
#'
#' @examples
#' # Specify survival curves and obtain contrasts.
#' # Percentile difference is calculated at S(t) = 0.8. Obtain contrasts.
#' results <- convert_contrast_ph(scale_trmt = 10, scale_ctrl = 6, tau = 4,
#'                                percentile = 80, plot_curves = TRUE)
#' print(results)
#'
#' # Specify scale in control group and hazard ratio.
#' # Obtain scale in treatment group and remaining contrasts.
#' results <- convert_contrast_ph(scale_ctrl = 6, tau = 4, HR = 0.6)
#' print(results)
#'
#' # Specify scale in treatment group, shape and median survival.
#' # Obtain scale in control group and remaining contrasts.
#' results <- convert_contrast_ph(scale_trmt = 10, shape = 1.5,
#'                                tau = 4, median_diff = 2)
#' print(results)
#'
convert_contrast_ph <- function(
    scale_trmt = NULL,
    scale_ctrl = NULL,
    shape = 1,
    parameterization = 1,
    RMSTD = NULL,
    RMSTR = NULL,
    tau = NULL,
    HR = NULL,
    median_diff = NULL,
    percentile_diff = NULL,
    percentile = 50,
    survival_diff = NULL,
    plot_curves = FALSE) {
  # error management -------------------------------------------------------------
  number_of_defined_scales <- sum(!is.null(scale_trmt), !is.null(scale_ctrl))
  number_of_defined_contrasts <- sum(
    !is.null(HR),
    !is.null(median_diff),
    !is.null(percentile_diff),
    !is.null(survival_diff),
    !is.null(RMSTR),
    !is.null(RMSTD)
  )

  stopifnot(
    "error: please specify either both scale parameters, or one scale
            parameter and one contrast" = xor(
      number_of_defined_scales == 2 &&
        number_of_defined_contrasts == 0,
      number_of_defined_scales == 1 &&
        number_of_defined_contrasts == 1
    )
  )

  stopifnot(
    # throw error when parameterization misspecified
    "parameterization must be defined as either 1, 2, or 3" = parameterization ==
      1 ||
      parameterization == 2 ||
      parameterization == 3
  )

  # convert parameters -----------------------------------------------------------

  # converts scale parameter wrt parameterization
  if (parameterization == 2) {
    if (!is.null(scale_trmt)) {
      scale_trmt <- 1 / (scale_trmt^(1 / shape))
    }
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / (scale_ctrl^(1 / shape))
  }
  # converts scale parameter wrt parameterization
  if (parameterization == 3) {
    if (!is.null(scale_trmt)) {
      scale_trmt <- 1 / scale_trmt
    }
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / scale_ctrl
  }

  # obtain missing scale parameter -----------------------------------------------

  # specifies previously undefined scale parameter if HR is chosen as input contrast
  if (!is.null(HR)) {
    if (!is.null(scale_trmt)) {
      scale_ctrl <- scale_trmt * HR
    } else {
      scale_trmt <- scale_ctrl / HR
    }
  }

  # if median difference is chosen as input parameter, median in one group is
  # calculated from scale and shape parameters, and input from other group is
  # calculated from median of first group minus median difference. Undefined
  # Scale and shape parameters are calculated from median.
  if (!is.null(median_diff)) {
    if (!is.null(scale_trmt)) {
      median_time_trmt <- (-log(0.5))^(1 / shape) * scale_trmt
      median_time_ctrl <- median_time_trmt - median_diff
      scale_ctrl <- median_time_ctrl * (-log(0.5))^(-1 / shape)
    } else {
      median_time_ctrl <- (-log(0.5))^(1 / shape) * scale_ctrl
      median_time_trmt <- median_time_ctrl + median_diff
      scale_trmt <- median_time_trmt * (-log(0.5))^(-1 / shape)
    }
  }

  # percentile diff. procedure, similar to median difference
  if (!is.null(percentile_diff)) {
    if (!is.null(scale_trmt)) {
      percentile_time_trmt <- (-log(percentile / 100))^(1 / shape) * scale_trmt
      percentile_time_ctrl <- percentile_time_trmt - percentile_diff
      scale_ctrl <- percentile_time_ctrl * (-log(percentile / 100))^(-1 / shape)
    } else {
      percentile_time_ctrl <- (-log(percentile / 100))^(1 / shape) * scale_ctrl
      percentile_time_trmt <- percentile_time_ctrl + percentile_diff
      scale_trmt <- percentile_time_trmt * (-log(percentile / 100))^(-1 / shape)
    }
  }

  # survival diff. procedure similar to median difference
  if (!is.null(survival_diff)) {
    if (!is.null(scale_trmt)) {
      survival_trmt <- stats::pweibull(tau, scale = scale_trmt, shape = shape, lower.tail = F)
      survival_ctrl <- survival_trmt - survival_diff
      if (survival_ctrl <= 0) {
        stop(
          "survival in treatment group minus survival difference is equal to or less than zero"
        )
      }
      scale_ctrl <- tau * (-log(survival_ctrl))^(-1 / shape)
    } else {
      survival_ctrl <- stats::pweibull(tau, scale = scale_ctrl, shape = shape, lower.tail = F)
      survival_trmt <- survival_ctrl + survival_diff
      if (survival_trmt <= 0) {
        stop(
          "survival in control group minus survival difference is equal to or less than zero"
        )
      }
      scale_trmt <- tau * (-log(survival_trmt))^(-1 / shape)
    }
  }

  # calculate RMST in unspecified group when RMSTR or RMSTD are given as input

  if (xor(!is.null(RMSTD), !is.null(RMSTR))) {
    # when either RMSTR or RMSTD defined
    if (!is.null(scale_trmt)) {
      # if scale_ctrl is undefined, calculate it from RMST of trmt and RMSTD/RMSTR
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
      scale_ctrl_temp <- stats::uniroot(
        find_root_weibull,
        shape = shape,
        tau = tau,
        RMST = RMST_ctrl,
        lower = 0.000001,
        upper = 100000,
        tol = 0.0001
      )$root
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

      scale_trmt <- stats::uniroot(
        find_root_weibull,
        shape = shape,
        tau = tau,
        RMST = RMST_trmt,
        lower = 0.000001,
        upper = 100000,
        tol = 0.0001
      )$root
    }
    if (is.null(scale_ctrl)) scale_ctrl <- scale_ctrl_temp
  }

  # calculate missing contrasts --------------------------------------------------

  # calculate HR
  if (is.null(HR)) {
    HR <- scale_ctrl / scale_trmt
  }

  # calculate median_diff
  if (is.null(median_diff)) {
    median_time_trmt <- (-log(0.5))^(1 / shape) * scale_trmt
    median_time_ctrl <- (-log(0.5))^(1 / shape) * scale_ctrl
    median_diff <- median_time_trmt - median_time_ctrl
  }

  # calculate percentile_diff
  if (is.null(percentile_diff)) {
    percentile_time_trmt <- (-log(percentile / 100))^(1 / shape) * scale_trmt
    percentile_time_ctrl <- (-log(percentile / 100))^(1 / shape) * scale_ctrl
    percentile_diff <- percentile_time_trmt - percentile_time_ctrl
  }

  # calculate survival_diff
  if (is.null(survival_diff)) {
    survival_trmt <- stats::pweibull(tau, shape = shape, scale = scale_trmt, lower.tail = F)
    survival_ctrl <- stats::pweibull(tau, shape = shape, scale = scale_ctrl, lower.tail = F)
    survival_diff <- survival_trmt - survival_ctrl
  }
  # calculate RMSTs, and RMSTR or RMSTD, if not defined already.
  if (!exists("RMST_trmt")) {
    RMST_trmt <- stats::integrate(
      stats::pweibull,
      shape = shape,
      scale = scale_trmt,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
  }
  if (!exists("RMST_ctrl")) {
    RMST_ctrl <- stats::integrate(
      stats::pweibull,
      shape = shape,
      scale = scale_ctrl,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
  }
  if (is.null(RMSTD)) {
    RMSTD <- RMST_trmt - RMST_ctrl
  }
  if (is.null(RMSTR)) {
    RMSTR <- RMST_trmt / RMST_ctrl
  }

  if (plot_curves) {
    x <- NULL
    graphics::curve(
      stats::pweibull(x, scale = scale_trmt, shape = shape, lower.tail = FALSE),
      col = "darkblue",
      xlab = "t",
      ylab = "S(t)",
      ylim = c(0, 1),
      xlim = c(0, 1.5 * tau),
      lwd = 2
    )
    graphics::curve(
      stats::pweibull(x, scale = scale_ctrl, shape = shape, lower.tail = FALSE),
      col = "red",
      add = TRUE,
      lwd = 2
    )
    graphics::abline(v = tau, col = "black", lwd = 2)
    graphics::text(
      x = tau,
      y = 0.1,
      pos = 4,
      labels = bquote("time horizon " * tau * " = " * .(tau)),
      cex = .8
    )
    graphics::legend(
      "bottomleft",
      legend = c(
        paste0(
          "Treatment group with \n",
          "scale = ",
          round(scale_trmt, 2),
          " and shape = ",
          round(shape, 2)
        ),
        paste0(
          "Control group with \n",
          "scale = ",
          round(scale_ctrl, 2),
          " and shape = ",
          round(shape, 2)
        )
      ),
      col = c("darkblue", "red"),
      lty = 1:1,
      y.intersp = 1.5,
      bty = "n",
      cex = 0.8
    )
  }

  # prepare list of all results and return
  results_df <- list(
    "scale trmt" = scale_trmt,
    "scale ctrl" = scale_ctrl,
    "shape" = shape,
    "hazard ratio" = HR,
    "median difference" = median_diff,
    "percentile difference" = percentile_diff,
    "percentile" = percentile,
    "survival difference at tau" = survival_diff,
    "RMST trmt" = RMST_trmt,
    "RMST ctrl" = RMST_ctrl,
    "RMSTD" = RMSTD,
    "RMSTR" = RMSTR
  )
  return(results_df)
}

find_root_weibull <- function(unknown_scale, shape, tau, RMST) {
  stats::integrate(
    stats::pweibull,
    shape = shape,
    scale = unknown_scale,
    lower = 0,
    upper = tau,
    lower.tail = F
  )$value -
    RMST
}
