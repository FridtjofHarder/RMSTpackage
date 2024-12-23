#' Calculates contrasts between fully specified survival functions
#'
#' Function for calculating contrasts between two fully specified Weibull survival
#' functions defined by scale and shape.
#'
#' Survival curves need to be defined by \dfn{scale} and \dfn{shape} parameter, in the standard parameterization
#' defined by \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}.
#' Possible contrasts are: \itemize{
#' \item Restricted mean survival time difference (RMSTD)
#' \item Restricted mean survival time ratio (RMSTR)
#' \item Hazard ratio (HR)
#' \item Median difference
#' \item Percentile difference at a certain percentile
#' \item Survival difference at time \eqn{t}}
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group.
#' Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group.
#' Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: specifies Weibull distributed survival as
#' \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: specifies Weibull distributed survival as
#' \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: specifies Weibull distributed survival as
#' \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param percentile A scalar specifying the percentile at which to evaluate the
#' percentile difference.
#' @param t A scalar specifying the time at which to evaluate the survival difference.
#' @param tau A scalar specifying the time horizon for evaluating RMSTs.
#' @param plot_curves Boolean. Creates a plot if \code{TRUE}.
#'
#' @return Returns the contrast between two fully specified survival curves.
#' Defaults to calculating all contrasts.
#'
#' @export
#'
#' @examples
#' calculate_contrast(scale_trmt = 2,
#' scale_ctrl = 1,
#' parameterization = 1,
#' percentile = NULL,
#' t = NULL,
#' tau = 1,
#' plot_curves = TRUE)
calculate_contrast <- function(scale_trmt = NULL,
                               shape_trmt = 1,
                               scale_ctrl = NULL,
                               shape_ctrl = 1,
                               parameterization = 1,
                               percentile = NULL,
                               t = NULL,
                               tau = NULL,
                               plot_curves = TRUE){
# error management -------------------------------------------------------------

  stopifnot( # throw error when parameterization misspecified
    "parameterization must be defined as either 1, 2, or 3" =
      parameterization == 1 || parameterization == 2 || parameterization == 3
  )
  # throw error when functions are misspecified
  if (is.null(scale_trmt) || is.null(scale_ctrl)) {
    stop("please specify scale parameters in both treatment and survival group")
  }

# adjust scale parameter if parameterization != 1-------------------------------

  if (parameterization == 2) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / (scale_trmt ^ (1 / shape_trmt))
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / (scale_ctrl ^ (1 / shape_ctrl))
  }

  if (parameterization == 3) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / scale_trmt
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / scale_ctrl
  }

# main function ----------------------------------------------------------------
  browser()
  # define weibull functions for treatment and control
  weibull_ctrl <- function(q){
    stats::pweibull(q, shape = shape_ctrl, scale = scale_ctrl,
             lower.tail = FALSE)
  }
  weibull_trmt <- function(q){
    stats::pweibull(q, shape = shape_trmt, scale = scale_trmt,
             lower.tail = FALSE)
  }

  # median survival time and median survival time difference
  median_surv_time_ctrl <- uniroot(function(q) {weibull_ctrl(q) - 0.5},
                                   lower = 0, upper = 1000)$root
  median_surv_time_trmt <- uniroot(function(q) {weibull_trmt(q) - 0.5},
                                   lower = 0, upper = 1000)$root
  median_surv_time_diff <- median_surv_time_trmt - median_surv_time_ctrl

  # hazard ratio
  if(shape_ctrl == 1) h_ctrl <- 1 / scale_ctrl else {
    h_ctrl <-  function(t) get_weibull_h(t, shape_ctrl, scale_ctrl)
  }
  if(shape_trmt == 1) h_trmt <- 1 / scale_trmt else {
    h_trmt <-  function(t) get_weibull_h(t, shape_trmt, scale_trmt)
  }

  h_trmt <-  function(t) get_weibull_h(t, shape_trmt, scale_trmt)
  h_ratio <- function(t) h_trmt(t) / h_ctrl(t)

  # calculate RMSTD and RMSTR for treatment and control curve
  if (!is.null(tau)){
    RMST_trmt <- stats::integrate(
      stats::pweibull,
      shape = shape_trmt,
      scale = scale_trmt,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMST_ctrl <- stats::integrate( # calculate RMST for control curve
      stats::pweibull,
      shape = shape_ctrl,
      scale = scale_ctrl,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMSTD <- RMST_trmt - RMST_ctrl
    RMSTR <- RMST_trmt / RMST_ctrl
  }

  # hazard ratio

  # percentile survival

  # pointwise survival difference

  if (plot_curves) {
    x = NULL # to silence check() notes
    # plots treatment curve
    graphics::curve(stats::pweibull(x, scale = scale_trmt, shape = shape_trmt,
                                    lower.tail = FALSE),
                    col = "darkblue", xlab = "t", ylab = "S(t)", lwd = 3,
                    ylim = c(0, 1), xlim = c(0, 1.5*tau), las = 1)
    # plots control curve
    graphics::curve(stats::pweibull(x, scale = scale_ctrl, shape = shape_ctrl,
                                    lower.tail = FALSE),
                    col = "red", lwd = 3, add = TRUE)
    graphics::abline(v = tau, col = "black", lwd = 3)
    # mark tau if tau is defined
    if (!is.null(tau)) {
    graphics::text(x = tau, y = 0.1, pos = 4,
                   labels = bquote("time horizon " * tau * " = " * .(tau)))
    }
    graphics::legend("bottomleft",
                     legend = c(paste0("treatment group with \n", "scale = ",
                                         round(scale_trmt, 2), " and shape = ",
                                         round(shape_trmt, 2)),
                                paste0("control group with \n",   "scale = ",
                                         round(scale_ctrl, 2), " and shape = ",
                                         round(shape_ctrl, 2))),
                     col=c("darkblue", "red"), lty = 1, y.intersp = 2,
                     bty = "n", cex = 1)
    graphics::title("Survival curves in control group and treatment group")
  }
  return(list("RMSTD" = RMSTD, "RMSTR" = RMSTR))
}

#' Hazard function for weibull distributed survival
#'
#' @param t time
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @returns
get_weibull_h <- function(t, shape, scale){
  return((shape / scale) * (t / scale)^(shape - 1))
}
