#' Calculates contrasts between fully specified survival functions
#'
#' Function for calculating contrasts between two fully specified Weibull survival
#' functions defined by scale and shape.
#'
#' Survival curves need to be defined by \dfn{scale} and \dfn{shape} parameter, in the standard parameterization
#' defined by \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}. Possible contrasts are: \itemize{
#' \item Restricted mean survival time difference (RMSTD)
#' \item Restricted mean survival time ratio (RMSTR)
#' \item Hazard ratio (HR)
#' \item Median difference
#' \item Percentile difference at a certain percentile
#' \item Survival difference at time \eqn{t}}
#'
#' @param scale_trmt A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_trmt A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_trmt} \eqn{=1}, simplifying to exponential survival.
#' @param scale_ctrl A scalar \eqn{>0} specifying the \dfn{scale parameter} in the treatment group.
#' @param shape_ctrl A scalar \eqn{>0} specifying the \dfn{shape parameter} in the treatment group. Defaults to \code{shape_ctrl} \eqn{=1}, simplifying to exponential survival.
#' @param parameterization One of: \itemize{
#' \item \code{parameterization = 1}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}},
#' \item \code{parameterization = 2}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterization = 3}: specifies Weibull distributed survival as \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}}.}
#' @param output A string specifying whether to calculate the difference in
#' RMSTs (\code{"RMSTD"}), the ratio between RMSTs (\code{"RMSTR"}), the HR
#' (\code{"HR"}), the median difference (\code{"median_diff"}), the difference
#' in time (\code{"percentile_diff"}) for reaching a survival percentile as
#' specified in (\code{percentile}), and the survival difference (\code{"survival_diff"})
#' at time specified in \code{t}.
#' (\code{"HR"}).
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
#' output = "RMSTD",
#' percentile = NULL,
#' t = NULL,
#' tau = 1,
#' plot_curves = TRUE)
calculate_contrast <- function(scale_trmt = NULL,
                               shape_trmt = 1,
                               scale_ctrl = NULL,
                               shape_ctrl = 1,
                               parameterization = 1,
                               output = "RMSTD",
                               percentile = NULL,
                               t = NULL,
                               tau = NULL,
                               plot_curves = TRUE){
# error management --------------------------------------------------------

  stopifnot( # throw error when output misspecified
    "output misspecified" = output == "RMSTD" ||
      output == "RMSTR" || output == "HR" || output == "median_diff" || output == "percentile_diff"
    || output == "survival_diff"
  )
  stopifnot( # throw error when parameterization misspecified
    "parameterization must be defined as either 1, 2, or 3" = parameterization == 1 ||
      parameterization == 2 || parameterization == 3
  )
  # throw error when functions misspecified
  if(is.null(scale_trmt)||is.null(scale_ctrl)){
    stop("please specify scale parameters for both treatment and survival group")
  }

  # throw error when tau should have been specified
  if((output == "all" || output == "RMSTR" || output == "RMSTD") &&  !is.numeric(tau)){
    stop("please specify valid time horizon tau")
  }

  # throw error when percentile should have been specified
  if(output == "all" || output == "percentine_diff")
    stop("please specify valid percentile at which to evaluate percentile difference"
         = !is.numeric(percentile))

  # throw error when t should have been specified
  if(output == "all" || output == "survival_diff")
    stop("please specify valid time t at which to evaluate survival difference"
         = !is.numeric(t))

# main function -----------------------------------------------------------

  # converts scale parameter wrt parameterization
  if (parameterization == 2) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / (scale_trmt ^ (1 / shape_trmt))
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / (scale_ctrl ^ (1 / shape_ctrl))
  }
  # converts scale parameter wrt parameterization
  if (parameterization == 3) {
    if (!is.null(scale_trmt)) scale_trmt <- 1 / scale_trmt
    if (!is.null(scale_ctrl)) scale_ctrl <- 1 / scale_ctrl
  }

  if(output == "all" || output == "RMSTD" || output == "RMSTR"){
    RMST_trmt <- integrate( # calculate RMST for treatment curve
      pweibull,
      shape = shape_trmt,
      scale = scale_trmt,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMST_ctrl <- integrate( # calculate RMST for control curve
      pweibull,
      shape = shape_ctrl,
      scale = scale_ctrl,
      lower = 0,
      upper = tau,
      lower.tail = F
    )$value
    RMSTD <- RMST_trmt - RMST_ctrl
    RMSTR <- RMST_trmt / RMST_ctrl
  }

  if(output == "all" || output == "RMSTD" || output == "RMSTR")

  if(plot_curves){
    curve(pweibull(x, scale = scale_trmt, shape = shape_trmt, lower.tail = FALSE),
          col = "green", xlab = "t", ylab = "S(t)", ylim = c(0, 1), xlim = c(0, 1.5*tau))
    curve(pweibull(x, scale = scale_ctrl, shape = shape_ctrl, lower.tail = FALSE),
          col = "red", add = TRUE)
    abline(v = tau, col = "blue")
    text(x = tau, y = 0.1, pos = 4, labels = paste("time horizon \U1D70F  =", tau))
    legend("bottomleft", legend=c(paste0("treatment group with \n", "scale =",
                                         round(scale_trmt, 2), " and shape =",
                                         round(shape_trmt, 2)),
                                  paste0("control group with \n", "scale =",
                                         round(scale_ctrl, 2), " and shape =",
                                         round(shape_ctrl, 2))),
           col=c("green", "red"), lty=1:1, y.intersp = 1.5, bty = "n", cex = 0.8)


  }
  return(list("RMSTD" = RMSTD, "RMSTR" = RMSTR))
}

