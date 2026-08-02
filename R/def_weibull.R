#' Defines a single global Weibull curve
#'
#' The Weibull curve is constructed from two pairs of \eqn{t} and \eqn{S(t)} on a survival curve. It has constant shape and scale over the full time range; it is not piecewise.
#' Piecewise definition from given points on a survival curve is currently only supported for \code{\link[=def_pexp]{piecewise exponential function}}.
#'
#' @param time A vector of two time points, both \eqn{> 0}
#' @param surv Two survival probabilities \eqn{S(t)} at the time points specified by \code{time}. Both survival probabilities must be \eqn{> 0} and \eqn{< 1}.
#' @param plot Logical. Specifies whether to plot the resulting function.
#'
#' @return Returns a named list with one shape parameter, and one scale parameter, defining a Weibull survival \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}.
#' @seealso \code{\link{def_pexp}}
#' @export
#'
#' @examples
#' def_weibull(time = c(1, 2), surv = c(0.9, 0.8))

def_weibull <- function(time, surv, plot = FALSE){
  # error management --------------------------------------------------------
  stopifnot("Provide exactly two time points" = length(time) == 2)
  stopifnot("Provide exactly two values for the survival curve" = length(surv) == 2)
  stopifnot("Time points must be > 0" = all(time > 0))
  stopifnot("Survival must be above 0 and below 1" = all(surv > 0 & surv < 1))
  stopifnot("Inputs of time must differ from each other" = !any(duplicated(time)))
  stopifnot("Inputs of surv must differ from each other" = !any(duplicated(surv)))
  stopifnot("Survival curve must be decreasing: smaller time must correspond to larger surv input" =
              order(time) == rev(order(surv)))
  shape <- log((-log(surv[1])) / (-log(surv[2]))) / log(time[1] / time[2])
  scale <- time[1] / (-log(surv[1]))^(1 / shape)
  if(plot){
    x <- NULL
    graphics::curve(
      ppweibull::ppweibull(
        x,
        rate = 1 / scale^shape,
        alpha = shape,
        lower.tail = FALSE
      ),
      col = "darkblue",
      xlab = "t",
      ylab = "S(t)",
      ylim = c(0, 1),
      main = paste0("Weibull survival with shape = ", round(shape, 2), " and scale = ", round(scale, 2)),
      yaxt = "n"
    )
    graphics::axis(
      2,
      at = seq(1, 0, by = -0.2),
      labels = paste0(seq(100, 0, by = -20), "%"),
      las = 1
    )
  }
  return(list(shape = shape, scale = scale))
}
