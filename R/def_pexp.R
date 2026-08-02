#' Defines a piecewise exponential curve
#'
#' The piecewise exponential curve is constructed from pairs of \eqn{t} and \eqn{S(t)}, where the time points \eqn{t} serve as the breakpoints between the piecewise functions.
#'
#' @param breakpoints A vector of breakpoints in increasing order between the piecewise functions. First element must be \code{> 0}.
#' @param surv Survival probabilities \eqn{S(t)} at each breakpoint. The final pair in \code{surv} and \code{breakpoints} will serve for extrapolation beyond the second to last pair.
#' @param plot Logical. Specifies whether to plot the resulting function.
#'
#' @return Returns a named list with breakpoints and hazards for each interval. Breakpoints will have length of hazard vector \eqn{-1}.
#' @seealso \code{\link{def_weibull}}
#' @export
#'
#' @examples
#' def_pexp(breakpoints = c(1, 2), surv = c(0.9, 0.8))

def_pexp <- function(breakpoints = NULL, surv, plot = FALSE){
  # error management --------------------------------------------------------
  stopifnot("Number of breakpoints must be equal to number of elements in surv vector" = length(breakpoints) == length(surv))
  stopifnot("Breakpoints must be in strictly increasing order without duplicates" = all(diff(breakpoints) > 0))
  stopifnot("First breakpoint must be > 0" = breakpoints[1] > 0)
  stopifnot("Survival must be strictly decreasing" = all(diff(surv) < 0))
  stopifnot("Survival must be above 0 and below 1" = all(surv < 1 & surv > 0))
  breakpoints_temp <- c(0, breakpoints)
  surv_temp <- c(1, surv)
  interval_hazards <- -log(surv_temp[-1] / surv_temp[-length(surv_temp)]) / diff(breakpoints_temp)
  hazards <- c(interval_hazards, interval_hazards[length(interval_hazards)])
  if(plot){
    x <- NULL
    graphics::curve(
      ppweibull::ppweibull(
        x,
        rate = hazards,
        alpha = rep(1, length(hazards)),
        t = breakpoints_temp,
        lower.tail = FALSE
      ),
      col = "darkblue",
      xlab = "t",
      ylab = "S(t)",
      ylim = c(0, 1),
      xlim = c(0, 1.5 * breakpoints[length(breakpoints)]),
      main = paste0("Piecewise exponential survival with hazards = c(", paste(round(hazards, 2), collapse = ", "), ") and breakpoints at c(", paste(round(breakpoints, 2), collapse = ", "), ")"),
      yaxt = "n"
    )
    graphics::axis(
      2,
      at = seq(1, 0, by = -0.2),
      labels = paste0(seq(100, 0, by = -20), "%"),
      las = 1
    )
    graphics::segments(x0 = breakpoints, y0 = 0, y1 = surv, col = "black", lwd = 2)
  }
  return(hazards)
}
