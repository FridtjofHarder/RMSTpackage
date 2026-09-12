#' Defines a piecewise exponential curve
#'
#' The piecewise exponential survival function is constructed from \eqn{n} pairs of breakpoints
#' \eqn{t_1, \dots, t_n} and survival probabilities \eqn{S(t_1), \dots, S(t_n)}.
#' Let \eqn{t_0 = 0} and \eqn{S(t_0) = 1}. For each interval \eqn{(t_{i-1}, t_i]}, \eqn{i = 1, \dots, n-1},
#' a constant hazard \eqn{h_i} is computed such that
#' \deqn{S(t_i) = S(t_{i-1}) \exp\{-h_i (t_i - t_{i-1})\}.}
#' The final hazard \eqn{h_n} is defined on \eqn{(t_{n-1}, \infty)} and is determined by the last pair
#' \eqn{(t_n, S(t_n))} via
#' \deqn{S(t_n) = S(t_{n-1}) \exp\{-h_n (t_n - t_{n-1})\},}
#' and then held constant for all \eqn{t > t_{n-1}}.
#'
#' @param breakpoints Specifies a vector of breakpoints in increasing order, with the first
#'   element \code{> 0}. These define the interval boundaries \eqn{t_1, \dots, t_n}.
#' @param surv Specifies a vector of survival probabilities \eqn{S(t_1), \dots, S(t_n)} at each breakpoint.
#'   Must have the same length as \code{breakpoints}, with values in \eqn{(0, 1]} and strictly decreasing.
#' @param plot Logical. If \code{TRUE}, plots the resulting piecewise exponential survival function.
#'
#' @return A named list with components:
#'   - breakpoints: The input vector of breakpoints t1, ..., t_(n-1).
#'   Note that the ultimate element has been removed from the input breakpoints vector.
#'   - hazards: A numeric vector of length n containing the constant hazards
#'     h1, ..., hn for the intervals (0, t1], (t1, t2], ..., (t_(n-1), Inf).
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
  hazards <- -log(surv_temp[-1] / surv_temp[-length(surv_temp)]) / diff(breakpoints_temp)
  breakpoints_red <- breakpoints[-length(breakpoints)]
  if(plot){
    x <- NULL
    graphics::curve(
      ppweibull::ppweibull(
        x,
        rate = hazards,
        alpha = rep(1, length(hazards)),
        t = c(0, breakpoints_red),
        lower.tail = FALSE
      ),
      col = "darkblue",
      lwd = 2,
      xlab = "t",
      ylab = "S(t)",
      ylim = c(0, 1),
      xlim = c(0, 1.5 * breakpoints_red[length(breakpoints_red)]),
      main = c(
        paste0(
          "Piecewise exponential survival with hazards = c(",
          paste(round(hazards, 2), collapse = ", "),
          ")"
        ),
        paste0(
          "and breakpoints at c(",
          paste(round(breakpoints_red, 2), collapse = ", "),
          ")"
        )
      ),
      yaxt = "n"
    )
    graphics::axis(
      2,
      at = seq(1, 0, by = -0.2),
      labels = paste0(seq(100, 0, by = -20), "%"),
      las = 1
    )
    graphics::segments(x0 = breakpoints_red, y0 = 0, y1 = surv, col = "black", lwd = 2)
  }
  return(list(breakpoints = breakpoints_red,
              hazards = hazards))
}
