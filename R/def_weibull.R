#' Defines a single global Weibull curve
#'
#' The Weibull curve is constructed from two time points \eqn{t_1} and \eqn{t_2}, and two survival probabilities \eqn{S_1(t_1)} and \eqn{S_2(t_2)} on a survival curve. It has constant shape and scale over the full time range; it is not piecewise.
#' Piecewise definition of a survival curve is currently only supported for the \code{\link[=def_pexp]{piecewise exponential function}}.
#'
#' @param time A vector of two time points \eqn{t_1} and \eqn{t_2}, both \eqn{> 0}
#' @param surv Two survival probabilities \eqn{S_1(t_1)} and \eqn{S_2(t_2)} at the time points specified by \code{time}. Both survival probabilities must be \eqn{> 0} and \eqn{< 1}.
#' @param plot Logical. Specifies whether to plot the resulting function.
#' @param parameterisation Define only if Weibull function is specified, not for piecewise exponential survival. One of: \itemize{
#' \item \code{parameterisation = 1}: Default. Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(\mathrm{scale} * t)^\mathrm{shape})}},
#' \item \code{parameterisation = 2}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-\mathrm{scale} * t^\mathrm{shape})}},
#' \item \code{parameterisation = 3}: Specifies Weibull distributed survival as \cr \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape})}}. This is the parameterisation used for the base \R{} function \code{stats::pweibull()}.}
#'
#'
#' @return Returns a named list with one shape parameter, and one scale parameter, defining a Weibull survival \eqn{S(t) = 1- F(t) = \exp{(-(t/\mathrm{scale})^\mathrm{shape}))}}.
#' @seealso \code{\link{def_pexp}}
#' @export
#'
#' @examples
#' def_weibull(time = c(1.5, 2), surv = c(0.9, 0.8))

def_weibull <- function(time, surv, plot = FALSE, parameterisation = 1){
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
  # scale follows parameterisation 1, scale_external follows user-specified parameterisation
  scale <- scale_external <- (-log(surv[1]))^(1 / shape) / time[1]
  if(parameterisation !=1){
    scale_external <- rereparameterise(parameterisation, scale, shape)
  }
  if(plot){
    curve_data <- data.frame(t = seq(0, 1, length.out = 500))
    curve_data$S <- stats::pweibull(
      curve_data$t,
      scale = scale,
      shape = shape,
      lower.tail = FALSE
    )

    plot_obj <- ggplot2::ggplot(curve_data, ggplot2::aes(x = t, y = S)) +
      ggplot2::geom_line(color = "darkblue", linewidth = 1) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),
        labels = paste0(seq(0, 100, by = 20), "%"),
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        x = "t",
        y = "S(t)",
        title = paste0("Weibull survival with shape = ", round(shape, 2), " and scale = ", round(scale_external, 2))
      ) +
      ggplot2::theme_bw()

    print(plot_obj)
  }
  return(list(shape = shape, scale = scale_external))
}
