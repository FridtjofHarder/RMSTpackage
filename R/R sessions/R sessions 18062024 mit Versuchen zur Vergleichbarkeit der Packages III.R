sample_size_closed_form_noninf <- function (arm0, arm1, test = list(test = "weighted logrank"), 
          power = 0.8, alpha = 0.025, sides = 1, margin = 0) 
{
  if (!inherits(test[[1]], "list")) {
    p1 <- arm1$size/(arm0$size + arm1$size)
    p0 <- 1 - p1
    design <- calc_design(arm0, arm1, test)
    out <- (sqrt(design$sigma2) * stats::qnorm(1 - alpha/sides) + 
              sqrt(design$tsigma2) * stats::qnorm(power))^2/design$(delta^2 * 
      c(p0, p1)
    if (sides == 2) {
      i <- 1
      arm0$size <- out[1]
      arm1$size <- out[2]
      cont <- T
      while (cont) {
        temp <- power_two_arm(arm0, arm1, test, alpha, 
                              sides)
        if (temp > power) {
          i <- i + 1
          arm0$size <- out[1] - i * p0
          arm1$size <- out[2] - i * p1
        }
        else {
          out <- out - (i - 1) * c(p0, p1)
          cont <- F
        }
      }
    }
    out <- c(out, sum(out), out[1] * prob_event(arm0), out[2] * 
               prob_event(arm1), out[1] * prob_event(arm0) + out[2] * 
               prob_event(arm1))
    names(out) <- c("n0", "n1", "n", "d0", "d1", "d")
    return(out)
  }
  else {
    out <- c()
    for (i in 1:length(test)) {
      label = ifelse("label" %in% names(test[[i]]), test[[i]]$label, 
                     i)
      out <- rbind(out, c(label, size_two_arm(arm0, arm1, 
                                              test[[i]], power, alpha, sides)))
    }
    out <- data.frame(out)
    names(out)[1] <- "test"
    return(out)
  }
}
