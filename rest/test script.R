# scale_trt <- 1.3
# shape_trt <- 1.1
# HR <- 1.05
# tau <- 1.4
# pweibull(median_time, scale = scale_trt, shape = shape_trt, lower.tail = F)
#
#
# median_time <- (-log(0.5))^(1/shape_trt)*scale_trt
# scale_ctrl <- median_time*(-log(0.5))^(-1/shape_trt)
#
# result <- integrate(pweibull, scale = scale_trt, shape = shape_trt, lower.tail = F, lower = 0, upper = tau)
#
# convert_contrast(scale_trt = scale_trt, HR = HR, output = "RMSTR", tau = tau)
#
#
# ?convert_contrast
#
# testit <- function() warning("testit")
# testit() ## shows call
# testit <- function() warning("problem in testit", call. = FALSE)
# testit() ## no call


