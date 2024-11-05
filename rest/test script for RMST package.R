scale_trt <- 1.2
shape_trt <- 1.1
HR <- 1.05
tau <- 4
survival <- pweibull(median_time, scale = scale_trt, shape = shape_trt, lower.tail = F)
median_time*(-log(survival))^(-1/shape_trt)


median_time <- (-log(0.5))^(1/shape_trt)*scale_trt
scale_ctrl <- median_time*(-log(0.5))^(-1/shape_trt)

result <- integrate(pweibull, scale = scale_trt, shape = shape_trt, lower.tail = F, lower = 0, upper = tau)

result <- convert_contrast(scale_trt = scale_trt, shape_trt = shape_trt, parameterization = 1, median_diff = 0.1, percentile_diff = 0.1, percentile = NULL,
                             survival_diff = NULL, t = NULL, output = "RMSTR", tau = tau)


?convert_contrast

testit <- function() warning("testit")
testit() ## shows call
testit <- function() warning("problem in testit", call. = FALSE)
testit() ## no call

a = NULL
b = NULL
ablist <- list(a,b)
sapply(ablist, is.null)
