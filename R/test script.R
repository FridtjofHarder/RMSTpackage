scale_trt <- 1.3
shape_trt <- 1.1
HR <- 1.05
tau <- 1.4
pweibull(1, scale = scale_trt, shape = shape_trt, lower.tail = F)
result <- integrate(pweibull, scale = scale_trt, shape = shape_trt, lower.tail = F, lower = 0, upper = tau)

convert_contrast(scale_trt = scale_trt, HR = HR, output = "RMSTR", tau = tau)
?convert_contrast


