scale <- 1.3
shape <- 1.2
tau <- 1.3
time <- 1.4
survival <- pweibull(q = time, shape = shape, scale = scale, lower.tail = F)
RMST <- integrate(
  pweibull,
  shape = shape,
  scale = scale,
  lower = 0,
  upper = tau,
  lower.tail = F
)$value



pweibull(1, shape = shape, scale = scale, lower.tail = F)
Integral <- integrate(pweibull, scale = scale, shape = shape, lower = 0, upper = tau, lower.tail = F)
function_to_find_root <- function()

fnToFindRoot = function(unknown_scale) {
  integrate(pweibull, shape = shape,  scale = unknwonscale,  lower = 0,  upper = tau,  lower.tail = F)$value-RMST
}

unknown_scale <- uniroot(fnToFindRoot, c(0, 10), tol = 0.0001)$root
unknown_scale + shape + survival + tau + time

find_root_weibull <- function(unknown_scale, shape_passed){
  integrate(pweibull, shape = shape,  scale = unknown_scale,  lower = 0,  upper = tau,  lower.tail = F)$value-RMST
}



start.time <- Sys.time()
for(i in 1:10000){
uniroot(find_root_weibull, shape_passed = shape, f.lower = -RMST, lower = 0.0001, upper = 1000, tol = 0.0001, trace = 4)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

