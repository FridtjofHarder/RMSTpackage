scale_trt <- 1.4
scale_ctrl <- 1
margin = 0.2
M <- 100

result <- compare_sample_size(scale_trt = scale_trt, scale_ctrl = scale_ctrl,
                              accrual_time = 1, follow_up_time = 10, tau = 1, sides = 1, alpha = 0.025, power = 0.8, margin = margin, simulation = TRUE, M = 100)

weibull_trt <- pweibull(t, scale = scale_trt, shape = 1, lower.tail = F)

x<-1:10
df<-data.frame(x)
ggplot(df,aes(x))+
  stat_function(fun=function(t) pweibull(t, scale = scale_trt, shape = 1, lower.tail = F),
                colour = "red") +
  xlim(c(0, 2)) +
  ylim(c(0,1)) +
  stat_function(fun=function(t) pweibull(t, scale = scale_ctrl, shape = 1, lower.tail = F),
                colour = "green") +
  labs(title = "Survival in trt (red, scale = 1.4) vs. ctrl (green, scale = 1)") +
  xlab("t") +
  ylab(("S(t)"))

inferiority_margin = seq(0, 0.2, 0.01)
sample_size_npsurvSS = rep(NA, length(margin))
simulated_power_SSRMST = rep(NA, length(margin))

for(i in 1:length(margin)){
  sample_size_npsurvSS[i] <- compare_sample_size(scale_trt = scale_trt, scale_ctrl = scale_ctrl,
                                    accrual_time = 1, follow_up_time = 10, tau = 1, sides = 1, alpha = 0.025,
                                    power = 0.8, margin = margin[i], simulation = TRUE, M = 100)$`n calculated via closed form`
  simulated_power_SSRMST[i] <- compare_sample_size(scale_trt = scale_trt, scale_ctrl = scale_ctrl,
                                                    accrual_time = 1, follow_up_time = 10, tau = 1, sides = 1, alpha = 0.025,
                                                    power = 0.8, margin = margin[i], simulation = TRUE, M = 1000)$`simulated test power`$power1
}

base::plot(x = inferiority_margin, y = sample_size_npsurvSS)
base::plot(x = inferiority_margin, y = simulated_power_SSRMST, ylim = c(0,1))

