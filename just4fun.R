## just example code. Can delete anytime.

###example calculate_sample_size()
calculate_sample_size(
scale_trmt = 1.4,
scale_ctrl = 1,
accrual_time = 1,
follow_up_time = 10,
tau = 1,
RMSTD_simulation = TRUE,
RMSTR_simulation = TRUE,
cox_ph_simulation = TRUE,
loss_scale = 1,
M = 100,
simulation_sample_size = 100)

### example simulate_data()
data1 <- simulate_data(scale = 1,
accrual_time = 1,
follow_up_time = 1,
loss_scale = 1,
sample_size = 100,
label = 0,
plot_curve = TRUE,
plot_data = FALSE,
plot_recruitment = FALSE)

data2 <- simulate_data(scale = 0.5,
                       accrual_time = 1,
                       follow_up_time = 1,
                       loss_scale = 1,
                       sample_size = 100,
                       label = "trmt",
                       plot_curve = TRUE,
                       plot_data = FALSE,
                       plot_recruitment = FALSE)

data <- rbind(data1, simulate_data(scale = 0.5,
                                         accrual_time = 1,
                                         follow_up_time = 1,
                                         loss_scale = 1,
                                         sample_size = 100,
                                         label = 1,
                                         plot_curve = TRUE,
                                         plot_data = FALSE,
                                         plot_recruitment = FALSE))


result <-  rmst2(data$observations, data$status, data$label, tau = 100,
                 alpha = 0.05)



