files.sources = list.files() # source all required functions from working directory
sapply(files.sources, source)

if (!require("pacman")) install.packages("pacman") # install all required packages
pacman::p_load(devtools, numDeriv, npsurvSS, SSRMST, survmixer, ggplot2, hesim, survival)

source("sample_size_functions.R")
compare_sample_size()


steps <- 10
final_horizon <- 5

df <- data.frame(matrix(ncol = 8, nrow = steps))
for(i in 1:steps){
  horizon <- i*final_horizon/steps
  df[i,] <- compare_sample_size(horizon = horizon)
}
