rm(list=ls())

source("SamplingFunction.R")

data = read.csv('pairs_beta_birds_environ_cookie_20250708.csv')

cookies = as.numeric(as.character(data$cookie))
year = data$year_val

x = data$mean_beta_shannon_9_env

#: Run models for each measure and save in RData object
y = data$mean_beta_C80
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_C80.RData')

y = data$mean_beta_sn_200
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_sn_200.RData')

y = data$mean_beta_q2
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_q2.RData')

y = data$mean_beta_q0
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_q0.RData')



# Remove the outlier cookie 5004
data = data[data$cookie != 5004,]

cookies = as.numeric(as.character(data$cookie))
year = data$year_val

x = data$mean_beta_shannon_9_env

#: Run models for each measure and save in RData object
y = data$mean_beta_C80
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_C80_outlier.RData')

y = data$mean_beta_sn_200
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_sn_200_outlier.RData')

y = data$mean_beta_q2
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_q2_outlier.RData')

y = data$mean_beta_q0
mod = MCMC(y,x,cookies,year,R=20000,burn=10000)
save(list=c('mod'),file='Results/Model_mean_beta_q0_outlier.RData')
