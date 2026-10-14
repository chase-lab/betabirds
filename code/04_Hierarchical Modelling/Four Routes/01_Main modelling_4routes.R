################################################################################
########                                                                ########
########              Fitting Bayesian Hierarchical Models              ########
########                                                                ########
################################################################################

# This code fits Bayesian hierarchical models
# to log-transformed values of beta diversity and environmental heterogeneity
# for three beta diversity metrics: q0, q1, and q2 derived from
# random sampling of 4 routes within each spatial polygon using the sampling function
# contained in 00_SamplingFunction.R.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

# Define file path
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


###############################################
###                                         ###
###      Data and Function Loading          ###
###                                         ###
###############################################

# Load the sampling function
source(paste(path, "04_Hierarchical Modelling/00_SamplingFunction.R", sep = ""))

# Define the path again
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"

# Load the data for the beta diversity and environmental heterogeneity
data <- read.csv(paste(
  path,
  'four_beta_birds_environ_cookie_20250709.csv',
  sep = ""
))


###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################

cookies <- as.numeric(as.character(data$cookie))
year <- data$year_val

x <- log(data$mean_beta_shannon_9_env)


###############################################
###                                         ###
###                Modelling                ###
###                                         ###
###############################################

# Q2
y <- log(data$mean_beta_2)
mod1 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod2 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod3 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod4 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
save(
  list = c('mod1', 'mod2', 'mod3', 'mod4'),
  file = paste(
    path,
    'Model Results/Four Routes/Model_mean_beta_q2.RData',
    sep = ""
  )
)

# Q0
y <- log(data$mean_beta_0)
mod1 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod2 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod3 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod4 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
save(
  list = c('mod1', 'mod2', 'mod3', 'mod4'),
  file = paste(
    path,
    'Model Results/Four Routes/Model_mean_beta_q0.RData',
    sep = ""
  )
)

# Q1
y <- log(data$mean_beta_1)
mod1 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod2 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod3 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
mod4 <- MCMC(y, x, cookies, year, R = 10000, burn = 5000)
save(
  list = c('mod1', 'mod2', 'mod3', 'mod4'),
  file = paste(
    path,
    'Model Results/Four Routes/Model_mean_beta_q1.RData',
    sep = ""
  )
)
