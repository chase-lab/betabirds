################################################################################
###                                                                          ###
###                       Trends of EH and Beta Div                          ###
###                                                                          ###
################################################################################

# This code extracts the posterior mean and 95% credible intervals for the
# trends in EH and the trends in beta diversity for q=0, q=1, and q=2

# EH = A
# Beta = B

# Clean workspace
rm(list = ls())
gc()

# Load required packages
library(tidyverse)

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


#####################################
###                               ###
###       Data Loading            ###
###                               ###
#####################################

# Load the raw data
data <- read.csv(paste(
  path,
  'pairs_beta_birds_environ_cookie_20250708.csv',
  sep = ""
))


#####################################
###                               ###
###           q=0                 ###
###                               ###
#####################################

# Load the model object that includes all spatial polygons
load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q0.RData',
  sep = ""
))

# Create a list of the models
mod.list <- list(mod1, mod2, mod3, mod4)

# Define cookies
cookies <- data$cookie

# Beta: 3 dimensional array with 10000 rows
# 2 columns (2 coefficients or parameters per row)
# 74 layers (1 per cookie)

# Make dataframe for all estimates
all.estimates_q0 <- data.frame()

for (j in seq_along(mod.list)) {
  mod <- mod.list[[j]]

  # Processing to extract model predictions
  n <- nrow(mod$beta[, 2, ])
  m <- ncol(mod$beta[, 2, ])
  B <- mod$eta[, 2, ]
  A <- mod$beta[, 2, ]

  df <- data.frame()
  for (i in 1:m) {
    df <- rbind(
      df,
      data.frame(
        A = A[, i],
        B = B[, i],
        Variable = paste0("Cookie", cookies[i]),
        Model = paste0("Mod", j)
      )
    )
  }

  # Bind results
  all.estimates_q0 <- rbind(all.estimates_q0, df)
}


# Compute means and 95% CIs by cookie
cookie.means <- all.estimates_q0 %>%
  group_by(Variable) %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Compute means and 95% CIs across all cookies
all.means <- all.estimates_q0 %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Save the results
write.csv(
  cookie.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q0_bycookie_Trend Means and CIs.csv",
  row.names = FALSE
)


write.csv(
  all.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q0_Trend Means and CIs.csv",
  row.names = FALSE
)

# Clean up workspace
rm(
  A,
  all.estimates_q0,
  all.means,
  B,
  cookie.means,
  df,
  mod,
  mod.list,
  mod1,
  mod2,
  mod3,
  mod4
)


#####################################
###                               ###
###           q=1                 ###
###                               ###
#####################################

# Load the model object that includes all spatial polygons
load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q1.RData',
  sep = ""
))

# Create a list of the models
mod.list <- list(mod1, mod2, mod3, mod4)

# Define cookies
cookies <- data$cookie

# Beta: 3 dimensional array with 10000 rows
# 2 columns (2 coefficients or parameters per row)
# 74 layers (1 per cookie)

# Make dataframe for all estimates
all.estimates_q1 <- data.frame()

for (j in seq_along(mod.list)) {
  mod <- mod.list[[j]]

  # Processing to extract model predictions
  n <- nrow(mod$beta[, 2, ])
  m <- ncol(mod$beta[, 2, ])
  B <- mod$eta[, 2, ]
  A <- mod$beta[, 2, ]

  df <- data.frame()
  for (i in 1:m) {
    df <- rbind(
      df,
      data.frame(
        A = A[, i],
        B = B[, i],
        Variable = paste0("Cookie", cookies[i]),
        Model = paste0("Mod", j)
      )
    )
  }

  # Bind results
  all.estimates_q1 <- rbind(all.estimates_q1, df)
}


# Compute means and 95% CIs by cookie
cookie.means <- all.estimates_q1 %>%
  group_by(Variable) %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Compute means and 95% CIs across all cookies
all.means <- all.estimates_q1 %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Save the results
write.csv(
  cookie.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q1_bycookie_Trend Means and CIs.csv",
  row.names = FALSE
)


write.csv(
  all.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q1_Trend Means and CIs.csv",
  row.names = FALSE
)

# Clean up workspace
rm(
  A,
  all.estimates_q1,
  all.means,
  B,
  cookie.means,
  df,
  mod,
  mod.list,
  mod1,
  mod2,
  mod3,
  mod4
)


#####################################
###                               ###
###           q=2                 ###
###                               ###
#####################################

# Load the model object that includes all spatial polygons
load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q2.RData',
  sep = ""
))

# Create a list of the models
mod.list <- list(mod1, mod2, mod3, mod4)

# Define cookies
cookies <- data$cookie

# Beta: 3 dimensional array with 10000 rows
# 2 columns (2 coefficients or parameters per row)
# 74 layers (1 per cookie)

# Make dataframe for all estimates
all.estimates_q2 <- data.frame()

for (j in seq_along(mod.list)) {
  mod <- mod.list[[j]]

  # Processing to extract model predictions
  n <- nrow(mod$beta[, 2, ])
  m <- ncol(mod$beta[, 2, ])
  B <- mod$eta[, 2, ]
  A <- mod$beta[, 2, ]

  df <- data.frame()
  for (i in 1:m) {
    df <- rbind(
      df,
      data.frame(
        A = A[, i],
        B = B[, i],
        Variable = paste0("Cookie", cookies[i]),
        Model = paste0("Mod", j)
      )
    )
  }

  # Bind results
  all.estimates_q2 <- rbind(all.estimates_q2, df)
}


# Compute means and 95% CIs by cookie
cookie.means <- all.estimates_q2 %>%
  group_by(Variable) %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Compute means and 95% CIs across all cookies
all.means <- all.estimates_q2 %>%
  summarise(
    mean_EH = mean(A),
    EH_lower_CI = quantile(A, 0.025),
    EH_upper_CI = quantile(A, 0.975),
    mean_beta = mean(B),
    beta_lower_CI = quantile(B, 0.025),
    beta_upper_CI = quantile(B, 0.975),
    .groups = "drop"
  )

# Save the results
write.csv(
  cookie.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q2_bycookie_Trend Means and CIs.csv",
  row.names = FALSE
)


write.csv(
  all.means,
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Model_Results/q2_Trend Means and CIs.csv",
  row.names = FALSE
)

# Clean up workspace
rm(
  A,
  all.estimates_q2,
  all.means,
  B,
  cookie.means,
  df,
  mod,
  mod.list,
  mod1,
  mod2,
  mod3,
  mod4
)
