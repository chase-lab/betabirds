rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggExtra)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(ggmagnify)
library(ggfx)
library(ggforce)
library(reshape2)
library(GGally)

data <- read.csv('four_beta_birds_environ_cookie_20250709.csv')
data <- data[data$cookie != 5004, ]

cookies <- data$cookie
year <- data$year_val

mat <- matrix(NA, 10000, 4)

load('Results/Model_mean_beta_sn_200_outlier.RData')
mat[, 1] <- mod$phi[, 2]

load('Results/Model_mean_beta_q0_outlier.RData')
mat[, 2] <- mod$phi[, 2]

load('Results/Model_mean_beta_q2_outlier.RData')
mat[, 3] <- mod$phi[, 2]

load('Results/Model_mean_beta_C80_outlier.RData')
mat[, 4] <- mod$phi[, 2]


colnames(mat) <- names(data)[4:7]
df_long <- melt(mat)
colnames(df_long) <- c("Sample", "Variable", "Value")
levels(df_long$Variable) <- c("beta_sn_200", "beta_q0", "beta_q2", "beta_C80")

ggplot(df_long, aes(x = Value, color = Variable, fill = Variable)) +
  geom_density(alpha = 0.5) +
  labs(
    fill = "Beta diversity (US Birds)",
    color = "Beta diversity (US Birds)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.box.background = element_rect(fill = "white", color = "black")
  ) +
  geom_vline(xintercept = 0, col = 'red3', linetype = 'dashed', lwd = 0.8) +
  labs(
    title = "Posterior distribution of the slopes",
    subtitle = "Delta Beta diversity (slope) ~ Delta Environment (slope)",
    x = '',
    y = ""
  )

ggsave("Figures/4_slopes_outlier.pdf", width = 12, height = 8)

mean(mat[, 1] > 0) # Probability of beta_sn_200 > 0
mean(mat[, 2] > 0) # Probability of beta_q0 > 0
mean(mat[, 3] < 0) # Probability of beta_q2 < 0
mean(mat[, 4] > 0) # Probability of beta_C80 > 0
