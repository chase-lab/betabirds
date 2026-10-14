################################################################################
###                                                                          ###
###                   Figure S20: Posterior Distributions                    ###
###                                                                          ###
################################################################################

# This code plots the posterior distributions of the models associating
# the rate of change in beta diversity (q0, q1, and q2) with
# the rate of change in environmental heterogeneity.

# With all points

# Supplemental Figure

# Clean workspace
rm(list = ls())
gc()

# Load required packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(viridis)

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


#####################################
###                               ###
###       Data Loading and        ###
###          Processing           ###
###                               ###
#####################################

# Set up matrix
mat <- matrix(NA, 20000, 3)

# Load and extract results from each model
load(paste(
  path,
  'Model Results/Four Routes/Model_mean_beta_q0.RData',
  sep = ""
))

mat[, 1] <- c(mod1$phi[, 2], mod2$phi[, 2], mod3$phi[, 2], mod4$phi[, 2])

load(paste(
  path,
  'Model Results/Four Routes/Model_mean_beta_q1.RData',
  sep = ""
))

mat[, 2] <- c(mod1$phi[, 2], mod2$phi[, 2], mod3$phi[, 2], mod4$phi[, 2])

load(paste(
  path,
  'Model Results/Four Routes/Model_mean_beta_q2.RData',
  sep = ""
))

mat[, 3] <- c(mod1$phi[, 2], mod2$phi[, 2], mod3$phi[, 2], mod4$phi[, 2])


# Adjust matrix to prepare for plotting
colnames(mat) <- c("q0", "q1", "q2")
df_long <- reshape2::melt(mat)
colnames(df_long) <- c("Sample", "Variable", "Value")
levels(df_long$Variable) <- c("beta_q0", "beta_q1", "beta_q2")


#####################################
###                               ###
###      Define Probabilities     ###
###                               ###
#####################################

# Probability of beta_q0 > 0
mean(mat[, 1] > 0)

# Probability of beta_q1 > 0
mean(mat[, 2] > 0)

# Probability of beta_q2 < 0
mean(mat[, 3] < 0)

#####################################
###                               ###
###      Define 95% CIs           ###
###        and means              ###
#####################################

mean(df_long$Value[df_long$Variable == "beta_q0"])
quantile(df_long$Value[df_long$Variable == "beta_q0"], 0.025)
quantile(df_long$Value[df_long$Variable == "beta_q0"], 0.975)

mean(df_long$Value[df_long$Variable == "beta_q1"])
quantile(df_long$Value[df_long$Variable == "beta_q1"], 0.025)
quantile(df_long$Value[df_long$Variable == "beta_q1"], 0.975)

mean(df_long$Value[df_long$Variable == "beta_q2"])
quantile(df_long$Value[df_long$Variable == "beta_q2"], 0.025)
quantile(df_long$Value[df_long$Variable == "beta_q2"], 0.975)


#####################################
###                               ###
###          Plotting             ###
###                               ###
#####################################

p <- ggplot(df_long, aes(x = Value, color = Variable, fill = Variable)) +
  geom_density(alpha = 0.9) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 22)
  ) +
  scale_x_continuous(
    name = expression(atop(" " * beta * " -diversity slope ~ EH slope")),
    limits = c(-2.3, 1),
    breaks = c(-2, -1, 0, 1)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "beta_q0" = "#4D9DC9",
      "beta_q1" = "#EEBC4D",
      "beta_q2" = "#4DBB9D"
    ),
    labels = c("beta_q0" = "q0", "beta_q1" = "q1", "beta_q2" = "q2")
  ) +
  scale_color_manual(
    name = "",
    values = c("beta_q0" = "black", "beta_q1" = "black", "beta_q2" = "black"),
    labels = c("beta_q0" = "q0", "beta_q1" = "q1", "beta_q2" = "q2")
  ) +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed', lwd = 0.8)


ggsave(
  paste(path, "Manuscript Figures/FigS20_dist_4routes.jpg", sep = ""),
  p,
  device = "jpeg",
  height = 8,
  width = 12,
  dpi = 300
)

ggsave(
  paste(path, "Manuscript Figures/FigS20_dist_4routes.pdf", sep = ""),
  p,
  device = "pdf",
  height = 8,
  width = 12,
  dpi = 300
)
