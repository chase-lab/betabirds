################################################################################
###                                                                          ###
###      Plotting posterior distribution of Trends of EH and Beta Div        ###
###                                                                          ###
################################################################################

# This code posterior distributions for the 74 metacommunities used in the
# main analysis that calculated beta diversity from all possible pairs of
# communities in the metacommunity
# trends in EH and the trends in beta diversity for q=0, q=1, and q=2

# EH = A
# Beta = B

# Clean workspace
rm(list = ls())
gc()

# Load required packages
library(tidyverse)
library(ggplot2)
library(viridis)
library(tidybayes)


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

for (j in 1:length(mod.list)) {
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

# Convert the log values to % change
all_q0_transformed <- all.estimates_q0 %>%
  mutate(
    across(
      c(-Variable, -Model),
      ~ (10^.x - 1) * 100
    )
  )


####################
###              ###
###   Plotting   ###
###              ###
####################

# EH
EH.plot <- ggplot(all_q0_transformed, aes(x = A, fill = factor(Variable))) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = "% Change in environmental heterogeneity\nper time interval"
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__EH trend posterior distribution.jpg",
    sep = ""
  ),
  EH.plot,
  device = "jpeg",
  height = 7,
  width = 9,
  dpi = 300
)

# Zoomed in EH
EH.plot_zoom <- ggplot(
  all_q0_transformed,
  aes(x = A, fill = factor(Variable))
) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = "% Change in \nenvironmental heterogeneity\nper time interval",
    limits = c(-1, 1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, color = "gray30")) +
  theme(axis.title.y = element_text(size = 20, color = "gray30")) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__EH zoomed trend posterior distribution.jpg",
    sep = ""
  ),
  EH.plot_zoom,
  device = "jpeg",
  height = 7,
  width = 7,
  dpi = 300
)


#q0
q0.plot <- ggplot(all_q0_transformed, aes(x = B, fill = factor(Variable))) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-11.5, 13.3)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q0 trend posterior distribution.jpg",
    sep = ""
  ),
  q0.plot,
  device = "jpeg",
  height = 7,
  width = 9,
  dpi = 300
)

# Zoomed in EH
q0.plot_zoom <- ggplot(
  all_q0_transformed,
  aes(x = A, fill = factor(Variable))
) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-1, 1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, color = "gray30")) +
  theme(axis.title.y = element_text(size = 20, color = "gray30")) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q0 zoomed trend posterior distribution.jpg",
    sep = ""
  ),
  q0.plot_zoom,
  device = "jpeg",
  height = 7,
  width = 7,
  dpi = 300
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

for (j in 1:length(mod.list)) {
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

# Convert the log values to % change
all_q1_transformed <- all.estimates_q1 %>%
  mutate(
    across(
      c(-Variable, -Model),
      ~ (10^.x - 1) * 100
    )
  )


####################
###              ###
###   Plotting   ###
###              ###
####################

# q1
q1.plot <- ggplot(all_q1_transformed, aes(x = B, fill = factor(Variable))) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-11.5, 13.3)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q1 trend posterior distribution.jpg",
    sep = ""
  ),
  q1.plot,
  device = "jpeg",
  height = 7,
  width = 9,
  dpi = 300
)

# Zoomed in EH
q1.plot_zoom <- ggplot(
  all_q1_transformed,
  aes(x = A, fill = factor(Variable))
) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-1, 1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, color = "gray30")) +
  theme(axis.title.y = element_text(size = 20, color = "gray30")) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q1 zoomed trend posterior distribution.jpg",
    sep = ""
  ),
  q1.plot_zoom,
  device = "jpeg",
  height = 7,
  width = 7,
  dpi = 300
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

for (j in 1:length(mod.list)) {
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

# Convert the log values to % change
all_q2_transformed <- all.estimates_q2 %>%
  mutate(
    across(
      c(-Variable, -Model),
      ~ (10^.x - 1) * 100
    )
  )


####################
###              ###
###   Plotting   ###
###              ###
####################

# q2
q2.plot <- ggplot(all_q2_transformed, aes(x = B, fill = factor(Variable))) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-11.5, 13.3)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q2 trend posterior distribution.jpg",
    sep = ""
  ),
  q2.plot,
  device = "jpeg",
  height = 7,
  width = 9,
  dpi = 300
)

# Zoomed in EH
q2.plot_zoom <- ggplot(
  all_q2_transformed,
  aes(x = A, fill = factor(Variable))
) +
  geom_density(alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = expression(atop(
      "% change in " * beta * " -diversity ",
      "per time interval"
    )),
    limits = c(-1, 1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_discrete(name = "", palette = turbo(74)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, color = "gray30")) +
  theme(axis.title.y = element_text(size = 20, color = "gray30")) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "Manuscript Figures/FigS__q2 zoomed trend posterior distribution.jpg",
    sep = ""
  ),
  q2.plot_zoom,
  device = "jpeg",
  height = 7,
  width = 7,
  dpi = 300
)
