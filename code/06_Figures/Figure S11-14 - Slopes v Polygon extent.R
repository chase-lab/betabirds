###############################################################################
########                                                               ########
########          Figure S11-14 : Investigating the effect of          ########
########                        the extent of spatial polygons         ########
###############################################################################

# Author: Kimberly Thompson

# This code investigates whether any systematic patterns exist in the extent
# of spatial polygons (i.e., metacommunities) and the trends
# (i.e. slopes) of EH and beta diversity metrics (for q=0, 1, and 2)

# Change in EH and beta derived from posterior means of Bayesian regression
# models

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(tidyverse) # Data cleaning
library(sf)
library(ggplot2)
library(forcats)

###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


# Load the cookies/spatial polygons
cookies_sf <- sf::st_read(paste(
  path,
  "cookies_routes_assembled_20252305.shp",
  sep = ""
))

# Load the raw data of EH and B
data <- read.csv(paste(
  path,
  "pairs_beta_birds_environ_cookie_20250708.csv",
  sep = ""
))

# Load the model output for each beta metric
load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q0.RData',
  sep = ""
))

# Define the results to plot from
modq0 <- mod1
rm(mod1, mod2, mod3, mod4)


load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q1.RData',
  sep = ""
))

# Define the results to plot from
modq1 <- mod1
rm(mod1, mod2, mod3, mod4)


load(paste(
  path,
  'Model Results/Route Pairs/Model_mean_beta_q2.RData',
  sep = ""
))

# Define the results to plot from
modq2 <- mod1
rm(mod1, mod2, mod3, mod4)


###############################################
###                                         ###
###             Data Processing             ###
###                                         ###
###############################################

# Create unique df with only ID__200 (cookie name) and Area
cookies_area <- cookies_sf %>%
  select(ID__200, are_km2) %>%
  distinct(.keep_all = TRUE)

cookies <- data$cookie


###############################################
###                                         ###
###     Calculating Change in EH and B      ###
###                                         ###
###############################################

# EH as the response will be calculated along with each metric (but it should
# remain roughly the same - i.e, only slight differences because we're
# drawing different samples from the posterior each time)

###################
###    B q0     ###
###################

m_q0 <- ncol(modq0$beta[, 2, ])
B_q0 <- modq0$eta[, 2, ]
A_q0 <- modq0$beta[, 2, ]

df_q0 <- data.frame()
for (i in 1:m_q0) {
  df_q0 <- rbind(
    df_q0,
    data.frame(
      A_q0 = A_q0[, i],
      B_q0 = B_q0[, i],
      Variable = paste0("Cookie", cookies[i])
    )
  )
}

# Compute means and 95% credible intervals
means_q0 <- df_q0 %>%
  group_by(Variable) %>%
  summarise(
    mean_A_q0 = mean(A_q0),
    A_q0_lower = quantile(A_q0, 0.025),
    A_q0_upper = quantile(A_q0, 0.975),
    mean_B_q0 = mean(B_q0),
    B_q0_lower = quantile(B_q0, 0.025),
    B_q0_upper = quantile(B_q0, 0.975),
    .groups = "drop"
  )

# Clean up
rm(A_q0, B_q0, df_q0, modq0)


###################
###    B q1     ###
###################

m_q1 <- ncol(modq1$beta[, 2, ])
B_q1 <- modq1$eta[, 2, ]
A_q1 <- modq1$beta[, 2, ]

df_q1 <- data.frame()
for (i in 1:m_q1) {
  df_q1 <- rbind(
    df_q1,
    data.frame(
      A_q1 = A_q1[, i],
      B_q1 = B_q1[, i],
      Variable = paste0("Cookie", cookies[i])
    )
  )
}

# Compute means and 95% credible intervals
means_q1 <- df_q1 %>%
  group_by(Variable) %>%
  summarise(
    mean_A_q1 = mean(A_q1),
    A_q1_lower = quantile(A_q1, 0.025),
    A_q1_upper = quantile(A_q1, 0.975),
    mean_B_q1 = mean(B_q1),
    B_q1_lower = quantile(B_q1, 0.025),
    B_q1_upper = quantile(B_q1, 0.975),
    .groups = "drop"
  )

# Clean up
rm(A_q1, B_q1, df_q1, modq1)


###################
###    B q2     ###
###################

m_q2 <- ncol(modq2$beta[, 2, ])
B_q2 <- modq2$eta[, 2, ]
A_q2 <- modq2$beta[, 2, ]

df_q2 <- data.frame()
for (i in 1:m_q2) {
  df_q2 <- rbind(
    df_q2,
    data.frame(
      A_q2 = A_q2[, i],
      B_q2 = B_q2[, i],
      Variable = paste0("Cookie", cookies[i])
    )
  )
}

# Compute means and 95% credible intervals
means_q2 <- df_q2 %>%
  group_by(Variable) %>%
  summarise(
    mean_A_q2 = mean(A_q2),
    A_q2_lower = quantile(A_q2, 0.025),
    A_q2_upper = quantile(A_q2, 0.975),
    mean_B_q2 = mean(B_q2),
    B_q2_lower = quantile(B_q2, 0.025),
    B_q2_upper = quantile(B_q2, 0.975),
    .groups = "drop"
  )

# Clean up
rm(A_q2, B_q2, df_q2, modq2)


###############################################
###                                         ###
###     Combining area and posterior        ###
###                 means                   ###
###############################################

# Prepare means dfs: extract numeric ID from "Cookie107" --> "107"
means_q0 <- means_q0 %>%
  mutate(ID = as.numeric(str_remove(Variable, "Cookie")))

means_q1 <- means_q1 %>%
  mutate(ID = as.numeric(str_remove(Variable, "Cookie")))

means_q2 <- means_q2 %>%
  mutate(ID = as.numeric(str_remove(Variable, "Cookie")))

# Prepare cookies_area: rename ID column
cookies_area_df <- cookies_area %>%
  st_drop_geometry() %>% # removes geometry, keeps as data.frame
  rename(ID = ID__200) %>%
  mutate(ID = as.numeric(ID)) %>%
  select(ID, are_km2)

# Merge: all rows from means dfs
means_q0 <- left_join(means_q0, cookies_area_df, by = "ID")

means_q1 <- left_join(means_q1, cookies_area_df, by = "ID")

means_q2 <- left_join(means_q2, cookies_area_df, by = "ID")

# In each dataframe A is change in environmental heterogeneity and
# B is change in beta metric

###############################################
###                                         ###
###               Plotting                  ###
###                                         ###
###############################################

###############
###   EH    ###
###############

# Create ordered factor (smallest are_km2 at bottom)
means_q0$y_label <- fct_reorder(factor(means_q0$ID), means_q0$are_km2)

# Create custom labels: every 10th position shows are_km2, others blank
label_positions <- seq(1, 74, by = 10) # 1, 11, 21, 31, 41, 51, 61, 71
y_labels <- rep("", 74)
y_labels[label_positions] <-
  number(
    means_q0$are_km2[order(means_q0$are_km2)][label_positions],
    accuracy = 1
  )

# Environmental Heterogeneity
p <- ggplot(
  means_q0,
  aes(x = mean_A_q0, xmin = A_q0_lower, xmax = A_q0_upper, y = y_label)
) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, color = "red") +
  scale_y_discrete(labels = y_labels) +
  labs(x = "Change in EH", y = expression("Spatial polygon extent" ~ (km^2))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.x = element_text(size = 22, face = "bold"))

# Save the plot
ggsave(
  paste(path, "Manuscript Figures/FigS11_EH v polygon extent.jpg", sep = ""),
  p,
  device = "jpeg",
  height = 7,
  width = 11,
  dpi = 300
)


###############
###   Q0    ###
###############

q0 <- ggplot(
  means_q0,
  aes(x = mean_B_q0, xmin = B_q0_lower, xmax = B_q0_upper, y = y_label)
) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, color = "red") +
  scale_y_discrete(labels = y_labels) +
  labs(
    x = "Change in Beta Diversity q=0",
    y = expression("Spatial polygon extent" ~ (km^2))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.x = element_text(size = 22, face = "bold"))

# Save the plot
ggsave(
  paste(
    path,
    "Manuscript Figures/FigS12_Beta q0 v polygon extent.jpg",
    sep = ""
  ),
  q0,
  device = "jpeg",
  height = 7,
  width = 11,
  dpi = 300
)


###############
###   q1    ###
###############

# Create ordered factor (smallest are_km2 at bottom)
means_q1$y_label <- fct_reorder(factor(means_q1$ID), means_q1$are_km2)

q1 <- ggplot(
  means_q1,
  aes(x = mean_B_q1, xmin = B_q1_lower, xmax = B_q1_upper, y = y_label)
) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, color = "red") +
  scale_y_discrete(labels = y_labels) +
  labs(
    x = "Change in Beta Diversity q=1",
    y = expression("Spatial polygon extent" ~ (km^2))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.x = element_text(size = 22, face = "bold"))

# Save the plot
ggsave(
  paste(
    path,
    "Manuscript Figures/FigS13_Beta q1 v polygon extent.jpg",
    sep = ""
  ),
  q1,
  device = "jpeg",
  height = 7,
  width = 11,
  dpi = 300
)


###############
###   q2    ###
###############

# Create ordered factor (smallest are_km2 at bottom)
means_q2$y_label <- fct_reorder(factor(means_q2$ID), means_q2$are_km2)

q2 <- ggplot(
  means_q2,
  aes(x = mean_B_q2, xmin = B_q2_lower, xmax = B_q2_upper, y = y_label)
) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, color = "red") +
  scale_y_discrete(labels = y_labels) +
  labs(
    x = "Change in Beta Diversity q=2",
    y = expression("Spatial polygon extent" ~ (km^2))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.y = element_text(size = 22, face = "bold")) +
  theme(axis.title.x = element_text(size = 22, face = "bold"))

# Save the plot
ggsave(
  paste(
    path,
    "Manuscript Figures/FigS14_Beta q2 v polygon extent.jpg",
    sep = ""
  ),
  q2,
  device = "jpeg",
  height = 7,
  width = 11,
  dpi = 300
)
