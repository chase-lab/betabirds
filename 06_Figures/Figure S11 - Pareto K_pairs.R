################################################################################
########                                                                ########
########              Plotting the results of PSIS test                 ########
########                                                                ########
################################################################################

# This code creates a plot of Pareto k values to see which spatial polygons
# are more influential and potentially outliers.

# Pareto k values (PSIS results) are taken from
# diagnostics run on Bayesian hierarchical models fit to
# to log-transformed values of beta diversity and environmental heterogeneity
# for three beta diversity metrics: q0, q1, and q2 derived from
# pairs of routes within each spatial polygon using the sampling function
# contained in 00_SamplingFunction.R.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(coda)
library(loo)
library(tidyverse)
library(ggplot2)

# Define file path
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


###############################################
###                                         ###
###               Data Loading              ###
###                                         ###
###############################################

# Load the diagnostic data
load(paste(path, 'Model Results/Route Pairs/Diagnostics.RData', sep = ""))

# Load the raw data
data <- read.csv(paste(
  path,
  'pairs_beta_birds_environ_cookie_20250708.csv',
  sep = ""
))


###############################################
###                                         ###
###            Data Preparation             ###
###                                         ###
###############################################

names <- c("q0", "q1", "q2")

M <- NULL
P <- NULL


for (d in 1:3) {
  M <- rbind(
    M,
    matrix(
      unlist(lapply(dgn[[d]], function(l) {
        if (!is.null(names(l))) c(max(l$rhat), min(l$ess))
      })),
      ncol = 2,
      byrow = T
    )
  )
  P <- cbind(P, dgn[[d]]$pk)
}

summary(M[, 1])
summary(M[, 2])

colnames(P) <- names

cookies <- unique(as.character(data$cookie))

# Assume P is your matrix and labnames is length 74
df <- as.data.frame(P)
df$id <- 1:nrow(df)
df$label <- cookies

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(
    cols = colnames(P),
    names_to = "variable",
    values_to = "value"
  )


###############################################
###                                         ###
###               Plotting                  ###
###                                         ###
###############################################

pareto <- ggplot(df_long, aes(x = id, y = value)) +
  geom_point() +
  geom_text(aes(label = label), vjust = -0.5, size = 3) +
  # facet_wrap(~ variable, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  geom_hline(
    yintercept = 1,
    linetype = 'dashed',
    col = 'red3',
    linewidth = 0.8
  ) +

  # Polygon (outlier) 5004
  geom_point(
    data = df_long[df_long$label == "5004", ],
    aes(x = id, y = value),
    size = 4,
    shape = 24,
    col = 'black',
    fill = "orange",
    alpha = 0.7
  ) +
  # Polygon (outlier) 5013
  geom_point(
    data = df_long[df_long$label == "5013", ],
    aes(x = id, y = value),
    size = 4,
    shape = 24,
    col = 'black',
    fill = "orange",
    alpha = 0.7
  ) +
  facet_wrap(~variable) +
  labs(x = "", y = "Pareto k-statistic")


setwd(
  "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Figures"
)
ggsave("Figure S_ pareto.pdf", pareto, width = 9, height = 7)
ggsave(
  "Figure S_ pareto.jpg",
  pareto,
  device = "jpeg",
  width = 9,
  height = 7,
  dpi = 300
)


idOut <- apply(P, 2, function(x) cbind(which(x > 1), x[which(x > 1)]))
lapply(idOut, function(l) cookies[l])

save(
  list = c('idOut'),
  file = paste(path, 'Model Results/Route Pairs/Outliers.RData', sep = "")
)
