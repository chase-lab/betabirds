################################################################################
###                                                                          ###
###                        Plot of habitat associations                      ###
###                                                                          ###
################################################################################

# This code creates a plot of the model predictions relating relative abundance
# at time interval 1: 1973-1976 and time interval 12: 2017-2020 to the number
# of habitat associations of each bird species detected.

# Posterior predictions of the change in relative abundance between later period
# and earlier period (t12 - t1) for groups of species with different numbers of
# habitat associations

# so negative value of change in relative abundance means that for that habitat
# association abundances were reduced in the later time period relative to the
# earlier time period.

# For the model that removed house sparrows as a sensitivity analysis

# Clean workspace
rm(list = ls())
gc()


# Load required packages
library(tidyverse)
library(ggplot2)
library(viridis)
library(brms)
library(tidybayes)
library(MetBrewer)
library(forcats)

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


#####################################
###                               ###
###       Data Loading            ###
###                               ###
#####################################

# Read in the model
model <- readRDS(paste(
  path,
  "01_Analysis/IUCN Habitats/Model Objects/Habitats_wo_housesparrow.rds",
  sep = ""
))

# Read in underlying data of the model
hab.freq_long <- read.csv(
  paste(
    path,
    "01_Analysis/IUCN Habitats/Abundance and Habitat data.csv",
    sep = ""
  ),
  header = TRUE
)


#####################################
###                               ###
###    Overall Time x Habitat     ###
###       Data Processing         ###
#####################################

# Define factors
hab.freq_long$ID_200 <- as.factor(hab.freq_long$ID_200)
hab.freq_long$Species1 <- as.factor(hab.freq_long$Species1)
hab.freq_long$Time <- as.factor(hab.freq_long$Time)
hab.freq_long$n_habitat_associations <- as.factor(
  hab.freq_long$n_habitat_associations
)


# Draw predicted means at each time x habitat combination on the response scale
# then subtract t1 from t12 within each habitat-association level

# Set up prediction grid
hab_levels <- levels(hab.freq_long$n_habitat_associations)
hab_levels <- hab_levels[-7]
time_levels <- levels(hab.freq_long$Time)

nd_diff <- expand_grid(
  Time = time_levels, # "t1","t12"
  n_habitat_associations = hab_levels
) %>%
  mutate(.row_id = row_number())

# Get fitted draws on the response scale (Rel abundance between 0 and 1)
diff_draws <- nd_diff %>%
  add_epred_draws(
    object = model,
    value = ".value",
    re_formula = NA # omit random effects
  )

# Compute t12 - t1 via a self-join
# split t1 and t12 into separate objects
diff_t1 <- diff_draws %>%
  filter(Time == "t1") %>%
  select(.draw, n_habitat_associations, value_t1 = .value)

diff_t12 <- diff_draws %>%
  filter(Time == "t12") %>%
  select(.draw, n_habitat_associations, value_t12 = .value)

# join by draw and habitat level, then compute difference
diff_draws_diff <- diff_t12 %>%
  inner_join(diff_t1, by = c(".draw", "n_habitat_associations")) %>%
  mutate(diff = value_t12 - value_t1)

# Sanity check
diff_draws_diff %>%
  summarise(
    n = n(),
    n_na = sum(is.na(diff)),
    n_nonna = sum(!is.na(diff))
  )

# Define factor order
diff_draws_diff <- diff_draws_diff %>%
  mutate(
    n_habitat_associations = factor(
      n_habitat_associations,
      levels = c("6", "5", "4", "3", "2", "1")
    )
  ) %>%
  arrange(n_habitat_associations) # 1–6 drawn in this row order

# Sort so that 1,2,3 come first and 6, 5, 4, last
diff_draws_diff <- diff_draws_diff %>%
  mutate(
    n_habitat_associations = factor(
      n_habitat_associations,
      levels = c("1", "2", "3", "4", "5", "6") # legend order
    )
  ) %>%
  arrange(
    n_habitat_associations %in% c("6", "5", "4"), # FALSE (1–3) first, TRUE (4–6) last
    n_habitat_associations # within those, 1<2<3<4<5<6
  )

# check
head(diff_draws_diff$n_habitat_associations) # should start with "1"
tail(diff_draws_diff$n_habitat_associations) # should end with "6","5","4"


# summarize credible intervals for the difference
summary <- diff_draws_diff %>%
  group_by(n_habitat_associations) %>%
  summarize(
    mean.change = mean(diff),
    lower = quantile(diff, 0.025),
    upper = quantile(diff, 0.975)
  )

# Save summary
write.csv(
  summary,
  paste(
    path,
    "01_Analysis/IUCN Habitats/Posterior Dist_Hab Assoc Model_woHousesparrow.csv",
    sep = ""
  ),
  row.names = FALSE
)


#####################################
###                               ###
###      Define Probabilities     ###
###                               ###
#####################################

# Probability of hab 1 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "1"] > 0)

# Probability of hab 2 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "2"] > 0)

# Probability of hab 3 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "3"] > 0)

# Probability of hab 4 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "4"] > 0)

# Probability of hab 5 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "5"] > 0)

# Probability of hab 6 > 0
mean(diff_draws_diff$diff[diff_draws_diff$n_habitat_associations == "6"] > 0)

#####################################
###                               ###
###    Overall Time x Habitat     ###
###            Plot               ###
#####################################

# t12 - t1
# Posterior mass mostly below 0 for hab assoc = 6 means that relative abundance
# was higher at t1 than at t12

# Posterior mass mostly above 0 for hab assoc = 4 means that relative abundance
# was higher at t12 than at t1

hab.plot <- ggplot(
  diff_draws_diff,
  aes(
    x = diff,
    fill = factor(n_habitat_associations),
    colour = factor(n_habitat_associations)
  )
) +
  geom_density(aes(alpha = factor(n_habitat_associations))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  # Define limits and breaks to match original plot
  scale_x_continuous(
    name = "Change in expected relative abundance",
    limits = c(-0.011, 0.0051),
    breaks = c(-0.010, -0.005, 0.000, 0.005),
    labels = c("-0.010", "-0.005", "0.000", "0.005")
  ) +
  scale_y_continuous(name = "Density") +
  scale_color_manual(
    name = "# of Habitat\n Associations",
    values = c(
      "1" = "black",
      "2" = "black",
      "3" = "black",
      "4" = "black",
      "5" = "black",
      "6" = "black"
    )
  ) +
  scale_fill_manual(
    name = "# of Habitat\n Associations",
    values = c(
      "1" = viridis(6)[6],
      "2" = viridis(6)[5],
      "3" = viridis(6)[4],
      "4" = viridis(6)[3],
      "5" = viridis(6)[2],
      "6" = viridis(6)[1]
    )
  ) +
  scale_alpha_manual(
    name = "# of Habitat\n Associations",
    values = c("1" = 0.3, "2" = 0.3, "3" = 0.3, "4" = 0.7, "5" = 0.7, "6" = 0.7)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    legend.title = element_text(size = 18, face = "bold", color = "gray30")
  ) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.position = "right")

# Save the plot
# ggsave(paste(path,
#              "01_Analysis/IUCN Habitats/Overall t12-t1 habitat associations_no_housesparrow.jpg",
#              sep = ""),
#        hab.plot,
#        device = "jpeg",
#        height = 7, width = 13, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for Publishing/Manuscript Figures/FigS21_Overall t12-t1 habitat associations_no_housesparrow.jpg",
    sep = ""
  ),
  hab.plot,
  device = "jpeg",
  height = 7,
  width = 13,
  dpi = 300
)
