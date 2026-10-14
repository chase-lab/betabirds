################################################################################
###                                                                          ###
###         Which species are in the 6-habitat association group             ###
###                                                                          ###
################################################################################

# This code creates a plot of the posterior predictions for how relative
# abundance changed from time interval 12: 2017-2020 to
# time interval 1: 1973-1976 for species with 6 habitat associations (extreme
# generalists)

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
  "01_Analysis/IUCN Habitats/Model Objects/Habitats.rds",
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

# Load dataframe with common names for birds
common <- read.csv(
  paste(
    path,
    "01_Analysis/IUCN Habitats/Summary_Abundance and Habitat.csv",
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


#####################################
###     Time x Habitat with       ###
###    random effect of species   ###
###       Data Processing         ###
#####################################

# Predictions that include the random effect of species. (still on the
# response scale = Relative Abundance)

# Build the prediction grid over ID_200, Time, and habitats
species_levels <- levels(hab.freq_long$Species1)
time_levels <- levels(hab.freq_long$Time)
hab_levels <- levels(hab.freq_long$n_habitat_associations)
hab_levels <- hab_levels[-7]

nd_species <- expand_grid(
  Species1 = species_levels,
  Time = time_levels,
  n_habitat_associations = hab_levels
)

# Get expected values with random effects included - for each species we get
# specific deviation in intercept and slope
# Variation among ID_200 conditional on fixed effects, averaged over species
species_epred_draws <- nd_species %>%
  add_epred_draws(
    object = model,
    value = ".value", # expected RelAbundance_beta on response scale
    re_formula = ~ (Time || Species1) # to average over all species
  )

# Same logic to find the difference between time periods
species_t1 <- species_epred_draws %>%
  filter(Time == "t1") %>%
  select(Species1, n_habitat_associations, .draw, value_t1 = .value)

species_t12 <- species_epred_draws %>%
  filter(Time == "t12") %>%
  select(Species1, n_habitat_associations, .draw, value_t12 = .value)

species_diff_draws <- species_t12 %>%
  inner_join(
    species_t1,
    by = c("Species1", "n_habitat_associations", ".draw")
  ) %>%
  mutate(diff = value_t12 - value_t1)

# To make a plot of all species summarize the mean and 95% CI
species_summary <- species_diff_draws %>%
  select(Species1, n_habitat_associations, diff) %>%
  group_by(Species1) %>%
  summarize(
    mean = mean(diff),
    lower = quantile(diff, 0.025),
    upper = quantile(diff, 0.975)
  )

# Add common names
common <- common %>%
  select(English_Common_Name, Species1) %>%
  distinct()

species_summary <- merge(species_summary, common, by = "Species1", all.x = TRUE)

# Add number of habitat associations
hab.assoc <- hab.freq_long %>%
  select(Species1, n_habitat_associations) %>%
  distinct()

species_summary <- merge(
  species_summary,
  hab.assoc,
  by = "Species1",
  all.x = TRUE
)


# who is in the 6 group
species_summary %>% filter(n_habitat_associations == "6")

six.group <- unique(hab.freq_long$Species1[
  hab.freq_long$n_habitat_associations == 6
])
# [1] Butorides virescens    Falco sparverius       Tyto alba
# [4] Camptostoma imberbe    Euphagus cyanocephalus Passer domesticus
# [7] Catharus ustulatus

# Now it's time to group by species with 6 habitat associations
six.group <- species_diff_draws %>%
  filter(
    Species1 == "Butorides virescens" |
      Species1 == "Falco sparverius" |
      Species1 == "Tyto alba" |
      Species1 == "Camptostoma imberbe" |
      Species1 == "Euphagus cyanocephalus" |
      Species1 == "Catharus ustulatus" |
      Species1 == "Passer domesticus"
  )

# Define the factor levels based on how I want it to plot in ggplot
# What gets plotted last will be in front
six.group$Species1 <- factor(
  six.group$Species1,
  levels = c(
    "Tyto alba",
    "Butorides virescens",
    "Falco sparverius",
    "Catharus ustulatus",
    "Camptostoma imberbe",
    "Euphagus cyanocephalus",
    "Passer domesticus"
  )
)

#####################################
###     Time x Habitat with       ###
###    random effect of species   ###
###           Plotting            ###
#####################################

# Plotting species with 6 habitat associations

sixgroup.plot <- ggplot(
  six.group,
  aes(x = diff, fill = factor(Species1), colour = factor(Species1))
) +
  geom_density(aes(alpha = factor(Species1))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(
    name = "Change in expected relative abundance",
    limits = c(-0.05, 0.01)
  ) +
  scale_y_continuous(name = "Density") +
  scale_color_manual(
    name = "",
    values = c(
      "Butorides virescens" = "black",
      "Camptostoma imberbe" = "black",
      "Catharus ustulatus" = "black",
      "Euphagus cyanocephalus" = "black",
      "Falco sparverius" = "black",
      "Passer domesticus" = "black",
      "Tyto alba" = "black"
    ),
    labels = c(
      "Butorides virescens" = "Green heron",
      "Camptostoma imberbe" = "Northern beardless tyrannulet",
      "Catharus ustulatus" = "Swainson's thrush",
      "Euphagus cyanocephalus" = "Brewer's blackbird",
      "Falco sparverius" = "American kestrel",
      "Passer domesticus" = "House sparrow",
      "Tyto alba" = "Barn owl"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "Butorides virescens" = turbo(100)[7],
      "Camptostoma imberbe" = cividis(7)[6],
      "Catharus ustulatus" = turbo(100)[50],
      "Euphagus cyanocephalus" = cividis(7)[4],
      "Falco sparverius" = turbo(100)[40],
      "Passer domesticus" = turbo(100)[100],
      "Tyto alba" = turbo(100)[20]
    ),
    labels = c(
      "Butorides virescens" = "Green heron",
      "Camptostoma imberbe" = "Northern beardless tyrannulet",
      "Catharus ustulatus" = "Swainson's thrush",
      "Euphagus cyanocephalus" = "Brewer's blackbird",
      "Falco sparverius" = "American kestrel",
      "Passer domesticus" = "House sparrow",
      "Tyto alba" = "Barn owl"
    )
  ) +
  scale_alpha_manual(
    name = "",
    values = c(
      "Butorides virescens" = 0.5,
      "Camptostoma imberbe" = 0.3,
      "Catharus ustulatus" = 0.3,
      "Euphagus cyanocephalus" = 0.3,
      "Falco sparverius" = 0.5,
      "Passer domesticus" = 0.7,
      "Tyto alba" = 0.3
    ),
    labels = c(
      "Butorides virescens" = "Green heron",
      "Camptostoma imberbe" = "Northern beardless tyrannulet",
      "Catharus ustulatus" = "Swainson's thrush",
      "Euphagus cyanocephalus" = "Brewer's blackbird",
      "Falco sparverius" = "American kestrel",
      "Passer domesticus" = "House sparrow",
      "Tyto alba" = "Barn owl"
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_text(size = 22, face = "bold")) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(legend.position = "right")

# Save this plot
# ggsave(paste(path,
#              "01_Analysis/IUCN Habitats/Six group species.jpg",
#              sep = ""),
#        sixgroup.plot,
#        device = "jpeg",
#        height = 7, width = 13, dpi = 300)

# Save this plot
ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/FigS21_Six group species.jpg",
    sep = ""
  ),
  sixgroup.plot,
  device = "jpeg",
  height = 7,
  width = 13,
  dpi = 300
)
