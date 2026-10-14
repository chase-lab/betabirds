################################################################################
###                                                                          ###
###        Change in relative abundance (2017-2020 minus 1973-1976)          ###
###         for each species grouped by # of habitat associations            ###
################################################################################

# This code creates a plot of the model predictions relating relative abundance
# at time interval 1: 1973-1976 and time interval 12: 2017-2020 to the number
# of habitat associations of each bird species detected.

# Predictions made including random effect of species.

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
###     Time x Habitat with       ###
###    random effect of species   ###
###       Data Processing         ###
#####################################

# Predictions that include the random effect of species. (still on the
# response scale = Relative Abundance)

# Define factors
hab.freq_long$ID_200 <- as.factor(hab.freq_long$ID_200)
hab.freq_long$Species1 <- as.factor(hab.freq_long$Species1)
hab.freq_long$Time <- as.factor(hab.freq_long$Time)
hab.freq_long$n_habitat_associations <- as.factor(
  hab.freq_long$n_habitat_associations
)


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

# Extract n for each level of habitat association (for table S3)
species_summary %>%
  filter(n_habitat_associations == "1") %>%
  n_distinct()

species_summary %>%
  filter(n_habitat_associations == "2") %>%
  n_distinct()

species_summary %>%
  filter(n_habitat_associations == "3") %>%
  n_distinct()

species_summary %>%
  filter(n_habitat_associations == "4") %>%
  n_distinct()

species_summary %>%
  filter(n_habitat_associations == "5") %>%
  n_distinct()

species_summary %>%
  filter(n_habitat_associations == "6") %>%
  n_distinct()

#####################################
###     Time x Habitat with       ###
###    random effect of species   ###
###           Plotting            ###
#####################################

# Define habitat association as an ordered factor, order the birds
# by this and then by their mean relative abundance change
dat <- species_summary |>
  mutate(
    n_habitat_associations = factor(
      n_habitat_associations,
      levels = sort(unique(n_habitat_associations))
    )
  ) |>
  arrange(n_habitat_associations, mean)

# Add spacer rows (one between each group) for plotting aestethics
spacers <- dat |>
  group_by(n_habitat_associations) |>
  slice(1) |>
  ungroup() |>
  mutate(
    English_Common_Name = paste0(" ", n_habitat_associations, " spacer"),
    mean = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )

# Combine
dat_plot <- bind_rows(dat, spacers) |>
  arrange(n_habitat_associations, mean) |>
  mutate(
    English_Common_Name = factor(
      English_Common_Name,
      levels = unique(English_Common_Name) # keep desired vertical order
    )
  )

# To create a subset of which lines could be labeled with a common name and/or
# silhouette, reduce the data to only species whose credible intervals do not
# overlap zero
signif <- dat %>%
  filter(!(lower < 0 & upper > 0))

# Write this table to aid in adding labels and silhouettes to plot
write.csv(
  signif,
  paste(
    path,
    "01_Analysis/IUCN Habitats/Significant Species Rel Ab changes.csv",
    sep = ""
  ),
  row.names = FALSE
)


# With a big legend
all.species.legend <- ggplot(
  dat_plot,
  aes(x = mean, y = English_Common_Name, color = n_habitat_associations)
) +
  geom_point(
    data = subset(dat_plot, !grepl("spacer$", English_Common_Name)),
    size = 1
  ) +
  geom_errorbar(
    data = subset(dat_plot, !grepl("spacer$", English_Common_Name)),
    aes(xmin = lower, xmax = upper, color = n_habitat_associations),
    width = 0.2 # thickness of the horizontal bar
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(name = "Change in expected relative abundance (%)") +
  scale_y_discrete(
    name = "",
    # keep all levels including spacers
    limits = levels(dat_plot$English_Common_Name),
    # hide spacer labels so gaps look blank
    labels = function(x) ifelse(grepl("spacer$", x), "", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) + # add extra space above/below
  scale_color_manual(
    name = "# of Habitat\n Associations",
    values = c(
      "1" = met.brewer("Klimt", 6)[6],
      "2" = met.brewer("Klimt", 6)[5],
      "3" = met.brewer("Klimt", 6)[4],
      "4" = met.brewer("Klimt", 6)[3],
      "5" = met.brewer("Klimt", 6)[2],
      "6" = met.brewer("Klimt", 6)[1]
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    axis.title.y = element_text(size = 22, face = "bold", color = "gray30")
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 22, face = "bold", color = "gray30"),
    legend.title = element_text(size = 24, face = "bold", color = "gray30"),
    legend.key.size = unit(1.5, "lines")
  ) + # enlarge key boxes
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4, # larger points in legend
        linetype = 1 # if you want thicker errorbar line type shown
      )
    )
  )


# ggsave(paste(path,
#              "01_Analysis/IUCN Habitats/All species rel abundance change_LEGEND.jpg",
#              sep = ""),
#        all.species.legend,
#        device = "jpeg",
#        height = 11, width = 10, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/FigS22_All species rel abundance change_LEGEND.jpg",
    sep = ""
  ),
  all.species.legend,
  device = "jpeg",
  height = 11,
  width = 10,
  dpi = 300
)


# For manual labeling find the top and bottom 2 species for each habitat group
species.data <- dat_plot %>%
  arrange(n_habitat_associations, mean)

write.csv(
  species.data,
  paste(
    path,
    "01_Analysis/IUCN Habitats/All species rel ab change.csv",
    sep = ""
  ),
  row.names = FALSE
)
