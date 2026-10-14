################################################################################
###                                                                          ###
###           Plotting: Shifts in Hab Breadth ~ Shifts in Land Cover         ###
###                                                                          ###
################################################################################

# This code creates a plot of the distance each metacommunity shifted in where
# it sat in habitat breadth space (from 1973-1976 and 2017-2020) as a function of
# the distance each metacommunity shifted in where it sat in habitat composition
# space.

# Compute total variation distance for each metacommunity - this will be
# comparable across habitat breadth and land cover.

# Values range from 0 to 1 and can be interpreted as the fraction of the
# composition that would need to be reallocated among categories to turn the
# t1 composition into the t12 composition.
#
# dtv = 0.1 means that 10% of total composition was reallocated

# Clean workspace
rm(list = ls())
gc()


# Load required packages
library(tidyverse)
library(ggplot2)
library(viridis)
library(MetBrewer)
library(forcats)
library(sf)
library(ggmagnify)
library(patchwork)

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


#####################################
###                               ###
###       Data Loading            ###
###                               ###
#####################################

# Read in underlying data of the model
hab.freq_long <- read.csv(
  paste(
    path,
    "01_Analysis/IUCN Habitats/Abundance and Habitat data.csv",
    sep = ""
  ),
  header = TRUE
)

# Read in proportional land cover data for each route in each time interval
lc.prop <- read.csv(
  paste(
    path,
    "01_Analysis/BBS and Landcover/Proportional Land cover for each BBS rt_1973 to 2020.csv",
    sep = ""
  ),
  header = TRUE
)

# Load sf object of metacommunities
cookies <- st_read(paste(
  path,
  "03_Code/Code for Publishing/cookies_routes_assembled_20252305.shp",
  sep = ""
))

# Load USA and Canada shapefile
usa.shape <- sf::st_read(paste(path, "00_Data/raw/USACANAB.shp", sep = ""))


#####################################
###                               ###
###        Data Processing        ###
###                               ###
#####################################

# For future sanity check, sum land cover columns to confirm they sum to 1
lc.prop$row.total <- rowSums(lc.prop[, c(4:18)])
# If not 1 they are all at 0.99

# Filter to only interval 1 and 12 and only the routes that exist
# in cookies
lc.prop_sub <- lc.prop %>%
  filter(
    (Interval == 1 | Interval == 12) &
      unique_route %in% cookies$uniq_rt
  )

# Change name of route column in lc.prop_sub to match cookies
names(lc.prop_sub)[3] <- "uniq_rt"

# Join cookie info with land cover proportions
lc.prop_sub <- lc.prop_sub %>%
  left_join(cookies[, c(1:2)], by = "uniq_rt")

# Combine land cover types based on what we aggregated in the analysis
lc.prop_sub <- lc.prop_sub %>%
  rowwise() %>%
  mutate(
    agriculture = sum(c_across(c(hay.pasture, cropland))),
    forest = sum(c_across(c(forest.dec, forest.evgrn, forest.mixed))),
    shrub.grassland = sum(c_across(c(shrubland, grassland))),
    wetland = sum(c_across(c(wetland.woody, wetland.herb)))
  ) %>%
  ungroup()

# Remove component classes that went into combination
lc.prop_sub <- lc.prop_sub %>%
  select(
    -hay.pasture,
    -cropland,
    -forest.dec,
    -forest.evgrn,
    -forest.mixed,
    -wetland.woody,
    -wetland.herb,
    -shrubland,
    -grassland
  )

# Remove also unclassified and ice.snow and mining as these have only zero values
lc.prop_sub <- lc.prop_sub %>%
  select(-unclassified, -ice.snow, -mining)

# Rearrange columns - remove row total column
lc.prop_sub <- lc.prop_sub[, c(1, 2, 8, 3, 4:6, 10:13, 9)]

# Recheck row sums to make sure sum is still 1 or close to 1
lc.prop_sub$row.total <- rowSums(lc.prop_sub[, c(5:11)])
# It passes the check

# Change some names
names(lc.prop_sub)[2] <- "Time"
names(lc.prop_sub)[3] <- "ID_200"


# Define factors
lc.prop_sub$ID_200 <- as.factor(lc.prop_sub$ID_200)
lc.prop_sub$Time <- as.factor(lc.prop_sub$Time)

# Convert lc.prop_sub to long form
lc.prop_long <- lc.prop_sub %>%
  pivot_longer(
    cols = c(
      open.water,
      urban.developed,
      barren,
      shrub.grassland,
      agriculture,
      forest,
      wetland
    ),
    names_to = "Land.cover.category",
    values_to = "Proportion"
  )

# Aggregate across unique routes to get "community totals"
# First aggregate by ID_200 x Time x breadth x draw to get the sum across
# species
lc.summarized <- lc.prop_long %>%
  group_by(ID_200, Time, Land.cover.category) %>%
  summarize(
    sum_prop = sum(Proportion),
    .groups = "drop"
  )

# Now sum_prop is the total predicted relative abundance across all routes
# observed in that metacommunity for that time and land cover category

# Normalize so that per metacommunity relative abundances sum to 1
lc.normalized <- lc.summarized %>%
  group_by(ID_200, Time) %>%
  mutate(
    total_sum = sum(sum_prop),
    prop = sum_prop / total_sum
  ) %>%
  ungroup()


#####################################
###                               ###
###        Habitat-association    ###
###                               ###
#####################################

# Aggregate species' relative abundances within each metacommunity, time, and
# n_hab_assoc, then normalize within each ID_200 x Time to ensure every
# composition sums to exactly 1.

# Use rel.ab not RelAbundance_beta because RelAbundance_beta has been adjusted
# away from zero for Beta regression and is not the raw composition.

hab_normalized <- hab.freq_long %>%
  group_by(ID_200, Time, n_habitat_associations) %>%
  summarise(
    sum_rel_ab = sum(rel.ab, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(ID_200, Time) %>%
  mutate(
    total_sum = sum(sum_rel_ab),
    prop = sum_rel_ab / total_sum
  ) %>%
  ungroup()

# Ensure every metacomm has all six hab association classes at both times.
# Missing classes should be treated as zero --- not dropped --- because an
# absent category represents zero composition.

hab_levels <- sort(unique(hab_normalized$n_habitat_associations))

hab_complete <- hab_normalized %>%
  mutate(
    n_habitat_associations = as.character(n_habitat_associations)
  ) %>%
  complete(
    ID_200,
    Time,
    n_habitat_associations = as.character(hab_levels),
    fill = list(sum_rel_ab = 0, prop = 0)
  )

# Pivot to one row per ID_200 x habitat class with the two time-specific
# proportions in separate columns
hab_tv_by_class <- hab_complete %>%
  select(ID_200, Time, n_habitat_associations, prop) %>%
  pivot_wider(
    names_from = Time,
    values_from = prop,
    values_fill = 0
  )

# Calculate total variation distance per metacommunity
hab_tv <- hab_tv_by_class %>%
  group_by(ID_200) %>%
  summarise(
    TVD_habitat = 0.5 * sum(abs(t12 - t1)),
    .groups = "drop"
  )


#####################################
###                               ###
###       Land cover              ###
###                               ###
#####################################

# Recode time, complete absent land cover classes as zero, pivot wider,
# calculate the distance
lc_levels <- sort(unique(lc.normalized$Land.cover.category))

lc_complete <- lc.normalized %>%
  mutate(
    Time = recode(
      as.character(Time),
      "1" = "t1",
      "12" = "t12"
    )
  ) %>%
  complete(
    ID_200,
    Time,
    Land.cover.category = lc_levels,
    fill = list(sum_prop = 0, total_sum = 0, prop = 0)
  )

lc_tv_by_class <- lc_complete %>%
  select(ID_200, Time, Land.cover.category, prop) %>%
  pivot_wider(
    names_from = Time,
    values_from = prop,
    values_fill = 0
  )

lc_tv <- lc_tv_by_class %>%
  group_by(ID_200) %>%
  summarise(
    TVD_land_cover = 0.5 * sum(abs(t12 - t1)),
    .groups = "drop"
  )

hab_tv$ID_200 <- as.factor(hab_tv$ID_200)

tv_comparison <- hab_tv %>%
  full_join(lc_tv, by = "ID_200") %>%
  mutate(
    difference = TVD_land_cover - TVD_habitat
  )


summary(tv_comparison$TVD_habitat)
summary(tv_comparison$TVD_land_cover)

tv_comparison %>%
  arrange(desc(difference)) %>%
  print(n = 74)

# Positive difference means the land-cover composition changes more than the
# habitat-association composition for that metacommunity, on this 0-1
# reallocation scale.

#####################################
###                               ###
###        Map Preparation        ###
###                               ###
#####################################

# Filter cookies to have only one row per ID_200
cookies <- cookies %>%
  select(ID__200, Freq, are_km2, geometry) %>%
  distinct()

# Change name of ID_200
names(cookies)[1] <- "ID_200"

# Define the coordinate system
albers <- sp::CRS("epsg:5070")

# Define coordintate system as albers:
usa.shape <- sf::st_set_crs(usa.shape, albers)

# Crop to USA
usa.shape2 <- usa.shape[
  usa.shape$STATE != "AK" &
    usa.shape$STATE != "NWT" &
    usa.shape$STATE != "YT" &
    usa.shape$STATE != "BC" &
    usa.shape$STATE != "QUE" &
    usa.shape$STATE != "LAB" &
    usa.shape$STATE != "SAS" &
    usa.shape$STATE != "xx" &
    usa.shape$STATE != "MAN" &
    usa.shape$STATE != "NFD" &
    usa.shape$STATE != "ONT" &
    usa.shape$STATE != "NS" &
    usa.shape$STATE != "NB" &
    usa.shape$STATE != "PEI" &
    usa.shape$STATE != "ALB",
]

# Join DTV with cookie geometry
cookies <- cookies %>%
  full_join(tv_comparison, by = "ID_200")

#####################################
###                               ###
###          Plotting             ###
###                               ###
#####################################

# Scatter plot
range(tv_comparison$difference, na.rm = TRUE)

max_abs_deviation <- max(
  abs(tv_comparison$difference),
  na.rm = TRUE
)

diff.plot <- ggplot(
  data = tv_comparison,
  aes(x = TVD_land_cover, y = TVD_habitat)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    aes(fill = difference),
    shape = 21,
    color = "black",
    stroke = 0.8,
    size = 6
  ) +
  scale_x_continuous(
    name = "Land-cover compositional shift",
    limits = c(0, 0.32),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name = "Habitat-breadth\ncompositional shift",
    limits = c(0, 0.32),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_equal() +
  scale_fill_gradient2(
    name = "Relative temporal shift\n",
    low = "#2166AC",
    mid = "grey90",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_deviation, max_abs_deviation),
    breaks = c(-max_abs_deviation, 0, max_abs_deviation),
    labels = c(
      "More habitat-breadth\nchange",
      "Equal change",
      "More land-cover\nchange"
    ),
    oob = scales::squish
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
    legend.title = element_text(size = 20, face = "bold", color = "gray30")
  ) +
  theme(legend.text = element_text(size = 20))

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig6_Hab shift v LC shift scatter.jpg",
    sep = ""
  ),
  diff.plot,
  device = "jpeg",
  height = 7,
  width = 11,
  dpi = 300
)


# Map plot
# version that colors polygons by distance
map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(data = cookies, aes(fill = difference), color = "black", size = 1) +
  scale_fill_gradient2(
    name = "Relative temporal shift\n",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_deviation, max_abs_deviation),
    breaks = c(-max_abs_deviation, 0, max_abs_deviation),
    labels = c(
      "More habitat-breadth\nchange",
      "Equal change",
      "More land-cover\nchange"
    ),
    oob = scales::squish
  ) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "none")


ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig6_Hab shift v LC shift map.jpg",
    sep = ""
  ),
  map.cookies,
  device = "jpeg",
  height = 5,
  width = 9,
  dpi = 300
)
