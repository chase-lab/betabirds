################################################################################
###                                                                          ###
###           Plotting: Habitat Associations and Metacommunities             ###
###                                                                          ###
################################################################################

# This code creates plots of the model predictions relating relative abundance
# at time interval 1: 1973-1976 and time interval 12: 2017-2020 to the number
# of habitat associations of each bird species detected and examines variation
# in metacommunities.

# Creates a multivariate space with each number of habitat associations as the
# vertex.

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

# Define factors
hab.freq_long$ID_200 <- as.factor(hab.freq_long$ID_200)
hab.freq_long$Species1 <- as.factor(hab.freq_long$Species1)
hab.freq_long$Time <- as.factor(hab.freq_long$Time)
hab.freq_long$n_habitat_associations <- as.factor(
  hab.freq_long$n_habitat_associations
)

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

#####################################
###                               ###
###       Model Predictions       ###
###                               ###
#####################################

# Build the metacommunity, time, breadth grid

# Use the observed species–metacommunity combinations
# so we only predict for species that actually occur in each metacommunity.

# Start from the long data; keep the columns needed for prediction
# This gives one row per ID_200 x Species1 x Time x Breadth
nd_species <- hab.freq_long %>%
  distinct(ID_200, Species1, Time, n_habitat_associations)


# Now get posterior draws of the expected response for each row,
# including both species and ID_200 random effects
species_epred_draws <- nd_species %>%
  add_epred_draws(
    object = model,
    value = "epred", # expected RelAbundance_beta on response scale
    re_formula = NULL # include all REs: Species1 and ID_200
  )


# Each row of nd_species is repeated for many posterior draws
# epred is on the 0-1 proportion scale for a single species in that
# metacommunity, time, and breadth

# Aggregate across species to get "community totals"
# First aggregate by ID_200 x Time x breadth x draw to get the sum across
# species
species_sum_draws <- species_epred_draws %>%
  group_by(.draw, ID_200, Time, n_habitat_associations) %>%
  summarize(
    sum_epred = sum(epred),
    .groups = "drop"
  )


# Now sum_epred is the total predicted relative abundance across all species
# observed in that metacommunity for that time and breadth, in that draw

# Within a given ID_200 x Time x .draw, summing over all breadth levels will
# give the total predicted abundance across all species, which is not
# constrained to 1 because the model treats species independently

# So we will normalize to get the composition of each metacommunity to sum to 1
comp_draws <- species_sum_draws %>%
  group_by(.draw, ID_200, Time) %>%
  mutate(
    total_sum = sum(sum_epred),
    prop = sum_epred / total_sum
  ) %>%
  ungroup()

# prop is the predicted share of total community relative abundance associated
# with each habitat-breadth class within each metacommunity and time.
# For each ID_200 x Time x .draw these values sum to 1.

# Summarize over draws
# comp_summary <- comp_draws %>%
#   group_by(ID_200, Time, n_habitat_associations) %>%
#   summarize(
#     mean_prop = mean(prop),
#     lower = quantile(prop, 0.025),
#     upper = quantile(prop, 0.975),
#     .groups = "drop"
#   )

#####################################
###                               ###
###        Geometry Setup         ###
###                               ###
#####################################

# Assign each habitat-association class to one vertex of a regular hexagon

# e.g. Using angles 90, 30, -30, -90, -150, 150, the six vertices are evenly
# spaced around a circle.

# So if we want 1 and 6 to be opposite one another, we just adjust the
# pairings of the habitat associations and angles accordingly

# Then compute the weighted centroid of those six vertices using the predicted
# composition weights which gives one point for t1 and one for t12 for
# each metacommunity

# with one and 6 next to each other
# vertex_df <- tibble(
#   n_habitat_associations = factor(c("1","2","3","4","5","6"),
#                                   levels = c("1","2","3","4","5","6")),
#   angle = c(90, 30, -30, -90, -150, 150) * pi / 180, # in radians
#   vx = cos(angle),
#   vy = sin(angle)
# )

# with one and 6 opposite
vertex_df <- tibble(
  n_habitat_associations = factor(
    c("1", "3", "5", "6", "4", "2"),
    levels = c("1", "3", "5", "6", "4", "2")
  ),
  angle = c(90, 30, -30, -90, -150, 150) * pi / 180, # in radians
  vx = cos(angle),
  vy = sin(angle)
)


##################################################
###   Finding weighted centroid in t1 and t12  ###
###     for each metacommunity                 ###
##################################################

# Compute per-draw centroids
# Computing for each draw will preserve posterior uncertainty
comp_draws_pts <- comp_draws %>%
  left_join(vertex_df, by = "n_habitat_associations") %>%
  group_by(.draw, ID_200, Time) %>%
  summarize(
    px_draw = sum(prop * vx),
    py_draw = sum(prop * vy),
    .groups = "drop"
  )

# summarize centroids (mean and credible intervals)
centroid_summary <- comp_draws_pts %>%
  group_by(ID_200, Time) %>%
  summarize(
    px_mean = mean(px_draw),
    px_lo = quantile(px_draw, 0.025),
    px_hi = quantile(px_draw, 0.975),
    py_mean = mean(py_draw),
    py_lo = quantile(py_draw, 0.025),
    py_hi = quantile(py_draw, 0.975),
    .groups = "drop"
  )

# Check the distribution of the per-draw centroid differences
# between t1 and t12
# per-draw centroids (comp_draws_pts from above)
t1_pts <- comp_draws_pts %>% filter(Time == "t1")
t12_pts <- comp_draws_pts %>% filter(Time == "t12")

diffs_draws <- inner_join(
  t12_pts,
  t1_pts,
  by = c(".draw", "ID_200"),
  suffix = c("_12", "_1")
) %>%
  mutate(dist = sqrt((px_draw_12 - px_draw_1)^2 + (py_draw_12 - py_draw_1)^2))

# summarize per metacommunity
diffs_summary <- diffs_draws %>%
  group_by(ID_200) %>%
  summarize(
    mean_dist = mean(dist),
    med_dist = median(dist),
    p95_dist = quantile(dist, 0.95),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_dist))

head(diffs_summary, 10)


#####################################
###                               ###
###    Single global hexagon      ###
###                               ###
#####################################

### This is the plot we're working with ###
# full posterior cloud
hex_pts_all <- comp_draws_pts # (.draw, ID_200, Time, px_draw, py_draw)

spokes_df <- vertex_df # just a semantic alias

# Build a label object for labeling vertices of the hexagon
label_df <- vertex_df %>%
  mutate(
    label_x = 1.12 * vx,
    label_y = 1.12 * vy,
    label = as.character(n_habitat_associations)
  )

# Inspect range of centroid means to inform inset boundaries
centroid_summary %>%
  summarise(
    xmin = min(px_mean),
    xmax = max(px_mean),
    ymin = min(py_mean),
    ymax = max(py_mean)
  )

# Hexagon for conceptual figure
hex_blank <- ggplot() +
  # Draw the hexagon
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA,
    colour = "grey50"
  ) +
  # Draw the radial spokes
  geom_segment(
    data = spokes_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey60",
    linewidth = 0.6
  ) +
  # Add labels to each vertex
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = label),
    size = 8,
    colour = "black",
    fontface = "plain"
  ) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), clip = "off") +
  theme_void()

# ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Blank Hexagon.jpg",
#              sep = ""),
#        hex_blank,
#        device = "jpeg",
#        height = 5, width = 8, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig4a_Blank Hexagon.jpg",
    sep = ""
  ),
  hex_blank,
  device = "jpeg",
  height = 5,
  width = 8,
  dpi = 300
)

# Hexagon with ellipses and metacommunity means
p_hex <- ggplot() +
  # Draw the hexagon
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA,
    colour = "grey50"
  ) +
  # Draw the radial spokes
  geom_segment(
    data = spokes_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey60",
    linewidth = 0.6
  ) +
  # Add labels to each vertex
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = label),
    size = 8,
    colour = "black",
    fontface = "plain"
  ) +
  # Ellipse for Time 1
  stat_ellipse(
    data = hex_pts_all %>% filter(Time == "t1"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.95,
    colour = "purple",
    alpha = 0.5,
    linewidth = 2
  ) +
  # Ellipse for time 12
  stat_ellipse(
    data = hex_pts_all %>% filter(Time == "t12"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.95,
    colour = "orange",
    alpha = 0.7,
    linewidth = 2
  ) +
  # Metacommunity mean values for each time interval
  geom_point(
    data = centroid_summary,
    aes(x = px_mean, y = py_mean, colour = Time, fill = Time),
    size = 1,
    shape = 21,
    stroke = 0.5,
    alpha = 0.4
  ) +
  geom_rect(
    aes(xmin = -0.15, xmax = 0.10, ymin = -0.052, ymax = 0.25),
    color = "black",
    fill = NA,
    linewidth = 0.5
  ) +
  scale_color_manual(
    name = "",
    values = c("t1" = "black", "t12" = "black"),
    labels = c("t1" = "1973 - 1976", "t12" = "2017 - 2020")
  ) +
  scale_fill_manual(
    name = "",
    values = c("t1" = "purple", "t12" = "orange"),
    labels = c("t1" = "1973 - 1976", "t12" = "2017 - 2020")
  ) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), clip = "off") +
  theme_void() +
  # theme(legend.position = "none")
  theme(legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 8)))

p_hex

# ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Hab assoc Hexagon_pts and ellipses.jpg",
#              sep = ""),
#        p_hex,
#        device = "jpeg",
#        height = 5, width = 8, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig4b_Hab assoc Hexagon_pts and ellipses.jpg",
    sep = ""
  ),
  p_hex,
  device = "jpeg",
  height = 5,
  width = 8,
  dpi = 300
)


# Zoomed in version
p_hex_zoom <- ggplot() +
  # Draw the radial spokes
  geom_segment(
    data = spokes_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey60",
    linewidth = 0.6
  ) +
  # Ellipse for Time 1
  stat_ellipse(
    data = hex_pts_all %>% filter(Time == "t1"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.95,
    colour = "purple",
    alpha = 0.5,
    linewidth = 2
  ) +
  # Ellipse for time 12
  stat_ellipse(
    data = hex_pts_all %>% filter(Time == "t12"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.95,
    colour = "orange",
    alpha = 0.7,
    linewidth = 2
  ) +
  # Metacommunity mean values for each time interval
  geom_point(
    data = centroid_summary,
    aes(x = px_mean, y = py_mean, colour = Time, fill = Time),
    size = 4,
    shape = 21,
    alpha = 0.4
  ) +
  scale_color_manual(
    name = "",
    values = c("t1" = "black", "t12" = "black"),
    labels = c("t1" = "1973 - 1976", "t12" = "2017 - 2020")
  ) +
  scale_fill_manual(
    name = "",
    values = c("t1" = "purple", "t12" = "orange"),
    labels = c("t1" = "1973 - 1976", "t12" = "2017 - 2020")
  ) +
  coord_equal(xlim = c(-0.15, 0.10), ylim = c(-0.052, 0.25)) +
  theme_void() +
  theme(legend.position = "none")


# ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Hab assoc Hexagon_Zoomed pts and ellipses.jpg",
#              sep = ""),
#        p_hex_zoom,
#        device = "jpeg",
#        height = 5, width = 8, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig4c_Hab assoc Hexagon_Zoomed pts and ellipses.jpg",
    sep = ""
  ),
  p_hex_zoom,
  device = "jpeg",
  height = 5,
  width = 8,
  dpi = 300
)


# To add arrows between each t1 and t12 points, build a wide dataframe
# with one row per metacommunity and columns for both time coordinates
centroids_wide <- centroid_summary %>%
  select(ID_200, Time, px_mean, py_mean) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = c(px_mean, py_mean),
    names_sep = "_"
  )

# To define the arrow direction, define a helper function that identifies the
# category of the first spoke reached from a point in a specified direction.
# s tells you how far you travel along the imaginary arrow extension before
# hitting each spoke.
# u tells you whether the intersection lies on the actual outward spoke from
# the origin, rather than on its opposite, unplotted extension
# the smallest positive s defines the first spoke crossed
first_spoke_crossed <- function(x0, y0, dx, dy, vertex_df, tolerance = 1e-10) {
  cross2d <- function(ax, ay, bx, by) ax * by - ay * bx

  denom <- cross2d(dx, dy, vertex_df$vx, vertex_df$vy)
  s <- cross2d(-x0, -y0, vertex_df$vx, vertex_df$vy) / denom
  u <- cross2d(-x0, -y0, dx, dy) / denom

  valid <- is.finite(s) &
    is.finite(u) &
    abs(denom) > tolerance &
    s > tolerance &
    u >= -tolerance

  if (!any(valid)) {
    return(NA_character_)
  }

  first_index <- which(valid)[which.min(s[valid])]
  as.character(vertex_df$n_habitat_associations[first_index])
}

# Create to.direction and from.direction in a single pipeline
hab_levels <- as.character(vertex_df$n_habitat_associations)

centroids_wide_dir <- centroids_wide %>%
  rowwise() %>%
  mutate(
    dx = px_mean_t12 - px_mean_t1,
    dy = py_mean_t12 - py_mean_t1,

    arrow_length = sqrt(dx^2 + dy^2),

    arrow_angle = atan2(dy, dx),

    # Extend forward from the arrowhead in the same direction
    to.spoke = if (arrow_length < 0.005) {
      "no substantial change"
    } else {
      first_spoke_crossed(
        x0 = px_mean_t12,
        y0 = py_mean_t12,
        dx = dx,
        dy = dy,
        vertex_df = vertex_df
      )
    },

    # Extend backward from the tail, opposite movement direction
    from.spoke = if (arrow_length < 0.005) {
      "no substantial change"
    } else {
      first_spoke_crossed(
        x0 = px_mean_t1,
        y0 = py_mean_t1,
        dx = -dx,
        dy = -dy,
        vertex_df = vertex_df
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    to.spoke = factor(
      to.spoke,
      levels = c(hab_levels, "no substantial change")
    ),
    from.spoke = factor(
      from.spoke,
      levels = c(hab_levels, "no substantial change")
    )
  )


# Validate the definition of to.direction and from.direction on the hexagon
# by plotting a long forward extension from each arrowhead
extension_scale <- 3

forward_extensions <- centroids_wide_dir %>%
  filter(arrow_length >= 0.005) %>%
  mutate(
    x_start = px_mean_t12,
    y_start = py_mean_t12,
    x_end = px_mean_t12 + extension_scale * dx,
    y_end = py_mean_t12 + extension_scale * dy
  )

ggplot() +
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(vx, vy),
    fill = NA,
    colour = "grey40"
  ) +
  geom_segment(
    data = vertex_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey60"
  ) +
  geom_segment(
    data = forward_extensions,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end,
      colour = to.spoke
    ),
    linetype = "dashed",
    alpha = 0.65
  ) +
  coord_equal(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) +
  theme_void()


# Diagnosing NA values in the to.spoke and from.spoke columns
# Isolate the rows
na_to <- centroids_wide_dir %>%
  filter(is.na(to.spoke)) %>%
  select(
    ID_200,
    px_mean_t1,
    py_mean_t1,
    px_mean_t12,
    py_mean_t12,
    dx,
    dy,
    arrow_length,
    arrow_angle,
    to.spoke,
    from.spoke
  )

na_from <- centroids_wide_dir %>%
  filter(is.na(from.spoke)) %>%
  select(
    ID_200,
    px_mean_t1,
    py_mean_t1,
    px_mean_t12,
    py_mean_t12,
    dx,
    dy,
    arrow_length,
    arrow_angle,
    to.spoke,
    from.spoke
  )

na_to
na_from

# Map all failed cases at once
na_extensions <- centroids_wide_dir %>%
  filter(is.na(to.spoke) | is.na(from.spoke)) %>%
  mutate(
    ext_scale = 10,

    # Forward extension from t12
    to_xend = px_mean_t12 + ext_scale * dx,
    to_yend = py_mean_t12 + ext_scale * dy,

    # Backward extension from t1
    from_xend = px_mean_t1 - ext_scale * dx,
    from_yend = py_mean_t1 - ext_scale * dy,

    id_label = as.character(ID_200)
  )

ggplot() +
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA,
    colour = "grey35",
    linewidth = 0.6
  ) +
  geom_segment(
    data = vertex_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey65"
  ) +
  geom_segment(
    data = na_extensions,
    aes(
      x = px_mean_t12,
      y = py_mean_t12,
      xend = to_xend,
      yend = to_yend
    ),
    colour = "#D55E00",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_segment(
    data = na_extensions,
    aes(
      x = px_mean_t1,
      y = py_mean_t1,
      xend = from_xend,
      yend = from_yend
    ),
    colour = "#0072B2",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_text(
    data = na_extensions,
    aes(
      x = (px_mean_t1 + px_mean_t12) / 2,
      y = (py_mean_t1 + py_mean_t12) / 2,
      label = id_label
    ),
    size = 3
  ) +
  coord_equal(xlim = c(-0.15, 0.10), ylim = c(-0.052, 0.25)) +
  theme_void()

# We can see the direction of the habitat association for these NA examples,
# which are NA because the trajectory moves towards the wide part of the pie
# slice rather than across a spoke.
# Overwrite these manually based on arrow plot and diagnostic plot
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "241"] <- "4"
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "272"] <- "2"
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "5036"] <- "3"
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "5046"] <- "1"
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "5048"] <- "3"


# Create another measure of direction from the arrow vector
# For each metacommunity, calculate the movement angle using arctan2 which returns
# the angle between the x-axis and the vector from the origin to (x,y)
# returns the direction of the arrow in radians, measured from the positive
# x axis
centroids_wide_dir <- centroids_wide_dir %>%
  rowwise() %>%
  mutate(
    dx = px_mean_t12 - px_mean_t1,
    dy = py_mean_t12 - py_mean_t1,

    arrow_length = sqrt(dx^2 + dy^2),

    arrow_angle = atan2(dy, dx),

    to.direction = {
      # Smallest circular angular separation between arrow and each vertex
      # This is important because it correctly treats angles near +pi and -pi
      # as adjacent rather than as far apart. e.g., Angle of -179 should be
      # considered very close to a vertex at 180
      angular_difference <- abs(
        atan2(
          sin(arrow_angle - vertex_df$angle),
          cos(arrow_angle - vertex_df$angle)
        )
      )

      as.character(
        vertex_df$n_habitat_associations[which.min(angular_difference)]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    to.direction = factor(
      to.direction,
      levels = as.character(vertex_df$n_habitat_associations)
    )
  )

# For metacommunities with no change, or an extremely tiny arrow, label those
# separately rather than forcing them into a vertex class
centroids_wide_dir <- centroids_wide_dir %>%
  mutate(
    to.direction = if_else(
      arrow_length < 0.005,
      "no substantial change",
      as.character(to.direction)
    ),
    to.direction = factor(
      to.direction,
      levels = c(
        as.character(vertex_df$n_habitat_associations),
        "no substantial change"
      )
    )
  )

# From direction
centroids_wide_dir <- centroids_wide_dir %>%
  rowwise() %>%
  mutate(
    # distances from this arrow tail to each vertex
    from.direction = {
      dists <- sqrt(
        (px_mean_t1 - vertex_df$vx)^2 +
          (py_mean_t1 - vertex_df$vy)^2
      )
      vertex_df$n_habitat_associations[which.min(dists)]
    }
  ) %>%
  ungroup()

# To color as the same colors in the map, merge centroids_wide with
# diffs summary
centroids_wide_withdist <- merge(
  centroids_wide_dir,
  diffs_summary,
  by = "ID_200",
  all.x = TRUE
)


# Zoomed in version with arrows
p_hex_zoom_arrows <- ggplot() +
  # Draw the radial spokes
  geom_segment(
    data = spokes_df,
    aes(x = 0, y = 0, xend = vx, yend = vy),
    linetype = "dotted",
    colour = "grey60",
    linewidth = 0.6
  ) +
  # arrows from t1 to t12
  geom_segment(
    data = centroids_wide_withdist,
    aes(
      x = px_mean_t1,
      y = py_mean_t1,
      xend = px_mean_t12,
      yend = py_mean_t12,
      color = mean_dist
    ),
    arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
    alpha = 0.8,
    linewidth = 0.9
  ) +
  scale_color_gradient2(
    name = "Distance\n",
    low = viridis(100)[10],
    mid = viridis(100)[50],
    high = viridis(100)[90],
    midpoint = mean(diffs_summary$mean_dist)
  ) +
  coord_equal(xlim = c(-0.15, 0.10), ylim = c(-0.052, 0.25)) +
  theme_void() +
  theme(
    legend.position = c(1.1, 0.5), # x, y in [0, 1] npc coordinates
    legend.justification = c("left", "center"),
    legend.box.margin = margin(0, 20, 0, 0)
  ) + # extra right margin
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


# ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Hab assoc Hexagon_Zoomed arrows.jpg",
#              sep = ""),
#        p_hex_zoom_arrows,
#        device = "jpeg",
#        height = 5, width = 8, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig4d_Hab assoc Hexagon_Zoomed arrows.jpg",
    sep = ""
  ),
  p_hex_zoom_arrows,
  device = "jpeg",
  height = 5,
  width = 8,
  dpi = 300
)

# Save distance dataframe
write.csv(
  centroids_wide_withdist,
  paste(
    path,
    "01_Analysis/IUCN Habitats/Multivariate distance_bird comms.csv",
    sep = ""
  ),
  row.names = FALSE
)


# Map of distances

# Join distances with cookie geometry
cookies <- cookies %>%
  full_join(centroids_wide_withdist, by = "ID_200")

# Create new factor based on the vertex the start of the arrow is closest to
# (i.e. the breadth composition in 1973-1976) and the vertex the arrow is
# moving towards (i.e. the vertex the arrow is pointing to in 2017-2020)
to_from_levels <- c(
  "1_1",
  "1_2",
  "2_1",
  "2_2",
  "3_1",
  "3_2",
  "4_1",
  "4_2",
  "4_3"
)

cookies <- cookies %>%
  mutate(
    to_from = paste(
      as.character(to.direction),
      as.character(from.direction),
      sep = "_"
    ),
    to_from = factor(
      to_from,
      levels = to_from_levels
    )
  )

# version that colors polygons by distance
map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(data = cookies, aes(fill = mean_dist)) +
  scale_fill_gradient2(
    name = "Distance\n",
    low = viridis(100)[10],
    mid = viridis(100)[50],
    high = viridis(100)[90],
    midpoint = mean(diffs_summary$mean_dist)
  ) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


### Version that colors polygon by direction it's moving towards

# This will be new version for the supplement

# Clean ordered factor for plotting - only 1, 2, 3, 4 are present
cookies <- cookies %>%
  mutate(to.spoke = factor(to.spoke, levels = c("1", "2", "3", "4")))

map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(
    data = cookies,
    aes(fill = to.spoke),
    # size = mean_dist),
    color = "black",
    alpha = 0.7
  ) +
  scale_fill_manual(
    name = "Trajectory moving\n towards",
    values = c(
      "1" = met.brewer("Austria", 7)[2],
      "2" = met.brewer("Austria", 7)[3],
      "3" = met.brewer("Austria", 7)[4],
      "4" = met.brewer("Austria", 7)[5]
    )
  ) +
  # scale_size_continuous(name = "Distance\n",
  #                       limits = c(0.017, 0.081),
  #                       breaks = c(0.02, 0.04, 0.06, 0.08),
  #                       range = c(0.2, 1.4)) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  theme(legend.position = "none")

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig S_ Trajectory Map.jpg",
    sep = ""
  ),
  map.cookies,
  device = "jpeg",
  height = 5,
  width = 9,
  dpi = 300
)


# Version that colors polygons by where it moved from and where it's going to

# Clean ordered factor for plotting
cookies <- cookies %>%
  mutate(
    to_from = factor(
      to_from,
      levels = c("1_1", "2_2", "2_1", "3_2", "4_3", "3_1", "4_2", "4_1")
    )
  )

map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(
    data = cookies,
    aes(fill = to_from, size = mean_dist),
    color = "black"
  ) +
  scale_fill_manual(
    name = "Habitat Breadth\nDirection\n",
    values = c(
      "1_1" = "grey50",
      "2_2" = "grey50",
      "2_1" = met.brewer("Austria", 7)[5],
      "3_2" = met.brewer("Austria", 7)[5],
      "4_3" = met.brewer("Austria", 7)[5],
      "3_1" = met.brewer("Austria", 7)[7],
      "4_2" = met.brewer("Austria", 7)[7],
      "4_1" = met.brewer("Austria", 7)[4]
    ),
    labels = c(
      "1_1" = "No change: specialist-dominated",
      "2_2" = "No change: generalist-dominated",
      "2_1" = "Toward generalists: from 1 to 2",
      "3_2" = "Toward generalists: from 2 to 3",
      "4_3" = "Toward generalists: from 3 to 4",
      "3_1" = "Toward generalists: from 1 to 3",
      "4_2" = "Toward generalists: from 2 to 4",
      "4_1" = "Toward generalists: from 1 to 4"
    )
  ) +
  scale_size_continuous(
    name = "Distance\n",
    limits = c(0.017, 0.081),
    breaks = c(0.02, 0.04, 0.06, 0.08),
    range = c(0.2, 1.4)
  ) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


# ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Hab Assoc Multivariate Distance Map.jpg",
#              sep = ""),
#        map.cookies,
#        device = "jpeg",
#        height = 5, width = 9, dpi = 300)

ggsave(
  paste(
    path,
    "03_Code/Code for publishing/Manuscript Figures/Fig4e_Hab Assoc Multivariate Distance Map.jpg",
    sep = ""
  ),
  map.cookies,
  device = "jpeg",
  height = 5,
  width = 9,
  dpi = 300
)
