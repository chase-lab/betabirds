################################################################################
###                                                                          ###
###                        Multivariate Land cover                           ###
###                                                                          ###
################################################################################

# This code creates plots of relating relative abundance of each land cover 
# category at time interval 1: 1973-1976 and time interval 12: 2017-2020 to 
# in metacommunities.

# Creates a multivariate space with each land cover category as the 
# vertex.


# Clean workspace
rm(list = ls() )
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

# Read in proportional land cover data for each route in each time interval
lc.prop <- read.csv(paste(path,
                       "01_Analysis/BBS and Landcover/Proportional Land cover for each BBS rt_1973 to 2020.csv",
                       sep = ""),
                    header = TRUE)

# Load sf object of metacommunities
cookies <- st_read(paste(path, "03_Code/Code for Publishing/cookies_routes_assembled_20252305.shp",
                         sep = ""))

# Load USA and Canada shapefile
usa.shape <- sf :: st_read(paste(path, "00_Data/raw/USACANAB.shp", sep = ""))


#####################################
###                               ###
###        Data Processing        ###
###                               ###
#####################################

# For future sanity check, sum land cover columns to confirm they sum to 1
lc.prop$row.total <- rowSums(lc.prop[ , c(4:18)])
# If not 1 they are all at 0.99

# Filter to only interval 1 and 12 and only the routes that exist 
# in cookies
lc.prop_sub <- lc.prop %>%
  filter((Interval == 1 | Interval == 12) &
           unique_route %in% cookies$uniq_rt)

# Change name of route column in lc.prop_sub to match cookies
names(lc.prop_sub)[3] <- "uniq_rt"

# Join cookie info with land cover proportions
lc.prop_sub <- lc.prop_sub %>%
  left_join(cookies[ , c(1:2)], by = "uniq_rt")

# Combine land cover types based on what we aggregated in the analysis
lc.prop_sub <- lc.prop_sub %>%
  rowwise() %>%
  mutate(agriculture = sum(c_across(c(hay.pasture, cropland))),
         forest = sum(c_across(c(forest.dec, forest.evgrn, forest.mixed))),
         shrub.grassland = sum(c_across(c(shrubland, grassland))),
         wetland = sum(c_across(c(wetland.woody, wetland.herb)))
  ) %>%
  ungroup()

# Remove component classes that went into combination
lc.prop_sub <- lc.prop_sub %>%
  select(-hay.pasture, -cropland, -forest.dec, -forest.evgrn, -forest.mixed,
         -wetland.woody, -wetland.herb, -shrubland, -grassland)

# Remove also unclassified and ice.snow and mining as these have only zero values
lc.prop_sub <- lc.prop_sub %>%
  select(-unclassified, -ice.snow, -mining)

# Rearrange columns - remove row total column
lc.prop_sub <- lc.prop_sub[ , c(1, 2, 8, 3, 4:6, 10:13, 9)]

# Recheck row sums to make sure sum is still 1 or close to 1
lc.prop_sub$row.total <- rowSums(lc.prop_sub[ , c(5:11)])
# It passes the check

# Change some names
names(lc.prop_sub)[2] <- "Time"
names(lc.prop_sub)[3] <- "ID_200"


# Define factors
lc.prop_sub$ID_200 <- as.factor(lc.prop_sub$ID_200)
lc.prop_sub$Time <- as.factor(lc.prop_sub$Time)


# Filter cookies to have only one row per ID_200
cookies <- cookies %>%
  select(ID__200, Freq, are_km2, geometry) %>%
  distinct()

# Change name of ID_200
names(cookies)[1] <- "ID_200"

# Define the coordinate system
albers = sp :: CRS("epsg:5070")

# Define coordintate system as albers:
usa.shape <- sf :: st_set_crs( usa.shape, albers )

# Crop to USA
usa.shape2 <- usa.shape[usa.shape$STATE != "AK" & usa.shape$STATE != "NWT" &
                          usa.shape$STATE != "YT" & usa.shape$STATE != "BC" &
                          usa.shape$STATE != "QUE" & usa.shape$STATE != "LAB" &
                          usa.shape$STATE != "SAS" & usa.shape$STATE != "xx" &
                          usa.shape$STATE != "MAN" & usa.shape$STATE != "NFD" &
                          usa.shape$STATE != "ONT" & usa.shape$STATE != "NS" &
                          usa.shape$STATE != "NB" & usa.shape$STATE != "PEI" &
                          usa.shape$STATE != "ALB", ]


#####################################
###                               ###
###       Further Processing      ###
###   to set up metacomm grid     ###
#####################################

# Convert lc.prop_sub to long form
lc.prop_long <- lc.prop_sub %>%
  pivot_longer(cols = c(open.water, urban.developed, barren, shrub.grassland,
                        agriculture, forest, wetland),
               names_to = "Land.cover.category",
               values_to = "Proportion")

# Aggregate across unique routes to get "community totals"
# First aggregate by ID_200 x Time x breadth x draw to get the sum across
# species
lc.summarized <- lc.prop_long %>%
  group_by(ID_200, Time, Land.cover.category) %>%
  summarize(
    sum_prop = sum(Proportion),
    .groups   = "drop"
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

# prop is the share of total relative abundance associated
# with each land cover class within each metacommunity and time.
# For each ID_200 x Time these values sum to 1.


#####################################
###                               ###
###        Geometry Setup         ###
###                               ###
#####################################

# Assign each land cover class to one vertex of a regular septagon

# 360 / 7 = 51.42857

# # e.g. Using angles 90, 38.57143, -12.85714, -64.2857, -115.7143, -167.1429,
#           -141.4286
#           the 7 vertices are evenly 
# spaced around a circle.



# So if we want forest to be at 90, with urban and ag roughly opposite,
# we just adjust the 
# pairings of the habitat associations and angles accordingly

# Then compute the weighted centroid of those seven vertices using the predicted
# composition weights which gives one point for t1 and one for t12 for
# each metacommunity


# with forest and agriculture opposite
vertex_df <- tibble(
  Land.cover.category = factor(c("forest", "shrub.grassland", 
                                 "barren", 
                                 "agriculture", "urban.developed",
                                 "open.water", "wetland"), 
                               levels = c("forest", "shrub.grassland", 
                                          "barren", 
                                          "agriculture", "urban.developed",
                                          "open.water", "wetland")),
  angle = c(90, 38.57143, -12.85714, -64.2857, -115.7143, -167.1429, 141.4286) * pi / 180, # in radians
  vx = cos(angle),
  vy = sin(angle)
)


##################################################
###   Finding weighted centroid in t1 and t12  ###
###     for each metacommunity                 ###
##################################################

# Compute centroids
comp_pts <- lc.normalized %>%
  left_join(vertex_df, by = "Land.cover.category") %>%
  group_by(ID_200, Time) %>%
  summarize(
    px_draw = sum(prop * vx),
    py_draw = sum(prop * vy),
    .groups = "drop"
  )

# Check the distribution of the per-draw centroid differences
# between t1 and t12
# per-draw centroids (comp_draws_pts from above)
t1_pts  <- comp_pts %>% filter(Time == "1")
t12_pts <- comp_pts %>% filter(Time == "12")

diffs_draws <- inner_join(t12_pts, t1_pts, by = c("ID_200"), 
                          suffix = c("_12","_1")) %>%
  mutate(dist = sqrt((px_draw_12 - px_draw_1)^2 + (py_draw_12 - py_draw_1)^2))

# summarize per metacommunity
diffs_summary <- diffs_draws %>%
  group_by(ID_200) %>%
  summarize(
    mean_dist = mean(dist),
    med_dist  = median(dist),
    p95_dist  = quantile(dist, 0.95),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_dist))

head(diffs_summary, 10)

# For easier plotting change Time values to t1 and t12
comp_pts <- comp_pts %>%
  mutate(
    Time = case_when(
      Time == "1"  ~ "t1",
      Time == "12" ~ "t12",
      TRUE         ~ as.character(Time)
    ),
    Time = factor(Time, levels = c("t1", "t12"))
  )

#####################################
###                               ###
###    Single global septagon     ###
###                               ###
#####################################

# summarizes where all metacommunities sit in "composition space" at t1 and t12

# full posterior cloud
sep_pts_all <- comp_pts  # (ID_200, Time, px_draw, py_draw)

spokes_df <- vertex_df  # just a semantic alias

# Define labels 
vertex_df$labels <- c("Forest", "Shrub/\nGrassland", "Barren",
                      "Agriculture", "Urban", "Water", "Wetland")

# Build a label object for labeling vertices of the octagon
label_df <- vertex_df %>%
  mutate(
    label_x = 1.12 * vx,
    label_y = 1.12 * vy,
    label = as.character(labels)
  )

# Shrub/grassland label should be a bit further to the right so it doesn't 
# overlap the plot space
# So change the x value manually
label_df$label_x[label_df$Land.cover.category == "shrub.grassland"] 
# 0.8756512 originally
label_df$label_x[label_df$Land.cover.category == "shrub.grassland"] <- 1.12

# Inspect range of centroid means to inform inset boundaries
comp_pts %>%
  summarise(
    xmin = min(px_draw),
    xmax = max(px_draw),
    ymin = min(py_draw),
    ymax = max(py_draw)
  )

# Septagon for conceptual figure
sep_blank <- ggplot() +
  # Draw the septagon
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA, colour = "grey50"
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

ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Blank Septagon.jpg",
             sep = ""),
       sep_blank,
       device = "jpeg",
       height = 5, width = 8, dpi = 300)



# Septagon with ellipses and metacommunity means
p_sep <- ggplot() +
  # Draw the octagon
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA, colour = "grey50"
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
    data = sep_pts_all %>% filter(Time == "t1"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.90,
    colour = "purple",
    alpha = 0.5,
    linewidth = 2
  ) +
  # Ellipse for time 12
  stat_ellipse(
    data = sep_pts_all %>% filter(Time == "t12"),
    aes(x = px_draw, y = py_draw),
    type = "norm",
    level = 0.90,
    colour = "orange",
    alpha = 0.5,
    linewidth = 2
  ) +
  # Metacommunity mean values for each time interval
  geom_point(
    data = comp_pts,
    aes(x = px_draw, y = py_draw, colour = Time,
        fill = Time),
    size = 3,
    shape = 21,
    stroke = 0.5,
    alpha = 0.6
  ) +
  scale_color_manual(name = "",
                     values = c("t1" = "black",
                                "t12" = "black"),
                     labels = c("t1" = "1973 - 1976",
                                "t12" = "2017 - 2020")) +
  scale_fill_manual(name = "",
                    values = c("t1" = "purple",
                               "t12" = "orange"),
                    labels = c("t1" = "1973 - 1976",
                               "t12" = "2017 - 2020")) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), clip = "off") +
  theme_void() +
  theme(legend.position = c(1.1, 0.5), # x, y in [0, 1] npc coordinates
        legend.justification = c("left", "center"),
        legend.box.margin = margin(0, 20, 0, 0)) + # extra right margin
  theme(legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 8)))

# Legend cuts off but we can use legend from other multivariate plot

ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Land cover Septagon_pts and ellipses.jpg",
             sep = ""),
       p_sep,
       device = "jpeg",
       height = 5, width = 8, dpi = 300)




# To add arrows between each t1 and t12 points, build a wide dataframe
# with one row per metacommunity and columns for both time coordinates
centroids_wide <- comp_pts %>%
  select(ID_200, Time, px_draw, py_draw) %>%
  tidyr::pivot_wider(
    names_from  = Time,
    values_from = c(px_draw, py_draw),
    names_sep   = "_"
  )

# To define the arrow direction, define a helper function that identifies the
# category of the first spoke reached from a point in a specified direction.
# s tells you how far you travel along the imaginary arrow extension before
# hitting each spoke.
# u tells you whether the intersection lies on the actual outward spoke from
# the origin, rather than on its opposite, unplotted extension
# the smallest positive s defines the first spoke crossed
spoke_diagnostics <- function(x0, y0, dx, dy, vertex_df,
                              tolerance = 1e-10) {
  
  cross2d <- function(ax, ay, bx, by) ax * by - ay * bx
  
  denom <- cross2d(dx, dy, vertex_df$vx, vertex_df$vy)
  
  s <- cross2d(-x0, -y0, vertex_df$vx, vertex_df$vy) / denom
  
  u <- cross2d(-x0, -y0, dx, dy) / denom
  
  tibble(
    Land.cover.category = as.character(vertex_df$Land.cover.category),
    vertex_angle_deg = vertex_df$angle * 180 / pi,
    vx = vertex_df$vx,
    vy = vertex_df$vy,
    denom = denom,
    s = s,
    u = u,
    nonparallel = abs(denom) > tolerance,
    ahead = s > tolerance,
    outward_spoke = u >= -tolerance,
    valid = nonparallel & ahead & outward_spoke
  )
}

# Create to.direction and from.direction in a single pipeline
landcover_levels <- as.character(vertex_df$Land.cover.category)

centroids_wide_dir <- centroids_wide %>%
  rowwise() %>%
  mutate(
    dx = px_draw_t12 - px_draw_t1,
    dy = py_draw_t12 - py_draw_t1,
    
    arrow_length = sqrt(dx^2 + dy^2),
    
    arrow_angle = atan2(dy, dx),
    
    # Extend forward from the arrowhead in the same direction
    to.spoke = if (arrow_length < 0.005) {
      "no substantial change"
    } else {
      spoke_diagnostics(
        x0 = px_draw_t12,
        y0 = py_draw_t12,
        dx = dx,
        dy = dy,
        vertex_df = vertex_df
      )
    },
    
    # Extend backward from the tail, opposite movement direction
    from.spoke = if (arrow_length < 0.005) {
      "no substantial change"
    } else {
      spoke_diagnostics(
        x0 = px_draw_t1,
        y0 = py_draw_t1,
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
      levels = c(landcover_levels, "no substantial change")
    ),
    from.spoke = factor(
      from.spoke,
      levels = c(landcover_levels, "no substantial change")
    )
  )


# Validate the definition of to.direction and from.direction on the septagon
# by plotting a long forward extension from each arrowhead
extension_scale <- 3

forward_extensions <- centroids_wide_dir %>%
  filter(arrow_length >= 0.005) %>%
  mutate(
    x_start = px_draw_t12,
    y_start = py_draw_t12,
    x_end = px_draw_t12 + extension_scale * dx,
    y_end = py_draw_t12 + extension_scale * dy
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
      x = x_start, y = y_start,
      xend = x_end, yend = y_end,
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
    px_draw_t1, py_draw_t1,
    px_draw_t12, py_draw_t12,
    dx, dy,
    arrow_length,
    arrow_angle,
    to.spoke,
    from.spoke
  )

na_from <- centroids_wide_dir %>%
  filter(is.na(from.spoke)) %>%
  select(
    ID_200,
    px_draw_t1, py_draw_t1,
    px_draw_t12, py_draw_t12,
    dx, dy,
    arrow_length,
    arrow_angle,
    to.spoke,
    from.spoke
  )

na_to
na_from

# Specifically compare arrow lengths
centroids_wide_dir %>%
  mutate(
    spoke_status = case_when(
      is.na(to.spoke) & is.na(from.spoke) ~ "both NA",
      is.na(to.spoke) ~ "to.spoke NA",
      is.na(from.spoke) ~ "from.spoke NA",
      TRUE ~ "classified"
    )
  ) %>%
  group_by(spoke_status) %>%
  summarise(
    n = n(),
    min_arrow = min(arrow_length),
    median_arrow = median(arrow_length),
    max_arrow = max(arrow_length),
    .groups = "drop"
  )

# Inspect a single problematic to.spoke case
id_check <- na_to$ID_200[1]

row_check <- centroids_wide_dir %>%
  filter(ID_200 == id_check) %>%
  slice(1)

spoke_diagnostics(
  x0 = row_check$px_draw_t12,
  y0 = row_check$py_draw_t12,
  dx = row_check$dx,
  dy = row_check$dy,
  vertex_df = vertex_df
)

# Inspect a single problematic from.spoke case
id_check <- na_from$ID_200[1]

row_check <- centroids_wide_dir %>%
  filter(ID_200 == id_check) %>%
  slice(1)

spoke_diagnostics(
  x0 = row_check$px_draw_t1,
  y0 = row_check$py_draw_t1,
  dx = -row_check$dx,
  dy = -row_check$dy,
  vertex_df = vertex_df
)

# Map all failed cases at once
na_extensions <- centroids_wide_dir %>%
  filter(is.na(to.spoke) | is.na(from.spoke)) %>%
  mutate(
    ext_scale = 3,
    
    # Forward extension from t12
    to_xend = px_draw_t12 + ext_scale * dx,
    to_yend = py_draw_t12 + ext_scale * dy,
    
    # Backward extension from t1
    from_xend = px_draw_t1 - ext_scale * dx,
    from_yend = py_draw_t1 - ext_scale * dy,
    
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
      x = px_draw_t12, y = py_draw_t12,
      xend = to_xend, yend = to_yend
    ),
    colour = "#D55E00",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_segment(
    data = na_extensions,
    aes(
      x = px_draw_t1, y = py_draw_t1,
      xend = from_xend, yend = from_yend
    ),
    colour = "#0072B2",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_text(
    data = na_extensions,
    aes(
      x = (px_draw_t1 + px_draw_t12) / 2,
      y = (py_draw_t1 + py_draw_t12) / 2,
      label = id_label
    ),
    size = 3
  ) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  theme_void()

# We can see that for the two to.spoke NA examples, it's just because the 
# head is too close to spoke but doesn't actually cross it.
# Overwrite these manually based on arrow plot and diagnostic plot
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "129"] <- "agriculture"
centroids_wide_dir$to.spoke[centroids_wide_dir$ID_200 == "5031"] <- "urban.developed"


# Create direction from the arrow vector as a movement angle
# For each metacommunity, calculate the movement angle using arctan2 which returns
# the angle between the x-axis and the vector from the origin to (x,y)
# returns the direction of the arrow in radians, measured from the positive
# x axis
centroids_wide_dir <- centroids_wide_dir %>%
  rowwise() %>%
  mutate(
    
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
        vertex_df$Land.cover.category[which.min(angular_difference)]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    to.direction = factor(
      to.direction,
      levels = as.character(vertex_df$Land.cover.category)
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
        as.character(vertex_df$Land.cover.category),
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
      dists <- sqrt((px_draw_t1 - vertex_df$vx)^2 +
                      (py_draw_t1 - vertex_df$vy)^2)
      vertex_df$Land.cover.category[which.min(dists)]
    }
  ) %>%
  ungroup()

# To color as the same colors in the map, merge centroids_wide with 
# diffs summary
centroids_wide_withdist <- merge(centroids_wide_dir,
                                 diffs_summary,
                                 by = "ID_200",
                                 all.x = TRUE)


# Zoomed in version with arrows
p_sep_arrows <- ggplot() +
  # Draw the octagon
  geom_polygon(
    data = vertex_df %>%
      arrange(angle) %>%
      bind_rows(slice(., 1)),
    aes(x = vx, y = vy),
    fill = NA, colour = "grey50"
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
  # arrows from t1 to t12
  geom_segment(
    data = centroids_wide_withdist,
    aes(
      x    = px_draw_t1,
      y    = py_draw_t1,
      xend = px_draw_t12,
      yend = py_draw_t12,
      color = mean_dist
    ),
    arrow   = arrow(length = unit(0.02, "npc"), type = "closed"),
    alpha   = 0.8,
    linewidth = 0.9
  ) +
  scale_color_gradient2(name = "Distance",
                        low = viridis(100)[10],
                        mid = viridis(100)[50],
                        high = viridis(100)[90],
                        midpoint = mean(diffs_summary$mean_dist)) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), clip = "off") +
  theme_void() +
  theme(legend.position = c(1.1, 0.5), # x, y in [0, 1] npc coordinates
        legend.justification = c("left", "center"),
        legend.box.margin = margin(0, 20, 0, 0)) + # extra right margin
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) 


ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Land cover Septagon_arrows.jpg",
             sep = ""),
       p_sep_arrows,
       device = "jpeg",
       height = 5, width = 8, dpi = 300)


# Save distance dataframe
write.csv(centroids_wide_withdist,
          paste(path, "01_Analysis/Land Cover/Multivariate distance_land cover.csv",
                sep = ""),
          row.names = FALSE)

# Map of distances

# Join distances with cookie geometry
cookies <- cookies %>%
  full_join(centroids_wide_withdist, by = "ID_200")


# Version with polygons colored by distance
map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(data = cookies, aes(fill = mean_dist)) +
  scale_fill_gradient2(name = "Distance",
                       low = viridis(100)[10],
                       mid = viridis(100)[50],
                       high = viridis(100)[90],
                       midpoint = mean(diffs_summary$mean_dist)) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size=18)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) 


ggsave(paste(path, "01_Analysis/Figures/Multivariate Plots/Land cover Multivariate Distance Map.jpg",
             sep = ""),
       map.cookies,
       device = "jpeg",
       height = 5, width = 9, dpi = 300)


# Version that colors polygon by spoke it crosses first

# Clean ordered factor for plotting
cookies <- cookies %>%
  mutate(to.spoke = factor(to.spoke,
                           levels = c("agriculture",
                                      "barren",
                                      "forest",
                                      "open.water",
                                      "shrub.grassland",
                                      "urban.developed",
                                      "wetland")))
map.cookies <- ggplot() +
  geom_sf(data = usa.shape2) +
  geom_sf(data = cookies, aes(fill = to.spoke),
          color = "black", alpha = 0.7) +
  scale_fill_manual(name = "Trajectory moving\n towards",
                    values = c("agriculture" = met.brewer("Lakota", 6)[2],
                               "barren" = met.brewer("Lakota", 6)[1],
                               "forest" = met.brewer("Lakota", 6)[5],
                               "shrub.grassland" = met.brewer("Lakota", 6)[4],
                               "urban.developed" = met.brewer("Lakota", 6)[3],
                               "wetland" = met.brewer("Lakota", 6)[6]),
                    labels = c("agriculture" = "Agriculture",
                               "barren" = "Barren",
                               "forest" = "Forest",
                               "shrub.grassland" = "Shrub/Grassland",
                               "urban.developed" = "Urban",
                               "wetland" = "Wetland")) +
  # scale_size_continuous(name = "Distance\n",
  #                       limits = c(0.01, 0.35),
  #                       breaks = c(0.1, 0.2, 0.3, 0.4),
  #                       range = c(0.2, 1.4)) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size=18)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))


ggsave(paste(path, "03_Code/Code for Publishing/Manuscript Figures/FigS21a_LC trajectory map.jpg",
             sep = ""),
       map.cookies,
       device = "jpeg",
       height = 5, width = 9, dpi = 300)
