###############################################################################
########                                                               ########
########             Figure S : Map of subsample of routes             ########
########                          and cookies                          ########
###############################################################################

# Author: Kimberly Thompson

# This code generated a map for supplementary figure S_ of routes
# that had been sampled at least once within each 4 year interval from
# 1973-2020
# and the cookies (spatial blocks) we generated

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(tidyverse) # Data cleaning
library(sf)
library(ggplot2)

###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"

# Load the shapefile of BBS routes
routes <- sf::st_read(paste(
  path,
  "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_Rtes_Linestring.shp",
  sep = ""
))

# Load the USA/Canada shapefile
usa <- sf::st_read(paste(path, "/00_Data/Raw/USACANAB.shp", sep = ""))

# Load the cookies/spatial blocks
cookies <- sf::st_read(paste(
  path,
  "/00_Data/Processed/cookies_routes_assembled_20252305.shp",
  sep = ""
))


###############################################
###                                         ###
###             Data Processing             ###
###                                         ###
###############################################

# Define the coordinate system
albers <- sp::CRS("epsg:5070")

# Define coordintate system as albers:
usa <- sf::st_set_crs(usa, albers)

# Crop to USA
usa <- usa[
  usa$STATE != "AK" &
    usa$STATE != "NWT" &
    usa$STATE != "YT" &
    usa$STATE != "BC" &
    usa$STATE != "QUE" &
    usa$STATE != "LAB" &
    usa$STATE != "SAS" &
    usa$STATE != "xx" &
    usa$STATE != "MAN" &
    usa$STATE != "NFD" &
    usa$STATE != "ONT" &
    usa$STATE != "NS" &
    usa$STATE != "NB" &
    usa$STATE != "PEI" &
    usa$STATE != "ALB",
]


###############################################
###                                         ###
###             Make the figure             ###
###                                         ###
###############################################

map.plot <- ggplot() +
  geom_sf(data = usa) +
  geom_sf(data = cookies, fill = "orange", alpha = 0.3) +
  geom_sf(data = routes, color = "navyblue", size = 1) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  # theme(axis.ticks.x = element_blank()) +
  # theme(axis.ticks.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


path2 <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Figures/Supplementary Materials/"
ggsave(
  paste(path2, "Map of routes and cookies.jpg", sep = ""),
  plot = map.plot,
  device = "jpg",
  width = 9,
  height = 5,
  dpi = 300
)


###############################################
###                                         ###
###   Histogram of polygon characteristics  ###
###                                         ###
###############################################

# Number of routes per polygon

num.routes <- cookies %>%
  group_by(ID__200) %>%
  summarize(Total = n_distinct(uniq_rt), .groups = "drop")

mean(num.routes$Total)
# 7.79

range(num.routes$Total)
# 3 to 20

route.hist <- ggplot() +
  geom_histogram(data = num.routes, aes(x = Total), bins = 15) +
  geom_vline(xintercept = 7.8, color = "red", linewidth = 2) +
  scale_x_continuous(name = "Number of Routes per Polygon") +
  scale_y_continuous(name = "Polygon count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20))


path2 <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Figures/Supplementary Materials/"
ggsave(
  paste(path2, "Num of Route per Polygon histogram.jpg", sep = ""),
  plot = route.hist,
  device = "jpg",
  width = 5,
  height = 5,
  dpi = 300
)


# Polygon area
pol.area <- cookies %>%
  select(ID__200, are_km2) %>%
  distinct()


mean(pol.area$are_km2)
# 15039.67

range(pol.area$are_km2)
# 10272.09 20781.98

area.hist <- ggplot() +
  geom_histogram(data = pol.area, aes(x = are_km2), bins = 25) +
  geom_vline(xintercept = 15040, color = "red", linewidth = 2) +
  scale_x_continuous(
    name = expression("Polygon area" ~ (km^2)),
    # limits = c(10000, 21000),
    breaks = c(10000, 15000, 20000)
  ) +
  scale_y_continuous(name = "Polygon count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20))

ggsave(
  paste(path2, "Polygon area histogram.jpg", sep = ""),
  plot = area.hist,
  device = "jpg",
  width = 5,
  height = 5,
  dpi = 300
)
