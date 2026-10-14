################################################################################
########                                                                ########
########    Calculating proportional land cover for each BBS route      ########
########       from the harmonized NLCD/Conus backcasted data           ########
########                                                                ########
################################################################################

# This code calculates the proportion of each land cover category within
# each BBS route
# the period of 1966 - 2023 using two USGS land cover products that have been
# harmonized so that they all match in resolution and land
# cover classes:
# National Land Cover Database land cover data for the conterminous USA:
#    resampled to 250 meters and reclassified to match CONUS backcasted data
# CONUS backcasted data for the conterminous USA was modeled for 1938-1992
#    using various input data. 250m resolution. The classes of NLCD
#    were reclassified to match this data.

# Categories of land cover are the following:
# 0 Unclassified
# 1 Open Water
# 17 Perennial Ice/Snow
# 2 Urban/Developed
# 7 Barren
# 8 Deciduous Forest
# 9 Evergreen Forest
# 10 Mixed Forest
# 12 Shrubland
# 11 Grassland
# 14 Hay/Pasture
# 13 Cultivated Cropland
# 16 Woody Wetland
# 15 Herbaceous Wetland
# 6 Mining

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory


# Load required packages
library(terra)
library(sf)
library(sp)
library(tidyverse)

print("libraries loaded")


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Read in the linestring shapefile of the BBS routes
# setwd("~/share/groups/MAS/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/BBS/Sample_1973to2021_by4")
setwd("/gpfs1/data/iupdate/Modeling_Inputs2/Beta_Birds_Project")
bbs.routes <- sf::st_read("BBS_Rtes_Linestring.shp")

# Pare down file to only be unique route, geometry, country and state columns
bbs.routes <- bbs.routes[, c("CntryNm", "StateNm", "uniq_rt", "geometry")]


# Define the coordinate system
albers <- sp::CRS(SRS_string = "EPSG:5070")

# Make sure bbs.routes are in albers
bbs.routes <- sf::st_transform(bbs.routes, crs = albers)

print("BBS shapefile loaded")


###############################################
###                                         ###
###         Set up blank dataframe          ###
###                                         ###
###############################################

lulc.predictors <- data.frame(
  Year = rep(seq(from = 1973, to = 2021, by = 1), length(bbs.routes$uniq_rt)),
  unique_route = rep(
    bbs.routes$uniq_rt,
    each = length(seq(from = 1973, to = 2021, by = 1))
  ),
  unclassified = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  open.water = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  ice.snow = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  urban.developed = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  barren = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.dec = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.evgrn = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.mixed = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  shrubland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  grassland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  hay.pasture = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  cropland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  wetland.woody = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  wetland.herb = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  mining = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  )
)


# Fill predictor columns with 0 values
lulc.predictors[, c(3:17)] <- 0.00

# Create a dataframe in which to store the area (m2)
lulc.area <- data.frame(
  Year = rep(seq(from = 1973, to = 2021, by = 1), length(bbs.routes$uniq_rt)),
  unique_route = rep(
    bbs.routes$uniq_rt,
    each = length(seq(from = 1973, to = 2021, by = 1))
  ),
  unclassified = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  open.water = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  ice.snow = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  urban.developed = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  barren = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.dec = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.evgrn = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  forest.mixed = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  shrubland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  grassland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  hay.pasture = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  cropland = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  wetland.woody = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  wetland.herb = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  ),
  mining = integer(
    length(bbs.routes$uniq_rt) *
      length(seq(from = 1973, to = 2021, by = 1))
  )
)

# Fill predictor columns with 0 values
lulc.area[, c(3:17)] <- 0.00


print("blank dataframes created")

###############################################
###                                         ###
###        Preparing the parallel job       ###
###                                         ###
###############################################

# Parallelizing over the 49 years
# harmonized conus and nlcd data
# lulc.path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/Land_cover/Harmonized Data NLCD_CONUSb"
lulc.path <- "/gpfs1/data/iupdate/Modeling_Inputs2/Land_Cover/Harmonized_NLCD_Conusb"

# list the land cover files
lulc.files <- list.files(lulc.path)

# Define the array (aka climate variable) ID
i <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

lulc.file <- lulc.files[i]

print("parallel job set up")

###############################################
###                                         ###
###       Extract Proportional Cover        ###
###           for each BBS Route            ###
###############################################

# Read in the land cover raster for the first year
setwd(lulc.path)
lulc.raster <- terra::rast(lulc.file)

# Extract the year of the raster file for use as an index further down
year <- gsub("\\D", "", lulc.file)

# For NLCD reclassified, this will also include the resolution so pare down to
# just the year
year <- substr(year, 1, 4)


# Project the raster to be in the same projection as the bbs routes
# This projection uses meters and is therefore appropriate for calculating planar distances
# Use method "near" (nearest neighbor) which is good for categorical variables

lulc.raster <- terra::project(lulc.raster, "epsg:5070", method = "near")


print("LULC raster projected to albers")

##################################################

#######################################
#### Convert raster to categorical ####

levels(lulc.raster) <- data.frame(
  ID = c(0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17),
  Class = c(
    "unclassified",
    "open.water",
    "urban.developed",
    "mining",
    "barren",
    "forest.dec",
    "forest.evgrn",
    "forest.mixed",
    "grassland",
    "shrubland",
    "cropland",
    "hay.pasture",
    "wetland.herb",
    "wetland.woody",
    "ice.snow"
  )
)

print("Raster converted to categorical")

#######################################

for (j in 1:length(unique(lulc.predictors$unique_route))) {
  # number of routes

  ###################################################
  ### Determine the spatial weights of each route ###

  # Create a buffer around the route (in meters)
  buffered.route <- sf::st_buffer(bbs.routes$geometry[j], dist = 400)

  # Create vector of the linestring (bbs.route)
  tmp.route <- terra::vect(sf::st_geometry(buffered.route))

  # Crop the spatRaster to the linestring (route)
  cropped.lyr <- terra::crop(lulc.raster, tmp.route, snap = "out")

  # Creates a raster of the linestring where the values of the cells is the
  # length in each cell
  x <- terra::rasterizeGeom(tmp.route, cropped.lyr, "count")

  # Calculate the total length of the route
  # (Sum function not working in terra, so have to use raster package here)
  total <- raster::cellStats(raster::raster(x), stat = 'sum')

  # Divide the values of the spatraster by the total length of the route to get a proportion for each cell
  x <- x / total

  ###################################################

  ########################################
  ### Calculate the proportional cover ###

  # Extract the land cover value of each raster cell the BBS route touches
  tmp.df <- terra::extract(
    cropped.lyr,
    tmp.route,
    method = "simple",
    factors = TRUE,
    touches = TRUE
  )

  # Add the proportional cover value of each raster cell to the tmp.df
  tmp.df$Proportion <- unique(terra::values(x))[unique(terra::values(x) > 0)]

  # Add the area covered by each cell, with 250 meter resolution, each raster
  # cell represents 62,500 meters squared
  tmp.df$Area <- 62500

  # Sum the proportions of each cover type
  tmp.df_sum <- tmp.df %>%
    dplyr::group_by(Class) %>%
    dplyr::summarise(Sum = sum(Proportion), Area_m2 = sum(Area))

  ########################################

  #############################################################################
  ### Add the proportional cover for each cover type to the blank dataframe ###

  for (k in 1:length(tmp.df_sum$Class)) {
    # for each cover type represented in the route

    for (m in 3:length(colnames(lulc.predictors))) {
      # each land cover column

      if (colnames(lulc.predictors)[m] == as.character(tmp.df_sum$Class[k])) {
        lulc.predictors[
          lulc.predictors$Year == as.integer(year) &
            lulc.predictors$unique_route ==
              unique(lulc.predictors$unique_route)[j],
          m
        ] <-
          tmp.df_sum$Sum[k]

        lulc.area[
          lulc.area$Year == as.integer(year) &
            lulc.area$unique_route == unique(lulc.area$unique_route)[j],
          m
        ] <-
          tmp.df_sum$Area_m2[k]
      } # end of if statement for filling lulc.predictors dataframe
    } # end of the m loop
  } # end of the k loop

  # Clean up workspace
  rm(tmp.route, cropped.lyr, x, total, tmp.df, tmp.df_sum)
  gc()

  print(paste("Year: ", year, "- Route: ", j, "- Completed"))
} # end of j loop


# Save the file
setwd("/gpfs1/data/iupdate/Analysis_Output2/Beta_Birds_Project/Land_Cover")

write.csv(
  lulc.predictors,
  file = paste("BBS_proportion_", as.integer(year), ".csv", sep = ""),
  row.names = FALSE
)

write.csv(
  lulc.area,
  file = paste("BBS_area_", as.integer(year), ".csv", sep = ""),
  row.names = FALSE
)
