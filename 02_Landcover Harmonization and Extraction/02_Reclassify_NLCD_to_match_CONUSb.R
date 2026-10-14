################################################################################
########                                                                ########
########    Reclassifying NLCD data to match CONUS backcasted data      ########
########                                                                ########
################################################################################

# This code reclassifies land cover data for the conterminous USA from
# the National Land Cover Database for years 1985 - 2023
# to match the class values in the land cover
# data in the CONUS backcasting products (also from
# USGS).

# Produces rasters:
# 'NLCD_2001_250m_reclassified.tif'
# and on for the years listed above

# land cover classes converted from (NLCD original) to (CONUS backcasted)
# 0	Unclassified                    -->       0 Unclassified
# 11	Open Water                    -->       1 Open Water
# 12	Perennial Ice / Snow          -->      17 Perennial Ice/Snow
# 21	Developed, Open Space         -->       2 Urban/Developed
# 22	Developed, Low Intensity      -->       2 Urban/Developed
# 23	Developed, Medium Intensity   -->       2 Urban/Developed
# 24	Developed, High Intensity     -->       2 Urban/Developed
# 31	Barren Land (Rock/Sand/Clay)  -->       7 Barren
# 41	Deciduous Forest              -->       8 Deciduous Forest
# 42	Evergreen Forest              -->       9 Evergreen Forest
# 43	Mixed Forest                  -->      10 Mixed Forest
# 52	Shrub/Scrub                   -->      12 Shrubland
# 71	Grassland/herbaceous          -->      11 Grassland
# 81	Pasture/Hay                   -->      14 Hay/Pasture
# 82	Cultivated Crops              -->      13 Cultivated Cropland
# 90	Woody Wetlands                -->      16 Woody Wetland
# 95	Emergent Herbaceous Wetlands  -->      15 Herbaceous Wetland

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory


# Load required packages
library(terra)


###############################################
###                                         ###
###        Preparing the loop               ###
###                                         ###
###############################################

# Define the path from which to read the rasters
lulc.path <- "~/share/groups/mas_data/00_data/processed/NLCD_USA_250m/"

# list the land cover files but we don't want files with .tif.aux pattern
lulc.files <- list.files(lulc.path, pattern = ".tif")

# Define the reclassification table
reclass_table <- cbind(
  c(0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
  c(0, 1, 17, 2, 2, 2, 2, 7, 8, 9, 10, 12, 11, 14, 13, 16, 15)
)

# Define path to write reclassified rasters
path2 <- "~/share/groups/mas_data/00_data/processed/NLCD_USA_250m/Reclassified to match CONUS_back/"


###############################################
###                                         ###
###        Reclassifying each raster        ###
###                                         ###
###############################################

start.time <- Sys.time()

for (i in 1:length(lulc.files)) {
  lc.rast <- terra::rast(paste(lulc.path, lulc.files[i], sep = ""))

  # Extract the year
  year <- substr(lulc.files[i], 6, 9)

  # Reclassify the raster based on the reclass table
  lc.rast_reclass <- terra::classify(lc.rast, reclass_table)

  # Write the reclassified rasters
  terra::writeRaster(
    lc.rast_reclass,
    paste(
      path2,
      paste("NLCD_", year, "_250m_reclassified.tif", sep = ""),
      sep = ""
    )
  )
  # Clean up workspace
  rm(lc.rast, lc.rast_reclass)
  gc()

  print(i)
} # end of i loop

end.time <- Sys.time()
total <- end.time - start.time
print(total)
