################################################################################
########                                                                ########
########        Resampling NLCD data to a broader extent                ######## 
########                                                                ########
################################################################################


# This code converts/resamples land cover data for the conterminous USA from
# the National Land Cover Database for 1985 - 2023 to match the extent and 
# resolution of land cover
# data in the CONUS backcasting (also from 
# USGS and spanning 1938 - 1992). 
# Resolution of the NLCD data is 30 meters. Resolution of the CONUS data is 
# 250 meters.
# All rasters are in the same projection: Albers Conical Equal Area.

# Produces rasters:
# 'NLCD_2001_250m.tif'
# with files made for each of the years listed above

# land cover classes for NLCD are the following:
# 0	Unclassified
# 11	Open Water
# 12	Perennial Ice / Snow
# 21	Developed, Open Space
# 22	Developed, Low Intensity
# 23	Developed, Medium Intensity
# 24	Developed, High Intensity
# 31	Barren Land (Rock/Sand/Clay)
# 41	Deciduous Forest
# 42	Evergreen Forest
# 43	Mixed Forest
# 52	Shrub/Scrub
# 71	Grassland/herbaceous
# 81	Pasture/Hay
# 82	Cultivated Crops
# 90	Woody Wetlands
# 95	Emergent Herbaceous Wetlands



########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory


# Load required packages
library(terra)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Read in the CONUS raster to which the NLCD data will be resampled to
# Conus raster already projected to match the CRS of the NLCD rasters
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/Land_cover/"

conus <- terra :: rast(paste(path, "CONUS_Backcasting_proj_y1992.tif", sep = ""))



###############################################
###                                         ###
###        Preparing the loop               ###
###                                         ###
###############################################

lulc.path <- "~/share/groups/mas_data/00_data/original/NLCD_USA/"

# list the land cover files 
lulc.files <- list.files(lulc.path, pattern = ".tif")


###############################################
###                                         ###
###        Resampling each raster           ###
###                                         ###
###############################################

start.time <- Sys.time()

for (i in 1:length(lulc.files)) {
  
  lc.rast <- terra :: rast(paste(lulc.path, lulc.files[i], sep = ""))
  
  # Extract the year
  year <- substr(lulc.files[i], 20, 23)
  
  # Resample with method near as it is the preferred method when values 
  # represent land cover classes
  lc.rast <- terra :: resample(lc.rast, conus, method = "near")
  
  # Write the resampled raster as a tif file (so format matches the 
  # Conus rasters)
  path2 <- "~/share/groups/mas_data/00_data/processed/NLCD_USA_250m/"
  terra :: writeRaster(lc.rast,
                       paste(path2, 
                             paste("NLCD_", year, "_250m.tif", sep = ""),
                             sep = ""))
  # Clean up workspace
  rm(lc.rast)
  gc()
  
  print(i)
  
} # end of i loop

end.time <- Sys.time()
total <- end.time - start.time
print(total)
