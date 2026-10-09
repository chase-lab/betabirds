################################################################################
########                                                                ########
########    Reclassifying NLCD and CONUS backcasted data based on       ######## 
########               the results of the accuracy assessment           ########
################################################################################


# This code reclassifies land cover data for the conterminous USA from
# the National Land Cover Database for years that it overlaps with Conus:
# 1985 - 1992
# and Conus backcasted products from the same period
# based on which classes should be aggregated (determined as a result of 
# the accuracy assessment)

# land cover classes unchanged
# 0  Unclassified         -->       0  Unclassified
# 1  Open Water           -->       1  Open Water
# 17 Perennial Ice/Snow   -->       17 Perennial Ice/Snow
# 2  Urban/Developed      -->       2  Urban/Developed
# 7  Barren               -->       7  Barren

# land cover classes converted
# 8  Deciduous Forest     -->       8  Forest
# 9  Evergreen Forest     -->       8  Forest
# 10 Mixed Forest         -->       8  Forest
# 12 Shrubland            -->       12 Shrub.Grassland
# 13 Grassland            -->       12 Shrub.Grassland
# 14 Hay/Pasture          -->       14 Agriculture
# 13 Cultivated Cropland  -->       14 Agriculture
# 16 Woody Wetland        -->       16 Wetland
# 15 Herbaceous Wetland   -->       16 Wetland




########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory


# Load required packages
library(terra)


###############################################
###                                         ###
###        Preparing the loop for NLCD      ###
###                                         ###
###############################################

# Define the path from which to read the NLCD rasters
lulc.path <- "~/share/groups/mas_data/00_data/processed/NLCD_USA_250m/Reclassified to match CONUS_back/"

# list the land cover files but we don't want files with .tif.aux pattern 
lulc.files <- list.files(lulc.path)

# Filter to only 1985-1992
lulc.files <- lulc.files[c(1:8)]

# Define the reclassification table
reclass_table <- cbind(
  c(0, 1, 17, 2,  7,  8,  9, 10, 12, 11, 14, 13, 16, 15),
  c(0, 1, 17, 2, 7, 8, 8, 8, 12, 12, 14, 14, 16, 16)
)

# Define path to write reclassified rasters
path2 <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/Land_cover/NLCD_reclass_based_on_accuracy/"


###############################################
###                NLCD                     ###
###        Reclassifying each raster        ###
###                                         ###
###############################################

start.time <- Sys.time()

for (i in 1:length(lulc.files)) {
  
  lc.rast <- terra :: rast(paste(lulc.path, lulc.files[i], sep = ""))
  
  # Extract the year
  year <- substr(lulc.files[i], 6, 9)
  
  # Reclassify the raster based on the reclass table
  lc.rast_reclass <- terra :: classify(lc.rast,
                                       reclass_table)
  
  # Write the reclassified rasters
  terra :: writeRaster(lc.rast_reclass,
                       paste(path2, 
                             paste("NLCD_", year, "_250m_reclassified_accuracy.tif", sep = ""),
                             sep = ""))
  # Clean up workspace
  rm(lc.rast, lc.rast_reclass)
  gc()
  
  print(i)
  
} # end of i loop

end.time <- Sys.time()
total <- end.time - start.time
print(total)

# clean up workspace
rm(reclass_table, end.time, i, lulc.files,lulc.path, path2, start.time,
   total, year)

###############################################
###                                         ###
###        Preparing the loop for CONUS     ###
###                                         ###
###############################################

# Define the path from which to read the Conus rasters
lulc.path <- "~/share/groups/mas_data/00_data/original/CONUS backcasting/"

# list the land cover files but we don't want files with .tif.aux pattern 
lulc.files <- list.files(lulc.path, pattern = ".tif")

# Filter to only 1985-1992
lulc.files <- lulc.files[c(48:55)]

# Define the reclassification table
reclass_table <- cbind(
  c(0, 1, 17, 2,  7,  8,  9, 10, 12, 11, 14, 13, 16, 15),
  c(0, 1, 17, 2, 7, 8, 8, 8, 12, 12, 14, 14, 16, 16)
)

# Define path to write reclassified rasters
path2 <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/Land_cover/Conus_reclass_based_on_accuracy/"


###############################################
###           Conus backcasted              ###
###        Reclassifying each raster        ###
###                                         ###
###############################################

start.time <- Sys.time()

for (i in 1:length(lulc.files)) {
  
  lc.rast <- terra :: rast(paste(lulc.path, lulc.files[i], sep = ""))
  
  # Extract the year
  year <- substr(lulc.files[i], 20, 23)
  
  # Reclassify the raster based on the reclass table
  lc.rast_reclass <- terra :: classify(lc.rast,
                                       reclass_table)
  
  # Write the reclassified rasters
  terra :: writeRaster(lc.rast_reclass,
                       paste(path2, 
                             paste("CONUSb_", year, "_250m_reclassified_accuracy.tif", sep = ""),
                             sep = ""))
  # Clean up workspace
  rm(lc.rast, lc.rast_reclass)
  gc()
  
  print(i)
  
} # end of i loop

end.time <- Sys.time()
total <- end.time - start.time
print(total)