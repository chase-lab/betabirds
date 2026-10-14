################################################################################
########                                                                ########
########        Projecting Conus backcasted raster to NLCD              ########
########                                                                ########
################################################################################

# Prior to resampling the NLCD rasters to 250 meters (to match the CONUS
# backcasted data), we need to make sure that
# the projections of the two rasters match one another.

# Produces rasters:
# 'CONUS_backcasted_1992_projected.tif'

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory


# Load required packages
library(terra)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Read in a CONUS raster to which the NLCD data will be resampled to
# It doesn't matter which one as they all have the same extent and resolution
path <- "~/share/groups/mas_data/00_data/original/CONUS backcasting/"
# path <- "/gpfs1/data/iupdate/Modeling_Inputs2/Land_Cover/"

conus <- terra::rast(paste(path, "CONUS_Backcasting_y1992.tif", sep = ""))


# Read in an NLCD raster
# Also doesn't matter which one as we will project Conus to match the NLCD
lulc.path <- "~/share/groups/mas_data/00_data/original/NLCD_USA/"

lc.rast <- terra::rast(paste(
  path,
  "Annual_NLCD_LndCov_1985_CU_C1V0.tif",
  sep = ""
))

print("Data Loaded")


###############################################
###                                         ###
###     Projecting the Conus raster to      ###
###         the same proj as NLCD           ###
###############################################

start.time <- Sys.time()
print(paste("projection starting at ", start.time, sep = ""))

# project conus to same CRS as lc.rast
target_crs <- crs(lc.rast)

conus.new <- terra::project(conus, target_crs, res = 250)

end.time <- Sys.time()
print(paste("Finished projecting at ", end.time, sep = ""))
print(paste("It took ", end.time - start.time, sep = ""))

# Clean up workspace
rm(conus, lc.rast)
gc()

###############################################
###                                         ###
###          Writing the file               ###
###                                         ###
###############################################

print("Starting to write the file")

setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/Land_cover"
)


terra::writeRaster(
  conus.new,
  "CONUS_Backcasting_proj_y1992.tif",
  overwrite = TRUE
)

print("Finished writing the file")

###############################################
###                                         ###
###           Session Diagnostics           ###
###                                         ###
###############################################

# Sorted list of memory storage
print(sort(sapply(ls(), function(x) {
  format(object.size(get(x)), unit = 'auto')
})))

# Total memory usage
print("Total Memory Used in R session: ")
print(object.size(x = lapply(ls(), get)), units = "Mb")
