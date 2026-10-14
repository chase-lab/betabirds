################################################################################
########                                                                ########
########     Merging output of proportional land cover calculations     ########
########       for each BBS route from NLCD & CONUS land cover data     ########
################################################################################

# This code takes output from a script run on the cluster which
# calculated the proportion and the area of each land cover category within
# each BBS route
# from the period of 1966-2019 using the NLCD and CONUS land cover datasets
# and merges them into a single csv file.

# There is not quite one file for every year = 46 files --> 1 file

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
library(tidyverse)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# List the files that show proportional land cover for each route
lulc.path <- "~/share/groups/MAS/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/NLCD_CONUS_BBS_Eve"

eve.files <- list.files(lulc.path, pattern = "proportion")

eve.files_area <- list.files(lulc.path, pattern = "area")


###############################################
###                                         ###
###         Set up blank dataframe          ###
###                                         ###
###############################################

lulc.predictors <- data.frame(
  Year = integer(),
  unique_route = character(),
  unclassified = numeric(),
  open.water = numeric(),
  ice.snow = numeric(),
  urban.developed = numeric(),
  barren = numeric(),
  forest.dec = numeric(),
  forest.evgrn = numeric(),
  forest.mixed = numeric(),
  shrubland = numeric(),
  grassland = numeric(),
  hay.pasture = numeric(),
  cropland = numeric(),
  wetland.woody = numeric(),
  wetland.herb = numeric(),
  mining = numeric()
)


lulc.area <- data.frame(
  Year = integer(),
  unique_route = character(),
  unclassified = numeric(),
  open.water = numeric(),
  ice.snow = numeric(),
  urban.developed = numeric(),
  barren = numeric(),
  forest.dec = numeric(),
  forest.evgrn = numeric(),
  forest.mixed = numeric(),
  shrubland = numeric(),
  grassland = numeric(),
  hay.pasture = numeric(),
  cropland = numeric(),
  wetland.woody = numeric(),
  wetland.herb = numeric(),
  mining = numeric()
)


###############################################
###                                         ###
###       Pare each csv and then merge      ###
###                                         ###
###############################################

for (i in 1:length(eve.files_area)) {
  # Read in the year file
  setwd(lulc.path)
  test <- read.csv(eve.files_area[i], header = TRUE)

  # Extract the year from the eve.file name
  file.year <- as.integer(gsub("\\D", "", eve.files_area[i]))

  # Pare the file down to the one year it contains
  test <- test[test$Year == file.year, ]

  # Add the pared file to the blank dataframe
  lulc.area <- rbind(lulc.area, test)

  print(i)
} # end of i loop


for (i in 1:length(eve.files)) {
  # Read in the year file
  setwd(lulc.path)
  test <- read.csv(eve.files[i], header = TRUE)

  # Extract the year from the eve.file name
  file.year <- as.integer(gsub("\\D", "", eve.files[i]))

  # Pare the file down to the one year it contains
  test <- test[test$Year == file.year, ]

  # Add the pared file to the blank dataframe
  lulc.predictors <- rbind(lulc.predictors, test)

  print(i)
} # end of i loop


###############################################
###                                         ###
###            Save the merged CSVs         ###
###                                         ###
###############################################

setwd(
  "~/share/groups/MAS/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/NLCD_CONUS_BBS_Merged"
)

write.csv(
  lulc.predictors,
  "LC Proportions for each BBS Route NLCD_CONUS.csv",
  row.names = FALSE
)
write.csv(
  lulc.area,
  "LC Area m2 for each BBS Route NLCD_CONUS.csv",
  row.names = FALSE
)
