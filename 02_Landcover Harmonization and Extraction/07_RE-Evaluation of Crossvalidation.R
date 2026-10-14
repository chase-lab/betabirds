################################################################################
########                                                                ########
########    Evaluating classification consistency between NLCD data     ########
########             and CONUS backcasted data with graphs              ########
########                 after reclassification based on accuracy       ########
################################################################################

# This code evaluates the misclassifications between NLCD and Conus for
# First for 1000 random points drawn for each cover class across the entire
# contiguous US.
# Second for 1000 random points drawn from each cookie/polygon.
# Based on the reclassification (aggregation of classes) that was based on
# where there were systematic misclassifications in between the 2 products.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory


# Load required packages
library(terra)
library(tidyverse)
library(ggplot2)
library(viridis)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Define path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/"

# List yearly classification files (per cookie)
yr.files <- list.files(paste(path, "Cookie Classification by year", sep = ""))

# List yearly classification files (for all USA)
yr.files_usa <- list.files(paste(path, "USA Classification by year", sep = ""))


# Load the shapefile of the cookies
setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/BBS/Sample_1973to2021_by4/Cookies"
)
cookies <- sf::st_read("cookies_routes_assembled_20252305.shp")

# Load the USA shapefile
setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/raw"
)
usa.shape <- sf::st_read("USACANAB.shp")

###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################

# Define the coordinate system
albers <- sp::CRS(SRS_string = "EPSG:5070")

# Set CRS of USA
usa.shape <- sf::st_set_crs(usa.shape, albers)

# Crop to USA
usa.shape <- usa.shape[
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

# Set CRS of cookies
cookies <- sf::st_set_crs(cookies, albers)

# Make a list of the classes
Class <- c(
  "Unclassified",
  "Open Water",
  "Urban/Developed",
  "Barren",
  "Forest",
  "Shrub/Grassland",
  "Agriculture",
  "Wetland",
  "Perennial Ice/Snow"
)

# Create dataframe to store results for cookies
overall.misclass <- data.frame()

# Create dataframe to store results for all USA
overall.misclass_usa <- data.frame()


###############################################
###                                         ###
###    Cross-tabulate misclassifications    ###
###               for cookies               ###
###############################################

for (i in 1:length(yr.files)) {
  # Read in the file
  setwd(paste(path, "Cookie Classification by year", sep = ""))
  yr <- read.csv(yr.files[i])

  # Extract the year
  year <- substr(yr.files[i], 23, 26)

  for (j in 1:length(Class)) {
    # Filter nlcd to each cover
    yr.filtered <- yr %>%
      filter(nlcd_class == Class[j])

    if (length(yr.filtered$nlcd_class) == 0) {
      next
    }

    # Proportion misclassified for each class and each cookie
    misclass <- yr.filtered %>%
      group_by(polygon_id) %>%
      count(conus_class) %>%
      mutate(proportion = n / sum(n))

    # Add the nlcd class
    misclass$nlcd_class <- Class[j]

    # Add the year
    misclass$Year <- year

    # Reorder columns
    misclass <- misclass[, c(6, 1, 5, 2:4)]

    # Add to overall dataframe
    overall.misclass <- rbind(overall.misclass, misclass)
  } # end of j loop

  print(i)
}


# Write the results
setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Cookie Agreement"
)
write.csv(
  overall.misclass,
  "Misclassification by class, cookie, and year.csv",
  row.names = FALSE
)

# Clean up workspace
rm(misclass, overall.misclass, yr, yr.filtered, i, j, year)


###############################################
###                                         ###
###    Cross-tabulate misclassifications    ###
###               for all USA               ###
###############################################

for (i in 1:length(yr.files_usa)) {
  # Read in the file
  setwd(paste(path, "USA Classification by year", sep = ""))
  yr <- read.csv(yr.files_usa[i])

  # Extract the year
  year <- substr(yr.files_usa[i], 20, 23)

  for (j in 1:length(Class)) {
    # Filter nlcd to each cover
    yr.filtered <- yr %>%
      filter(nlcd_class == Class[j])

    if (length(yr.filtered$nlcd_class) == 0) {
      next
    }

    # Proportion misclassified for each class and each cookie
    misclass <- yr.filtered %>%
      count(conus_class) %>%
      mutate(proportion = n / sum(n))

    # Add the nlcd class
    misclass$nlcd_class <- Class[j]

    # Add the year
    misclass$Year <- year

    # Reorder columns
    misclass <- misclass[, c(5, 4, 1:3)]

    # Add to overall dataframe
    overall.misclass_usa <- rbind(overall.misclass_usa, misclass)
  } # end of j loop

  print(i)
}


# Write the results
setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation"
)
write.csv(
  overall.misclass_usa,
  "Misclassification by class and year_USA.csv",
  row.names = FALSE
)
