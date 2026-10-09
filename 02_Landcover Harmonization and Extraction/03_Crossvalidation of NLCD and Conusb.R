################################################################################
########                                                                ########
########    Comparing classification consistency between NLCD data      ######## 
########                 and CONUS backcasted data                      ########
################################################################################


# This code compares classification accuracy between National Land Cover 
# Database data (1985-2023) and CONUS backcasted data (1938-1992) for their 
# overlapping period of 1985 - 1992 using stratified sampling.

# Exectuted first for the entire contiguous USA and then again only using
# the cookies (i.e., the areas for which we are interested in the land cover).


# Random sampling can overestimate agreement because nearby points are likely
# to be similar. Spatial blocking ensures that test data are spatially 
# independent from training data, leading to a more honest assessment of 
# inconsistencies.(1)
# However, While using spatial blocks reduces bias due to spatial 
# autocorrelation, some research argues that spatial cross-validation is not
# theoretically justified for map accuracy assessment, as it may not directly
# reflect the accuracy of the map itself but rather the predictive performance
# of a model (5). For land cover map agreement (rather than model prediction),
# traditional stratified random sampling is still widely used and recommended
# (5).
# 1: Stock (2025) Choosing blocks for spatial cross-validation: lessons from a
#         marine remote sensing case study. Frontiers in remote sensing.
# 5: Wadoux et al. (2021) Spatial cross-validation is not the right way to 
#         evaluate map accuracy. Ecological Modelling: 457.



########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory


# Load required packages
library(terra)
library(caret)
library(sf)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Read in resampled (and reclassified) NLCD and Conus backcasted data
path <- "~/share/groups/mas_data/00_data/processed/NLCD_USA_250m/Reclassified to match CONUS_back/"
nlcd.files <- list.files(path)

# Limit to only 1985 - 1992, i.e., the first 8 files
nlcd.files <- nlcd.files[1:8]

# Read in the Conus data
path2 <- "~/share/groups/mas_data/00_data/original/CONUS backcasting"
conus.files <- list.files(path2)

# Limit to only the last 8 files
conus.files <- conus.files[49:56]

# Read in the shapefile of the cookies
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/BBS/Sample_1973to2021_by4/Cookies")
cookies <- sf :: st_read("cookies_routes_assembled_20252305.shp")

###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################

# Create a blank dataframe to save overall results
overall <- data.frame()

# Create a blank dataframe to save by class results
by.class <- data.frame()

# Same for cookie-specific analysis
overall.cookies <- data.frame()
by.class.cookies <- data.frame()

# Remove unique routes from cookies
cookies <- cookies[ , -1]
cookies <- cookies[!duplicated(cookies$ID__200), ]

# Define the coordinate system
albers = sp :: CRS(SRS_string = "EPSG:5070")

# Make sure cookies are in albers
cookies <- sf :: st_transform(cookies, crs = albers)


###############################################
###                                         ###
###     Loop to assess the agreement of     ###
###  Conus to NLCD for each year in all USA ###
###############################################

# Path to save intermediate output
path3 <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/Cookie Crossvalidation/USA Classification by year/"

for (i in 1:length(conus.files)) {
  
  # Read in the year's files
  setwd(path)
  nlcd <- terra :: rast(nlcd.files[i])
  
  setwd(path2)
  conus <- terra :: rast(conus.files[i])
  
  # Extract the year
  year <- substr(nlcd.files[i], 6, 9)
  
  ######################################
  ###   Project the conus raster to  ###
  ###    the same proj as NLCD       ###
  ######################################
  
  nlcd <- terra :: project(nlcd, "epsg:5070", method = "near")
  conus.proj <- terra :: project(conus, "epsg:5070", method = "near")
  
  ######################################
  ###      Define each raster as     ###
  ###           categorical          ###
  ######################################
  
  levels(conus.proj) <- data.frame(ID = c(0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                          15, 16, 17),
                                   Class = c("Unclassified",
                                             "Open Water",
                                             "Urban/Developed",
                                             "Unknown",
                                             "Barren",
                                             "Deciduous Forest",
                                             "Evergreen Forest",
                                             "Mixed Forest",
                                             "Grassland",
                                             "Shrubland",
                                             "Cultivated Cropland",
                                             "Hay/Pasture",
                                             "Herbaceuous Wetland",
                                             "Woody Wetland",
                                             "Perennial Ice/Snow"))
  
  levels(nlcd) <- data.frame(ID = c(0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                    15, 16, 17),
                             Class = c("Unclassified",
                                       "Open Water",
                                       "Urban/Developed",
                                       "Unknown",
                                       "Barren",
                                       "Deciduous Forest",
                                       "Evergreen Forest",
                                       "Mixed Forest",
                                       "Grassland",
                                       "Shrubland",
                                       "Cultivated Cropland",
                                       "Hay/Pasture",
                                       "Herbaceuous Wetland",
                                       "Woody Wetland",
                                       "Perennial Ice/Snow"))
  
  ######################################
  ###      Extract sample points     ###
  ###       to use in assessment     ###
  ######################################
  
  # Define sample points
  set.seed(456)
  
  # The goal is to asses how well CONUS matches NLCD so we sample from NLCD
  # and then extract the Conus values
  # Size represents the number of points for each class
  sample_points <- spatSample(nlcd,
                              size = 1000,
                              na.rm = TRUE, 
                              as.points = TRUE, method = "stratified",
                              each = TRUE)
  
  # Extract values from both datasets
  conus_values <- terra :: extract(conus.proj, sample_points)
  nlcd_values <- terra :: extract(nlcd, sample_points)
  
  # Combine extracted values into one dataframe
  comparison <- data.frame(nlcd_class = nlcd_values$Class,
                           conus_class = conus_values$Class)
  
  # Remove rows with NA or unwanted classes in either column
  comparison <- comparison[!is.na(comparison$nlcd_class) & 
                             !is.na(comparison$conus_class) &
                             !(comparison$nlcd_class %in% c("Unclassified",
                                                            "Unknown")) &
                             !(comparison$conus_class %in% c("Unclassified",
                                                             "Unknown")), ]
  
  # Write the comparison dataframe so that we can assess systematic
  # misclassification
  write.csv(comparison, file = paste(path3, "USA Classification_",
                                     year,
                                     ".csv",
                                     sep = ""),
            row.names = FALSE)
  
  ######################################
  ###      Generate the confusion    ###
  ###              matrix            ###
  ######################################
  
  # Get all class levels present in either raster
  all_levels <- union(levels(factor(comparison$nlcd_class)),
                      levels(factor(comparison$conus_class)))
  
  # confusionMatrix(prediction, reference)
  conf_matrix <- confusionMatrix(data = factor(comparison$conus_class,
                                               levels = all_levels),
                                 reference = factor(comparison$nlcd_class,
                                                    levels = all_levels))
  
  ######################################
  ###          Store results         ###
  ###                                ###
  ######################################
  
  overall.tmp <- data.frame(conf_matrix$overall)
  
  # add Year
  overall.tmp$Year <- year
  
  # Make rownames a column
  overall.tmp <- cbind(rownames(overall.tmp), overall.tmp)
  
  # Remove rownames
  rownames(overall.tmp) <- NULL
  
  # change names
  names(overall.tmp) <- c("Metric", "Conus.match.to.NLCD", "Year")
  
  # add to overall dataframe
  overall <- rbind(overall, overall.tmp)
  
  
  by.class.tmp <- data.frame(conf_matrix$byClass)
  
  # add Year
  by.class.tmp$Year <- year
  
  # Make rownames a column
  by.class.tmp <- cbind(rownames(by.class.tmp), by.class.tmp)
  
  # Remove rownames
  rownames(by.class.tmp) <- NULL
  
  # change names
  names(by.class.tmp)[1] <- "Class"
  
  # add to by.class dataframe
  by.class <- rbind(by.class, by.class.tmp)
  
  # Clean up workspace
  rm(by.class.tmp, comparison, conf_matrix, conus, conus_values,
     conus.proj, nlcd, nlcd_values, overall.tmp, sample_points,
     all_levels, year)

  print(i)
  
}


###############################################
###                                         ###
###           Write the results             ###
###                                         ###
###############################################

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover")
write.csv(overall, "Overall agreement of Conus to NLCD.csv", row.names = FALSE)
write.csv(by.class, "By Class agreement of Conus to NLCD.csv", row.names = FALSE)


###############################################
###                                         ###
###     Loop to assess the agreement of     ###
###     Conus to NLCD for each year in      ###
###            only the cookies             ###
###############################################

# Initialize storage of polygon-level results
polygon_overall <- data.frame()
polygon_byclass <- data.frame()

for (i in 1:length(conus.files)) {
  
  # Read in the year's files
  setwd(path)
  nlcd <- terra :: rast(nlcd.files[i])
  
  setwd(path2)
  conus <- terra :: rast(conus.files[i])
  
  # Extract the year
  year <- substr(nlcd.files[i], 6, 9)
  
  ######################################
  ###   Project the conus raster to  ###
  ###    the same proj as NLCD       ###
  ######################################
  
  nlcd <- terra :: project(nlcd, "epsg:5070", method = "near")
  conus.proj <- terra :: project(conus, "epsg:5070", method = "near")
  
  
  ######################################
  ###      Define each raster as     ###
  ###           categorical          ###
  ######################################
  
  levels(conus.proj) <- data.frame(ID = c(0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                          15, 16, 17),
                                   Class = c("Unclassified",
                                             "Open Water",
                                             "Urban/Developed",
                                             "Unknown",
                                             "Barren",
                                             "Deciduous Forest",
                                             "Evergreen Forest",
                                             "Mixed Forest",
                                             "Grassland",
                                             "Shrubland",
                                             "Cultivated Cropland",
                                             "Hay/Pasture",
                                             "Herbaceuous Wetland",
                                             "Woody Wetland",
                                             "Perennial Ice/Snow"))
  
  levels(nlcd) <- data.frame(ID = c(0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                    15, 16, 17),
                             Class = c("Unclassified",
                                       "Open Water",
                                       "Urban/Developed",
                                       "Unknown",
                                       "Barren",
                                       "Deciduous Forest",
                                       "Evergreen Forest",
                                       "Mixed Forest",
                                       "Grassland",
                                       "Shrubland",
                                       "Cultivated Cropland",
                                       "Hay/Pasture",
                                       "Herbaceuous Wetland",
                                       "Woody Wetland",
                                       "Perennial Ice/Snow"))
  
  ######################################
  ###      Extract sample points     ###
  ###       to use in assessment     ###
  ######################################
  
  # Define sample points
  set.seed(654)
  
  # Convert cookies to spatvector
  cookies.proj.vect <- terra :: vect(cookies)
  
  # The goal is to asses how well CONUS matches NLCD so we sample from NLCD
  # and then extract the Conus values
  # Size represents the number of points for each class
  sample_points <- spatSample(cookies.proj.vect,
                              size = c(rep(1000, 74)),
                              method = "random")
  
  # Extract values from both datasets
  conus_values <- terra :: extract(conus.proj, sample_points)
  nlcd_values <- terra :: extract(nlcd, sample_points)
  
  # Combine extracted values into one dataframe, including polygon ID
  comparison <- data.frame(nlcd_class = nlcd_values$Class,
                           conus_class = conus_values$Class,
                           polygon_id = sample_points$ID__200)
  
  # Save this dataframe for each year
  setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/Cookie Crossvalidation/Classification by year")
  write.csv(comparison, paste("nlcd_conus_comparison_", year, ".csv", sep = ""),
            row.names = FALSE)
  
  # Remove rows with NA or unwanted classes in either column
  # comparison <- comparison[!is.na(comparison$nlcd_class) & 
  #                            !is.na(comparison$conus_class) &
  #                            !(comparison$nlcd_class %in% c("Unclassified",
  #                                                           "Unknown")) &
  #                            !(comparison$conus_class %in% c("Unclassified",
  #                                                            "Unknown")), ]
  
  ######################################
  ###      Generate the confusion    ###
  ###              matrix            ###
  ######################################
  
  # Get all class levels present in either raster
  all_levels <- union(levels(factor(comparison$nlcd_class)),
                      levels(factor(comparison$conus_class)))
  
  # Split data by polygon_ID
  polygon_list <- split(comparison, comparison$polygon_id)
  
  # Calculate metrics for each polygon
  for(poly_id in names(polygon_list)) {
    
    poly_data <- polygon_list[[poly_id]]
    
    # Skip polygons with insufficient data
    if(nrow(poly_data) < 10) next  # Adjust threshold as needed
    
    # Generate confusion matrix
    conf_matrix <- confusionMatrix(
      factor(poly_data$conus_class, levels = all_levels),
      factor(poly_data$nlcd_class, levels = all_levels)
    )
    
    # Store overall metrics
    overall_temp <- data.frame(
      Metric = names(conf_matrix$overall),
      Value = conf_matrix$overall,
      Year = year,
      Polygon_ID = poly_id,
      stringsAsFactors = FALSE
    )
    
    # Store class-specific metrics
    byclass_temp <- data.frame(
      Class = rownames(conf_matrix$byClass),
      conf_matrix$byClass,
      Year = year,
      Polygon_ID = poly_id,
      stringsAsFactors = FALSE
    )
    
    # Combine results
    polygon_overall <- rbind(polygon_overall, overall_temp)
    polygon_byclass <- rbind(polygon_byclass, byclass_temp)
    
    print(poly_id)
  }
  
  
  # confusionMatrix(prediction, reference) all cookies together
  conf_matrix <- confusionMatrix(data = factor(comparison$conus_class,
                                               levels = all_levels),
                                 reference = factor(comparison$nlcd_class,
                                                    levels = all_levels))
  
  ######################################
  ###          Store results         ###
  ###                                ###
  ######################################
  
  overall.tmp <- data.frame(conf_matrix$overall)
  
  # add Year
  overall.tmp$Year <- year
  
  # Make rownames a column
  overall.tmp <- cbind(rownames(overall.tmp), overall.tmp)
  
  # Remove rownames
  rownames(overall.tmp) <- NULL
  
  # change names
  names(overall.tmp) <- c("Metric", "Conus.match.to.NLCD", "Year")
  
  # add to overall dataframe
  overall.cookies <- rbind(overall.cookies, overall.tmp)
  
  
  by.class.tmp <- data.frame(conf_matrix$byClass)
  
  # add Year
  by.class.tmp$Year <- year
  
  # Make rownames a column
  by.class.tmp <- cbind(rownames(by.class.tmp), by.class.tmp)
  
  # Remove rownames
  rownames(by.class.tmp) <- NULL
  
  # change names
  names(by.class.tmp)[1] <- "Class"
  
  # add to by.class dataframe
  by.class.cookies <- rbind(by.class.cookies, by.class.tmp)
  
  # Clean up workspace
  rm(by.class.tmp, comparison, conf_matrix, conus, conus_values,
     conus.proj, nlcd, nlcd_values, overall.tmp, sample_points,
     all_levels, year)
  
  print(i)
  
}


###############################################
###                                         ###
###           Write the results             ###
###                                         ###
###############################################

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/Cookie Crossvalidation")
write.csv(overall.cookies, "Overall agreement of Conus to NLCD.csv", row.names = FALSE)
write.csv(by.class.cookies, "By Class agreement of Conus to NLCD.csv", row.names = FALSE)
write.csv(polygon_overall, "Overall by polygon agreement of Conus to NLCD.csv", row.names = FALSE)
write.csv(polygon_byclass, "By Class by polygon agreement of Conus to NLCD.csv", row.names = FALSE)


