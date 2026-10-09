################################################################################
###                                                                          ###
###            Calculating Abundance and defining bird species               ###
###               detected in first and last time interval                   ###
################################################################################

# This code calculates the abundance of bird species detected in the first
# (1973-1976) and last (2017-2020) time interval and also generates a 
# species list for which to extract habitat associations.


# Clean workspace
rm(list = ls() )
gc() 

# Load required packages
library(tidyverse)
library(sf)

# Define file path
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


#####################################
###                               ###
###         Data Loading          ###
###                               ###
#####################################

# Site x species matrices for Interval 1 (1973-1976) and Interval 12 (2017-2020)
ss.mat.int1 <- 
  read.csv(paste(path,
                 "03_Code/Code for publishing/Site x Species Matrix_1973_1976.csv",
                 sep = ""),
           header = TRUE)

ss.mat.int12 <- 
  read.csv(paste(path,
                 "03_Code/Code for publishing/Site x Species Matrix_2017_2020.csv",
                 sep = ""),
           header = TRUE)


# BBS species list (AOU and species name)
species.list <- 
  read.csv(paste(path,
                 "03_Code/Code for publishing/SpeciesList.csv",
                               sep = ""), header = TRUE)

# Load the cookie information
cookies <- 
  st_read(paste(path,
                "03_Code/Code for publishing/cookies_routes_assembled_20252305.shp",
                sep = ""))


#####################################
###                               ###
###       Data Processing         ###
###                               ###
#####################################

# First reduce rows of the site x species matrix to only those routes used 
# in the analyses

# Create a vector of unique routes
routes <- cookies$uniq_rt

# Filter site x species matrices to only contain those routes
ss.mat.int1 <- ss.mat.int1 %>%
  filter(unique_route %in% routes)

ss.mat.int12 <- ss.mat.int12 %>%
  filter(unique_route %in% routes)

# Change column name to match site x species matrices
names(cookies)[1] <- "unique_route"

# Add cookie column to the site x species matrices
ss.mat.int1 <- merge(ss.mat.int1, as.data.frame(cookies[ , c(1:2)]),
                     by = "unique_route",
                     all.x = TRUE)

ss.mat.int12 <- merge(ss.mat.int12, as.data.frame(cookies[ , c(1:2)]),
                     by = "unique_route",
                     all.x = TRUE)

# Remove geometry column
ss.mat.int1$geometry <- NULL
ss.mat.int12$geometry <- NULL

# Remove Route column and Reorder columns
ss.mat.int1 <- ss.mat.int1[ , c(488, 2:487)]
ss.mat.int12 <- ss.mat.int12[ , c(526, 2:525)]


# Calculate abundance of each species across all routes
abundance1 <- ss.mat.int1 %>%
  group_by(ID__200) %>%
  summarise(
    across(
      .cols = where(is.numeric),   # sum all numeric columns
      .fns  = ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


abundance12 <- ss.mat.int12 %>%
  group_by(ID__200) %>%
  summarise(
    across(
      .cols = where(is.numeric),   # sum all numeric columns
      .fns  = ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Pivot longer
abundance_long1 <- abundance1 %>%
  pivot_longer(
    cols = -ID__200,
    names_to = "Species.code",
    values_to = "abundance"
  )

abundance_long12 <- abundance12 %>%
  pivot_longer(
    cols = -ID__200,
    names_to = "Species.code",
    values_to = "abundance"
  )

# Remove rows where abundance = 0
abundance_long1 <- abundance_long1[abundance_long1$abundance != 0, ]
abundance_long12 <- abundance_long12[abundance_long12$abundance != 0, ]

# Save the abundance dataframes
write.csv(abundance_long1,
          paste(path, "01_Analysis/IUCN Habitats/abundance_1973_1976.csv",
                sep = ""),
          row.names = FALSE)

write.csv(abundance_long12,
          paste(path, "01_Analysis/IUCN Habitats/abundance_2017_2020.csv",
                sep = ""),
          row.names = FALSE)



#####################################
###                               ###
###       Filter Species List     ###
###                               ###
#####################################

# Create vector of AOUs in both abundance dfs
aou1 <- abundance_long1$Species.code
aou12 <- abundance_long12$Species.code

# Remove the leading X
aou1 <- gsub("X", "", aou1)
aou12 <- gsub("X", "", aou12)

# Combine aou vectors and then subset to unique values
aou_all <- unique(c(aou1, aou12))

# Subset the BBS list to the AOUs that are present in abundance 1 and abundance12
species.list_subset <- species.list %>%
  filter(AOU %in% aou_all)

# 505 in aou_all vs. 504 in species list subset
# what are we missing
`%!in%` <- Negate(`%in%`)
which(aou_all %!in% species.list_subset$AOU )
# 244
aou_all[244]
# AOU: 4890
# 4890 is only listed in interval 1 and has an abundance of 13

# Create species column in species list
species.list_subset$Species1 <- paste(species.list_subset$Genus, "",
                                      species.list_subset$Species)

# function that collapses any run of whitespace to a single space and trims
# ends
species.list_subset$Species1 <- str_squish(species.list_subset$Species1)


# Write as csv 
species.list_subset %>%
  distinct(Species1) %>%
  rename(scientific_name = Species1) %>%
  write.csv(paste(path,
                  "01_Analysis/IUCN Habitats/species_list_base.csv",
                  sep = ""),
                  row.names = FALSE)



