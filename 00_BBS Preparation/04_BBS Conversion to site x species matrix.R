#############################################################################################
########                                                                             ########
########                  Creating Site x Species Matrix with BBS Data               ########
########                                                                             ########
#############################################################################################

# This code converts the sample size selected of the North American Breeding Bird
# Survey Data to site x species matrices, one matrix per year.
# Sample includes routes that have been surveyed at least once within 4-year
# intervals from 1973-2021.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory


library(tidyverse) # Data cleaning
library(R.utils) # time intervals


###############################################
###                                         ###
###      Load the observation data          ###
###                                         ###
###############################################

path <- "~/share/groups/MAS/01_projects/iUpdate/Beta Diversity with Pao_paper/"

# Load the survey data
raw <- read.csv(
  paste(
    path,
    "00_Data/processed/BBS/Sample_1973to2021_by4/BBS_StateLev.csv",
    sep = ""
  ),
  header = TRUE
)


###############################################
###                                         ###
###           Data Manipulation             ###
###                                         ###
###############################################

# Pare df to only unique route, AOU (species), and Species Total, and Year
raw <- raw[, c("unique_route", "AOU", "SpeciesTotal", "Year")]


###############################################
###                                         ###
###      Create a site x species matrix     ###
###         for each year interval          ###
###############################################

# Function to round up to nearest integer
round_up <- function(x) {
  ceiling(x) # ceiling() always rounds up
}

# Define year intervals
year_intervals <- list(
  c(1973, 1976),
  c(1977, 1980),
  c(1981, 1984),
  c(1985, 1988),
  c(1989, 1992),
  c(1993, 1996),
  c(1997, 2000),
  c(2001, 2004),
  c(2005, 2008),
  c(2009, 2012),
  c(2013, 2016),
  c(2017, 2020)
)

# Loop through each interval
for (i in seq_along(year_intervals)) {
  # Get current interval
  start_year <- year_intervals[[i]][1]
  end_year <- year_intervals[[i]][2]

  # Filter by year range
  sub.years <- raw %>%
    filter(Year >= start_year & Year <= end_year)

  # Calculate average abundance across years, then round up
  aggregated_data <- sub.years %>%
    group_by(unique_route, AOU) %>%
    summarize(SpeciesTotal = round_up(mean(SpeciesTotal, na.rm = TRUE))) %>%
    ungroup()

  # Transpose the species into columns
  sub.matrix <- reshape2::dcast(
    aggregated_data,
    unique_route ~ AOU,
    value.var = 'SpeciesTotal',
    fun.aggregate = NULL
  )

  # Replace NA values with 0s
  sub.matrix[is.na(sub.matrix)] <- 0

  # Add X prefix to species ID columns (starting from column 2)
  colnames(sub.matrix)[-1] <- paste0("X", colnames(sub.matrix)[-1])

  # Save this dataframe
  setwd(
    "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/BBS/Sample_1973to2021_by4/Site x Species Matrices"
  )

  write.csv(
    sub.matrix,
    paste(
      "Site x Species Matrix_",
      start_year,
      "_",
      end_year,
      ".csv",
      sep = ""
    ),
    row.names = FALSE
  )

  print(paste("Completed interval ", start_year, "-", end_year, sep = ""))
}


ssmatrix <- list.files(getwd())

for (i in seq_along(ssmatrix)) {
  test <- read.csv(ssmatrix[i])

  print(i)

  print(length(unique(test$unique_route)))
}
