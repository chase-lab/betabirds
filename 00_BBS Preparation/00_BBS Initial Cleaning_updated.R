###############################################################################
########                                                               ########
########               Downloading and Cleaning BBS Data               ########
########                                                               ########
###############################################################################

# This code does initial
# cleaning and preprocessing of the North American Breeding Bird
# Survey Dataset including: stop level data, migrants (stop-level format),
# routes, and weather.

# Main preprocessing step is to add a unique route ID to each dataset, which
# simply combines the state ID with the route ID, so that each dataset can
# talk to one another

# Website for data download
# https://www.sciencebase.gov/catalog/item/5ea04e9a82cefae35a129d65

# state level: 1966-2023, 50 stops along each route aggregated into 10-stop
#              bins
# stop level: 1997-2023, data for each of 50 stops recorded

# Download includes:
# States - 62 individual files, 1 for each state and province - 1966-2023
# Stop_Data - 10 files, with 5-9 states represented in each file, data for
#             all 50 stops for each route, officially 1997-2023 but some other
#             older years sneak in
# Migrants - counts of non-breeding birds, stop data format, data for
#             all 50 stops for each route, officially 1997-2023 but some other
#             older years sneak in
# Migrant Summary - counts for non-breeding birds, state data format,
#                   1966-2023
# Routes - Route information plus the lat long coordinate of the route's
#          starting point
# Species List - Links AOU (4-digit species ID code) with order, family, genus
#                and species
# Vehicle Data - vehicle counts and excessive noise ??? no info in metadata
# Weather - route quality information and weather conditions - more info in
#           metadata

# Produces datasets:

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(tidyverse) # Data cleaning


###############################################
###                                         ###
###          Stop-Level Data                ###
###             Processing                  ###
###############################################

# Combine data into 1 CSV file, add a unique identifier for each route

setwd(
  "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/raw/BBS/States"
)
state.list <- list.files(getwd())


# Create blank dataframe
master <- data.frame()

for (i in 1:length(state.list)) {
  tmp.file <- read.csv(state.list[i], header = TRUE)

  master <- rbind(master, tmp.file)

  print(i)
}

# Column RPID represents 3-digit run protocol ID
# 101: Standard BBS Protocol: 3-minute counts, 1 observer, single run per year
# 102: Standard BBS Protocol, Replicate 1: 3-minute counts, 1 observer
#      (same or different person), second run in year
# 103: Standard BBS Protocol, Replicate 1: 3-minute counts, 1 observer
#      (same or different person), second run in year
# 104: Standard BBS Protocol, Replicate 3: 3-minute counts, 1 observer
#      (same or different person), fourth run in year.
# 203: Experimental Protocol: Double-Observer Independent: 3-minute counts,
#      2 observer method, single run per year, Secondary Observer
# 501: Experimental Protocol: Distance and time-of-detection methods: Three
#      1-minute sample periods; two distance bands (0-50 m, >50-400m)
# 502: Experimental Protocol: Distance and time-of-detection methods, Rep1:
#      Three 1-minute sample periods; two distance bands (0-50 m, >50-400m),
#      second run in year.

# Frequency of each protocol type
master %>%
  count(RPID) %>%
  mutate(prop = prop.table(n))

# We will keep only protocol 101 which represents 99.4% of the data
master <- master %>%
  filter(RPID == 101)

# Define both statenum and Route as character columns so that trailing zeros
# are retained
master$StateNum <- as.character(master$StateNum)
master$Route <- as.character(master$Route)

# Create a unique route identifier
master$unique_route <- paste(master$StateNum, master$Route, sep = "_")

# Remove Alaska and Canada
master <- master[master$CountryNum != 124 & master$StateNum != 3, ]


# Reorder columns
master <- master[, c(15, 6, 1:5, 7, 14)]

# clean up workspace
rm(tmp.df, tmp.file, total.years, i, path, route.num, state.list, years)


###############################################
###                                         ###
###                 Routes                  ###
###               Processing                ###
###############################################

path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"

routes <- read.csv(
  paste(path, "/00_Data/Raw/BBS/Routes.csv", sep = ""),
  header = TRUE
)

# Active: active (1) or discontinued (0)

# RouteTypeID: indicates if a route was established along a roadside or
# on a body of water (1 = Roadside, 2 = water)

# RouteTypeDetailID: Indicates route length and selection criteria
# 1: Location of route randomly established, 50 point counts
# 2: Location of route non-randomly established, 50 point counts
# 3: Location of route non-randomly established, < 50 point counts

# Frequency of DetailID types
routes %>%
  count(RouteTypeDetailID) %>%
  mutate(prop = prop.table(n))

# Remove type 3 which accounts for 0.003% of routes
routes <- routes %>%
  filter(RouteTypeDetailID != 3)

# Define both statenum and Route as character columns so that trailing zeros
# are retained
routes$StateNum <- as.character(routes$StateNum)
routes$Route <- as.character(routes$Route)

# Create a unique route identifier
routes$unique_route <- paste(routes$StateNum, routes$Route, sep = "_")

# Remove Alaska and Canada
routes <- routes %>%
  filter(CountryNum != 124 & StateNum != 3)

###############################################
###                                         ###
###            Weather/Quality              ###
###               Processing                ###
###############################################

# RPID, unique_route, year range, routes-to-retain, Total abundance

weather <- read.csv(
  paste(path, "/00_Data/Raw/BBS/Weather.csv", sep = ""),
  header = TRUE
)

# keep only protocol 101
weather <- weather %>%
  filter(RPID == 101)

# Define both statenum and Route as character columns so that trailing zeros
# are retained
weather$StateNum <- as.character(weather$StateNum)
weather$Route <- as.character(weather$Route)

# Create a unique route identifier
weather$unique_route <- paste(weather$StateNum, weather$Route, sep = "_")

# Quality Current ID: indicates if the sampling took place under suitable
# weather conditions and met suitable time, data, and route completion criteria
# 1 - Data meet quality control criteria
# 0 - Data do not meet one or more quality control criteria

# Remove samplings with a 0 for quality control
weather <- weather %>%
  filter(QualityCurrentID == 1)

# Runtype: identifies which data were collected consistent with all standard
# BBS criteria
# 1: consistent
# 0: not consistent

# Remove samplings with a 0 for RunType
weather <- weather %>%
  filter(RunType == 1)

# Remove Alaska and Canada
weather <- weather %>%
  filter(CountryNum != 124 & StateNum != 3)


###############################################
###                                         ###
###      Last step in Initial Cleaning      ###
###                                         ###
###############################################

# Each dataset has been pared to retain only routes with high quality and
# adherence to protocols

length(unique(master$unique_route))
# 4228

length(unique(routes$unique_route))
# 4300

length(unique(weather$unique_route))
# 3942

# The routes in weather should be used to further pare the master and routes
# dataframes

final.routes <- unique(weather$unique_route)

# Filter the weather df according to final.routes
master <- master %>%
  filter(unique_route %in% final.routes)

# Filter the routes df according to final.routes
routes <- routes %>%
  filter(unique_route %in% final.routes)


# Write the resulting files
write.csv(
  master,
  paste(path, "/00_Data/Processed/BBS/1st Cleaning/BBS_StateLev.csv", sep = ""),
  row.names = FALSE
)

write.csv(
  routes,
  paste(path, "/00_Data/Processed/BBS/1st Cleaning/BBS_Routes.csv", sep = ""),
  row.names = FALSE
)

write.csv(
  weather,
  paste(
    path,
    "/00_Data/Processed/BBS/1st Cleaning/BBS_Weather_Quality.csv",
    sep = ""
  ),
  row.names = FALSE
)
