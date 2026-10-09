###############################################################################
########                                                               ########
########                  Determining Ideal Sample                     ######## 
########                                                               ########
###############################################################################


# This code determines the best sample size of BBS routes by checking how many
# routes were consecutively surveyed within 3, 4, and 5 year intervals.


########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() )
gc() #releases memory

library(tidyverse) # Data cleaning
library(sf)



###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Set path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"

# Load State level breeding data
bbs.data <- read.csv(paste(path, 
                           "/00_Data/Processed/BBS/Final_Cleaning/BBS_StateLev.csv",
                           sep = ""),
                     header = TRUE)

routes <- read.csv(paste(path,
                            "/00_Data/Processed/BBS/Final_Cleaning/BBS_Routes.csv",
                            sep = ""), header = TRUE)

weather <- read.csv(paste(path,
                              "/00_Data/Processed/BBS/Final_Cleaning/BBS_Weather_Quality.csv",
                              sep = ""), header = TRUE)


# Load the spatial point file
spatial.data_pts <- sf :: st_read(paste(path,
                                    "/00_Data/Processed/BBS/BBS_All_Rtes_Point.shp",
                                    sep = ""))


# Load the spatial linestring file
spatial.data_ls <- sf :: st_read(paste(path,
                                        "/00_Data/Processed/BBS/BBS_All_Rtes_Linestring.shp",
                                        sep = ""))


###############################################
###                                         ###
###         Determining which routes        ###
###   have been consecutively surveyed      ###
###############################################


# Function to create intermediate data frames and filter routes sampled in all intervals
find_consecutive_routes <- function(start_year, end_year, interval_length, df) {
  # Generate all intervals within the range
  intervals <- seq(start_year, end_year - interval_length + 1, by = interval_length)
  
  # Create an empty data frame to store routes for each interval
  intermediate_df <- data.frame(matrix(NA, nrow = length(unique(df$unique_route)), 
                                       ncol = length(intervals)))
  
  # Iterate through intervals and populate columns with routes sampled in all years
  for (i in seq_along(intervals)) {
    interval_start <- intervals[i]
    interval_end <- interval_start + interval_length - 1
    col_name <- paste0(interval_start, "-", interval_end)
    
    routes_in_interval <- unique(df$unique_route[df$Year >= interval_start & df$Year <= interval_end])
    
    # Add routes to the intermediate data frame
    names(intermediate_df)[i] <- col_name
    intermediate_df[, i] <- c(routes_in_interval,
                              rep(NA, length(unique(df$unique_route)) - length(routes_in_interval)))
  }
  
  # Get a list of unique, non-NA values from each column
  transect_lists <- lapply(intermediate_df, function(col) unique(na.omit(col)))
  
  # Find the intersection across all columns
  common_transects <- Reduce(intersect, transect_lists)
  
  # Apply a function to each column that keeps only common_transects, others become NA
  filtered_df <- as.data.frame(
    lapply(intermediate_df, function(col) ifelse(col %in% common_transects, col, NA))
  )
  
  # Remove rows that are all NA
  filtered_df <- filtered_df[rowSums(!is.na(filtered_df)) > 0, ]
  
  # Vector of final routes
  final_routes <- Reduce(intersect, lapply(intermediate_df, unique))
  
  return(list(intermediate_df = intermediate_df, filtered_df = filtered_df,
              final_routes = final_routes))
}

###########################
#### 3-year intervals   ###
###########################

# Define parameters
start_years <- seq(1967, 1980)         # Starting years
end_year <- max(bbs.data$Year)          # Last year in the dataset
interval_length <- 3                  # Interval length (e.g., 5 years)

# Initialize a list to store results for each start year
results_list_3 <- list()

# Iterate through start years and apply the function
for (start_year in start_years) {
  result <- find_consecutive_routes(start_year, end_year, interval_length, bbs.data)
  
  # Store intermediate data frames and final routes for each start year
  results_list_3[[paste0(start_year, "-", end_year)]] <- result
}


# Extract lengths of final_routes for each element in results_list
final_routes_lengths <- lapply(results_list_3, function(x) length(x$final_routes))

# Convert to a named vector for easier interpretation
final_routes_lengths <- unlist(final_routes_lengths)

# Print the lengths
print(final_routes_lengths)

# Create a dataframe to store results
all.intervals <- data.frame()

all.intervals <- rbind(all.intervals, final_routes_lengths)
names(all.intervals) <- names(final_routes_lengths)




###########################
#### 4-year intervals   ###
###########################

# Define parameters
start_years <- seq(1967, 1980)         # Starting years
end_year <- max(bbs.data$Year)          # Last year in the dataset
interval_length <- 4                  # Interval length (e.g., 5 years)

# Initialize a list to store results for each start year
results_list_4 <- list()

# Iterate through start years and apply the function
for (start_year in start_years) {
  result <- find_consecutive_routes(start_year, end_year, interval_length, bbs.data)
  
  # Store intermediate data frames and final routes for each start year
  results_list_4[[paste0(start_year, "-", end_year)]] <- result
}


# Extract lengths of final_routes for each element in results_list
final_routes_lengths <- lapply(results_list_4, function(x) length(x$final_routes))

# Convert to a named vector for easier interpretation
final_routes_lengths <- unlist(final_routes_lengths)

# Print the lengths
print(final_routes_lengths)

all.intervals <- rbind(all.intervals, final_routes_lengths)



###########################
####  5-year intervals   ###
###########################

# Define parameters
start_years <- seq(1967, 1980)         # Starting years
end_year <- max(bbs.data$Year)          # Last year in the dataset
interval_length <- 5                  # Interval length (e.g., 5 years)

# Initialize a list to store results for each start year
results_list_5 <- list()

# Iterate through start years and apply the function
for (start_year in start_years) {
  result <- find_consecutive_routes(start_year, end_year, interval_length, bbs.data)
  
  # Store intermediate data frames and final routes for each start year
  results_list_5[[paste0(start_year, "-", end_year)]] <- result
}


# Extract lengths of final_routes for each element in results_list
final_routes_lengths <- lapply(results_list_5, function(x) length(x$final_routes))

# Convert to a named vector for easier interpretation
final_routes_lengths <- unlist(final_routes_lengths)

# Print the lengths
print(final_routes_lengths)

all.intervals <- rbind(all.intervals, final_routes_lengths)


# Add interval length column to summary
all.intervals$Interval.Length <- c(3, 4, 5)

# Save the all.intervals dataframe
write.csv(all.intervals, paste(path, "01_Analysis/BBS_Preparation/Consec Sampled Routes.csv",
                               sep = ""),
            row.names = FALSE)




###############################################
###                                         ###
###       Choosing sample size and          ###
###         filtering dataframes            ###
###############################################

# We are looking for a balance between sample size and temporal coverage
# 1973 - 2023 in 4 year intervals
# 786 routes

# Check the geographic coverage
routes.to.check <- results_list_4$`1973-2023`$final_routes

# Sample coverage
sample.cov <- spatial.data_pts %>%
  filter(uniq_rt %in% routes.to.check)

plot(sample.cov$geometry)

# The coverage is pretty good  - let's move forward with this sample
# 1973 - 2023 (which is really 2021) in 4 year intervals

# Filter the bbs.data
bbs.data <- bbs.data %>%
  filter(unique_route %in% routes.to.check)

# filter the route data
routes <- routes %>%
  filter(unique_route %in% routes.to.check)

# Filter the weather data
weather <- weather %>%
  filter(unique_route %in% routes.to.check)

# Filter the spatial points
spatial.data_pts <- spatial.data_pts %>%
  filter(uniq_rt %in% routes.to.check)

# Filter the linestrings
spatial.data_ls <- spatial.data_ls %>%
  filter(uniq_rt %in% routes.to.check)

# Make sure that bbs.data and weather are subset to the correct years for the
# sample size
bbs.data <- bbs.data %>%
  filter(bbs.data$Year >= 1973 & bbs.data$Year <= 2021)

weather <- weather %>%
  filter(weather$Year >= 1973 & weather$Year <= 2021)



###############################################
###                                         ###
###           Write the files               ###
###                                         ###
###############################################

write.csv(bbs.data, paste(path,
                          "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_StateLev.csv", 
                          sep = ""), row.names = FALSE)

write.csv(routes, paste(path,
                            "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_Routes.csv",
                            sep = ""), row.names = FALSE)

write.csv(weather, paste(path,
                              "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_Weather_Quality.csv",
                              sep = ""), row.names = FALSE)

st_write(spatial.data_pts, 
         paste(path, "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_Rtes_Point.shp",
               sep = ""))

st_write(spatial.data_ls, 
         paste(path, "/00_Data/Processed/BBS/Sample_1973to2021_by4/BBS_Rtes_Linestring.shp",
               sep = ""))




