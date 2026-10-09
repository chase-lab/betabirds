###############################################################################
########                                                               ########
########             Aggregating land cover data extracted for         ######## 
########               each BBS route into 4-year intervals            ########
###############################################################################


# This code finds the average land cover (areal cover of each land cover
# category) across 4-year intervals from 1973 - 2020.



########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() )
gc() #releases memory

library(tidyverse) # Data cleaning


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# List the yearly land cover files 
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/LC extracted for BBS routes"
lc.files <- list.files(path, pattern = "area")


###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################


# Define the start of the interval years
intervals <- seq(1973, 2020, by = 4)

# Define all the years
years <- seq(1973, 2020, by = 1)

# Make a key of interval (i) year combos
key <- data.frame(year = years, interval = rep(1:12, each = 4),
                  file = lc.files[1:48])

# Dataframe to store overall results
overall <- data.frame()


###############################################
###                                         ###
###       Calculate Average Cover per       ###
###               Interval                  ###
###############################################

for (i in 1:length(intervals)) {
  
  # set working directory for the files
  setwd(path)
  
  # Create blank dataframe to combine read in files
  tmp.file <- data.frame()
  
  for(j in 1:length(key$interval[key$interval == i])) {
    
    # Read in the four files in the interval
    file.info <- read.csv(key$file[key$interval == i][j])
    
    # Reduce to only the relevant year, year is column 1
    file.info <- file.info %>%
      filter(Year == key$year[key$interval == i][j])
    
    tmp.file <- rbind(tmp.file, file.info)
    
  }

  # Find the average land cover for each class, by route
  int.summary <- tmp.file %>%
    group_by(unique_route) %>%
    summarise(across(where(is.numeric), mean)) %>%
    mutate(Year = paste(key$year[key$interval == i][1], "-",
                        key$year[key$interval == i][4], sep = ""),
           Interval = i) %>%
    select(Year, Interval, unique_route, everything())
  
  # Add results to overall dataframe
  overall <- rbind(overall, int.summary)
  
  print(i)
  
}


# Write the results
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/BBS and Landcover")
write.csv(overall, "Areal Land cover for each BBS rt_1973 to 2020.csv", row.names = FALSE)
