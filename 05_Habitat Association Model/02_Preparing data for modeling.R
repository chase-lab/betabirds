###############################################################################
########                                                               ########
########                  Preparing Data for modeling                  ######## 
########                                                               ########
###############################################################################

# This code wrangles data on abundance and habitat associations so it can be
# used as the underlying data to model the relative abundance of each species
# as a function of 
# time and the number of habitat associations of each species, with 
# random effects for metacommunity and species.


########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() )
gc() #releases memory

library(tidyverse) # Data cleaning
library(brms)

path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Load habitat and abundance summary
final.summary <- read.csv(paste(path, "01_Analysis/IUCN Habitats/Summary_Abundance and Habitat.csv",
                                sep = ""),
                          header = TRUE)

# Load individual abundance dataframes
int1 <- read.csv(paste(path, "01_Analysis/IUCN Habitats/abundance_1973_1976_details.csv",
                       sep = ""),
                 header = TRUE)

int12 <- read.csv(paste(path, "01_Analysis/IUCN Habitats/abundance_2017_2020_details.csv",
                       sep = ""),
                 header = TRUE)


###############################################
###                                         ###
###           Data Preparation              ###
###                                         ###
###############################################

# Summarize the number of habitat associations of each species
hab.freq <- final.summary %>%
  group_by(Species1) %>%
  summarise(
    n_habitat_associations = n_distinct(Habitat_association, na.rm = TRUE)
  )

# Create relative abundance columns for both intervals
int1 <- int1 %>%
  group_by(ID__200) %>%
  mutate(
    relab1 = abundance / sum(abundance, na.rm = TRUE)
  ) %>%
  ungroup()

int12 <- int12 %>%
  group_by(ID__200) %>%
  mutate(
    relab12 = abundance / sum(abundance, na.rm = TRUE)
  ) %>%
  ungroup()

# Change name of abundance column to be unique in each df
names(int1)[3] <- "abundance1"
names(int12)[3] <- "abundance12"

# Combine int1 and int12
int_both <- int1 %>%
  select(ID__200, Species1, abundance1, relab1) %>%
  full_join(
    int12 %>% select(ID__200, Species1, abundance12, relab12),
    by = c("ID__200", "Species1")
  )

# add habitat breadth
hab.freq_full <- int_both %>%
  left_join(hab.freq, by = "Species1")


# We will remove the species that have 0 habitat associations (these would be
# for those listed as hybrids or only identified to the genus level or for 
# which no taxonomic match was found)

hab.freq_full <- hab.freq_full %>%
  filter(n_habitat_associations != 0)

# Now any NA values in either relative abundance column can be filled with 0 
# because it means that that species was not detected in a particular interval
hab.freq_full$relab1[is.na(hab.freq_full$relab1)] <- 0
hab.freq_full$relab12[is.na(hab.freq_full$relab12)] <- 0
hab.freq_full$abundance1[is.na(hab.freq_full$abundance1)] <- 0
hab.freq_full$abundance12[is.na(hab.freq_full$abundance12)] <- 0

# Pivot longer to add a time variable
hab.freq_long <- hab.freq_full %>%
  # pivot only the abundance columns to define Time
  pivot_longer(
    cols = c(abundance1, abundance12),
    names_to  = "Time",
    values_to = "abundance"
  ) %>%
  mutate(
    Time = if_else(Time == "abundance1", "t1", "t12"),
    rel.ab = if_else(Time == "t1", relab1, relab12)
  ) %>%
  select(ID__200, Species1, Time,
         abundance, rel.ab, n_habitat_associations)

# Define factors
hab.freq_long$ID__200 <- as.factor(hab.freq_long$ID__200)
hab.freq_long$Species1 <- as.factor(hab.freq_long$Species1)
hab.freq_long$Time <- as.factor(hab.freq_long$Time)
hab.freq_long$n_habitat_associations <- as.factor(hab.freq_long$n_habitat_associations) 

# Prepare for beta distribution by nudging any values of zero or 1 away from
# that min or max
eps <- 1e-6
hab.freq_long <- hab.freq_long %>%
  mutate(
    RelAbundance_beta = pmin(pmax(rel.ab, eps), 1 - eps)
  )

# There is one species with 7 habitat associations, but we will group it with
# the species with 6 habitat associations to prevent one factor level from 
# containing only 1 species.
# Change tyto alba to 6 associations rather than 7
hab.freq_long$n_habitat_associations[hab.freq_long$Species1 == "Tyto alba"] <- 6

# In brms, variable names can't contain double underscores so change name 
names(hab.freq_long)[1] <- "ID_200"

# Save the dataframe
write.csv(hab.freq_long,
          paste(path, "01_Analysis/IUCN Habitats/Abundance and Habitat data.csv",
                sep = ""),
          row.names = FALSE)


