###############################################################################
########                                                               ########
########                 Combining habitat associations                ########
########                          and abundances                       ########
###############################################################################

# We compiled habitat associations for bird species detected in the first time
# interval (1973-1976) or the last time interval (2017-2020) from the IUCN
# red list.

# This code wrangles that data and combines it with abundance estimates
# of the species for the first and last time interval.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(tidyverse) # Data cleaning

path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Load the habitat associations
species.habitats <- read.csv(
  paste(
    path,
    "00_Data/raw/IUCN_Habitat/species_list_with_habitats.csv",
    sep = ""
  ),
  header = TRUE
)

# Load the reclassifcation scheme
reclass.key <-
  read.csv(
    paste(
      path,
      "01_Analysis/IUCN Habitats/Reclassification Scheme.csv",
      sep = ""
    ),
    header = TRUE
  )

# Load the BBS species key (AOU plus species info)
species.key <-
  read.csv(
    paste(path, "03_Code/Code for publishing/SpeciesList.csv", sep = ""),
    header = TRUE
  )


# Load the abundance data for the first and last time interval
int1 <- read.csv(
  paste(path, "01_Analysis/IUCN Habitats/abundance_1973_1976.csv", sep = ""),
  header = TRUE
)


int12 <- read.csv(
  paste(path, "01_Analysis/IUCN Habitats/abundance_2017_2020.csv", sep = ""),
  header = TRUE
)


#####################################
###                               ###
###       Data Processing         ###
###   Abundances & Species Key    ###
#####################################

# First concatenate genus and species in the species.key
# Create species column in species list
species.key$Species1 <- paste(species.key$Genus, "", species.key$Species)

# function that collapses any run of whitespace to a single space and trims
# ends
species.key$Species1 <- str_squish(species.key$Species1)

# Remove the leading X from the AOU code in the abundance dataframe
int1$Species.code <- gsub("X", "", int1$Species.code)
int12$Species.code <- gsub("X", "", int12$Species.code)

# Change the name of the column to match what's in the species key
names(int1)[2] <- "AOU"
names(int12)[2] <- "AOU"

# Merge the abundance dataframes with the AOU, scientific name, and
# common name in the species key
int1_merged <- merge(
  int1,
  species.key[, c("AOU", "English_Common_Name", "Species1")],
  by = "AOU",
  all.x = TRUE
)

int12_merged <- merge(
  int12,
  species.key[, c("AOU", "English_Common_Name", "Species1")],
  by = "AOU",
  all.x = TRUE
)

# Based on previous script there should only be one NA value where the
# AOU did not match any species in the BBS species list: AOU 4890
which(is.na(int1_merged$English_Common_Name))
which(is.na(int12_merged$English_Common_Name))
# Confirmed

# Save abundance dataframes
write.csv(
  int1_merged,
  paste(
    path,
    "01_Analysis/IUCN Habitats/abundance_1973_1976_details.csv",
    sep = ""
  ),
  row.names = FALSE
)

write.csv(
  int12_merged,
  paste(
    path,
    "01_Analysis/IUCN Habitats/abundance_2017_2020_details.csv",
    sep = ""
  ),
  row.names = FALSE
)


#####################################
###                               ###
###       Data Processing         ###
###      Species Habitats         ###
#####################################

# Make sure that every cell in the column IUCN Subcategory has an entry
# as this is what we will use to reclassify the habitat associations
# according to the land cover types used in our study

# Turn truly empty strings into NA values in both columns
species.habitats$Habitats[species.habitats$Habitats == ""] <- NA
species.habitats$Habitats_detail[species.habitats$Habitats_detail == ""] <- NA


# For any blank cells in habitats_detail
# we will fill with either the broad column value

species.habitats$Habitats_detail <- ifelse(
  is.na(species.habitats$Habitats_detail),
  species.habitats$Habitats,
  species.habitats$Habitats_detail
)

# Remove the broad habitat category
species.habitats$Habitats <- NULL

# Rename reclass.key column to match col name of species habitats
names(reclass.key)[2] <- "Habitats_detail"

# Reclass.key$Habitats_detail and species.habitats$Habitats_detail
# differ in their capitalization
# Create lowercased vrsions of the key columns to enable merging
species.habitats$Habitats_detail <-
  tolower(species.habitats$Habitats_detail)

reclass.key$Habitats_detail <-
  tolower(reclass.key$Habitats_detail)

# Merge the species.habitats with reclassification scheme
species.habitats_merged <- merge(
  species.habitats,
  reclass.key[, c("Habitat_association", "Habitats_detail")],
  by = "Habitats_detail",
  all.x = TRUE
)

# There were still some discrepancies between how the detail was described in
# the reclass key and how it was described in the species list that require
# troubleshooting.
which(is.na(species.habitats_merged$Habitat_association))

# Reclassify these categories

# 10.1, 10.2 --> Open Water
ix <- startsWith(species.habitats_merged$Habitats_detail, "10.1") |
  startsWith(species.habitats_merged$Habitats_detail, "10.2")

species.habitats_merged$Habitat_association[ix] <- "Open Water"


# 12.1, 12.2, 12.3, 12.4, 12.5, 12.6 --> Wetland
ix2 <- startsWith(species.habitats_merged$Habitats_detail, "12.1") |
  startsWith(species.habitats_merged$Habitats_detail, "12.2") |
  startsWith(species.habitats_merged$Habitats_detail, "12.3") |
  startsWith(species.habitats_merged$Habitats_detail, "12.4") |
  startsWith(species.habitats_merged$Habitats_detail, "12.5") |
  startsWith(species.habitats_merged$Habitats_detail, "12.6")

species.habitats_merged$Habitat_association[ix2] <- "Wetland"


# 13.1, 13.2, 13.3, 13.5, 13.5 --> Wetland
ix3 <- startsWith(species.habitats_merged$Habitats_detail, "13.1") |
  startsWith(species.habitats_merged$Habitats_detail, "13.2") |
  startsWith(species.habitats_merged$Habitats_detail, "13.3") |
  startsWith(species.habitats_merged$Habitats_detail, "13.4") |
  startsWith(species.habitats_merged$Habitats_detail, "13.5")

species.habitats_merged$Habitat_association[ix3] <- "Wetland"


# 14.1, 14.2, 14.3 --> Agriculture
ix4 <- startsWith(species.habitats_merged$Habitats_detail, "14.1") |
  startsWith(species.habitats_merged$Habitats_detail, "14.2") |
  startsWith(species.habitats_merged$Habitats_detail, "14.3")

species.habitats_merged$Habitat_association[ix4] <- "Agriculture"


# 14.4, 14.5 --> Urban/Developed
ix5 <- startsWith(species.habitats_merged$Habitats_detail, "14.4") |
  startsWith(species.habitats_merged$Habitats_detail, "14.5")

species.habitats_merged$Habitat_association[ix5] <- "Urban/Developed"


# 14.6 --> Forest
ix6 <- startsWith(species.habitats_merged$Habitats_detail, "14.6")

species.habitats_merged$Habitat_association[ix6] <- "Forest"


# 15.1, 15.2, 15.3 --> Open Water
ix7 <- startsWith(species.habitats_merged$Habitats_detail, "15.1") |
  startsWith(species.habitats_merged$Habitats_detail, "15.2") |
  startsWith(species.habitats_merged$Habitats_detail, "15.3")

species.habitats_merged$Habitat_association[ix7] <- "Open Water"


# 15.4, 15.5 --> Barren
ix8 <- startsWith(species.habitats_merged$Habitats_detail, "15.4") |
  startsWith(species.habitats_merged$Habitats_detail, "15.5")

species.habitats_merged$Habitat_association[ix8] <- "Barren"


# 15.6, 15.9 --> Urban Developed
ix9 <- startsWith(species.habitats_merged$Habitats_detail, "15.6") |
  startsWith(species.habitats_merged$Habitats_detail, "15.9")

species.habitats_merged$Habitat_association[ix9] <- "Urban/Developed"


# 15.7, 15.8 --> Agriculture
ix10 <- startsWith(species.habitats_merged$Habitats_detail, "15.7") |
  startsWith(species.habitats_merged$Habitats_detail, "15.8")

species.habitats_merged$Habitat_association[ix10] <- "Agriculture"


# 3.4 --> Shrub/Grassland
ix11 <- startsWith(species.habitats_merged$Habitats_detail, "3.4")

species.habitats_merged$Habitat_association[ix11] <- "Shrub/Grassland"


# 5.5, 5.6, 5.7, 5.8 --> Wetland
ix12 <- startsWith(species.habitats_merged$Habitats_detail, "5.5") |
  startsWith(species.habitats_merged$Habitats_detail, "5.6") |
  startsWith(species.habitats_merged$Habitats_detail, "5.7") |
  startsWith(species.habitats_merged$Habitats_detail, "5.8")

species.habitats_merged$Habitat_association[ix12] <- "Wetland"


# 6. --> Barren
ix13 <- startsWith(species.habitats_merged$Habitats_detail, "6.")

species.habitats_merged$Habitat_association[ix13] <- "Barren"


# 9.10, 9.9 --> Wetland
ix14 <- startsWith(species.habitats_merged$Habitats_detail, "9.10") |
  startsWith(species.habitats_merged$Habitats_detail, "9.9")

species.habitats_merged$Habitat_association[ix14] <- "Wetland"


# Last one that was missed 5.10 --> Wetland
ix15 <- startsWith(species.habitats_merged$Habitats_detail, "5.10")

species.habitats_merged$Habitat_association[ix15] <- "Wetland"


# Now the only remaining NAs are for birds identified only to the genus level
# or birds with 2 species listed that had different habitat associations

### Summarize habitat associations by removing the detailed column and
### keeping only distinct rows

species.habitats_summarized <- species.habitats_merged %>%
  select(-Habitats_detail) %>%
  distinct()


#####################################
###                               ###
###  Putting everything together  ###
###                               ###
#####################################

# Change column name in summarized habitat associations to enable
# merging
names(species.habitats_summarized)[1] <- "Species1"

# Combine the summarized habitat associations and abundances from interval 1
final.summary <- merge(
  species.habitats_summarized,
  int1_merged[, c(2:5)],
  by = "Species1",
  all.x = TRUE
)

# Change name of abundance column to indicate interval1
names(final.summary)[7] <- "Int1_abundance"


# Combine the final summary with abundances from interval 12
final.summary <- merge(
  final.summary,
  int12_merged[, c(2:5)],
  by = "Species1",
  all.x = TRUE
)

# Change name of abundance column to indicate interval12
names(final.summary)[10] <- "Int12_abundance"

# Combine the two english common name columns
final.summary <- final.summary %>%
  mutate(
    English_Common_Name = dplyr::coalesce(
      English_Common_Name.x,
      English_Common_Name.y
    ),
    ID__200 = dplyr::coalesce(ID__200.x, ID__200.y)
  )

# Drop the x and y common name columns
final.summary <- final.summary %>%
  select(-English_Common_Name.x, -English_Common_Name.y, -ID__200.x, -ID__200.y)

# Adjust order of the columns
final.summary <- final.summary[, c(8, 1, 9, 2:7)]


# Save the summary
write.csv(
  final.summary,
  paste(
    path,
    "01_Analysis/IUCN Habitats/Summary_Abundance and Habitat.csv",
    sep = ""
  ),
  row.names = FALSE
)


#####################################
###                               ###
###    Supplemental Table S4      ###
###                               ###
#####################################

# Filter IUCN reclassification scheme to only the habitats that actually
# appear for the bird species in the study

unique(species.habitats$Habitats_detail)
# there are no NA values but there are blank values which means that the broad
# category should be used

which(species.habitats$Habitats_detail == "")

# Subset to those rows just to make sure filtering the reclass scheme by
# the unique values in Habitats_detail won't miss anything
blank.detail <- species.habitats[
  which(species.habitats$Habitats_detail == ""),
]
# these are the 23 species for which no taxonomic name match could be found
# the remaining species are found in rocky areas which is reclassified as barren

reclass.key <- reclass.key %>%
  mutate(Number.key = sub(" .*", "", IUCN.Subcategory))

species.habitats <- species.habitats %>%
  mutate(Number.key_hab = sub(" .*", "", Habitats_detail))
