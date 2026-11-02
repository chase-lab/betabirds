# betabirds
Code to reproduce analyses associated with ‘Bird community change only partially matches landscape habitat change across five decades in North America’

Code authors: Kimberly L. Thompson, Martha Paola Barajas Barbosa, Nicolas Bianco

Barajas Barbosa, M.P.; Thompson, K.L; Biano, N.; Klein, N.; Chase, J.M. 2026. Bird community change only partially matches landscape habitat change across five decades in North America. EcoEvoRxiv.
DOI: https://doi.org/10.32942/X2XD51

## Data
Data publicly available to reproduce results include:
North American Breeding Bird Survey data available from https://doi.org/10.5066/P9GS9K64.

North American Breeding Bird Survey route coordinates available from https://databasin.org/datasets/02fe0ebbb1b04111b0ba1579b89b7420/.

Annual National Land-cover Database data available from https://doi.org/10.5066/P94UXNTS.

Historical land-cover data are available from https://www.sciencebase.gov/catalog/item/59d3c73de4b05fe04cc3d1d1.

Habitat-association data are openly available from The IUCN Red List of Threatened Species Version 2026-1 at www.iucnredlist.org/en.

### Intermediate output data provided includes:

**betabirds/Output Data/cookies_routes_assembled_20252305.shp** Metacommunity spatial data, produced by script folders 01_Generating spatial polygons. 

**betabirds/Output Data/ four_beta_birds_environ_cookie_20250709.csv** Data showing the results of calculating average alpha- and gamma-diversity from 99 random draws of four BBS routes from each metacommunity. 

**betabirds/Output Data/ pairs_beta_birds_environ_cookie_20250709.csv** Data showing the results of calculating average alpha- and gamma-diversity for all pairs of BBS routes in each metacommunity. 

**betabirds/Output Data/Site x Species Matrix_1973_1976.csv** Data showing the combined abundance of each species in each BBS route for the period 1973-1976. The first column is unique_route, the ID of the BBS route, with the remaining column names corresponding to BBS AOU species codes combined with a leading “X”. For example, abundances of AOU 5220 are found in column ‘X5220’. Used with betabirds/Output Data/ cookies_routes_assembled_20252305.shp to allocate the abundances for each route to the appropriate metacommunity. Used with betabirds/Output Data/SpeciesList.csv to identify the species name from the AOU code.

**betabirds/Output Data/Site x Species Matrix_2017_2020.csv** Data showing the combined abundance of each species in each BBS route for the period 2017_2020. The first column is unique_route, the ID of the BBS route, with the remaining column names corresponding to BBS AOU species codes combined with a leading “X”. For example, abundances of AOU 5220 are found in column ‘X5220’. Used with betabirds/Output Data/ cookies_routes_assembled_20252305.shp to allocate the abundances for each route to the appropriate metacommunity. Used with betabirds/Output Data/SpeciesList.csv to identify the species name from the AOU code.

**betabirds/Output Data/SpeciesList.csv** Taxonomic data. Columns include: 
-	*Seq*: Arbitrary ID
-	*AOU*: BBS species code
-	*English_Common_Name*
-	*French_Common_Name*
-	*Order*
-	*Family*
-	*Genus*
-	*Species*

## R Scripts
R Script folders, which are listed below, are numbered and listed in the order they should be used. Within each folder, scripts are also numbered in the order they should be used. Several scripts within the 02_Landcover Harmonization and Extraction folder are run on a high-performance computing cluster to reduce computation time, which means that this workflow can not be continuously run on a local machine. However, the code from the scripts that use the HPC cluster can be used pedagogically or adapted to run on a local machine. Note however, that they are computationally intensive.

**00_BBS_Preparation** The scripts in this folder cleans the North American Breeding Bird Survey raw and spatial data, determines the temporal interval of the analysis, and creates route-level site  by species abundance matrices. Downloaded data and the data resulting from these processing scripts are not included, with the exception of the site by species matrix for the first and last time interval (1973-1976 and 2017-2020).

**01_Generating spatial polygons** The script in this folder groups BBS routes into metacommunities. The data resulting from this processing script is provided as ‘cookies_routes_assembled_20252305.shp’.

**02_Landcover Harmonization and Extraction** The scripts in this folder harmonize Annual National Land Cover Database maps and historical land-cover reconstruction maps, assess the accuracy and misclassification of the harmonization, aggregate some land-cover classes to improve accuracy, reassess the accuracy and misclassification based on the aggregation, and finally, extract the proportion and area of each land-cover category associated with each BBS route. These land-cover proportions and areas are then combined according by metacommunity. Original data associated with this code not provided, but are publicly available via the links provided above.

**03_Calculating beta diversity and EH** The scripts in this folder calculate average alpha-, gamma- and beta-diversity for Hill numbers of order 0, 1, and 2 averaged across all pairs of routes within a metacommunity and across 99 random draws of 4 routes per metacommunity. Original data associated with this code not provided.

**04_Hierarchical Modelling** 

**05_Habitat Association Model** 

**06_Figures** The scripts in this folder create main manuscript figures (Figures 2-5) and supplemental figures (S1, S10-21). The data output is not provided.

