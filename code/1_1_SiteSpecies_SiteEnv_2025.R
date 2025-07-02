# This R script calculates alpha, gamma, and beta diversity metrics for pairs of 
# biodiversity samples from the North american bird breeding survey since 1976. 
# The data is aggregated into 12 year-range and 74 cookies, 
# which is a spatially aggregated unite composed by routes building areas between 10k-20k km2

# Data wrangling 
library(tidyverse)
library(sf)
library(terra)
library(readr)
# Biodiversity metrics
# install.packages("remotes")
# remotes::install_github("T-Engel/betaC")
library(betaC)
library(iNEXT) # https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.pdf
setwd("G:/My Drive/_Projects/_Biodiversity_drivers/code_data")

#
# # # # Loading data ----
#

routes <- read_sf('Site x Species Matrices_20250505_by4/BBS_Rtes_Point.shp') # Edited by kim May 2025

load("cookies_routes_assembled_20252305.Rdata")   # Saved may 2025
cookies_routes_shp <- read_sf(dsn = "shp/cookies_routes_assembled_20252305.shp")

load("function_beta_C_ed.RData")

length(unique(routes$uniq_rt))                         # 786
length(unique(cookies_routes_assembled$uniq_rt))       # 577
length(unique(cookies_routes_assembled$ID_cookie_200)) # 74

# Species by site matrices
# Load each CSV file into a list of data frames
directory_path  <- "Site x Species Matrices_20250505_by4/Site x Species Matrices" # Updated data May 2025
csv_files       <- list.files(directory_path, pattern = "\\.csv$")
data_list       <- purrr::map(csv_files, ~read_csv(file.path(directory_path, .x)))
names(data_list)<- csv_files

cookie_ID <- cookies_routes_assembled %>% select(uniq_rt, ID_cookie_200) %>%
  rename(unique_route = uniq_rt) %>% as.data.frame() %>% select(-geometry)
year_range <- 1:length(data_list)

add_cookie_year <- function(df, year) {
  df <- df %>%
    left_join(cookie_ID, by = "unique_route") %>%
    mutate(Year = year)
  return(df)
}

list_SpecSite <- mapply(add_cookie_year, data_list, year_range, SIMPLIFY = FALSE)
list_SpecSite <- lapply(list_SpecSite, na.omit)
# 
# save(list_SpecSite, file = "list_SpecSite_20250527.RData")

load("list_SpecSite_20250527.RData")

# # subset test
list_SpecSite <- list_SpecSite[c(1, 2)]
target_values <- c("60", "68")
list_SpecSite <- lapply(list_SpecSite, function(df) {
  df %>% filter(ID_cookie_200 %in% target_values)
})

#
# Computing alpha and gamma for pairs of routes ----
#

# Preallocating lists 
beta_pairs_alphagamma <- list()
beta_pairs_bc <- list()

# Counters for lists
idx_alpha_gamma <- 1
idx_bc <- 1

for (i in seq_along(list_SpecSite)) {
  
  o <- list_SpecSite[[i]]
  print(o)
  
  for (cookie in unique(o$ID_cookie_200))  {
    
    tempo <- filter(o, ID_cookie_200 == cookie)
    
    # Creates a list of sites per unique_route
    site_list <- split(tempo, tempo$unique_route)
    route_ids <- names(site_list)
    
    # Generates pairs 
    route_pairs <- combn(route_ids, 2, simplify=FALSE)
    
    for (pair in route_pairs) {
      site1 <- site_list[[pair[1]]]
      site2 <- site_list[[pair[2]]]
      
      pair0 <- rbind(site1, site2)
      pair_data <- pair0 %>% select(-unique_route, -ID_cookie_200, -Year)
      
      # Metadata
      pair_id       <- paste(pair0$unique_route[1], "-", pair0$unique_route[2])
      ID_cookie_val <- unique(pair0$ID_cookie_200)
      year_val      <- unique(pair0$Year)
      
      # Alpha 
      tmp_alpha   <- estimateD(t(pair_data), q = c(0, 1, 2), datatype = "abundance", base="coverage", level=0.95)
      tmp_alpha$metric  <- "alpha"
      tmp_alpha$pair_id <- pair_id
      tmp_alpha$cookie  <- cookie
      tmp_alpha$year_val<- year_val
      
      # Gamma 
      tmp_gamma   <- estimateD(rowSums(t(pair_data)), q = c(0, 1, 2), datatype = "abundance", base="coverage", level=0.95)
      tmp_gamma$metric  <- "gamma"
      tmp_gamma$pair_id <- pair_id
      tmp_gamma$cookie  <- cookie
      tmp_gamma$year_val<- year_val
      
      # Collects results in list
      beta_pairs_alphagamma[[idx_alpha_gamma]] <- rbind(tmp_alpha, tmp_gamma)
      idx_alpha_gamma <- idx_alpha_gamma + 1
      
      # Beta C
      beta_65 <- beta_C(pair_data, 0.65, extrapolation = TRUE)
      beta_80 <- beta_C(pair_data, 0.80, extrapolation = TRUE)
      
      y <- as.data.frame(cbind(beta_65, beta_80))
      y$pair_id <- pair_id
      y$cookie  <- cookie
      y$year_val <- year_val
      
      beta_pairs_bc[[idx_bc]] <- y
      idx_bc <- idx_bc + 1
      
    }
  }
}

beta_pairs_alphagamma_df <- do.call(rbind, beta_pairs_alphagamma)
beta_pairs_bc_df         <- do.call(rbind, beta_pairs_bc)

# save(beta_pairs_alphagamma_df,  beta_pairs_bc_df,   # started 10:48 am -- 3 mins
     file = "results_beta_pairs_20250526.RData") 

#  Estimating Beta diversity -  https://methodsblog.com/2015/05/27/beta_diversity/

# Checking results ----
load("I:/synthesis/Paola_tmp/_Biodiversity_drivers/code_data/results_beta_pairs_20250526.RData") 

unique(beta_pairs_alphagamma_df$year_val)       # 12 
length(unique(beta_pairs_alphagamma_df$cookie)) # 74 cookies

alpha <- beta_pairs_alphagamma_df %>%
  filter(metric == "alpha") %>%
  group_by(year_val, cookie, pair_id, Order.q) %>%
  summarise(
    mean_qD = mean(qD, na.rm=TRUE),  # mean_qD is alpha
    mean_qD.LCL = mean(qD.LCL, na.rm=TRUE),
    mean_qD.UCL = mean(qD.UCL, na.rm=TRUE),
    .groups = "drop" ) %>%
  mutate(id = paste(Order.q, pair_id, cookie, year_val, sep = ";"))

gamma <- beta_pairs_alphagamma_df %>%
  filter(metric == "gamma") %>%
  mutate(id = paste(Order.q, pair_id, cookie, year_val, sep = ";"))  %>%
  select(id, qD,                    # qD is gamma
         qD.LCL, qD.UCL)

alpha_beta_gamma <- alpha %>% left_join(gamma, by = "id", keep = FALSE) %>%
  mutate(beta = qD/mean_qD)





