# This R script calculates environmental heterogeneity
# using the North american bird breeding survey.
# The data is spatially and temporally aggregated
# First into 12 year-range and second into 66 spatially aggregated
# spatial polygons (i.e., cookies) with areas between 10k-20k km2.

library(tidyverse)
library(data.table)
library(vegan)

setwd("G:/My Drive/_Projects/_Biodiversity_drivers/code_data")

# ----------------------- Loading and organizing data --------------------------

# LULC times series data: downloaded and harmonized
site_env_15 <- read.csv(
  "Site x LULC Matrices_by4/Proportional Land cover for each BBS rt_1973 to 2020.csv"
)
colnames(site_env_15)[3] <- "uniq_rt"
colnames(site_env_15)[2] <- "year_range"

length(unique(site_env_15$uniq_rt)) # 786

# site_env <- filter (site_env, Year == "1966") #  | Year == "1968")

# Spatially aggregated routes into cookies: spatial analysis
load("cookies_routes_assembled_20252305.Rdata")
cookies_routes_assembled1 <- cookies_routes_assembled %>%
  select(uniq_rt, ID_cookie_200)

# Iding routes with cookies ID
site_env_15_cookie <- site_env_15 %>%
  left_join(cookies_routes_assembled1, by = "uniq_rt") %>%
  drop_na()

length(unique(site_env_15_cookie$uniq_rt)) # 577
length(unique(site_env_15_cookie$ID_cookie_200)) # 74
length(unique(site_env_15_cookie$year_range)) # 12

# # Converts data to numeric and in percentages coverage
# site_env <- site_env %>%
# mutate(across(where(is.numeric), function(x) x * 100))

names(site_env_15_cookie)

site_env_conus_9 <- data.frame(
  # 9 classes
  uniq_rt = site_env_15_cookie$uniq_rt,
  ID_cookie_200 = site_env_15_cookie$ID_cookie_200,
  year_range = site_env_15_cookie$year_range,

  unclassified = site_env_15_cookie$unclassified,

  water = site_env_15_cookie$open.water,
  ice = site_env_15_cookie$ice.snow,
  urban = site_env_15_cookie$urban.developed,
  mining = site_env_15_cookie$mining,
  barren = site_env_15_cookie$barren,
  forest = site_env_15_cookie$forest.dec +
    site_env_15_cookie$forest.evgrn +
    site_env_15_cookie$forest.mixed,
  shrub.grass = site_env_15_cookie$shrubland + site_env_15_cookie$grassland,
  cropland = site_env_15_cookie$cropland + site_env_15_cookie$hay.pasture,
  wetland = site_env_15_cookie$wetland.woody + site_env_15_cookie$wetland.herb
)


#___________________________

# Beta pairs-routes computation ----

#___________________________

# Computing Beta per pairs per cookie
# Testing with 2 years and 2 cookies
# site_env_conus_9_ <- site_env_conus_9 %>% filter (year_range == "1" | year_range == "2")
# site_env_conus_9_ <- site_env_conus_9_ %>% filter (ID_cookie_200 == "93" | ID_cookie_200 == "116")

# Using Conus with 8 classes

beta_pairs_list <- list()
counter <- 1

for (cookie in unique(site_env_conus_9$ID_cookie_200)) {
  tempo <- filter(site_env_conus_9, ID_cookie_200 == cookie)

  for (year in unique(tempo$year_range)) {
    tempo_year <- filter(tempo, year_range == year)

    # Create a list of sites per unique_route
    site_list <- split(tempo_year, tempo_year$uniq_rt)
    route_ids <- names(site_list)

    # Generate all unique pairs of routes
    route_pairs <- combn(route_ids, 2, simplify = FALSE)

    for (pair in route_pairs) {
      site1 <- site_list[[pair[1]]]
      site2 <- site_list[[pair[2]]]

      pair0 <- rbind(site1, site2)
      pair_data <- pair0 %>% select(-(1:3)) # Assuming columns 1:3 = meta columns

      # Alpha, Gamma, Beta calculations
      alpha <- mean(specnumber(pair_data))
      gamma <- specnumber(colSums(pair_data))

      alpha_simpson <- mean(diversity(pair_data, index = "simpson"))
      gamma_simpson <- diversity(colSums(pair_data), index = "simpson")

      alpha_shannon <- mean(diversity(pair_data, index = "shannon"))
      gamma_shannon <- diversity(colSums(pair_data), index = "shannon")

      # Metadata
      pair_id <- paste(pair[1], "-", pair[2])
      ID_cookie_val <- unique(pair0$ID_cookie_200)
      year_val <- unique(pair0$year_range)

      result_row <- tibble(
        alpha = alpha,
        alpha_simpson = alpha_simpson,
        alpha_shannon = alpha_shannon,
        gamma = gamma,
        gamma_simpson = gamma_simpson,
        gamma_shannon = gamma_shannon,
        pair_id = pair_id,
        ID_cookie = ID_cookie_val,
        year = year_val
      )

      beta_pairs_list[[counter]] <- result_row
      counter <- counter + 1
    }
  }
}

## End data  pairs ----

beta_pairs_env_9 <- bind_rows(beta_pairs_list) # Final dataset

unique(beta_pairs_env_9$year) # Check results that all year and cookies were calculated
length(unique(beta_pairs_env_9$ID_cookie))

#___________________________

# Beta four-routes computation ----

#___________________________

# removing cookies with less than 4 routes
cookies_routes_assembled2 <- cookies_routes_assembled %>%
  as.data.frame() %>%
  select(uniq_rt, ID_cookie_200, Freq)

site_env_conus_9_tmp <- site_env_conus_9 %>%
  left_join(cookies_routes_assembled2, by = "uniq_rt") %>%
  filter(Freq > 4) %>%
  select(-Freq, -ID_cookie_200.y) %>%
  rename(ID_cookie_200 = 2)

length(unique(site_env_conus_9_tmp$ID_cookie_200)) # 66 cookies

beta_fours_list <- list()
counter <- 1
set.seed(123)

for (cookie in unique(site_env_conus_9_tmp$ID_cookie_200)) {
  tempo <- filter(site_env_conus_9_tmp, ID_cookie_200 == cookie)
  print(cookie)

  for (year in unique(tempo$year_range)) {
    tempo_year <- filter(tempo, year_range == year)
    print(year)

    # Create a list of sites per unique_route
    site_list <- split(tempo_year, tempo_year$uniq_rt)
    route_ids <- names(site_list)

    for (rep in 1:99) {
      # Repeat 99 times

      # Sample 4 unique routes
      sampled_routes <- sample(route_ids, 4, replace = FALSE)

      # Combine observations from those routes
      sampled_sites <- bind_rows(site_list[sampled_routes])

      # Remove metadata columns (adjust column exclusions as needed)
      sampled_data <- sampled_sites %>%
        select(-ID_cookie_200, -year_range) %>%
        column_to_rownames(var = "uniq_rt")

      # Compute metrics
      alpha <- mean(specnumber(sampled_data))
      gamma <- specnumber(colSums(sampled_data))

      alpha_simpson <- mean(diversity(sampled_data, index = "simpson"))
      gamma_simpson <- diversity(colSums(sampled_data), index = "simpson")

      alpha_shannon <- mean(diversity(sampled_data, index = "shannon"))
      gamma_shannon <- diversity(colSums(sampled_data), index = "shannon")

      # Metadata
      sampled_ids <- paste(sampled_routes, collapse = "-")
      ID_cookie_val <- unique(sampled_sites$ID_cookie_200)
      year_val <- unique(sampled_sites$year_range)

      result_row <- tibble(
        alpha = alpha,
        alpha_simpson = alpha_simpson,
        alpha_shannon = alpha_shannon,
        gamma = gamma,
        gamma_simpson = gamma_simpson,
        gamma_shannon = gamma_shannon,
        sampled_routes = sampled_ids,
        ID_cookie = ID_cookie_val,
        year = year_val,
        repetition = rep
      )

      beta_fours_list[[counter]] <- result_row
      counter <- counter + 1
    }
  }
}
beta_fours_env_9 <- bind_rows(beta_fours_list) # Final dataset

# Check results that all year and cookies were calculated
unique(beta_fours_env_9$year)
length(unique(beta_fours_env_9$ID_cookie))

## End data 4 routes ----
# both 8 classes takes about 10 minutes to run
save(beta_fours_env_9, file = "results_beta_four_env_9LULC_20250709.Rdata")
