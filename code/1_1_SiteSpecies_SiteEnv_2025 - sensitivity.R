# This R script calculates alpha, gamma, and beta diversity metrics for 4 routes
# using the North american bird breeding survey since 1976.
# The data is aggregated into 12 year-range and 74 cookies,
# which is also spatially aggregated into unites (aka cookies)
# composed by routes building areas between 10k-20k km2
# Written by M. Paola Barajas Barbosa

# # # # Loading packages ----
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

# # # # Loading data ----

routes <- read_sf('Site x Species Matrices_20250505_by4/BBS_Rtes_Point.shp') # Edited by kim May 2025

load("cookies_routes_assembled_20252305.Rdata") # Saved may 2025
cookies_routes_shp <- read_sf(dsn = "shp/cookies_routes_assembled_20252305.shp")
cookies_routes_assembled <- filter(cookies_routes_assembled, Freq > 4)

length(unique(routes$uniq_rt)) # 786
length(unique(cookies_routes_assembled$uniq_rt)) # 546
length(unique(cookies_routes_assembled$ID_cookie_200)) # 66

# Species by site matrices
# Load each CSV file into a list of data frames
directory_path <- "Site x Species Matrices_20250505_by4/Site x Species Matrices" # Updated data May 2025
csv_files <- list.files(directory_path, pattern = "\\.csv$")
data_list <- purrr::map(csv_files, ~ read_csv(file.path(directory_path, .x)))
names(data_list) <- csv_files

cookie_ID <- cookies_routes_assembled %>%
  select(uniq_rt, ID_cookie_200) %>%
  rename(unique_route = uniq_rt) %>%
  as.data.frame() %>%
  select(-geometry)
year_range <- 1:length(data_list)

add_cookie_year <- function(df, year) {
  df <- df %>%
    left_join(cookie_ID, by = "unique_route") %>%
    mutate(Year = year)
  return(df)
}

list_SpecSite <- mapply(
  add_cookie_year,
  data_list,
  year_range,
  SIMPLIFY = FALSE
)
list_SpecSite <- lapply(list_SpecSite, na.omit)

save(list_SpecSite, file = "list_SpecSite_sensitiv_20250527.RData")

load("list_SpecSite_sensitiv_20250527.RData")

# # subset test
list_SpecSite <- list_SpecSite[c(1, 2)]
target_values <- c("60", "68")
list_SpecSite <- lapply(list_SpecSite, function(df) {
  df %>% filter(ID_cookie_200 %in% target_values)
})

# ALPHA AND GAMMA COMPUTATION PAIRS ----

# Preallocate lists
beta_samples_alphagamma <- list()
beta_samples_bc <- list()

# Counters for lists
idx_alpha_gamma <- 1
idx_bc <- 1

set.seed(123) # For reproducibility

for (i in seq_along(list_SpecSite)) {
  o <- list_SpecSite[[i]]
  print(unique(o$Year))

  for (cookie in unique(o$ID_cookie_200)) {
    tempo <- filter(o, ID_cookie_200 == cookie)
    route_ids <- unique(tempo$unique_route)

    for (rep in 1:10) {
      # Repeat 10 times
      sampled_routes <- sample(route_ids, 4, replace = FALSE)
      sampled_data <- filter(tempo, unique_route %in% sampled_routes)
      sample_data <- sampled_data %>%
        select(-ID_cookie_200, -Year) %>%
        column_to_rownames(var = "unique_route")

      # Metadata
      sample_id <- paste(sampled_routes, collapse = "-")
      ID_cookie_val <- unique(sampled_data$ID_cookie_200)
      year_val <- unique(sampled_data$Year)

      # Alpha
      tmp_alpha <- estimateD(
        t(sample_data),
        q = c(0, 1, 2),
        datatype = "abundance",
        base = "coverage",
        level = 0.95
      )
      tmp_alpha$metric <- "alpha"
      tmp_alpha$sample_id <- sample_id
      tmp_alpha$cookie <- cookie
      tmp_alpha$year_val <- year_val
      tmp_alpha$replicate <- rep # Track the random sample

      # Gamma
      tmp_gamma <- estimateD(
        rowSums(t(sample_data)),
        q = c(0, 1, 2),
        datatype = "abundance",
        base = "coverage",
        level = 0.95
      )
      tmp_gamma$metric <- "gamma"
      tmp_gamma$sample_id <- sample_id
      tmp_gamma$cookie <- cookie
      tmp_gamma$year_val <- year_val
      tmp_gamma$replicate <- rep # Track the random sample

      # Collect results in list
      beta_samples_alphagamma[[idx_alpha_gamma]] <- rbind(tmp_alpha, tmp_gamma)
      idx_alpha_gamma <- idx_alpha_gamma + 1

      # Beta
      beta_65 <- beta_C(sample_data, 0.65, extrapolation = TRUE)
      beta_80 <- beta_C(sample_data, 0.80, extrapolation = TRUE)

      y <- as.data.frame(cbind(beta_65, beta_80))
      y$sample_id <- sample_id
      y$cookie <- cookie
      y$year_val <- year_val
      y$replicate <- rep # Track the random sample

      beta_samples_bc[[idx_bc]] <- y
      idx_bc <- idx_bc + 1
    }
  }
}

# Bind results into dataframes
beta_samples_alphagamma_df <- do.call(rbind, beta_samples_alphagamma)
beta_samples_bc_df <- do.call(rbind, beta_samples_bc)

# Save results
# save(beta_samples_alphagamma_df, beta_samples_bc_df, file = "results_beta_samples_sensitiv_20250526.RData")
