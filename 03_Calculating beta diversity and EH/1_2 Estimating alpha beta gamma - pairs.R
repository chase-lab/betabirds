# This R script calculates alpha, gamma, and beta diversity metrics for 2 routes
# using the North american bird breeding survey.
# The data is spatially and temporally aggregated
# First into 12 year-range and second into 66 spatially aggregated
# spatial polygons (i.e., cookies) with areas between 10k-20k km2.

library(tidyverse)
library(sf)
library(terra)
library(readr)
library(iNEXT) # https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.pdf


setwd("G:/My Drive/_Projects/_Biodiversity_drivers/code_data")
#___________________________

# 0. Organizing data ----

#___________________________

# load data
routes <- read_sf('Site x Species Matrices_20250505_by4/BBS_Rtes_Point.shp') # Edited May 2025
load("cookies_routes_assembled_20252305.Rdata") # Edited may 2025
cookies_routes_shp <- read_sf(dsn = "shp/cookies_routes_assembled_20252305.shp")

# organize data
cookies_routes_assembled <- filter(cookies_routes_assembled, Freq > 4)

length(unique(routes$uniq_rt)) # 786
length(unique(cookies_routes_assembled$uniq_rt)) # 546
length(unique(cookies_routes_assembled$ID_cookie_200)) # 66

# Species by site matrices  -  Load each CSV file into a list of data frames
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

# # subset for running tests
list_SpecSite <- list_SpecSite[c(1, 2)]
target_values <- c("60", "68")
list_SpecSite <- lapply(list_SpecSite, function(df) {
  df %>% filter(ID_cookie_200 %in% target_values)
})


#___________________________

# 1. Computation (Pairs) diversity of order q0,  q1 and q2  ----

#___________________________

# Preallocating lists
beta_pairs_alphagamma <- list()

# Counters for lists
idx_alpha_gamma <- 1


for (i in seq_along(list_SpecSite)) {
  o <- list_SpecSite[[i]]
  print(i)

  for (cookie in unique(o$ID_cookie_200)) {
    print(cookie)
    tempo <- filter(o, ID_cookie_200 == cookie)

    # Create a list of sites per unique_route
    site_list <- split(tempo, tempo$unique_route)
    route_ids <- names(site_list)

    # Generate pairs using combn
    route_pairs <- combn(route_ids, 2, simplify = FALSE)

    for (pair in route_pairs) {
      site1 <- site_list[[pair[1]]]
      site2 <- site_list[[pair[2]]]

      pair0 <- rbind(site1, site2)
      pair_data <- pair0 %>% select(-unique_route, -ID_cookie_200, -Year)

      # Metadata
      pair_id <- paste(pair0$unique_route[1], "-", pair0$unique_route[2])
      ID_cookie_val <- unique(pair0$ID_cookie_200)
      year_val <- unique(pair0$Year)

      # Alpha diversity estimation, sensu Chao et al q0, q1, q2
      tmp_alpha <- estimateD(
        t(pair_data),
        q = c(0, 1, 2),
        datatype = "abundance",
        base = "coverage",
        level = 0.95
      )
      tmp_alpha$metric <- "alpha"
      tmp_alpha$pair_id <- pair_id
      tmp_alpha$cookie <- cookie
      tmp_alpha$year_val <- year_val

      # Gamma diversity estimation, sensu Chao et al q0, q1, q2
      tmp_gamma <- estimateD(
        rowSums(t(pair_data)),
        q = c(0, 1, 2),
        datatype = "abundance",
        base = "coverage",
        level = 0.95
      )
      tmp_gamma$metric <- "gamma"
      tmp_gamma$pair_id <- pair_id
      tmp_gamma$cookie <- cookie
      tmp_gamma$year_val <- year_val

      # beta_pairs_alphagamma - - Collects results in list
      beta_pairs_alphagamma[[idx_alpha_gamma]] <- rbind(tmp_alpha, tmp_gamma)
      idx_alpha_gamma <- idx_alpha_gamma + 1
    }
  }
}

beta_pairs_alphagamma_df <- do.call(rbind, beta_pairs_alphagamma)


#___________________________

# 2. Checks of alpha and gamma  ----

#___________________________

# save(beta_pairs_alphagamma_df,
#      file = "results_beta_pairs_20250526.RData")

# load ("results_beta_pairs_20250526.RData")

unique(beta_pairs_alphagamma_df$year_val) # 12 year range
length(unique(beta_pairs_alphagamma_df$cookie)) # 74 cookies

beta_pairs_alphagamma_df %>% # all cookies were included across all year ranges
  group_by(year_val) %>%
  summarise(num_unique_cookies = n_distinct(cookie))

beta_pairs_alphagamma_df %>% # the amount of pairs vary per cookie because cookies have diff number of routes
  group_by(year_val, cookie) %>%
  summarise(num_unique_pairs = n_distinct(pair_id))

#  Computing Beta diversity using  alpha and gamma diversity (sensu Chao et al) -  https://methodsblog.com/2015/05/27/beta_diversity/
alpha <- beta_pairs_alphagamma_df %>%
  filter(metric == "alpha") %>%
  group_by(year_val, cookie, pair_id, Order.q) %>%
  summarise(
    mean_qD = mean(qD, na.rm = TRUE), # mean_qD is alpha
    mean_qD.LCL = mean(qD.LCL, na.rm = TRUE),
    mean_qD.UCL = mean(qD.UCL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(id = paste(Order.q, pair_id, cookie, year_val, sep = ";"))

gamma <- beta_pairs_alphagamma_df %>%
  filter(metric == "gamma") %>%
  mutate(id = paste(Order.q, pair_id, cookie, year_val, sep = ";")) %>%
  select(
    id,
    qD, # qD is gamma
    qD.LCL,
    qD.UCL
  )

alpha_beta_gamma <- alpha %>%
  left_join(gamma, by = "id", keep = FALSE) %>%
  mutate(beta = qD / mean_qD)

#___________________________

# End data  ----

#___________________________

# save(alpha_beta_gamma,
#      file = "results_alpha_beta_gamma.RData") # last checked Oct 2026
