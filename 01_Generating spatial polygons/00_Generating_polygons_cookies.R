library(sf)
library(terra)
library(tidyverse)
library(stats)
library(units)

setwd("G:/My Drive/_Projects/_Biodiversity_drivers/code_data")

#_______________________________________________________________________________

#                    BUILDING SPATIAL POLYGONS (COOKIES)                    ----

#_______________________________________________________________________________


# We create spatial polygons (i.e., cookies) containing minimum 4 routes and similar area, i.e., 
# between 10k - 20k km2


# Load shapefile with routes aggregated them into cookies
routes <- terra::vect('Site x Species Matrices_20250505_by4/BBS_Rtes_Point.shp')
crs(routes)

# Create an extent object for the contiguous United States
ext  <- c(-2353753.33649286, 2159145.22837876, 349706.638549438, 3150135.00833635)
res2 <- c(200000, 200000)  # 200 km or 40.000 km2

rast_template <- terra::rast(extent = ext,  resolution = res2)
rast2          <- rast_template
rast2[]        <- 1:ncell(rast2) # it creates a categorical raster with unique IDs
crs(rast2)     <- crs(routes)
crs(rast2)

plot(rast_template)
plot(rast2)
plot(routes, add = TRUE, col = "red")

#___________________________

# 1. Building the spatial polygons by extracting the cell ID from the raster ----

#___________________________

cookies_tmp <- rast2 %>%
  terra::extract(routes, touches = FALSE) %>%
  mutate(uniq_rt       = routes$uniq_rt,     # adds uniq_rt and factor columns
         ID_cookie_200 = as.factor(lyr.1),
         uniq_rt       = as.factor(uniq_rt)) %>%
  select(-lyr.1, -ID)

cookies_tmp <- merge(routes, cookies_tmp, all.x=TRUE, by.x = c('uniq_rt') ) # IDies routes with cookie
# writeVector(cookies_tmp, filename = "Site x Species Matrices_20250505_by4/cookies_routes_tmp1_20250506.shp") 


# Calculate area and number of routes per cookie
calculate_area_freq <- function(data, cookie, area) {
  cookie_sym <- sym(cookie) # converts string into a symbol to be used in the group_by function
  area_sym <- sym(area)
  
  data %>%
    group_by(!!cookie_sym) %>%
    mutate(Freq = n()) %>%             # Counts how many routes per cookie
    summarize(Freq = first(Freq)) %>%  # Keeps Freq in summarized result
    mutate(geometry = st_convex_hull(geometry)) %>%
    mutate(!!area_sym := st_area(geometry))  
}

cookies_tmp_0 <- cookies_tmp %>% st_as_sf() %>% select(uniq_rt, ID_cookie_200)

cookies_tmp1 <- calculate_area_freq(cookies_tmp_0, "ID_cookie_200", "area_m2") %>%
  mutate(area_km2 = as.numeric(area_m2) / 1e6) %>%  # Converts to km
  select(-area_m2) %>% 
  filter(Freq > 2) # removes spatial polygons with less than 3 routes 

length(unique(cookies_tmp1$ID_cookie_200)) # 103 cookies
hist(cookies_tmp1$Freq) 
hist(cookies_tmp1$area_km2) 
plot(cookies_tmp1["ID_cookie_200"])

checkcookiesnumber <- cookies_tmp1 %>% filter(area_km2 > 10000) 
length(unique(checkcookiesnumber$ID_cookie_200)) # cookies_tmp1 make 49 spatial polyons

# st_write(cookies_tmp1,'Site x Species Matrices_20250505_by4/cookies_tmp1_20250506.shp', delete_layer = TRUE)


#___________________________

# 2. Building the additional spatial polygons using cluster analysis ----
# Rescuing routes belonging to polygons smaller than 10k km2 AND 

#___________________________


# Work on routes belonging to spatial polygons smaller than 10k and single/pair routes
tmp <- st_join(cookies_tmp_0, # all 786 routes with a cookie ID
               cookies_tmp1,  # cookies with min 3 routes
               by = "ID_cookie_200") 

# Filters routes with NAs, i.e., single (an isolated point) or pair (line) routes that could not build a spatial polygons
# and cookies < 10.000 km2
tmp_1 <- tmp %>% filter(is.na(area_km2) | area_km2 < 10000)

coords   <- st_coordinates(tmp_1)
dist_mat <- dist(coords) # Euclidean distance matrix
hc       <- hclust(dist_mat, method = "complete") # Hierarchical clustering 

clusters <- cutree(hc, k = 60) # Cut tree into clusters
tmp_1$clusters <- clusters # Add cluster ID to the data

cookies_tmp2 <- tmp_1 %>% 
  calculate_area_freq("clusters", "area_m2") %>% 
  mutate(area_km2 = as.numeric(area_m2/1e6)) %>% # Converts to km
  filter(area_km2 > 10000 & area_km2 < 20000)  %>% select(-area_m2)

plot(cookies_tmp2["area_km2"])
# st_write(cookies_tmp2, dsn = "Site x Species Matrices_20250505_by4/cookies_tmp2_20250607.shp", delete_layer = TRUE)

# Matching spatial polygons (i.e., cookies) to routes ID
cookies_routes_tmp1 <- cookies_tmp1 %>%
  st_join(cookies_tmp_0, by = "ID_cookie_200") %>%
  filter(!is.na(area_km2) & area_km2 > 10000) %>%
  select(-ID_cookie_200.y)

length(unique(cookies_routes_tmp1$ID_cookie_200.x)) # 49 cookies? correct! 54 cookies where removed bc they were too small 

cookies_routes_tmp2 <- cookies_tmp2 %>%
  st_join(tmp_1, by = "ID_cookie_200") %>%
  select(-ID_cookie_200.x, -ID_cookie_200.y, -Freq.y, -area_km2.y, -clusters.y)

length(unique(cookies_routes_tmp2$clusters.x)) # 20 cookies? correct!
length(intersect(cookies_routes_tmp1$uniq_rt, cookies_routes_tmp2$uniq_rt)) # should be zero!

length(union(cookies_routes_tmp1$uniq_rt, cookies_routes_tmp2$uniq_rt)) # 561 routes already belong to a cookie

# HCLUST AGAIN! Filter out routes that are already in a cookie, run clustering again.  
routes_alreadyinacookie <- union(cookies_routes_tmp1$uniq_rt, cookies_routes_tmp2$uniq_rt)

routes_remaining <- cookies_tmp_0 %>%  
  filter(!uniq_rt %in% routes_alreadyinacookie) # what can we do with the 225 left our routes? Run again hclust

coords   <- st_coordinates(routes_remaining)
dist_mat <- dist(coords) # Euclidean distance matrix
hc       <- hclust(dist_mat, method = "complete") # Hierarchical clustering 

clusters <- cutree(hc, k = 50) # Cut tree into clusters
routes_remaining$clusters <- clusters # Add cluster ID to the data

cookies_tmp3 <- routes_remaining %>% 
  calculate_area_freq("clusters", "area_m2") %>% 
  mutate(area_km2 = as.numeric(area_m2/1e6)) %>% 
  filter(area_km2 > 10000 & area_km2 < 21000)%>% 
  select(-area_m2)

plot(cookies_tmp3["area_km2"]) # we made 5 cookies extra

cookies_routes_tmp3 <- cookies_tmp3 %>%
  st_join(routes_remaining, by = "clusters")%>%
  select(-ID_cookie_200, -clusters.y) 

# now we have 594 routes already belonging to a cookie, this is the sum of:
length(cookies_routes_tmp3$uniq_rt)+ 
  length(cookies_routes_tmp2$uniq_rt)+
  length(cookies_routes_tmp1$uniq_rt)

# But we have 192 routes reminding after building clusters twice that do not belong to any cookie
length(routes_remaining$uniq_rt) - length(cookies_routes_tmp3$uniq_rt)

#___________________________

# 3. Trimming spatial polygons (i.e., cookies) lager than 20k km2 ----

#___________________________


# Routes that make > 20.000 km2 cookies
tmp <- cookies_routes_tmp1 %>% filter(area_km2 > 20000) # filter routes of large cookies

length(unique(tmp$ID_cookie_200.x))# 9 cookies are > 20.000 km2
length(unique(tmp$uniq_rt))        # 135 routes 
plot(tmp["area_km2"])

A = as.data.frame(tmp) # change to data frame (this is routes of large cookies)
A$geometry= NULL 

# Getting routes coordinates belonging to big cookies
tmp1 <- left_join(cookies_tmp_0, # This is all 786 routes with cookie ID
                 A,              # routes belonging to >20k cookies
                 by = "uniq_rt") 

tmp1 <- tmp1 %>% filter(!is.na(ID_cookie_200.x)) %>%
  select(uniq_rt, ID_cookie_200, area_km2)

calculate_area <- function(data) {
  data %>%
    group_by(ID_cookie_200) %>%
    mutate(Freq = n()) %>%  
    summarize(Freq = first(Freq),
              geometry = st_combine(geometry)) %>%
    mutate(geometry = st_convex_hull(geometry)) %>%
    st_as_sf() %>%
    mutate(area_km2 = units::set_units(st_area(geometry), "km^2"))
}

area_before <- calculate_area(tmp1) # Calculates area

# Drop most distant point from each cookie
tmp_modified_1 <- tmp1 %>%
  group_by(ID_cookie_200) %>%
  group_modify(~ {
    # Combine and get centroid
    combined <- st_combine(.x$geometry)
    centroid <- st_centroid(st_convex_hull(combined))
    # Calculate distances
    distances <- st_distance(.x$geometry, centroid)
    # Drop most distant point
    .x[-which.max(distances), ]
  }) %>%
  ungroup()

area_after_1 <- calculate_area(tmp_modified_1) # Recalculate area
plot(area_after_1["ID_cookie_200"])



# This is the data needed for analyses 

setdiff(tmp1$uniq_rt, tmp_modified_1$uniq_rt) # 9 routes removed bc make cookies too large

cookies_routes_tmp4 <- left_join (tmp_modified_1, area_after_1, by = "ID_cookie_200") %>%
  filter (area_km2.y < set_units(20000, "km^2")) %>%
  select (ID_cookie_200, Freq, geometry.y,  area_km2.y, uniq_rt)

unique(cookies_routes_tmp4$ID_cookie_200) # 4 cookies trimmed. Great!


# ................


# Filters out cookies that remain > 20.000 km2
large_cookies_1 <- area_after_1 %>% 
  filter(area_km2 > set_units(20000, "km^2")) %>%
  select(ID_cookie_200)

tmp2 <- tmp1 %>% # Getting routes coordinates belonging to remaining > 20.000 km2 cookies
  st_join(large_cookies_1, by = "ID_cookie_200") %>%
  drop_na() %>%
  mutate(ID_cookie_200 = ID_cookie_200.x) %>% select(-ID_cookie_200.x, -ID_cookie_200.y)

area_before <- calculate_area(tmp2)

tmp_modified_2 <- tmp2 %>%
  group_by(ID_cookie_200) %>%
  group_modify(~ {
    combined <- st_combine(.x$geometry)
    centroid <- st_centroid(st_convex_hull(combined))
    distances <- st_distance(.x$geometry, centroid)
    .x[-which.max(distances), ]
  }) %>%
  ungroup()

area_after_2 <- calculate_area(tmp_modified_2) 
plot(area_after_2["ID_cookie_200"])

cookies_routes_tmp5 <- left_join (tmp_modified_2, area_after_2, by = "ID_cookie_200") %>%
  filter (area_km2.y < set_units(20000, "km^2")) %>%
  select (ID_cookie_200, Freq, geometry.y,  area_km2.y, uniq_rt)


# ................


large_cookies_2 <- area_after_2 %>%
  filter(area_km2 > set_units(20000, "km^2")) %>%
  select(ID_cookie_200)

tmp3 <- tmp1 %>% # Getting routes coordinates belonging to remaining > 20.000 km2 cookies
  st_join(large_cookies_2, by = "ID_cookie_200") %>%
  drop_na() %>%
  mutate(ID_cookie_200 = ID_cookie_200.x) %>% select(-ID_cookie_200.x, -ID_cookie_200.y)

area_before <- calculate_area(tmp3)

tmp_modified_3 <- tmp3 %>%
  group_by(ID_cookie_200) %>%
  group_modify(~ {
    combined <- st_combine(.x$geometry)
    centroid <- st_centroid(st_convex_hull(combined))
    distances <- st_distance(.x$geometry, centroid)
    .x[-which.max(distances), ]
  }) %>%
  ungroup()

area_after_3 <- calculate_area(tmp_modified_3) 
plot(area_after_3["ID_cookie_200"]) 
area_after_3$area_km2 # do not save it a cookies_routes_tmp6 because it is only one cookie. it needs to be trimmed again


# ................


large_cookies_3 <- area_after_3 %>%
  filter(area_km2 > set_units(20000, "km^2")) %>%
  select(ID_cookie_200)

tmp4 <- tmp1 %>% # Getting routes coordinates belonging to remaining > 20.000 km2 cookies
  st_join(large_cookies_3, by = "ID_cookie_200") %>%
  drop_na() %>%
  mutate(ID_cookie_200 = ID_cookie_200.x) %>% select(-ID_cookie_200.x, -ID_cookie_200.y)

area_before <- calculate_area(tmp4)

tmp_modified_4 <- tmp4 %>%
  group_by(ID_cookie_200) %>%
  group_modify(~ {
    combined <- st_combine(.x$geometry)
    centroid <- st_centroid(st_convex_hull(combined))
    distances <- st_distance(.x$geometry, centroid)
    .x[-which.max(distances), ]
  }) %>%
  ungroup()

area_after_4 <- calculate_area(tmp_modified_4) 
plot(area_after_4["ID_cookie_200"])
area_after_4$area_km2 # do not save it a cookies_routes_tmp6 because it is only one cookie. it needs to be trimmed again


# ................


large_cookies_4 <- area_after_4 %>%
  filter(area_km2 > set_units(20000, "km^2")) %>%
  select(ID_cookie_200)

tmp5 <- tmp1 %>%
  st_join(large_cookies_4, by = "ID_cookie_200") %>%
  drop_na() %>%
  mutate(ID_cookie_200 = ID_cookie_200.x) %>% select(-ID_cookie_200.x, -ID_cookie_200.y)

area_before <- calculate_area(tmp5)

tmp_modified_5 <- tmp5 %>%
  group_by(ID_cookie_200) %>%
  group_modify(~ {
    combined <- st_combine(.x$geometry)
    centroid <- st_centroid(st_convex_hull(combined))
    distances <- st_distance(.x$geometry, centroid)
    .x[-which.max(distances), ]
  }) %>%
  ungroup()

area_after_5 <- calculate_area(tmp_modified_5) 
plot(area_after_5["ID_cookie_200"])
area_after_5$area_km2 

cookies_routes_tmp6  <- left_join (tmp_modified_5, area_after_5, by = "ID_cookie_200") %>%
  filter (area_km2.y <= set_units(21000, "km^2")) %>%
  select (ID_cookie_200, Freq, geometry.y,  area_km2.y, uniq_rt)

cookies_routes_tmp4_trimming <- rbind(cookies_routes_tmp4, cookies_routes_tmp5, cookies_routes_tmp6)



# Checking that the large trimmed cookies are fine after trimming

sort(unique(tmp$ID_cookie_200.x))
sort(unique(cookies_routes_tmp4_trimming$ID_cookie_200)) # Alles OK! I have the same cookies. just trimmed :)


#___________________________

# 4. Merging edited cookies ----

#___________________________


# Initial ID using the grid 200 km2
unique(cookies_tmp$ID_cookie_200) ; length(unique(cookies_tmp$ID_cookie_200)) # 186
 
# This cookies ID are not modified from the grid IDs
unique(cookies_routes_tmp1$ID_cookie_200.x)
length(unique(cookies_routes_tmp1$ID_cookie_200.x)) # 49

cookies_routes_tmp1_vf <- cookies_routes_tmp1 %>%
  transmute(
    uniq_rt = as.factor(uniq_rt),
    ID_cookie_200 = as.factor(ID_cookie_200.x),
    Freq = as.integer(Freq),
    area_km2 = as.numeric(area_km2),
    geometry = geometry ) %>% 
  filter (area_km2 < 20000)  # I had not removed old big cookies here

length(unique(cookies_routes_tmp1_vf$ID_cookie_200)) # 40 only!

unique(cookies_routes_tmp4_trimming$ID_cookie_200)
length(unique(cookies_routes_tmp4_trimming$ID_cookie_200)) # 9

cookies_routes_tmp4_vf <- cookies_routes_tmp4_trimming %>%
  transmute(
    uniq_rt = as.factor(uniq_rt),
    ID_cookie_200 = as.factor(ID_cookie_200),
    Freq = as.integer(Freq),
    area_km2 = as.numeric(area_km2.y),
    geometry = geometry.y)

# We have to ID the clusters
unique(cookies_routes_tmp2$clusters.x)
unique(cookies_routes_tmp3$clusters.x)

cookies_routes_tmp2$ID_cookie_200 <- as.factor(cookies_routes_tmp2$clusters.x + 5000) # 20 cookies
cookies_routes_tmp3$ID_cookie_200 <- as.factor(cookies_routes_tmp3$clusters.x + 6000) # 5 cookies

cookies_routes_tmp2_vf <- cookies_routes_tmp2 %>%
  transmute(
    uniq_rt = as.factor(uniq_rt),
    ID_cookie_200 = as.factor(ID_cookie_200),
    Freq = as.integer(Freq.x),
    area_km2 = as.numeric(area_km2.x),
    geometry = geometry)

cookies_routes_tmp3_vf <- cookies_routes_tmp3 %>%
  transmute(
    uniq_rt = as.factor(uniq_rt),
    ID_cookie_200 = as.factor(ID_cookie_200),
    Freq = as.integer(Freq),
    area_km2 = as.numeric(area_km2),
    geometry = geometry)

cookies_routes_tmp4_vf <- st_as_sf(cookies_routes_tmp4_vf)

cookies_routes_assembled <- rbind(cookies_routes_tmp1_vf, cookies_routes_tmp2_vf, cookies_routes_tmp3_vf, cookies_routes_tmp4_vf)

intersect(cookies_routes_tmp1_vf$ID_cookie_200, cookies_routes_tmp4_vf$ID_cookie_200)

length(unique(cookies_routes_assembled$ID_cookie_200))  ## 74 cookies in TOTAL

plot(cookies_routes_assembled["ID_cookie_200"])


#___________________________

# 5.Saving assembled cookies ---- 

#___________________________


# save(cookies_routes_assembled, file = "cookies_routes_assembled_20252305.Rdata")  
# st_write(cookies_routes_assembled,'shp/cookies_routes_assembled_20252305.shp', delete_layer = TRUE)


#___________________________

# 6. Assessment of assembled cookies ---- 

#___________________________


load("cookies_routes_assembled_20252305.Rdata") 

df_summary <- cookies_routes_assembled %>%
  group_by(ID_cookie_200) %>%
  summarise(
    Freq = unique(Freq, na.rm = TRUE),         # or use mean(), median(), etc.
    area_km2 = unique(area_km2, na.rm = TRUE)
  )

a = ggplot(df_summary, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "gray50", color = "gray53")  +
  labs(title = "Histogram of Aggregated Freq", x = "Number of routes", y = "Count") +
  theme_minimal()
a

min_val <- min(df_summary$area_km2, na.rm = TRUE)
max_val <- max(df_summary$area_km2, na.rm = TRUE)
b = ggplot(df_summary, aes(x = area_km2)) +
  geom_histogram(fill = "gray50", color = "gray53") +
  scale_x_continuous(
    breaks = pretty(c(min_val, max_val), n = 5) , # Add more breaks for better spacing
    limits = c(min_val, max_val)
  ) +
  labs(title = "Histogram cookies area", x = "Cookies area", y = "Count") +
  theme_minimal()
b
ggpubr::ggarrange(a, b) 

pdf(file = "results_2025/Figures/Histogram_cookies.pdf",
    width =  5, height=  2,  paper = "a4")
ggpubr::ggarrange(a, b) 
dev.off()



