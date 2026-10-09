################################################################################
########                                                                ########
########    Evaluating classification consistency between NLCD data     ######## 
########             and CONUS backcasted data with graphs              ########
################################################################################

# This code visualizes classification accuracy between NLCD and Conus.
# First for 1000 random points drawn for each cover class across the entire
# contiguous US.
# Second for 1000 random points drawn from each cookie/polygon.

# Based on the reclassification/aggregating classes that were systematically
# misclassified when we ran the first crossvalidation (with all the land
# cover classes)


# Notes on the metrics:

### Overall results ###

# Accuracy: the proportion of correctly classified pixels (both datasets agree)
#           formula = correct pixels / total pixels
#           can be misleading if classes are imbalanced

# Accuracy lower and upper: 95% confidence interval bounds

# Cohen's kappa : Measures agreement between conus and nlcd, adjusting for 
#                 chance agreement
#                formula = (observed accuracy - expected accuracy) / 
#                          (1 - expected accuracy)
#                Interpretation: <0, less agreement than random chance
#                                0-0.2: Slight agreement
#                                0.2 - 0.4: Fair agreement 
#                                0.4-0.6: Moderate agreement
#                                0.6 - 0.8: Substantial agreement
#                                0.8-1: Almost perfect agreement
#                More robust than accuracy for imbalanced classes

# Accuracy Null: the accuracy achievable by always predicting the most 
#                frequent class in the reference data (NLCD)
#                e.g. if the most common class is agriculture at 80% then 
#                accuracy null = 0.8
#                Use as a baseline to compare against your observed accuracy
#                If your accuracy is 0.7 and accuracy null is 0.4, the 
#                harmonizaton adds value

# Accuracy pvalue: pvalue for a one-sided test comparing the observed accuracy
#                  to accuracy Null
#                  p < 0.05: observed accuracy is significantly better than just
#                  guessing the most common class


### By Class results ###

# Sensitivity: Proportion of actual pixels of a class that were correctly 
#              identified as that class (true positive rate)
#              formula = true positives / (true positives + false negative)
#              High specificity means most of the reference pixels were 
#              captured by Conus.

# Specificity: Proportion of pixels not belonging to a class that were
#              correctly identified as not being in that class (true negative
#               rate)
#              formula = true negatives / (true negatives + false positives)
#              High specificity means conus rarely mislabels other classes
#              as this class

# Pos.Pred.Value: Precision, Proportion of pixels conus labeled as a class
#                 that are actually that class in NLCD.
#                 formula = true positives / (true positives + false positives)
#                 High precision means most predicted "forest" pixels really
#                 are "forest".

# Neg.Pred.Value: Proportion of pixels not predicted as a class by conus that 
#                 are not that class in nlcd.
#                 formula = true negatives / (true negatives + false negatives)
#                 High value means few actual pixels were missed by conus.

# Precision: synonym for Pos.Pred.Value

# Recall: synonym for Sensitivity

# F1: harmonic mean of precision and recall
#     formula: 2 * [ (Precision x Recall) / (Precision + Recall) ]
#     Balances false positives and false negatives, high F1 means both
#     precision and recall are high

# Prevalence: Proportion of all reference pixels (NLCD) that belong to this 
#             class
#             shows how common each class is in the reference data

# Detection rate: Proportion of all pixels that are both predicted and actually
#                 are this class
#                 formula = true positives / total pixels

# Detection prevalence: Proportion of all pixels predicted by conus to be this 
#                       class, regardless of correctness
#                       Show how often conus assigns this class label

# Balanced Accuracy: Average of sensitivity and specificity
#                    formula = sensitivity + specificity / 2
#                    Useful for imbalanced classes, high value means both true
#                    positives and true negatives are well captured.

# Example result - I can adjust as necessary
# "CONUS and NLCD show substantial agreement (κ = 0.65), with 78% of pixels 
#  matching. However, systematic biases exist (McNemar’s p = 0.02), likely due 
#  to CONUS underestimating urban areas."


########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory


# Load required packages
library(terra)
library(tidyverse)
library(ggplot2)
library(viridis)
library(forcats)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Entire contiguous USA
# Read in the overall results and the results by cover class
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation")
overall <- read.csv("Overall agreement of Conus to NLCD.csv", header = TRUE)
by.class <- read.csv("By Class agreement of Conus to NLCD.csv", header = TRUE)
misclass.usa <- read.csv("Misclassification by class and year_USA.csv", header = TRUE)

# Cookies only
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Cookie Agreement")
overall.cookie <- read.csv("Overall agreement of Conus to NLCD.csv", header = TRUE)
by.class.cookie <- read.csv("By Class agreement of Conus to NLCD.csv", header = TRUE)
overall.by.polygon <- read.csv("Overall by polygon agreement of Conus to NLCD.csv", header = TRUE)
class.by.polygon <- read.csv("By Class by polygon agreement of Conus to NLCD.csv", header = TRUE)
misclass.cookie <- read.csv("Misclassification by class, cookie, and year.csv", header = TRUE)

# Load the shapefile of the cookies
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/processed/BBS/Sample_1973to2021_by4/Cookies")
cookies <- sf :: st_read("cookies_routes_assembled_20252305.shp")

# Load the USA shapefile
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/00_Data/raw")
usa.shape <- sf :: st_read("USACANAB.shp")

###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################

overall.wide <- pivot_wider(overall, names_from = Metric,
                            values_from = Conus.match.to.NLCD)

# Remove the beginning text from the class column
by.class$Class <- gsub(pattern = "Class: ", "", by.class$Class)


overall.cookie.wide <- pivot_wider(overall.cookie, names_from = Metric,
                            values_from = Conus.match.to.NLCD)

# Remove the beginning text from the class column
by.class.cookie$Class <- gsub(pattern = "Class: ", "", by.class.cookie$Class)

# Define the coordinate system
albers = sp :: CRS(SRS_string = "EPSG:5070")

# Set CRS of USA
usa.shape <- sf :: st_set_crs( usa.shape, albers )

# Crop to USA
usa.shape <- usa.shape[usa.shape$STATE != "AK" & usa.shape$STATE != "NWT" &
                          usa.shape$STATE != "YT" & usa.shape$STATE != "BC" &
                          usa.shape$STATE != "QUE" & usa.shape$STATE != "LAB" &
                          usa.shape$STATE != "SAS" & usa.shape$STATE != "xx" &
                          usa.shape$STATE != "MAN" & usa.shape$STATE != "NFD" &
                          usa.shape$STATE != "ONT" & usa.shape$STATE != "NS" &
                          usa.shape$STATE != "NB" & usa.shape$STATE != "PEI" &
                          usa.shape$STATE != "ALB", ]

# Set CRS of cookies
cookies <- sf :: st_set_crs( cookies, albers )



###############################################
###                                         ###
###          Summary Statistics             ###
###                                         ###
###############################################

summary.overall_usa <- overall %>%
  filter(Metric %in% c("Accuracy", "Kappa")) %>%
  group_by(Metric) %>%
  summarize(mean = mean(Conus.match.to.NLCD),
            min = min(Conus.match.to.NLCD),
            max = max(Conus.match.to.NLCD))

summary.overall_cookies <- overall.by.polygon %>%
  filter(Metric %in% c("Accuracy", "Kappa")) %>%
  group_by(Metric) %>%
  summarize(mean = mean(Value),
            min = min(Value),
            max = max(Value))

# Save dataframes
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation")
write.csv(summary.overall_usa, "Summary stats_USA.csv", row.names = FALSE)
write.csv(summary.overall_cookies, "Summary stats_cookies.csv", row.names = FALSE)


###############################################
###                                         ###
###         Visualizing Results             ###
###             Contiguous USA              ###
###############################################

# First, determine the range of x-axis
x_min <- min(overall$Year)
x_max <- max(overall$Year)
x_center <- (x_min + x_max) / 2

# Overall accuracy and Kappa
overall.plot <- ggplot(data = overall, aes(x = Year, y = Value, color = Metric)) +
  annotate("rect", xmin = x_min, xmax = x_max, 
           ymin = 0.0, ymax = 0.2, alpha = 0.6, fill = turbo(100)[81]) +
  annotate("rect", xmin = x_min, xmax = x_max, 
           ymin = 0.2, ymax = 0.4, alpha = 0.6, fill = turbo(100)[61]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.4, ymax = 0.6, alpha = 0.6, fill = turbo(100)[41]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.6, ymax = 0.8, alpha = 0.6, fill = turbo(100)[21]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.8, ymax = 1.0, alpha = 0.4, fill = turbo(100)[1]) +
  annotate("text", x = x_center,
           y = 0.05, label = "Slight Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.25, label = "Fair Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.45, label = "Moderate Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.65, label = "Substantial Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.85, label = "Almost Perfect Agreement", fontface = "italic", size = 5) +
  # Add dashed horizontal lines
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8),
             linetype = "dashed", color = "black") +
  geom_errorbar(
    data = overall.wide, 
    aes(x = Year, y = Accuracy, ymin = AccuracyLower, ymax = AccuracyUpper,
        color = "Accuracy"), 
    width = 0.2, inherit.aes = FALSE) +
  geom_point(data = overall.wide, aes(x = Year, y = Accuracy,
                                      color = "Accuracy"),
             size = 2) +
  geom_point(data = overall.wide, aes(x = Year, y = Kappa,
                                      color = "Kappa"),
             size = 2) +
  scale_color_manual(
    values = c("Accuracy" = "black", "Kappa" = "red"),
    labels = c("Overall Accuracy", "Cohen's Kappa")) +
  scale_y_continuous(name = "Accuracy", lim = c(0, 1.0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=18)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_text(size=18)) +
  theme(axis.title.y = element_text(size=18)) +
  theme(legend.title = element_blank())

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Graphs")
ggsave("Overall accuracy and Kappa_USA.jpg", overall.plot,
       width = 11, height = 9, dpi = 300,
       device = "jpeg")

# By class results
by.class.long <- pivot_longer(
  by.class,
  cols = c("Sensitivity", "Specificity", "Balanced.Accuracy"),
  names_to = "Metric",
  values_to = "Value"
)

by.class.plot <- ggplot(by.class.long, 
                        aes(x = Class, y = Value,
                            fill = Metric, color = Metric)) +
  geom_boxplot(width = 0.7, alpha = 0.5, position = "identity") +
  scale_color_manual(
    values = c("Sensitivity" = viridis(100)[20],
               "Specificity" = viridis(100)[80],
               "Balanced.Accuracy" = viridis(100)[50]),
    labels = c("Balanced Accuracy",
               "Sensitivity\nTrue Positive Rate",
               "Specificity\nTrue Negative Rate")) +
  scale_fill_manual(
    values = c("Sensitivity" = viridis(100)[20],
               "Specificity" = viridis(100)[80],
               "Balanced.Accuracy" = viridis(100)[50]),
    labels = c("Balanced Accuracy",
               "Sensitivity\nTrue Positive Rate",
               "Specificity\nTrue Negative Rate")) +
  scale_y_continuous(name = "Rate", lim = c(0.2, 1.0),
                     breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size=18, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_text(size=18)) +
  theme(axis.title.y = element_text(size=18)) +
  theme(legend.title = element_blank()) 

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Graphs")
ggsave("Cover class sensitivity and specificity_USA.jpg", by.class.plot,
       width = 13, height = 9, dpi = 300,
       device = "jpeg")



### Systematic misclassifications
setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Graphs/Misclass Plots")

# Fixed order for x-axis
nlcd_order <- c(
  "Agriculture", "Forest",
  "Shrub/Grassland", "Wetland",
  "Barren", "Urban/Developed", "Open Water", "Perennial Ice/Snow"
)

# Fixed order for fill categories/legend
conus_order <- c(
  "Agriculture", "Forest",
  "Shrub/Grassland", "Wetland",
  "Barren", "Urban/Developed", "Open Water", "Perennial Ice/Snow"
)

# One permanent colour assignment, named by class
conus_cols <- setNames(
  grDevices::hcl.colors(
    n = length(conus_order),
    palette = "Dynamic"
  ),
  conus_order
)

for (i in 1:length(unique(misclass.usa$Year))) {
  
  # Filter df to particular Year
  tmp.df <- misclass.usa %>%
    filter(Year == unique(Year)[i])
  
  tmp.df <- tmp.df %>%
    # Relevel nlcd_class to custom order
    mutate(nlcd_class = factor(nlcd_class, levels = nlcd_order),
           conus_class_ordered = factor(conus_class, levels = conus_order)) %>%
    # Within each nlcd_class, order conus_class by descending proportion
    arrange(nlcd_class, desc(proportion))
  
  # Grouped bar plot showing proportion of historical land cover classes (fill)
  # identified as the NLCD class (x axis)
  misclass.plot <- ggplot() +
    geom_bar(data = tmp.df, aes(x = nlcd_class, y = proportion,
                                fill = conus_class_ordered),
             color = "black",
             stat = "identity", position = position_dodge()) +
    scale_x_discrete(name = "NLCD Land Cover", drop = FALSE) +
    scale_y_continuous(name = "Proportion") +
    scale_fill_manual(values = conus_cols,
                      limits = conus_order,
                      drop = FALSE,
                      name = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(size=18, angle = 45, hjust = 1)) +
    theme(axis.text.y = element_text(size=18)) +
    theme(axis.title.x = element_text(size=18)) +
    theme(axis.title.y = element_text(size=18)) +
    theme(legend.title = element_blank()) +
    theme(plot.margin = unit(c(1,1,1,3), "cm")) # top, right, bottom, left
  
  ggsave(paste("USA_Systemtic Misclassifications_",
               unique(misclass.usa$Year)[i],
               ".jpg", sep = ""),
         misclass.plot,
         device = "jpeg",
         width = 11, height = 7, dpi = 300)
  
  print(i)
  
}


###############################################
###                                         ###
###         Visualizing Results             ###
###          Cookies as a group             ###
###############################################

# First, determine the range of your x-axis
x_min <- min(overall.cookie$Year)
x_max <- max(overall.cookie$Year)
x_center <- (x_min + x_max) / 2

# Overall accuracy and Kappa
overall.cookie.plot <- ggplot(data = overall.cookie,
                              aes(x = Year, y = Value, color = Metric)) +
  annotate("rect", xmin = x_min, xmax = x_max, 
           ymin = 0.0, ymax = 0.2, alpha = 0.6, fill = turbo(100)[81]) +
  annotate("rect", xmin = x_min, xmax = x_max, 
           ymin = 0.2, ymax = 0.4, alpha = 0.6, fill = turbo(100)[61]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.4, ymax = 0.6, alpha = 0.6, fill = turbo(100)[41]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.6, ymax = 0.8, alpha = 0.6, fill = turbo(100)[21]) +
  annotate("rect", xmin = x_min, xmax = x_max,
           ymin = 0.8, ymax = 1.0, alpha = 0.4, fill = turbo(100)[1]) +
  annotate("text", x = x_center,
           y = 0.05, label = "Slight Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.25, label = "Fair Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.45, label = "Moderate Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.65, label = "Substantial Agreement", fontface = "italic", size = 5) +
  annotate("text", x = x_center,
           y = 0.85, label = "Almost Perfect Agreement", fontface = "italic", size = 5) +
  # Add dashed horizontal lines
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8),
             linetype = "dashed", color = "black") +
  geom_errorbar(
    data = overall.cookie.wide, 
    aes(x = Year, y = Accuracy, ymin = AccuracyLower, ymax = AccuracyUpper,
        color = "Accuracy"), 
    width = 0.2, inherit.aes = FALSE) +
  geom_point(data = overall.cookie.wide, aes(x = Year, y = Accuracy,
                                      color = "Accuracy"),
             size = 2) +
  geom_point(data = overall.cookie.wide, aes(x = Year, y = Kappa,
                                      color = "Kappa"),
             size = 2) +
  scale_color_manual(
    values = c("Accuracy" = "black", "Kappa" = "red"),
    labels = c("Overall Accuracy", "Cohen's Kappa")) +
  scale_y_continuous(name = "Accuracy", lim = c(0, 1.0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=18)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_text(size=18)) +
  theme(axis.title.y = element_text(size=18)) +
  theme(legend.title = element_blank())

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Graphs")
ggsave("Overall accuracy and Kappa for cookies.jpg", overall.cookie.plot,
       width = 11, height = 9, dpi = 300,
       device = "jpeg")

# By class results
by.class.cookie.long <- pivot_longer(
  by.class.cookie,
  cols = c("Sensitivity", "Specificity", "Balanced.Accuracy"),
  names_to = "Metric",
  values_to = "Value"
)

# Remove unclassified/unknown
by.class.cookie.long <- by.class.cookie.long %>%
  filter(Class != "Unclassified" & Class != "Unknown")

by.class.cookie.plot <- ggplot(by.class.cookie.long, 
                        aes(x = Class, y = Value,
                            fill = Metric, color = Metric)) +
  geom_boxplot(width = 0.7, alpha = 0.5, position = "identity") +
  scale_color_manual(
    values = c("Sensitivity" = viridis(100)[20],
               "Specificity" = viridis(100)[80],
               "Balanced.Accuracy" = viridis(100)[50]),
    labels = c("Balanced Accuracy",
               "Sensitivity\nTrue Positive Rate",
               "Specificity\nTrue Negative Rate")) +
  scale_fill_manual(
    values = c("Sensitivity" = viridis(100)[20],
               "Specificity" = viridis(100)[80],
               "Balanced.Accuracy" = viridis(100)[50]),
    labels = c("Balanced Accuracy",
               "Sensitivity\nTrue Positive Rate",
               "Specificity\nTrue Negative Rate")) +
  scale_y_continuous(name = "Rate", lim = c(0.2, 1.0),
                     breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size=18, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.x = element_text(size=18)) +
  theme(axis.title.y = element_text(size=18)) +
  theme(legend.title = element_blank()) 

setwd("~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/01_Analysis/Land Cover/2nd Crossvalidation/Graphs")
ggsave("Cover class sensitivity and specificity_cookies.jpg", by.class.plot,
       width = 13, height = 9, dpi = 300,
       device = "jpeg")


###############################################
###                                         ###
###         Visualizing Results             ###
###         Cookies Individually            ###
###############################################

# Check the range of the values
names(overall.by.polygon)
range(overall.by.polygon$Value[overall.by.polygon$Metric == "Kappa"])

# Change the name of the polygon ID column to match the cookies shapefile
names(overall.by.polygon)[4] <- "ID__200"

# Remove the unique route column from cookies
cookies <- cookies[, -1]

# Remove duplicates from cookies so it's just cookie IDs
cookies <- unique(cookies)

# merge the shapefile with the accuracy information
polygon.accuracy <- merge(overall.by.polygon, cookies[ , c(1, 4)], by = "ID__200", all.x = TRUE)
polygon.accuracy <- st_as_sf(polygon.accuracy, crs = crs(cookies))


# MAP showing Cohen's Kappa by cookie #
kappa <- polygon.accuracy %>%
  filter(Metric == "Kappa")

# Find average kappa by year
kappa <- kappa %>%
  group_by(ID__200) %>%
  mutate(mean = mean(Value))

# Create the map
kappa.map <- ggplot() +
  geom_sf(data = usa.shape) +
  geom_sf(data = kappa, aes(fill = mean)) +
  theme_bw() +
  coord_sf(datum = st_crs(4735)) +
  theme(axis.text.x = element_text(size=18)) +
  theme(axis.text.y = element_text(size=18)) +
  # theme(axis.ticks.x = element_blank()) +
  # theme(axis.ticks.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  # theme(legend.position = "none") +
  # scale_x_continuous(name="Longitude") +
  # scale_y_continuous(name="Latitude") +
  scale_fill_continuous(name="Cohen's\nKappa",
                        type = "viridis") +
  # scale_color_gradient(name = "Total Years", low = "#440154", high = "#fde725",
  #                      breaks = c(20, 30, 40, 50)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=5)))
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       panel.background = element_blank())

library(plotly)

ggplotly(kappa.map, tooltip = "colour")
