################################################################################
###                                                                          ###
###               Converting Means and CIs of EH and Beta Div                ###
###                                to % change                               ###
################################################################################

# This code converts the posterior mean and 95% credible intervals for the
# trends in EH and the trends in beta diversity for q=0, q=1, and q=2 to 
# % change per time interval

# Clean workspace
rm(list = ls() )
gc() 

# Load required packages
library(tidyverse)

# Define file path
path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


#####################################
###                               ###
###       Data Loading            ###
###                               ###
#####################################

# Load means and CIs by cookie and overall
q0_all <- read.csv(paste(path, "01_Analysis/Model_Results/q0_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)

q1_all <- read.csv(paste(path, "01_Analysis/Model_Results/q1_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)

q2_all <- read.csv(paste(path, "01_Analysis/Model_Results/q2_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)


q0_cookie <- read.csv(paste(path, "01_Analysis/Model_Results/q0_bycookie_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)

q1_cookie <- read.csv(paste(path, "01_Analysis/Model_Results/q1_bycookie_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)

q2_cookie <- read.csv(paste(path, "01_Analysis/Model_Results/q2_bycookie_Trend Means and CIs.csv",
                         sep = ""),
                   header = TRUE)


#####################################
###                               ###
###       Data Processing         ###
###                               ###
#####################################

# Add identifying column and combine the overall data and the cookie data
q0_all$Type <- "q0"
q1_all$Type <- "q1"
q2_all$Type <- "q2"

q0_cookie$Type <- "q0"
q1_cookie$Type <- "q1"
q2_cookie$Type <- "q2"

all <- rbind(q0_all, q1_all, q2_all)
bycookie <- rbind(q0_cookie, q1_cookie, q2_cookie)

# Clean up workspace
rm(q0_all, q1_all, q2_all, q0_cookie, q1_cookie, q2_cookie)


#####################################
###                               ###
###       Converting values       ###
###         to % change           ###
#####################################

# this is the formula
# Percentage change = (10^(logarithmic change) - 1) x 100

all_transformed <- all %>%
  mutate(
    across(
      -Type,
      ~ (10^.x - 1) * 100
    )
  )

cookie_transformed <- bycookie %>%
  group_by(Variable) %>%
  mutate(
    across(
      -Type,
      ~ (10^.x - 1) * 100
    )
  )


# Save the results
write.csv(all_transformed,
          paste(paste(path, "01_Analysis/Model_Results/Percent Change Mean and CI_all.csv",
                      sep = "")),
          row.names = FALSE)

write.csv(cookie_transformed,
          paste(paste(path, "01_Analysis/Model_Results/Percent Change Mean and CI_bycookie.csv",
                      sep = "")),
          row.names = FALSE)

          