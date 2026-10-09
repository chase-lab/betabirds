################################################################################
###                                                                          ###
###                   Figure: Beta diversity q0 v. EH                        ###
###                                                                          ###
################################################################################

# This code plots the results of the model associating the rate of change in
# beta diversity (q=0) to the rate of change in environmental
# heterogeneity.

# Supplemental Figure 


# Clean workspace
rm(list = ls() )
gc() 

# Load required packages
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(viridis)

# Define file path
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


#####################################
###                               ###
###       Data Loading            ###
###                               ###
#####################################

# Load the raw data
data = read.csv(paste(path, 'four_beta_birds_environ_cookie_20250709.csv',
                      sep = ""))

# Load the model object that includes all spatial polygons
load(paste(path, 'Model Results/Four Routes/Model_mean_beta_q0_bothoutliers.RData',
           sep = ""))


#####################################
###                               ###
###       Data Processing         ###
###                               ###
#####################################

# Define the results to plot from
mod = mod1
rm(mod1, mod2, mod3, mod4)

# Beta: 3 dimensional array with 10000 rows
# 2 columns (2 coefficients or parameters per row)
# 74 layers (1 per cookie)

# Remove the influential points 
data <- data[data$cookie != "116" & data$cookie != "5004", ]

cookies = data$cookie
year = data$year_val

y = data$mean_beta_q0
x = data$mean_beta_shannon_9_env

# Processing for plotting from model object
n = nrow(mod$beta[,2,])
m = ncol(mod$beta[,2,])
B = mod$eta[,2,]
A = mod$beta[,2,]

df = data.frame()
for (i in 1:m) {
  df = rbind(df, data.frame(A = A[,i],
                            B = B[,i],
                            Variable = paste0("Cookie", cookies[i])))
}

# Compute means
means = df %>%
  group_by(Variable) %>%
  summarise(mean_A = mean(A), mean_B = mean(B), .groups = "drop")


#####################################
###                               ###
###          Plotting             ###
###                               ###
#####################################

####################
###              ###
###   Figure S_  ###
###              ###
####################

p = ggplot(df, aes(x = A, y = B, color = Variable)) + 
  stat_ellipse(aes(fill = Variable), type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.15, color = NA) +
  geom_point(data = means, aes(x = mean_A, y = mean_B, color = Variable),
             size = 3, shape = 16, col=1, alpha = 0.5) +
  theme_minimal() +
  # geom_abline(slope = 1, intercept = 0, color = "red", lwd = 0.5, alpha = 0.3) +
  geom_abline(slope = mean(mod$phi[,2]), intercept = mean(mod$phi[,1]),
              color = "darkorchid4", lwd = 1.5) + # previously "#1F2933"
  geom_abline(slope = quantile(mod$phi[,2],0.95), 
              intercept = quantile(mean(mod$phi[,1]), 0.95),
              linetype = "dashed", color = "darkorchid4", lwd=1.2) +
  geom_abline(slope = quantile(mod$phi[,2],0.05),
              intercept = quantile(mean(mod$phi[,1]), 0.05),
              linetype = "dashed", color = "darkorchid4",lwd=1.2) +
  geom_hline(yintercept = 0,lwd=0.8) + geom_vline(xintercept = 0,lwd=0.8) +
  scale_x_continuous("% Change in environmental heterogeneity per \ntime interval",
                     limits = c(-0.055, 0.055), 
                     breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     labels = c("-11.9", "-5.6", "0.0", "5.9", "12.2")) +
  scale_y_continuous(
    name = expression(atop("% Change in " * beta [q==0] * " per ",
    "time interval")),
    limits = c(-0.015, 0.02),
    breaks = c(-0.01, -0.00, 0.01, 0.02),
    labels = c("-2.3", "0.0", "2.3", "4.7")) +
  theme(legend.position = 'none',
        text=element_text(size=18),
        axis.title.x = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y = element_text(size = 20,  color = "black")) 

p = ggMarginal(p, type="density", size=5, lwd=1)

# Save plot
ggsave(paste(path, "Manuscript Figures/FigS__q0_fourroutes_no_outliers.jpg",
             sep = ""),
       p,
       device = "jpeg", 
       height = 7, width = 9, dpi = 300)

ggsave(paste(path, "Manuscript Figures/FigS__q0_fourroutes_no_outliers.pdf",
             sep = ""),
       p,
       device = "pdf", 
       height = 7, width = 9, dpi = 300)


