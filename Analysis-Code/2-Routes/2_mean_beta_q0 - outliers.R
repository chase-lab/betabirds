rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggExtra)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(ggmagnify)
library(ggfx)
library(ggforce)

data = read.csv('pairs_beta_birds_environ_cookie_20250708.csv')
data = data[data$cookie != 5004,]


cookies = as.numeric(as.character(data$cookie))
year = data$year_val

y = data$mean_beta_q0
x = data$mean_beta_shannon_9_env

load('Results/Model_mean_beta_q0_outlier.RData')

# Processing for plotting
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

# Plot
p = ggplot(df, aes(x = A, y = B, color = Variable)) + 
  stat_ellipse(aes(fill = Variable), type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.15, color = NA) +
  geom_point(data = means, aes(x = mean_A, y = mean_B, color = Variable),
             size = 3, shape = 16, col=1) +
  theme_minimal() +
  labs(title = "mean_beta_q0 ~ mean_beta_shannon_9_env",
       x = "Environmental heterogeneity - slope", y = "US Birds Beta diversity - slope", color = "Variable") +
  geom_abline(slope = mean(mod$phi[,2]), intercept = mean(mod$phi[,1]), color = "red3",lwd=1.2) +
  geom_abline(slope = quantile(mod$phi[,2],0.95), intercept = quantile(mean(mod$phi[,1]),0.95), linetype = "dashed", color = "red3",lwd=1) +
  geom_abline(slope = quantile(mod$phi[,2],0.05), intercept = quantile(mean(mod$phi[,1]),0.05), linetype = "dashed", color = "red3",lwd=1) +
  geom_hline(yintercept = 0,lwd=0.8) + geom_vline(xintercept = 0,lwd=0.8) +
  geom_label(data = means[(means$mean_A>0.027)&(means$mean_B>0.01),],
             aes(x = mean_A+0.014, y = mean_B, label = Variable), color='black',fill = "white") +
  geom_label(data = means[(means$mean_A>0.025)&(means$mean_A<0.027)&(means$mean_B>0.01),],
             aes(x = mean_A+0.006, y = mean_B+0.006, label = Variable), color='black',fill = "white") +
  theme(legend.position = 'none',
        text=element_text(size=18)) 

p = ggMarginal(p, type="density", size=5, lwd=1)

#Display it
p

# Save
ggsave("Figures/2_mean_beta_q0_outlier.pdf", plot = p, width = 12, height = 8)
