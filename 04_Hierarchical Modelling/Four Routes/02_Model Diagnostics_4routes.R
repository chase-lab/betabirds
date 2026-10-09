################################################################################
########                                                                ########
########                      Evaluating Model Fit                      ######## 
########                                                                ########
################################################################################


# This code runs model diagnostics on Bayesian hierarchical models fit to 
# to log-transformed values of beta diversity and environmental heterogeneity
# for three beta diversity metrics: q0, q1, and q2 derived from
# random sampling of 4 routes within each spatial polygon using the sampling function 
# contained in 00_SamplingFunction.R.



########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory

library(coda)
library(loo)

# Define file path
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"


###############################################
###                                         ###
###      Data and Function Loading          ###
###                                         ###
###############################################

# Load the sampling function
source(paste(path, "04_Hierarchical Modelling/00_SamplingFunction.R",
             sep = ""))


###############################################
###                                         ###
###             Data Preparation            ###
###                                         ###
###############################################

dgn = list()

names = c("q0","q1", "q2")

# Define file path again
path <- "I:/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/03_Code/Code for publishing/"



###############################################
###                                         ###
###                Diagnostics              ###
###                                         ###
###############################################

for (d in 1:3) {
  
  load(paste0(path, 'Model Results/Four Routes/Model_mean_beta_',
              names[d], '.RData'))
  
  # Rhat and effective sample size
  rhatLs = list()
  
  rhat = ess = rep(NA,66)
  
  for (i in 1:66) {
    mlist = list(mcmc(mod1$beta[,2,i]),
                 mcmc(mod2$beta[,2,i]),
                 mcmc(mod3$beta[,2,i]),
                 mcmc(mod4$beta[,2,i]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  
  rhatLs$beta2 = list(rhat=rhat,ess=ess)
  
  rhat = ess = rep(NA,66)
  for (i in 1:66) {
    mlist = list(mcmc(mod1$beta[,1,i]),
                 mcmc(mod2$beta[,1,i]),
                 mcmc(mod3$beta[,1,i]),
                 mcmc(mod4$beta[,1,i]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  rhatLs$beta1 = list(rhat=rhat,ess=ess)
  
  
  rhat = ess = rep(NA,66)
  for (i in 1:66) {
    mlist = list(mcmc(mod1$eta[,2,i]),
                 mcmc(mod2$eta[,2,i]),
                 mcmc(mod3$eta[,2,i]),
                 mcmc(mod4$eta[,2,i]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  rhatLs$eta2 = list(rhat=rhat,ess=ess)
  
  
  rhat = ess = rep(NA,66)
  for (i in 1:66) {
    mlist = list(mcmc(mod1$eta[,1,i]),
                 mcmc(mod2$eta[,1,i]),
                 mcmc(mod3$eta[,1,i]),
                 mcmc(mod4$eta[,1,i]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  rhatLs$eta1 = list(rhat=rhat,ess=ess)
  
  
  rhat = ess = rep(NA,66)
  for (i in 1:66) {
    mlist = list(mcmc(mod1$sigma2_Y[i,]),
                 mcmc(mod2$sigma2_Y[i,]),
                 mcmc(mod3$sigma2_Y[i,]),
                 mcmc(mod4$sigma2_Y[i,]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  rhatLs$sY = list(rhat=rhat,ess=ess)
  
  
  rhat = ess = rep(NA,66)
  for (i in 1:66) {
    mlist = list(mcmc(mod1$sigma2_X[i,]),
                 mcmc(mod2$sigma2_X[i,]),
                 mcmc(mod3$sigma2_X[i,]),
                 mcmc(mod4$sigma2_X[i,]))
    rhat[i] = gelman.diag(mlist)[1]$psrf[1]
    ess[i] = effectiveSize(mlist)
  }
  rhatLs$sX = list(rhat=rhat,ess=ess)
  
  
  
  mlist = list(mcmc(mod1$phi[,1]),
               mcmc(mod2$phi[,1]),
               mcmc(mod3$phi[,1]),
               mcmc(mod4$phi[,1]))
  rhatLs$phi1 = list(rhat=gelman.diag(mlist)[1]$psrf[1],ess=effectiveSize(mlist))
  
  
  mlist = list(mcmc(mod1$phi[,2]),
               mcmc(mod2$phi[,2]),
               mcmc(mod3$phi[,2]),
               mcmc(mod4$phi[,2]))
  rhatLs$phi2 = list(rhat=gelman.diag(mlist)[1]$psrf[1],ess=effectiveSize(mlist))
  
  mlist = list(mcmc(mod1$sigma2_v),
               mcmc(mod2$sigma2_v),
               mcmc(mod3$sigma2_v),
               mcmc(mod4$sigma2_v))
  rhatLs$sXY = list(rhat=gelman.diag(mlist)[1]$psrf[1],ess=effectiveSize(mlist))
  
  # Pareto k
  loo1 = loo(mod1$LLxy,save_psis = TRUE)
  rhatLs$pk = loo1$psis_object$diagnostics$pareto_k
  
  dgn[[d]] = rhatLs
  
  print(names[d])
} # end of diagnostics loop


###############################################
###                                         ###
###              Save Results               ###
###                                         ###
###############################################

names(dgn) = names

save(list = c('dgn'),
     file = paste(path, 'Model Results/Four Routes/Diagnostics.RData',
                  sep = ""))


