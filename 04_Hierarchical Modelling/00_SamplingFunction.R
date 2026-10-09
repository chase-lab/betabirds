################################################################################
########                                                                ########
########          Function for Bayesian Hierarchical Model              ######## 
########                                                                ########
################################################################################


# This code contains the function to fit the Bayesian hierarchical model
# to log-transformed values of beta diversity and environmental heterogeneity.



########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls() ) 
gc() #releases memory


# Load required packages
library(mvnfast)


MCMC = function(y0,x0,cookies,year,R,burn) {
  require(mvnfast)
  
  # y=Birds , x=Env
  y = (y0-mean(y0))/sd(y0)
  x = (x0-mean(x0))/sd(x0)
  
  n_cook = length(unique(cookies))
  n_time = length(unique(year))
  time = 1:n_time
  
  W = cbind(1,time)
  WtW = t(W)%*%W
  
  Y = matrix(y,n_cook,n_time)
  tWY = t(apply(Y, 1, function(y) t(W)%*%y))
  X = matrix(x,n_cook,n_time)
  tWX = t(apply(X, 1, function(x) t(W)%*%x))
  
  # Slope Env
  mu_beta0 = c(0,0)
  Omega_beta0 = diag(1,2)
  
  # Slope Birds
  mu_eta0 = c(0,0)
  Omega_eta0 = diag(1,2)
  
  # Slope slopes
  mu_phi0 = c(0,0)
  Omega_phi0 = diag(1,2)
  
  beta = eta = array(0,dim=c(R,2,n_cook))
  sigma2_Y = sigma2_X = matrix(1,n_cook,R)
  phi = matrix(1,R,2)
  sigma2_v = rep(1,R)
  
  LLx = LLy = LLxy = matrix(NA,R,n_cook)
  
  for (r in 2:R) {
    beta[r,,] = beta[r-1,,] 
    sigma2_Y[,r] = sigma2_Y[,r-1]
    eta[r,,] = eta[r-1,,] 
    sigma2_X[,r] = sigma2_X[,r-1]
    phi[r,] = phi[r-1,]
    sigma2_v[r] = sigma2_v[r-1]
    
    for (i in 1:n_cook) {
      Sigma_beta = solve(WtW/sigma2_X[i,r] + Omega_beta0)
      mu_beta = Sigma_beta %*% (tWX[i,]/sigma2_X[i,r] + Omega_beta0%*%mu_beta0)
      beta[r,,i] = rmvn(1,mu_beta,Sigma_beta)
      
      sigma2_X[i,r] = 1/rgamma(1, 0.01 + n_time/2, 0.01 + 1/2*sum((X[i,]-W%*%beta[r,,i])^2))
      
      Sigma_eta = solve(WtW/sigma2_Y[i,r] + Omega_eta0)
      mu_eta = Sigma_eta %*% (tWY[i,]/sigma2_Y[i,r] + Omega_eta0%*%mu_eta0)
      eta[r,,i] = rmvn(1,mu_eta,Sigma_eta)
      
      sigma2_Y[i,r] = 1/rgamma(1, 0.01 + n_time/2, 0.01 + 1/2*sum((Y[i,]-W%*%eta[r,,i])^2))
      
      LLx[r,i] = sum(dnorm(X[i,],beta[r,1,i] + beta[r,2,i]*time, sqrt(sigma2_X[i,r]),log=T))
      LLy[r,i] = sum(dnorm(Y[i,],eta[r,1,i] + eta[r,2,i]*time, sqrt(sigma2_Y[i,r]),log=T))
    }
    
    Z = cbind(1,beta[r,2,])
    Sigma_phi = solve(t(Z)%*%Z/sigma2_v[r] + Omega_phi0)
    mu_phi = Sigma_phi %*% (t(Z)%*%eta[r,2,]/sigma2_v[r] + Omega_phi0%*%mu_phi0)
    phi[r,] = rmvn(1,mu_phi,Sigma_phi)
    
    sigma2_v[r] = 1/rgamma(1, 0.01 + n_cook/2, 0.01 + 1/2*sum((eta[r,2,]-Z%*%phi[r,])^2))
    
    LLxy[r,] = dnorm(eta[r,2,],phi[r,1] + phi[r,2]*beta[r,2,], sqrt(sigma2_v[r]),log=T)
  }
  
  beta[,2,] = beta[,2,]*sd(x0)
  eta[,2,] = eta[,2,]*sd(y0)
  phi[,1] = phi[,1]*sd(y0)
  phi[,2] = phi[,2]*sd(y0)/sd(x0)
  
  list(beta = beta[(burn+1):R,,],
       sigma2_Y = sigma2_Y[,(burn+1):R],
       eta = eta[(burn+1):R,,],
       sigma2_X = sigma2_X[,(burn+1):R],
       phi = phi[(burn+1):R,],
       sigma2_v = sigma2_v[(burn+1):R],
       LLx = LLx[(burn+1):R,],
       LLy = LLy[(burn+1):R,],
       LLxy = LLxy[(burn+1):R,]
  )
}
