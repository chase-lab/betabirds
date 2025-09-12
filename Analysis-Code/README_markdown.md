\[2,4]-Routes folders contain the R scripts to reproduce the figures in the paper, for both the 4 routes and pairs datasets.



In each folder:

&nbsp;	- a **.csv** file with the data



&nbsp;	- **SamplingFunction.R** : the MCMC sampler used to estimate the Bayesian model.



&nbsp;	- **Main.R** : This script uses **SamplingFunction.R** to draw samples from the posterior. Results for all cookies (*Model\_mean\_beta\_\[C80,q0,q2,sn200].RData)* and when outlier cookie 5004 has been removed (*Model\_mean\_beta\_\[C80,q0,q2,sn200]\_outlier.RData)* are stored in **Results** folder as .RData files



&nbsp;	- **\[2,4]-\_mean\_beta\_\[C80,q0,q2,sn200].R** : R script to crete figures for 4 routes or pairs, and for each variable C80,q0,q2,sn200.



&nbsp;	- **\[2,4]-\_mean\_beta\_\[C80,q0,q2,sn200]-outliers.R** : R script to crete figures for 4 routes or pairs, and for each variable C80,q0,q2,sn200 *when outlier cookie 5004 has been removed*.



&nbsp;	- **\[2,4]-\_slopes.R** : R script to crete figures for 4 routes or pairs, with the posterior densities of the slope between C80,q0,q2,sn200 and "mean\_beta\_shannon\_9\_env". Posterior probabilities of being >,< 0 are computed in this script.



&nbsp;	- **\[2,4]-\_slopes-outliers.R** : R script to crete figures for 4 routes or pairs, with the posterior densities of the slope between C80,q0,q2,sn200 and "mean\_beta\_shannon\_9\_env" *when outlier cookie 5004 has been removed*. Posterior probabilities of being >,< 0 are computed in this script.



&nbsp;	- **Figures**: folder with the figures .pdf





