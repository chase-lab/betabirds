###############################################################################
########                                                               ########
########                  Modelling Habitat Associations               ########
########                                                               ########
###############################################################################

# This code models relative abundance of each species as a function of
# time and the number of habitat associations of each species, with
# random effects for metacommunity and species.

########## clean workspace and load required packages ####################

# clean workspace to improve efficiency: #
rm(list = ls())
gc() #releases memory

library(tidyverse) # Data cleaning
library(brms)

path <- "~/share/groups/mas/01_projects/iUpdate/Beta Diversity with Pao_paper/"


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Load the data with relative abundance and habitat associations
hab.freq_long <- read.csv(
  paste(
    path,
    "01_Analysis/IUCN Habitats/Abundance and Habitat data.csv",
    sep = ""
  ),
  header = TRUE
)


###############################################
###                                         ###
###               Modeling                  ###
###                                         ###
###############################################

# Define formula
formula_brm <- bf(
  RelAbundance_beta ~ Time *
    n_habitat_associations +
    (Time || Species1) +
    (Time || ID_200)
)

family_brm <- Beta()

priors_brm <- c(
  # Fixed effects on logit scale
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 5), class = "Intercept"),

  # Random effect SDs (Species1, ID__200)
  prior(exponential(1), class = "sd"), # generic for all group-level SDs

  # Beta precision (phi) on log scale
  prior(exponential(1), class = "phi")
)

# run the model
start.time <- Sys.time()
m2_brm <- brm(
  formula = formula_brm,
  data = hab.freq_long,
  family = family_brm,
  prior = priors_brm,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
end.time <- Sys.time()
end.time - start.time
# Took about 2 hours on the rstudio server

# Traceplots
plot(m2_brm)

# Posterior predictive checks
pp_check(m2_brm)

# Save the model object
saveRDS(
  m2_brm,
  file = paste(
    path,
    "01_Analysis/IUCN Habitats/Model Objects/Habitats.rds",
    sep = ""
  )
)
