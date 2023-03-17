


## Bayes model  pd.obs ####


# packages ####
library (tidyverse)
library (tidybayes)
library (brms)

# 0A Specify the pd.obs model 0####
## get default  priors  #### 
prior <- get_prior ( AMF ~ pd.obs + (1|gr(PlaSpe, cov = A)) 
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     , data = dataNL, data2 = list(A = A))

# model 0 


PD_fit0 <- brm(
  AMF ~ pd.obs  + (1|gr(PlaSpe, cov = A)) ,
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(PD_fit0)


# 0B check sampling quality of model 0 ####

# Look for rhat, ESS, Sd etc
summary (PD_fit0)

# check convergence #
launch_shinystan(PD_fit0)

# posterior predictive checks #

pp_check (PD_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")



# 1A Update to model 01 ####

#and then fit it again - adding more variables

PD_fit01 <- update(
  PD_fit0, formula = ~ . + DW_above ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 1B check sampling quality of model 01 ####

# Look for rhat, ESS, Sd etc

summary (PD_fit01)

# check convergence #
launch_shinystan(PD_fit01)

# posterior predictive checks #

pp_check (PD_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
PD_fit0 <- add_criterion(PD_fit0, "loo")

PD_fit01 <- add_criterion(PD_fit01, "loo")
loo_compare (PD_fit0, PD_fit01)
# best performing model will be named at top
#



# 2A Update to model 02 ####

#and then fit it again - adding more variables

PD_fit02 <- update(
  PD_fit01, formula = ~ . + DW_roots - DW_above ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 2B check sampling quality of model 02 ####

# Look for rhat, ESS, Sd etc
summary (PD_fit02)

# check convergence #
launch_shinystan(PD_fit02)

# posterior predictive checks #

pp_check (PD_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

PD_fit02 <- add_criterion(PD_fit02, "loo")
loo_compare (PD_fit0, PD_fit02)



# 3A Update to model 03 ####

#and then fit it again - adding more variables

PD_fit03 <- update(
  PD_fit02, formula = ~ .  - (1 | gr(PlaSpe, cov = A)) + (1 + DW_roots| gr(PlaSpe, cov = A)) ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 3B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (PD_fit03)

# check convergence #
launch_shinystan(PD_fit03)

# posterior predictive checks #

pp_check (PD_fit03, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

PD_fit03 <- add_criterion(PD_fit03, "loo")
loo_compare (PD_fit02, PD_fit03)
### Model 02 is the best ###


#### Using all properties ---- they are not independent observations I guess
# cannot be used in same model ###

## Bayes model  observed ####


# packages ####
library (tidyverse)
library (tidybayes)
library (brms)

# 0A Specify the Shannon  model 0####
## get default  priors  #### 
prior <- get_prior ( AMF ~  Shannon +  (1|gr(PlaSpe, cov = A)) 
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     , data = dataNL, data2 = list(A = A))

# model 0 


ALL_fit0 <- brm(
  AMF ~  Shannon + (1|gr(PlaSpe, cov = A)) ,
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(ALL_fit0)


# 0B check sampling quality of model 0 ####

# Look for rhat, ESS, Sd etc
summary (ALL_fit0)

# check convergence #
launch_shinystan(ALL_fit0)

# posterior predictive checks #

pp_check (ALL_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")



# 1A Update to model 01 ####

#and then fit it again - adding more variables

Sh_fit01 <- update(
  ALL_fit0, formula = ~ . + DW_above ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 1B check sampling quality of model 01 ####

# Look for rhat, ESS, Sd etc
summary (Sh_fit01)

# check convergence #
launch_shinystan(ALL_fit01)

# posterior predictive checks #

pp_check (Obs_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
ALL_fit0 <- add_criterion(ALL_fit0, "loo")

Sh_fit01 <- add_criterion(Sh_fit01, "loo")
loo_compare (ALL_fit0, Sh_fit01)
# best performing model will be named at top
#



# 2A Update to model 02 ####

#and then fit it again - adding more variables

Sh_fit02 <- update(
  Sh_fit01, formula = ~ . + DW_roots - DW_above ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 2B check sampling quality of model 02 ####

# Look for rhat, ESS, Sd etc
summary (Obs_fit02)

# check convergence #
launch_shinystan(Obs_fit02)

# posterior predictive checks #

pp_check (Obs_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

Obs_fit02 <- add_criterion(Obs_fit02, "loo")
loo_compare (Obs_fit0, Obs_fit02)



# 3A Update to model 03 ####

#and then fit it again - adding more variables

Obs_fit03 <- update(
  Obs_fit02, formula = ~ .  - (1 | gr(PlaSpe, cov = A)) + (1 + DW_roots| gr(PlaSpe, cov = A)) ,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 3B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (Obs_fit03)

# check convergence #
launch_shinystan(Obs_fit03)

# posterior predictive checks #

pp_check (Obs_fit03, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

Obs_fit03 <- add_criterion(Obs_fit03, "loo")
loo_compare (Obs_fit02, Obs_fit03)
### Model 02 is the best ###






