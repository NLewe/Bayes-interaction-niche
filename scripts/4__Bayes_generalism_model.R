






# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)

# data 
PCA_metric_data_sample <- readRDS ("results/PCA_metric_data_sample.rds")

dataNL_sample <- dataNL %>%  left_join (PCA_metric_data_sample) %>% 
  left_join(RelGen_E1_E2_sample %>%  select (sampleID, RelGenSpec) ) %>% 
  filter (!is.na (AMF),  !is.na (RelGenSpec))



# 0A Specify the RelGenSpecialsist model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ RelGenSpec * DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),#+ (1|PlantSpeciesfull),
                      #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                      #+ OTHEr than phylogenetic (environmental factors, niches)
                      data = dataNL_sample, data2 = list(A = A))

# model 0 


Gen_samples_fit0 <- brm(
  AMF ~ RelGenSpec * DW_roots  + (1 |gr(PlaSpe, cov = A))+ (1|PlantSpeciesfull),
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior0,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

#control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(Gen_samples_fit0)

# estimate of phylogenetic signal!!

hyp0 <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"

hyp0 <- hypothesis(Gen_samples_fit0, hyp0, class = NULL)

hyp0
#Note that the phylogenetic signal is just a synonym of the intra-class correlation (ICC) used in the context phylogenetic analysis.
#the intraclass correlation, or the intraclass correlation coefficient (ICC),[1] is a 
#descriptive statistic that can be used when quantitative measurements are made on units that are organized into groups. 
#It describes how strongly units in the same group resemble each other. 
#While it is viewed as a type of correlation, unlike most other correlation measures it operates 
#on data structured as groups, rather than data structured as paired observations. 


### means we ahave a phylogenetc signal 
plot (hyp0)

# 0B check sampling quality of model 0 ####

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit0)

# check convergence #
launch_shinystan(Gen_samples_fit0)

# posterior predictive checks #

pp_check (Gen_samples_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")



# 1A Update to model 01 ####

#and then fit it again - adding more variables

Gen_samples_fit01 <- update(
  Gen_samples_fit0, formula = ~ . - DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)




# 1B check sampling quality of model 01 ####

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit01)

# check convergence #
launch_shinystan(Gen_samples_fit01)

# posterior predictive checks #

pp_check (Gen_samples_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
Gen_samples_fit0 <- add_criterion(Gen_samples_fit0, "loo")

Gen_samples_fit01 <- add_criterion(Gen_samples_fit01, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit01)
# best performing model will be named at top
#



# 2A Update to model 02 ####

#and then fit it again - adding more variables

Gen_samples_fit02 <- update(
  Gen_samples_fit0, formula = ~ .  -  RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 2B check sampling quality of model 02 ####

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit02)

# check convergence #
launch_shinystan(Gen_samples_fit02)

# posterior predictive checks #

pp_check (Gen_samples_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

Gen_samples_fit02 <- add_criterion(Gen_samples_fit02, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit01, Gen_samples_fit02)





# 4 Update to model 04 ####
Gen_samples_fit03 <- update(
  Gen_samples_fit0, formula = ~ .  -  RelGenSpec:DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit03 <- add_criterion(Gen_samples_fit03, "loo")
loo_compare (Gen_samples_fit0,  Gen_samples_fit01, 
             Gen_samples_fit02, Gen_samples_fit03)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit03)
#and then fit it again - adding more variables

Gen_samples_fit04 <- update(
  Gen_samples_fit0, formula = ~ .  -  DW_roots - RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit04 <- add_criterion(Gen_samples_fit04, "loo")
loo_compare (Gen_samples_fit0,  Gen_samples_fit01, 
             Gen_samples_fit02, Gen_samples_fit03, 
             Gen_samples_fit04)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit04)

# 5A Update to model 05 ####

#and then fit it again - adding more variables

Gen_samples_fit05 <- update(
  Gen_samples_fit04, formula = ~ .  - (1|PlantSpeciesfull),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit05 <- add_criterion(Gen_samples_fit05, "loo")
loo_compare (Gen_samples_fit0,  Gen_samples_fit05, Gen_samples_fit04, 
             Gen_samples_fit03)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05)

# 6 Update to model 06 ####

#and then fit it again - adding more variables

Gen_samples_fit06 <- update(
  Gen_samples_fit04, formula = ~ .   + DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit06 <- add_criterion(Gen_samples_fit06, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit05, Gen_samples_fit06, Gen_samples_fit04)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit06)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

Gen_samples_fit07 <- update(
  Gen_samples_fit04, formula = ~ .   + DW_roots:DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit07 <- add_criterion(Gen_samples_fit07, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit07, Gen_samples_fit06, Gen_samples_fit04)

# 3B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit07)


# 5A Update to model 05 ####

#and then fit it again - adding more variables

Gen_samples_fit08 <- update(
  Gen_samples_fit04, formula = ~ .  + RelGenSpec:DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit08 <- add_criterion(Gen_samples_fit08, "loo")
loo_compare (Gen_samples_fit06, Gen_samples_fit0, Gen_samples_fit08, Gen_samples_fit04)



# 5A Update to model 05 ####

# #and then fit it again - adding more variables
# 
# Gen_samples_fit09 <- update(
#   Gen_samples_fit08, formula = ~ .  - DW_roots:RelGenSpec ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit09 <- add_criterion(Gen_samples_fit09, "loo")
# loo_compare (Gen_samples_fit0, Gen_samples_fit06, Gen_samples_fit09, Gen_samples_fit08)
# 
# 
# # 5A Update to model 05 ####
# 
# #and then fit it again - adding more variables
# 
# Gen_samples_fit10 <- update(
#   Gen_samples_fit08, formula = ~ .  - DW_above:DW_roots ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit10 <- add_criterion(Gen_samples_fit10, "loo")
# loo_compare (Gen_samples_fit08, Gen_samples_fit10, Gen_samples_fit09, Gen_samples_fit0)
# 
# 
# 
# Gen_samples_fit11 <- update(
#   Gen_samples_fit10, formula = ~ .  + RelGenSpec ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit11 <- add_criterion(Gen_samples_fit11, "loo")
# loo_compare (Gen_samples_fit08, Gen_samples_fit10, Gen_samples_fit09, Gen_samples_fit11)
# 
# 
# Gen_samples_fit12 <- update(
#   Gen_samples_fit10, formula = ~ .  - RelGenSpec:DW_above ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit12 <- add_criterion(Gen_samples_fit12, "loo")
# loo_compare (Gen_samples_fit08, Gen_samples_fit10, Gen_samples_fit11, Gen_samples_fit12)
# 
# Gen_samples_fit13 <- update(
#   Gen_samples_fit12, formula = ~ .  + DW_roots ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit13 <- add_criterion(Gen_samples_fit13, "loo")
# loo_compare (Gen_samples_fit08, Gen_samples_fit13, Gen_samples_fit09, Gen_samples_fit12)
# 
# 
# 
# Gen_samples_fit14 <- update(
#   Gen_samples_fit12, formula = ~ .  + RelGenSpec ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit14 <- add_criterion(Gen_samples_fit14, "loo")
# loo_compare (Gen_samples_fit14, Gen_samples_fit10, Gen_samples_fit12, Gen_samples_fit13)
# 
# Gen_samples_fit15 <- update(
#   Gen_samples_fit14, formula = ~ .  - RelGenSpec:DW_roots ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit15 <- add_criterion(Gen_samples_fit15, "loo")
# loo_compare (Gen_samples_fit13, Gen_samples_fit12, Gen_samples_fit14, Gen_samples_fit15)
# 
# 
# 
# Gen_samples_fit16 <- update(
#   Gen_samples_fit12, formula = ~ .  - DW_roots:RelGenSpec + RelGenSpec:DW_roots ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# Gen_samples_fit16 <- add_criterion(Gen_samples_fit16, "loo")
# loo_compare (Gen_samples_fit16, Gen_samples_fit15, Gen_samples_fit0, Gen_samples_fit12)
# 
# 
# 
# 





# check convergence #
launch_shinystan(Gen_samples_fit04)

# posterior predictive checks #

pp_check (Gen_samples_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

Gen_samples_fit03 <- add_criterion(Gen_samples_fit03, "loo")
loo_compare (Gen_samples_fit02, Gen_samples_fit03)
### Model 02 is the best ###


sum_mpd02_B <-  tidy (Gen_samples_fit04, effects = c("fixed"))





### get posteriors from the model####

# 
dataNL_sample %>%
  #group_by(PlaSpe) %>%
  #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
  add_epred_draws(Gen_samples_fit05) %>%
  ggplot(aes(x = RelGenSpec, y = AMF, color = ordered(PlaSpe))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = dataNL_sample) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~PlaSpe, scales = "free_x")


#conditions = make_conditions (Gen_samples_fit04, "RelGenSpec")


plot (conditional_effects(Gen_samples_fit05, effects = "DW_roots",  ndraws = 10000, spaghetti = F,  
                    prob = 0.9, conditions = conditions), points =F)

#plot (conditional_effects(Gen_samples_fit04, effects = "DW_roots:RelGenSpec",  ndraws = 10000, spaghetti = F,  
#                          prob = 0.9), points =T)

#conditions2 = make_conditions (Gen_samples_fit04, "DW_roots")


#plot (conditional_effects(Gen_samples_fit04, effects = "RelGenSpec",  ndraws = 10000, spaghetti = F,  prob = 0.9, conditions = conditions2), points = T)

plot (conditional_effects(Gen_samples_fit05, effects = "RelGenSpec:DW_roots",  ndraws = 10000, 
                          spaghetti = F,  prob = 0.5))


#conditions3= data.frame (DW_roots = c(0.25, 0.5, 0.75))
conditions2 = make_conditions (Gen_samples_fit04, "DW_roots")

conditional_effects(Gen_samples_fit04, effects = "RelGenSpec",  ndraws = 10000, spaghetti = F,  prob = 0.8, conditions = conditions2)

launch_shinystan(Gen_samples_fit04
                 )
#Run with better iterations
Gen_samples_fit04 <- update(
  Gen_samples_fit0, formula = ~ .  -  DW_roots - RelGenSpec ,
  newdata = dataNL_sample, chains = 8, cores = 8,
  iter = 10000, warmup = 4000
)

