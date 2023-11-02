




# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)

# data 
PCA_metric_data_sample <- readRDS ("results/PCA_metric_data_sample.rds")

dataNL_sample <- dataNL %>%  left_join (PCA_metric_data_sample) %>% 
  left_join(RelGen_E1_E2_sample ) %>% 
  filter (!is.na (AMF),  !is.na (Observed))



)# 0A Specify the MPD model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ Dim.1 + Dim.2 + Dim.3 + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),#+ (1|PlantSpeciesfull),
                      #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                      #+ OTHEr than phylogenetic (environmental factors, niches)
                      data = dataNL_sample, data2 = list(A = A))

# model 0 


PC_samples_fit0 <- brm(
  AMF ~ Dim.1 + Dim.2 + Dim.3  + (1 |gr(PlaSpe, cov = A))+ (1|PlantSpeciesfull),
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
summary(PC_samples_fit0)

# estimate of phylogenetic signal!!

hyp0 <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"

hyp0 <- hypothesis(PC_samples_fit0, hyp0, class = NULL)

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
summary (PC_samples_fit0)

# check convergence #
launch_shinystan(PC_samples_fit0)

# posterior predictive checks #

pp_check (PC_samples_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")



# 1A Update to model 01 ####

#and then fit it again - adding more variables

PC_samples_fit01 <- update(
  PC_samples_fit0, formula = ~ . + DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)




# 1B check sampling quality of model 01 ####

# Look for rhat, ESS, Sd etc
summary (PC_samples_fit01)

# check convergence #
launch_shinystan(PC_samples_fit01)

# posterior predictive checks #

pp_check (PC_samples_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
PC_samples_fit0 <- add_criterion(PC_samples_fit0, "loo")

PC_samples_fit01 <- add_criterion(PC_samples_fit01, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit01)
# best performing model will be named at top
#



# 2A Update to model 02 ####

#and then fit it again - adding more variables

PC_samples_fit02 <- update(
  PC_samples_fit01, formula = ~ .  - DW_above + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 2B check sampling quality of model 02 ####

# Look for rhat, ESS, Sd etc
summary (PC_samples_fit02)

# check convergence #
launch_shinystan(PC_samples_fit02)

# posterior predictive checks #

pp_check (PC_samples_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

PC_samples_fit02 <- add_criterion(PC_samples_fit02, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit02)



# 3A Update to model 03 ####

#and then fit it again - adding more variables

PC_samples_fit03 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit03 <- add_criterion(PC_samples_fit03, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit03, PC_samples_fit02)


# 4 Update to model 04 ####

#and then fit it again - adding more variables

PC_samples_fit04 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 + DW_roots - Dim.2,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit04 <- add_criterion(PC_samples_fit04, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit04, PC_samples_fit03, PC_samples_fit02)

summary (PC_samples_fit04)

launch_shinystan(PC_samples_fit04)

PC_samples_fit04B <- update(
  PC_samples_fit04, formula = ~ .   + DW_roots * Dim.3,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit04B <- add_criterion(PC_samples_fit04B, "loo")
loo_compare (PC_samples_fit04B, PC_samples_fit04, PC_samples_fit03, PC_samples_fit02)

summary (PC_samples_fit04B)

launch_shinystan(PC_samples_fit04B)

# Interaction term does not improve the model


# 5A Update to model 05 ####

#and then fit it again - adding more variables

PC_samples_fit05 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 - Dim.3 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit05 <- add_criterion(PC_samples_fit05, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit03, PC_samples_fit05, PC_samples_fit04)


# Look for rhat, ESS, Sd etc
summary (PC_samples_fit05)

# 6 Update to model 06 ####

#and then fit it again - adding more variables

PC_samples_fit06 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 - Dim.3 - Dim.2 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit06 <- add_criterion(PC_samples_fit06, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit03, PC_samples_fit06, PC_samples_fit04)


# Look for rhat, ESS, Sd etc
summary (PC_samples_fit06)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

PC_samples_fit07 <- update(
  PC_samples_fit06, formula = ~ .   + DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit07 <- add_criterion(PC_samples_fit07, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit07, PC_samples_fit06, PC_samples_fit04)

# 3B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (PC_samples_fit07)


# 5A Update to model 05 ####

#and then fit it again - adding more variables

PC_samples_fit08 <- update(
  PC_samples_fit06, formula = ~ .  + Dim.1 ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit08 <- add_criterion(PC_samples_fit08, "loo")
loo_compare (PC_samples_fit06, PC_samples_fit04, PC_samples_fit08, PC_samples_fit07)



# 5A Update to model 05 ####

#and then fit it again - adding more variables

PC_samples_fit09 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 - Dim.3 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit09 <- add_criterion(PC_samples_fit09, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit06, PC_samples_fit09, PC_samples_fit04)


# 5A Update to model 05 ####

#and then fit it again - adding more variables

PC_samples_fit10 <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 - Dim.2 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_samples_fit10 <- add_criterion(PC_samples_fit1, "loo")
loo_compare (PC_samples_fit0, PC_samples_fit10, PC_samples_fit06, PC_samples_fit04)










# check convergence #
launch_shinystan(PC_samples_fit04)

# posterior predictive checks #

pp_check (PC_samples_fit03, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

PC_samples_fit03 <- add_criterion(PC_samples_fit03, "loo")
loo_compare (PC_samples_fit02, PC_samples_fit03)
### Model 02 is the best ###


sum_mpd02_B <-  tidy (PC_samples_fit02, effects = c("fixed"))


PC_test  <- update(
  PC_samples_fit0, formula = ~ .  - Dim.1 - Dim.2 -Dim.3 + DW_roots:Dim.3 ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PC_test <- add_criterion(PC_test, "loo")
loo_compare (PC_test, PC_samples_fit0)



### get posteriors from the model

# 
dataNL_sample %>%
  #group_by(PlaSpe) %>%
  #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
  add_epred_draws(PC_samples_fit06) %>%
  ggplot(aes(x = DW_roots, y = AMF, color = ordered(PlaSpe))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = dataNL_sample) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~PlaSpe, scales = "free_x")


conditional_effects(PC_samples_fit02, ndraws = 100, spaghetti = F, mean = T, prob = 0.5)





