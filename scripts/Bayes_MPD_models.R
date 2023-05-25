

# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)


# 
data_metrics<-  
  dataNL_sample %>% 
  left_join (All_Metrics_E2_sample) %>%  
  filter (!is.na (AMF))
  

library (broom.mixed)
# 0A Specify the MPD model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ MPD + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     data = data_metrics, data2 = list(A = A))

# model 0 


mpd_fit0 <- brm(
  AMF ~ MPD  + (1 |gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
  data = data_metrics,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior0,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.8) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(mpd_fit0)

# estimate of phylogenetic signal!!

hyp0 <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"

hyp0 <- hypothesis(mpd_fit0, hyp0, class = NULL)

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
summary (mpd_fit0)

# check convergence #
launch_shinystan(mpd_fit0)

# posterior predictive checks #

pp_check (mpd_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")



# 1A Update to model 01 ####

#and then fit it again - adding more variables

mpd_fit01 <- update(
  mpd_fit0, formula = ~ . -MPD + MPD:DW_roots ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (mpd_fit01)



# 1B check sampling quality of model 01 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit01)

# check convergence #
launch_shinystan(mpd_fit01)

# posterior predictive checks #

pp_check (mpd_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
mpd_fit0 <- add_criterion(mpd_fit0, "loo")

mpd_fit01 <- add_criterion(mpd_fit01, "loo")
loo_compare (Gen_samples_fit01, mpd_fit01)
# best performing model will be named at top
#



# 2A Update to model 02 ####

#and then fit it again - adding more variables

mpd_fit02 <- update(
  mpd_fit01, formula = ~ . + DW_roots - DW_above ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 2B check sampling quality of model 02 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit02)

# check convergence #
launch_shinystan(mpd_fit02)

# posterior predictive checks #

pp_check (mpd_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

mpd_fit02 <- add_criterion(mpd_fit02, "loo")
loo_compare (mpd_fit0, mpd_fit02)



# 3A Update to model 03 ####

#and then fit it again - adding more variables

mpd_fit03 <- update(
  mpd_fit02, formula = ~ .  - (1 | gr(PlaSpe, cov = A)) + (1 + DW_roots| gr(PlaSpe, cov = A)) ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)





# 3B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit03)

# check convergence #
launch_shinystan(mpd_fit03)

# posterior predictive checks #

pp_check (mpd_fit03, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####

mpd_fit03 <- add_criterion(mpd_fit03, "loo")
loo_compare (mpd_fit02, mpd_fit0, mpd_fit03)
### Model 02 is the best ###


sum_mpd02_B <-  tidy (mpd_fit02, effects = c("fixed"))





### get posteriors from the model

# 
data_metrics %>%
  #group_by(PlaSpe) %>%
  #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
  add_epred_draws(mpd_fit01) %>%
  ggplot(aes(x = MPD, y = AMF, color = ordered(PlaSpe))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = data_metrics) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~PlaSpe, scales = "free_x")


plot (conditional_effects(mpd_fit01, effects = "MPD:DW_roots", ndraws = 1000, spaghetti = F, mean = T, prob = 0.9), points = T)





