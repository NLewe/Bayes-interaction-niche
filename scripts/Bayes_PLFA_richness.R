


## Bayes model  PLFA total  ####

# PLFA total = only Bacteria!!!! ## 

# packages ####
library (tidyverse)
library (tidybayes)
library (brms)

dataNLPL  <- 
  dataNLPL %>%  left_join(All_Metrics_E2_sample) %>%  filter (!is.na (totalPLFA)) %>% 
  left_join (RelGen_E1_E2_sample %>%  select (sampleID, RelGenSpec)) %>%  
  mutate (DWRA = DW_roots / DW_above)


# 0A Specify the Generalsim model 0####
## get default  priors  #### 
prior <- get_prior ( totalBact ~ RelGenSpec:DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull)
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     , data = dataNLPL, data2 = list(A = A))

# model 0 


PLFAtot_relGfit_0 <- brm(
  totalBact ~ RelGenSpec:DW_roots  + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) ,
  data = dataNLPL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

#control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(PLFAtot_relGfit_0)


# check convergence #
launch_shinystan(PLFAtot_relGfit_0)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_0, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")



# 1 Update to model 01 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_01 <- update(
  PLFAtot_relGfit_0, formula = ~ . - RelGenSpec:DW_roots + RelGenSpec ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_01)

PLFAtot_relGfit_01 %>%  bayes_R2()

# check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc


# check convergence #
launch_shinystan(PLFAtot_relGfit_01)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_01, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###
PLFAtot_relGfit_0 <- add_criterion(PLFAtot_relGfit_0, "loo")

PLFAtot_relGfit_01 <- add_criterion(PLFAtot_relGfit_01, "loo")
loo_compare (PLFAtot_relGfit_0, PLFAtot_relGfit_01)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_02 <- update(
  PLFAtot_relGfit_01, formula = ~ . + DW_roots  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

# check sampling quality of model 02 ##

# Look for rhat, ESS, Sd etc
summary (PLFAtot_relGfit_02)

# check convergence #
launch_shinystan(PLFAtot_relGfit_02)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_02, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

# Compare models ###

PLFAtot_relGfit_02 <- add_criterion(PLFAtot_relGfit_02, "loo")
loo_compare (PLFAtot_relGfit_0, PLFAtot_relGfit_02, PLFAtot_relGfit_01)



# 3 Update to model 03 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_03 <- update(
  PLFAtot_relGfit_01, formula = ~ .  + DW_above,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# check sampling quality of model 03 ###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_relGfit_03)

# check convergence #
launch_shinystan(PLFAtot_relGfit_03)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_03, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLFAtot_relGfit_03 <- add_criterion(PLFAtot_relGfit_03, "loo")
loo_compare (PLFAtot_relGfit_02, PLFAtot_relGfit_03, PLFAtot_relGfit_01, PLFAtot_relGfit_0)
### Model 02 is the best ###


# 4 Update to model 04 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_04 <- update(
  PLFAtot_relGfit_02, formula = ~ .  - RelGenSpec  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_04)

PLFAtot_relGfit_04 <- add_criterion(PLFAtot_relGfit_04, "loo")

loo_compare(PLFAtot_relGfit_04, PLFAtot_relGfit_0)



# 5 Update to model 05 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_05 <- update(
  PLFAtot_relGfit_04, formula = ~ .  -DW_roots + DW_above:DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
bayes_R2(PLFAtot_relGfit_04)




# check sampling quality of model 05###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_relGfit_05)

# check convergence #
launch_shinystan(PLFAtot_relGfit_05)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_05, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ##

PLFAtot_relGfit_05 <- add_criterion(PLFAtot_relGfit_05, "loo")
loo_compare (PLFAtot_relGfit_02, PLFAtot_relGfit_03, PLFAtot_relGfit_05, PLFAtot_relGfit_0)



# 6 Update to model 06 ####

#and then fit it again - adding more variables

PLFAtot_relGfit_06 <- update(
  PLFAtot_relGfit_05, formula = ~ .   - RelGenSpec ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# check sampling quality of model 06 ###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_relGfit_06)

# check convergence #
launch_shinystan(PLFAtot_relGfit_06)

# posterior predictive checks #

pp_check (PLFAtot_relGfit_06, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLFAtot_relGfit_06 <- add_criterion(PLFAtot_relGfit_06,"loo")
loo_compare (PLFAtot_relGfit_02, PLFAtot_relGfit_03, PLFAtot_relGfit_06, PLFAtot_relGfit_05)

#  7 Update to model 7 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_07 <- update(
  PLFAtot_relGfit_06, formula = ~ .   + RelGenSpec ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_07)


PLFAtot_relGfit_07 <- add_criterion(PLFAtot_relGfit_07, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06)

#  8 Update to model 8 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_08 <- update(
  PLFAtot_relGfit_07, formula = ~ .   + DW_roots - DW_roots:DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_08)


PLFAtot_relGfit_08 <- add_criterion(PLFAtot_relGfit_08, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_08)



#  9 Update to model 9 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_09 <- update(
  PLFAtot_relGfit_08, formula = ~ .   + DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_09)


PLFAtot_relGfit_09 <- add_criterion(PLFAtot_relGfit_09, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09)

#  9 Update to model 9B ####
#and then fit it again - adding more variables

PLFAtot_relGfit_09B <- update(
  PLFAtot_relGfit_09, formula = ~ .    + (1|gr(PlaSpe, cov = A))  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_09B)


PLFAtot_relGfit_09B <- add_criterion(PLFAtot_relGfit_09B, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_09B, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09)


#  9 Update to model 9C ####
#and then fit it again - adding more variables

PLFAtot_relGfit_09C <- update(
  PLFAtot_relGfit_09B, formula = ~ .    - (1 | PlantSpeciesfull)  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_09C)


PLFAtot_relGfit_09C <- add_criterion(PLFAtot_relGfit_09C, "loo")

loo_compare(PLFAtot_relGfit_09C, PLFAtot_relGfit_02, PLFAtot_relGfit_09B, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09)


#  9 Update to model 10 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_10 <- update(
  PLFAtot_relGfit_08, formula = ~ .   -DW_roots + DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_10)


PLFAtot_relGfit_10 <- add_criterion(PLFAtot_relGfit_10, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09, PLFAtot_relGfit_10)


#  11 Update to model 11 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_11 <- update(
  PLFAtot_relGfit_08, formula = ~ .   - DW_roots + DWRA ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PLFAtot_relGfit_11 %>% bayes_R2()


PLFAtot_relGfit_11 <- add_criterion(PLFAtot_relGfit_11, "loo")

loo_compare(PLFAtot_relGfit_05, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09, PLFAtot_relGfit_11)

#  9 Update to model 9 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_12 <- update(
  PLFAtot_relGfit_11, formula = ~ .   - RelGenSpec ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_12)


PLFAtot_relGfit_12 <- add_criterion(PLFAtot_relGfit_12, "loo")

loo_compare(PLFAtot_relGfit_12, PLFAtot_relGfit_02, PLFAtot_relGfit_07, 
            PLFAtot_relGfit_06, PLFAtot_relGfit_09)



### Plot Results of best model ####
dataNLPL %>%
  #group_by(PlaSpe) %>%
  #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
  add_epred_draws(PLFAtot_relGfit_09) %>%
  ggplot(aes(x = RelGenSpec             , y = totalPLFA, color = ordered(PlaSpe))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = dataNLPL) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
#+
 # facet_wrap(~PlaSpe, scales = "free_x")



plot (conditional_effects(PLFAtot_relGfit_09B, effects = "RelGenSpec",  ndraws = 1000, spaghetti = F,  
                          prob = 0.5, conditions = NULL), points =T)

plot (conditional_effects(PLFAtot_relGfit_09B, effects = "RelGenSpec:DW_above",  ndraws = 10000, spaghetti = F,  
                         prob = 0.9), points =T)

conditions2 = make_conditions (PLFAtot_relGfit_09B, "DW_above")


plot (conditional_effects(PLFAtot_relGfit_09B, effects = "RelGenSpec:DW_roots",  ndraws = 10000, 
                          spaghetti = F,  prob = 0.9, conditions = conditions2), points = T)

plot (conditional_effects(PLFAtot_relGfit_09, effects = "RelGenSpec",  ndraws = 10000, 
                          spaghetti = F,  prob = 0.9), points = T)

#conditions3= data.frame (DW_roots = c(0.25, 0.5, 0.75))
conditions2 = make_conditions (PLFAtot_relGfit_02, "DW_above")


plot (conditional_effects(PLFAtot_relGfit_02, effects = "RelGenSpec:DW_roots",  ndraws = 10000, spaghetti = F,  
                    prob = 0.5, conditions = conditions2), points = F)


#  

get_variables(PLFAtot_relGfit_02)
PLFAtot_relGfit_02 %>%  spread_draws(b_RelGenSpec, sigma) %>%  median_qi()


summary (lm (totalBact ~ DW_roots+RelGenSpec , data = dataNLPL))


