






# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library(GGally)

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# https://mc-stan.org/loo/articles/online-only/faq.html


# data 
PCA_metric_data_sample <- readRDS ("results/PCA_metric_data_sample.rds")

dataNL_sample <- dataNL %>%  left_join (PCA_metric_data_sample) %>% 
  left_join(RelGen_E1_E2_sample %>%  dplyr::select (sampleID, RelGenSpec  )) %>% 
  filter (!is.na (AMF),  !is.na (RelGenSpec)) %>% 
  filter (sampleID != "R18")  %>%   ## plant lost shoot biomass as it was almost dead
    dplyr::select (!starts_with ("Dim")) %>% 
  mutate (DWRA = DW_roots /DW_above )#%>% 
 # mutate (AMFlog = log (AMF), RGlog = log (RelGenSpec), DWR = log (DW_roots), DWA = log (DW_above))  
### these made the mdel terrible

# waht is the distribution of the data?
(hist <- ggplot(dataNL_sample, aes(x = DWRA )) +
    geom_histogram(bins = 40) +
    theme_classic())

plot (dataNL_sample$RelGenSpec, dataNL_sample$DWRA)

mod <- lm (AMF ~ RelGenSpec + DWRA , data = dataNL_sample)

plot (mod)





# right skewed for RelGenSpec, but left skewed for DW_roots- use diff family??

# 0A Specify  model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ RelGenSpec * DWRA + (1|gr(PlaSpe, cov = A)) ,
                       ### is only needed to account for diff between the species 
                      #+ OTHER than phylogenetic (environmental factors, niches)
                       family = skew_normal (),
                      data = dataNL_sample, data2 = list(A = A))

# model 0 
control = list(adapt_delta = 0.8) ## then rerun!!


Gen_samples_fit0 <- brm(
  AMF ~ RelGenSpec * DWRA + (1 |gr(PlaSpe, cov = A)),
  data = dataNL_sample,
  family = skew_normal(),
  data2 = list(A = A),
  prior = prior0,  chains = 4, cores = 8,
  iter = 2000, warmup = 500
)

##
Gen_samples_fit0 <- add_criterion(Gen_samples_fit0, c("loo", "waic"))


# model 0 
control = list(adapt_delta = 0.93) ## then rerun!!



##

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(Gen_samples_fit0)

# estimate of phylogenetic signal!!

# check convergence #
launch_shinystan(Gen_samples_fit0)

# posterior predictive checks #

pp_check (Gen_samples_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")

bayes_R2(Gen_samples_fit0)

# 1 Update to model 01 ####

#and then fit it again - adding more variables

Gen_samples_fit01 <- update(
  Gen_samples_fit0, formula = ~ . - RelGenSpec:DWRA,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit01 <- add_criterion(Gen_samples_fit01, c("loo", "waic"))


loo_compare (Gen_samples_fit0, Gen_samples_fit01)

Gen_samples_fit01B <- update(
  Gen_samples_fit01, formula = ~ . -(1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


#  check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit01)

# check convergence #
launch_shinystan(Gen_samples_fit01)

# posterior predictive checks #

pp_check (Gen_samples_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models ###
Gen_samples_fit0 <- add_criterion(Gen_samples_fit0, "loo")

Gen_samples_fit01 <- add_criterion(Gen_samples_fit01, "loo")
Gen_samples_fit01B <- add_criterion(Gen_samples_fit01B, "loo")

loo_compare (Gen_samples_fit0, Gen_samples_fit01, Gen_samples_fit01B)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

Gen_samples_fit02 <- update(
  Gen_samples_fit01B, formula = ~ .  -  DWRA 
)


#  check sampling quality of model 02 ###

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit02)

# check convergence #
launch_shinystan(Gen_samples_fit02)

# posterior predictive checks #

pp_check (Gen_samples_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models ###

Gen_samples_fit02 <- add_criterion(Gen_samples_fit02, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit01, Gen_samples_fit02)




# 3 Update to model 03 ####
Gen_samples_fit03 <- update(
  Gen_samples_fit0, formula = ~ .  -  RelGenSpec:DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit03B <- update(
  Gen_samples_fit03, formula = ~ .  -  RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
Gen_samples_fit03 <- add_criterion(Gen_samples_fit03, "loo")
Gen_samples_fit03B <- add_criterion(Gen_samples_fit03B, "loo")

loo_compare (Gen_samples_fit05Cskew,  Gen_samples_fit01, 
             Gen_samples_fit03B)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit03)

# 4 Update to model 04 ####

#and then fit it again - adding more variables

Gen_samples_fit04 <- update(
  Gen_samples_fit0, formula = ~ .  -  DW_roots - RelGenSpec ,
  newdata = dataNL_sample, chains = 6, cores = 8,
  iter = 8000, warmup = 4000
)

Gen_samples_fit04 <- add_criterion(Gen_samples_fit04, "loo")
loo_compare (Gen_samples_fit0,  Gen_samples_fit01, 
             Gen_samples_fit02, Gen_samples_fit03, 
             Gen_samples_fit04)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit04)
pp_check (Gen_samples_fit04, ndraws = 100)




# 5 Update to model 05 ####

#and then fit it again - adding more variables

Gen_samples_fit05 <- update(
  Gen_samples_fit04, formula = ~ .  - (1|PlantSpeciesfull),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


Gen_samples_fit05 <- add_criterion(Gen_samples_fit05, "loo")
loo_compare (Gen_samples_fit05C,  Gen_samples_fit05, Gen_samples_fit05skew, 
             Gen_samples_fit03)
pp_check (Gen_samples_fit05, ndraws = 50)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05)

# 5B Update to model 05B ####

#and then fit it again - adding more variables

Gen_samples_fit05B <- update(
  Gen_samples_fit05, formula = ~ . -  RelGenSpec*DW_roots + RelGenSpec,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit05B <- add_criterion(Gen_samples_fit05B, "loo")
loo_compare (Gen_samples_fit05B,  Gen_samples_fit05,fit2, 
             Gen_samples_fit03)
fit2 <- add_criterion(fit2, "loo")


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05B)


# 5c Update to model 05C ####

#and then fit it again - adding more variables

Gen_samples_fit05C <- update(
  Gen_samples_fit05B, formula = ~ .  + DW_roots,
  family = gaussian (),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit05Cskew <- update(
  Gen_samples_fit05C, family = skew_normal(),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
) 
Gen_samples_fit05Cskew <- add_criterion(Gen_samples_fit05Cskew, "loo")


Gen_samples_fit05C <- add_criterion(Gen_samples_fit05C, "loo")
loo_compare (Gen_samples_fit0,  Gen_samples_fit05, Gen_samples_fit05B, 
             Gen_samples_fit05C, Gen_samples_fit03)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05C)
launch_shinystan(Gen_samples_fit05C)
# 5c Update to model 05C ####
pp_check(Gen_samples_fit05C, ndraws = 100)
#and then fit it again - adding more variables

Gen_samples_fit05D <- update(
  Gen_samples_fit05B, formula = ~ .  + (1|PlantSpeciesfull),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit05E <- update(
  Gen_samples_fit05C, formula = ~ .  - RelGenSpec,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


Gen_samples_fit05D <- add_criterion(Gen_samples_fit05D, "loo")
Gen_samples_fit05E <- add_criterion(Gen_samples_fit05E, "loo")

Gen_samples_fit05B <- add_criterion(Gen_samples_fit05B, "waic")
Gen_samples_fit05D <- add_criterion(Gen_samples_fit05D, "waic")

Gen_samples_fit05C <- add_criterion(Gen_samples_fit05C, "waic")
loo_compare (Gen_samples_fit05D,  Gen_samples_fit05, Gen_samples_fit05B, 
             Gen_samples_fit05C, Gen_samples_fit05E)

loo_compare (Gen_samples_fit05D,  Gen_samples_fit05B, 
             Gen_samples_fit05C,criterion = "waic")


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05D)



# 6 Update to model 06 ####

#and then fit it again - adding more variables

Gen_samples_fit06 <- update(
  Gen_samples_fit05, formula = ~ .   + DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit06 <- add_criterion(Gen_samples_fit06, "loo")
loo_compare (Gen_samples_fit02, Gen_samples_fit05, Gen_samples_fit06, Gen_samples_fit05B)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit06)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

Gen_samples_fit07 <- update(
  Gen_samples_fit05B, formula = ~ .   + RelGenSpec:DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit07 <- add_criterion(Gen_samples_fit07, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit07, Gen_samples_fit06, Gen_samples_fit05B)

#  check sampling quality of model  ###

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit07)


# 8 Update to model 08 ####

#and then fit it again - adding more variables

Gen_samples_fit08 <- update(
  Gen_samples_fit05B, formula = ~ . -RelGenSpec + RelGenSpec:DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit08 <- add_criterion(Gen_samples_fit08, "loo")
loo_compare (Gen_samples_fit06, Gen_samples_fit04, Gen_samples_fit08, Gen_samples_fit05B, Gen_samples_fit07)
launch_shinystan(Gen_samples_fit08)


# 8 Update to model 08B ####

#and then fit it again - adding more variables

Gen_samples_fit08B <- update(
  Gen_samples_fit08, formula = ~ . +DW_roots + RelGenSpec:DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit08B <- add_criterion(Gen_samples_fit08B, "loo")
loo_compare (Gen_samples_fit04, Gen_samples_fit08, Gen_samples_fit05B, Gen_samples_fit08B)


launch_shinystan(Gen_samples_fit08B)


# 9 Update to model 09 ####

#and then fit it again - adding more variables

Gen_samples_fit09 <- update(
  Gen_samples_fit07, formula = ~ .  + DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit09 <- add_criterion(Gen_samples_fit09, "loo")
loo_compare (Gen_samples_fit04, Gen_samples_fit05B, Gen_samples_fit05C
             , Gen_samples_fit08, Gen_samples_fit05Cskew)

# 
# 10 Update to model 10 ####

#and then fit it again - adding more variables

Gen_samples_fit10 <- update(
  Gen_samples_fit08, formula = ~ .  + (1|PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit10 <- add_criterion(Gen_samples_fit10, "loo")
loo_compare (Gen_samples_fit08, Gen_samples_fit10, Gen_samples_fit05C, Gen_samples_fit05B)


# 11 Update to model 11 ####

Gen_samples_fit11 <- update(
  Gen_samples_fit10, formula = ~ .  + RelGenSpec + DW_above - RelGenSpec:DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit11 <- add_criterion(Gen_samples_fit11, "loo")
loo_compare (Gen_samples_fit05B, Gen_samples_fit10, Gen_samples_fit05C, Gen_samples_fit11)
# 
# 12 Update to model 12 ####
Gen_samples_fit12 <- update(
  Gen_samples_fit03, formula = ~ .  - RelGenSpec - DW_roots + RelGenSpec:DW_roots:DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit12 <- add_criterion(Gen_samples_fit12, "loo")
loo_compare (Gen_samples_fit05C, Gen_samples_fit10, Gen_samples_fit05B, Gen_samples_fit12)
# 
 Gen_samples_fit13 <- update(
   Gen_samples_fit0, formula = ~ .  - DW_roots - RelGenSpec + RelGenSpec:DW_above,
   newdata = dataNL_sample, chains = 4, cores = 8,
   iter = 5000, warmup = 2000
 )
# 
 Gen_samples_fit13 <- add_criterion(Gen_samples_fit13, "loo")
 loo_compare (Gen_samples_fit08, Gen_samples_fit13, Gen_samples_fit09, Gen_samples_fit12)
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

### Model  is the best ###

loo_compare (Gen_samples_fit0,  Gen_samples_fit01, 
             Gen_samples_fit02, Gen_samples_fit03, 
             Gen_samples_fit04, Gen_samples_fit08, 
             Gen_samples_fit05B, Gen_samples_fit05C, 
             Gen_samples_fit01B, Gen_samples_fit05, 
             Gen_samples_fit07, Gen_samples_fit09, 
             Gen_samples_fit10, Gen_samples_fit06)

# While model 5B only including RelGenSpec is the best based on loo croiterion, it is calculated with large amout of divergent transitions (222)
# The next best model is 05C RelGenSpec + DW_roots 

 tidy (Gen_samples_fit04, effects = c("fixed"))





### get posteriors from the model####

# 
dataNL_sample %>%
  #group_by(PlaSpe) %>%
  #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
  add_epred_draws(Gen_samples_fit02) %>%
  ggplot(aes(x = RelGenSpec, y = AMF, color = ordered(PlaSpe))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = dataNL_sample) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~PlaSpe, scales = "free_x")

# conditional effects ####
#conditions = make_conditions (Gen_samples_fit04, "RelGenSpec")


plot (conditional_effects(Gen_samples_fit01, effects = "DW_roots:RelGenSpec",  ndraws = 10000, spaghetti = F,  
                    prob = 0.9, conditions = conditions), points =T)
#plot (conditional_effects(Gen_samples_fit04, effects = "DW_roots:RelGenSpec",  ndraws = 10000, spaghetti = F,  
#                          prob = 0.9), points =T)

#conditions2 = make_conditions (Gen_samples_fit04, "DW_roots")


#plot (conditional_effects(Gen_samples_fit04, effects = "RelGenSpec",  ndraws = 10000, spaghetti = F,  prob = 0.9, conditions = conditions2), points = T)

 plot (conditional_effects(Gen_samples_fit01, effects = "RelGenSpec:DW_roots",  ndraws = 10000, 
                          spaghetti = F,  prob = 0.5
                          ), points = T)


#conditions3= data.frame (DW_roots = c(0.25, 0.5, 0.75))
conditions2 = make_conditions (Gen_samples_fit0, "DW_above")

plot (conditional_effects(Gen_samples_fit0, effects = "RelGenSpec:DW_roots",  ndraws = 10000, spaghetti = F,  
                          prob = 0.5, conditions = conditions2), points = T)

launch_shinystan(Gen_samples_fit05C )


#Run with better iterations
Gen_samples_fit04 <- update(
  Gen_samples_fit0, formula = ~ .  -  DW_roots - RelGenSpec ,
  newdata = dataNL_sample, chains = 8, cores = 8,
  iter = 10000, warmup = 4000
)

Gen_samples_fit04 <- add_criterion(Gen_samples_fit04, "loo")
