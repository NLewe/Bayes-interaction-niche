
# Part 8 Bayesian models for interaction generalism ####

# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library(GGally)

# Tutorials can be found here: 

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# https://mc-stan.org/loo/articles/online-only/faq.html


# data upload # 

dataNL_sample <- dataNL %>%  
  left_join(RelGen_E1_E2_sample %>%  dplyr::select (sampleID, RelGenSpec  )) %>% 
  filter (!is.na (AMF),  !is.na (RelGenSpec)) %>% 
  filter (sampleID != "R18")  %>%   ## plant lost shoot biomass as it was almost dead at harvest, removed from analysis
  dplyr::select (!starts_with ("Dim")) %>% 
  mutate (DW_roots = DW_roots /DW_above )

# AMF is the AMF biomass in the soil, calculated using the fatty acid biomarker for AMF
# AMFroot is the AMF biomass in the plant roots


# some data checks #
# what is the distribution of the data?
(hist <- ggplot(dataNL_sample, aes(x = DW_roots )) +
    geom_histogram(bins = 40) +
    theme_classic())

plot (dataNL_sample$RelGenSpec, dataNL_sample$DW_roots)

mod <- lm (AMF ~ RelGenSpec + DW_roots , data = dataNL_sample)

plot (mod)


# Note that all combinations of the variables have been run and the results compared.
# Model selection was done by different criteria (see publication)

# 0A Specify  model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ RelGenSpec * DW_roots + (1|gr(PlaSpe, cov = A)) ,
                      data = dataNL_sample, data2 = list(A = A))

# model 0 
#control = list(adapt_delta = 0.8) ## then rerun!!


Gen_samples_fit0 <- brm(
  AMF ~ RelGenSpec * DW_roots + (1 |gr(PlaSpe, cov = A)),
  data = dataNL_sample,
  data2 = list(A = A),
  prior = prior0,  chains = 4, cores = 8,
  iter = 2000, warmup = 500
)

##
Gen_samples_fit0 <- add_criterion(Gen_samples_fit0, c("loo", "waic"))


# model 0 
#control = list(adapt_delta = 0.93) ## then rerun!!



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
  Gen_samples_fit0, formula = ~ . - RelGenSpec:DW_roots,
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
bayes_R2 (Gen_samples_fit01)
bayes_R2 (Gen_samples_fit01B)

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
  Gen_samples_fit01B, formula = ~ .  -  DW_roots 
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

bayes_R2 (Gen_samples_fit02)



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

loo_compare (Gen_samples_fit03,  Gen_samples_fit01, 
             Gen_samples_fit03B)


bayes_R2(Gen_samples_fit03)
bayes_R2(Gen_samples_fit03B)

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit03)

# 4 Update to model 04 ####

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
pp_check (Gen_samples_fit04, ndraws = 100)

bayes_R2(Gen_samples_fit04)



# 5 Update to model 05 ####

#and then fit it again - adding more variables

Gen_samples_fit05 <- update(
  Gen_samples_fit01, formula = ~ .  + (1|PlantSpeciesfull),
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


Gen_samples_fit05 <- add_criterion(Gen_samples_fit05, "loo")
loo_compare (Gen_samples_fit05C,  Gen_samples_fit05, Gen_samples_fit05skew, 
             Gen_samples_fit03)
pp_check (Gen_samples_fit05, ndraws = 50)


# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit05)
bayes_R2(Gen_samples_fit05)





# 6 Update to model 06 ####

#and then fit it again - adding more variables

Gen_samples_fit06 <- update(
  Gen_samples_fit05, formula = ~ .   + DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit06 <- add_criterion(Gen_samples_fit06, "loo")
loo_compare (Gen_samples_fit02, Gen_samples_fit05, Gen_samples_fit06, Gen_samples_fit05B)

bayes_R2(Gen_samples_fit06)

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit06)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

Gen_samples_fit07 <- update(
  Gen_samples_fit05, formula = ~ .   + RelGenSpec:DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit07 <- add_criterion(Gen_samples_fit07, "loo")
loo_compare (Gen_samples_fit0, Gen_samples_fit07, Gen_samples_fit06, Gen_samples_fit05B)
bayes_R2(Gen_samples_fit07)

#  check sampling quality of model  ###

# Look for rhat, ESS, Sd etc
summary (Gen_samples_fit07)

bayes_R2(Gen_samples_fit07)

# 8 Update to model 08 ####

#and then fit it again - adding more variables

Gen_samples_fit08 <- update(
  Gen_samples_fit04, formula = ~ . + (1|PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit08 <- add_criterion(Gen_samples_fit08, "loo")
loo_compare (Gen_samples_fit06, Gen_samples_fit04, Gen_samples_fit08, Gen_samples_fit05B, Gen_samples_fit07)
launch_shinystan(Gen_samples_fit08)

bayes_R2(Gen_samples_fit08)




# 9 Update to model 09 ####

#and then fit it again - adding more variables

Gen_samples_fit09 <- update(
  Gen_samples_fit01, formula = ~ .  + DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit09 <- add_criterion(Gen_samples_fit09, "loo")
loo_compare (Gen_samples_fit04, Gen_samples_fit0, Gen_samples_fit05
             , Gen_samples_fit08, Gen_samples_fit02)

bayes_R2(Gen_samples_fit09)

# 10 Update to model 10 ####

#and then fit it again - adding more variables

Gen_samples_fit10 <- update(
  Gen_samples_fit05, formula = ~ .  -RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit10 <- add_criterion(Gen_samples_fit10, "loo")
loo_compare (Gen_samples_fit08, Gen_samples_fit10, Gen_samples_fit05C, Gen_samples_fit05B)

bayes_R2(Gen_samples_fit10)


# 11 Update to model 11 ####

Gen_samples_fit11 <- update(
  Gen_samples_fit04, formula = ~ .   + DW_above + (1|PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit11 <- add_criterion(Gen_samples_fit11, "loo")
loo_compare (Gen_samples_fit05, Gen_samples_fit10, Gen_samples_fit04, Gen_samples_fit11, Gen_samples_fit03B, Gen_samples_fit0, Gen_samples_fit03)
bayes_R2(Gen_samples_fit11)


# 12 Update to model 12 ####
Gen_samples_fit12 <- update(
  Gen_samples_fit03, formula = ~ .  -(1|gr(PlaSpe, cov = A)) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit12 <- add_criterion(Gen_samples_fit12, "loo")
loo_compare (Gen_samples_fit05C, Gen_samples_fit10, Gen_samples_fit05B, Gen_samples_fit12)
# 
Gen_samples_fit13 <- update(
  Gen_samples_fit12, formula = ~ .   - RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
# 
Gen_samples_fit13 <- add_criterion(Gen_samples_fit13, "loo")
loo_compare (Gen_samples_fit04, Gen_samples_fit13, Gen_samples_fit03B, Gen_samples_fit12)

bayes_R2 (Gen_samples_fit13)
bayes_R2 (Gen_samples_fit12)


# Model 14 ####

Gen_samples_fit14 <- update(
  Gen_samples_fit04, formula = ~ .   -(1|gr(PlaSpe, cov = A)) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

Gen_samples_fit14 <- add_criterion(Gen_samples_fit14, "loo")
loo_compare (Gen_samples_fit14, Gen_samples_fit10, Gen_samples_fit04, Gen_samples_fit13)

bayes_R2 (Gen_samples_fit14)



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

# While model 5B only including RelGenSpec is the best based on loo criterion, it is calculated with large amount of divergent transitions (222)
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


plot (conditional_effects(Gen_samples_fit04, effects = "DW_roots:RelGenSpec",  ndraws = 10000, spaghetti = F,  
                          prob = 0.9, conditions = conditions), points =T)
#plot (conditional_effects(Gen_samples_fit04, effects = "DW_roots:RelGenSpec",  ndraws = 10000, spaghetti = F,  
#                          prob = 0.9), points =T)

#conditions2 = make_conditions (Gen_samples_fit04, "DW_roots")


#plot (conditional_effects(Gen_samples_fit04, effects = "RelGenSpec",  ndraws = 10000, spaghetti = F,  prob = 0.9, conditions = conditions2), points = T)

plot (conditional_effects(Gen_samples_fit04, effects = "RelGenSpec:DW_roots",  ndraws = 10000, 
                          spaghetti = F,  prob = 0.5
), points = T)


#conditions3= data.frame (DW_roots = c(0.25, 0.5, 0.75))
conditions2 = make_conditions (Gen_samples_fit0, "DW_above")

plot (conditional_effects(Gen_samples_fit0, effects = "RelGenSpec:DW_roots",  ndraws = 10000, spaghetti = F,  
                          prob = 0.5, conditions = conditions2), points = T)

launch_shinystan(Gen_samples_fit05C )

##
#Run with better iterations
Gen_samples_fit04 <- update(
  Gen_samples_fit0, formula = ~ .  -  DW_roots - RelGenSpec ,
  newdata = dataNL_sample, chains = 8, cores = 8,
  iter = 10000, warmup = 4000
)

Gen_samples_fit04 <- add_criterion(Gen_samples_fit04, "loo")



### BEST MODEL and FIGURE ######

## model 2 additional
prior2 <- get_prior (   AMF  ~ RelGenSpec: DW_roots + (1|gr(PlaSpe, cov = A)),
                        ### is only needed to account for diff between the species 
                        #+ OTHER than phylogenetic (environmental factors, niches)
                        data = dataNL_sample, data2 = list(A = A))

fit2 <- brm(
  AMF  ~ RelGenSpec: DW_roots + (1|gr(PlaSpe, cov = A)),
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior2, 
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)


fit2 <- add_criterion(fit2, "loo")




## plot conditional effects model fit2 ####

plot_1 <- conditional_effects(fit2, effects = "RelGenSpec:DW_roots",  ndraws = 10000, spaghetti = F,  prob = 0.5)


##make  ggplot #####

Plot_RG_AMF <- 
  ggplot (plot_1$`RelGenSpec:DW_roots`, aes (x = RelGenSpec, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_bw ( ) + 
  theme (axis.title = element_text(size =12), 
         legend.position = "top") +
  xlab ( "Relative interaction generalism") +
  ylab(expression("NL 16:1"*omega*"5 in nmol/g soil")) +
  scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
  scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 

launch_shinystan(fit2)


