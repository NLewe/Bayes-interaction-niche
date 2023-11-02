
## Code for Lisa ####

# Part 8 Bayesian models for paper 
#Lewe, N, Deslippe, J "AMF biomass" ####

# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
# library(modelr)
# library (broom.mixed)
# library(GGally)
# library (lmerTest)
# library (rstatix)

# Tutorials can be found here: 

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# https://mc-stan.org/loo/articles/online-only/faq.html


# data  ####

dataNL_sample <- read_rds("dataNL_sample_LW.rds")
# covaraince structure (based on a phylogenetic tree of the used plant species) This was included to remove the effect of phylogeny - we had 8 plant species, from three families. 
#They have very different traits (we did not measure those traits), for example differences in roots or biochemical differences that might effect our dependant variables.
A <- read_rds ("covariance_str_tree_concat_plants.rds")

# PlaSpe and PlantSpeciesfull  categorical, they are the plant species names. BOTH are needed because the covariance structure (random factor) uses" PlaSpe" 
# whereas the plant species identity is coded by "PlantSpeciesfull"
# DW_above and DW_roots are the biomass (dry weights) of each plant individual's shoots (=aboveground biomass) and roots

# totalBact is the biomass of bacteria in the soil
# AMF is the  biomass of specifc fungi in the soil of each plant individual, 
# AMFroot is the biomass of those fungi in the  roots of each tested plant individual

# Dim.1 Dim.2 and Dim.3 are  principal componetns
# I determined/measured several plant traits for each plant individual. Because several of those plant traits are somewhat correlated, I ran Principal component analysis.
# I used the three PC as variables in the models (they describe more than 80% of the variation between the plant regarding the trats I measured)

### We hypothesised that the traits we measured had a positive effect (well, were positively related to) the biomass of the AMF (fungi) and bacteria in the soil.
## This code here is only for the fungal biomass #####



### BAYES ####

# I made notes on all these models by hand, comparing their convergence, the loo measures, bayes R2 and decided if the variables have any effect by checking 
#their errors and if the spread of the estimates passed Zero.

# 0A Specify  model 0####
## get default  priors  #### 
prior0 <- get_prior ( AMF ~ Dim.1 + Dim.2 +Dim.3 + DW_roots + DW_above + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
                      data = dataNL_sample, data2 = list(A = A))

# model 0 
#control = list(adapt_delta = 0.8) ## then rerun!!


PCA_axes_fit0 <- brm(
  AMF ~ Dim.1 + Dim.2 +Dim.3 + DW_roots + DW_above +     
    (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
  data = dataNL_sample,
  data2 = list(A = A),
  prior = prior0,  chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

##
PCA_axes_fit0 <- add_criterion(PCA_axes_fit0, c("loo", "waic"))





##

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries 
summary(PCA_axes_fit0)

bayes_R2 (PCA_axes_fit0)
# check convergence #
launch_shinystan(PCA_axes_fit0)

# posterior predictive checks #

pp_check (PCA_axes_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")
PCA_axes_fit0 <- add_criterion(PCA_axes_fit0, c("loo", "waic"))


# 1 Update to model 01 ####

#and then fit it again - adding more variables

PCA_axes_fit01 <- update(
  PCA_axes_fit0, formula = ~ . - Dim.1 - Dim.2 -DW_above,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit01 <- add_criterion(PCA_axes_fit01, c("loo", "waic"))

bayes_R2 (PCA_axes_fit01)


loo_compare (PCA_axes_fit0, PCA_axes_fit01)
#  check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit01)

# check convergence #
launch_shinystan(PCA_axes_fit01)

# posterior predictive checks #

pp_check (PCA_axes_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")



PCA_axes_fit01B <- update(
  PCA_axes_fit01, formula = ~ . -(1|gr(PlaSpe, cov = A))  ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)




# 3 Compare models ###
PCA_axes_fit0 <- add_criterion(PCA_axes_fit0, "loo")

PCA_axes_fit01 <- add_criterion(PCA_axes_fit01, "loo")
PCA_axes_fit01B <- add_criterion(PCA_axes_fit01B, "loo")

loo_compare (PCA_axes_fit0, PCA_axes_fit01, PCA_axes_fit01B)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

PCA_axes_fit02 <- update(
  PCA_axes_fit01B, formula = ~ .  - DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit02 <- add_criterion(PCA_axes_fit02, "loo")


#  check sampling quality of model 02 ###

# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit02)

bayes_R2 (PCA_axes_fit02)
# check convergence #
launch_shinystan(PCA_axes_fit02)

# posterior predictive checks #

pp_check (PCA_axes_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models ###

PCA_axes_fit02 <- add_criterion(PCA_axes_fit02, "loo")
loo_compare (PCA_axes_fit0, PCA_axes_fit01, PCA_axes_fit02)




# 3 Update to model 03 ####
PCA_axes_fit03 <- update(
  PCA_axes_fit0, formula = ~ .  -Dim.1 - Dim.2 - (1|PlantSpeciesfull)  ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit03 <- add_criterion(PCA_axes_fit03, "loo")


loo_compare (PCA_axes_fit02,  PCA_axes_fit01B, 
             PCA_axes_fit01, PCA_axes_fit03)


# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit03)

# 4 Update to model 04 ####

#and then fit it again - adding more variables

PCA_axes_fit04 <- update(
  PCA_axes_fit03, formula = ~ .  -  DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit04 <- add_criterion(PCA_axes_fit04, "loo")
loo_compare (PCA_axes_fit0,  PCA_axes_fit01, 
             PCA_axes_fit02, PCA_axes_fit03, 
             PCA_axes_fit04)


# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit04)
pp_check (PCA_axes_fit04, ndraws = 100)


bayes_R2 (PCA_axes_fit04)

# 5 Update to model 05 ####

#and then fit it again - adding more variables

PCA_axes_fit05 <- update(
  PCA_axes_fit01, formula = ~ .  - Dim.3 - DW_roots + Dim.3:DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


PCA_axes_fit05 <- add_criterion(PCA_axes_fit05, "loo")
loo_compare (PCA_axes_fit05,  PCA_axes_fit01, PCA_axes_fit02, 
             PCA_axes_fit04)
pp_check (PCA_axes_fit05, ndraws = 50)
bayes_R2 (PCA_axes_fit05)


# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit05)

# 5B Update to model 05B ####

#and then fit it again - adding more variables

PCA_axes_fit05B <- update(
  PCA_axes_fit01, formula = ~ . -Dim.3,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit05B <- add_criterion(PCA_axes_fit05B, "loo")
loo_compare (PCA_axes_fit05B,  PCA_axes_fit05, 
             PCA_axes_fit03)


# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit05B)
bayes_R2 (PCA_axes_fit05B)




# 6 Update to model 06 ####

#and then fit it again - adding more variables

PCA_axes_fit06 <- update(
  PCA_axes_fit01B, formula = ~ .   - (1|PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit06 <- add_criterion(PCA_axes_fit06, "loo")
loo_compare (PCA_axes_fit02, PCA_axes_fit05, PCA_axes_fit06, PCA_axes_fit05B)
bayes_R2 (PCA_axes_fit06)


# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit06)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

PCA_axes_fit07 <- update(
  PCA_axes_fit01, formula = ~ .   -DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit07 <- add_criterion(PCA_axes_fit07, c("loo", "waic"))
loo_compare (PCA_axes_fit0, PCA_axes_fit07, PCA_axes_fit01, PCA_axes_fit05)

#  check sampling quality of model  ###

# Look for rhat, ESS, Sd etc
summary (PCA_axes_fit07)
bayes_R2 (PCA_axes_fit07)


# 8 Update to model 08 ####

#and then fit it again - adding more variables

PCA_axes_fit08 <- update(
  PCA_axes_fit02, formula = ~ . -Dim.3 + Dim.3:DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit08 <- add_criterion(PCA_axes_fit08, "loo")
loo_compare (PCA_axes_fit06, PCA_axes_fit0, PCA_axes_fit02, PCA_axes_fit03, PCA_axes_fit05, PCA_axes_fit04, PCA_axes_fit08, PCA_axes_fit05B, PCA_axes_fit07, PCA_axes_fit01
)
launch_shinystan(PCA_axes_fit08)


# 8 Update to model 08B ####

#and then fit it again - adding more variables

PCA_axes_fit08B <- update(
  PCA_axes_fit08, formula = ~ . +DW_roots + RelGenSpec:DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit08B <- add_criterion(PCA_axes_fit08B, "loo")
loo_compare (PCA_axes_fit04, PCA_axes_fit08, PCA_axes_fit05B, PCA_axes_fit08B)


launch_shinystan(PCA_axes_fit08B)


# 9 Update to model 09 ####

#and then fit it again - adding more variables

PCA_axes_fit09 <- update(
  PCA_axes_fit06, formula = ~ .  - DW_roots - Dim.3 + Dim.3:DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit09 <- add_criterion(PCA_axes_fit09, "loo")
loo_compare (PCA_axes_fit06, PCA_axes_fit09, PCA_axes_fit05
             , PCA_axes_fit08, PCA_axes_fit04)


PCA_axes_fit09 %>%  bayes_R2


# 10 Update to model 10 ####

#and then fit it again - adding more variables

PCA_axes_fit10 <- update(
  PCA_axes_fit02, formula = ~ .  + totRootAMF ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit10 <- add_criterion(PCA_axes_fit10, "loo")
loo_compare (PCA_axes_fit08, PCA_axes_fit10, PCA_axes_fit05, PCA_axes_fit06)


# 11 Update to model 11 ####
# 
PCA_axes_fit11 <- update(
  PCA_axes_fit06, formula = ~ .  -DW_roots + totRootAMF,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PCA_axes_fit11 <- add_criterion(PCA_axes_fit11, "loo")
loo_compare (PCA_axes_fit05, PCA_axes_fit08, PCA_axes_fit06, PCA_axes_fit11)


bayes_R2 (PCA_axes_fit11)


###### A lot of more tests like that..... BEST model is at bottom ######
# # 
# # 12 Update to model 12 ####
# PCA_axes_fit12 <- update(
#   PCA_axes_fit03, formula = ~ .  - RelGenSpec - DW_roots + RelGenSpec:DW_roots:DW_above ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# PCA_axes_fit12 <- add_criterion(PCA_axes_fit12, "loo")
# loo_compare (PCA_axes_fit05C, PCA_axes_fit10, PCA_axes_fit05B, PCA_axes_fit12)
# # 
# PCA_axes_fit13 <- update(
#   PCA_axes_fit0, formula = ~ .  - DW_roots - RelGenSpec + RelGenSpec:DW_above,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# # 
# PCA_axes_fit13 <- add_criterion(PCA_axes_fit13, "loo")
# loo_compare (PCA_axes_fit08, PCA_axes_fit13, PCA_axes_fit09, PCA_axes_fit12)
# # 
# 
# 
# PCA_axes_fit14 <- update(
#   PCA_axes_fit12, formula = ~ .  + RelGenSpec ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# # 
# # PCA_axes_fit14 <- add_criterion(PCA_axes_fit14, "loo")
# # loo_compare (PCA_axes_fit14, PCA_axes_fit10, PCA_axes_fit12, PCA_axes_fit13)
# # 
# # PCA_axes_fit15 <- update(
# #   PCA_axes_fit14, formula = ~ .  - RelGenSpec:DW_roots ,
# #   newdata = dataNL_sample, chains = 4, cores = 8,
# #   iter = 5000, warmup = 2000
# # )
# # 
# # PCA_axes_fit15 <- add_criterion(PCA_axes_fit15, "loo")
# # loo_compare (PCA_axes_fit13, PCA_axes_fit12, PCA_axes_fit14, PCA_axes_fit15)
# # 
# # 
# # 
# # PCA_axes_fit16 <- update(
# #   PCA_axes_fit12, formula = ~ .  - DW_roots:RelGenSpec + RelGenSpec:DW_roots ,
# #   newdata = dataNL_sample, chains = 4, cores = 8,
# #   iter = 5000, warmup = 2000
# # )
# # 
# # PCA_axes_fit16 <- add_criterion(PCA_axes_fit16, "loo")
# # loo_compare (PCA_axes_fit16, PCA_axes_fit15, PCA_axes_fit0, PCA_axes_fit12)
# # 
# # 
# # 
# # 
# 
# # ### Model  is the best ###
# # 
# # loo_compare (PCA_axes_fit0,  PCA_axes_fit01, 
# #              PCA_axes_fit02, PCA_axes_fit03, 
# #              PCA_axes_fit04, PCA_axes_fit08, 
# #              PCA_axes_fit05B, PCA_axes_fit05C, 
# #              PCA_axes_fit01B, PCA_axes_fit05, 
# #              PCA_axes_fit07, PCA_axes_fit09, 
# #              PCA_axes_fit10, PCA_axes_fit06)
# 
# # While model 5B only including RelGenSpec is the best based on loo criterion, it is calculated with large amount of divergent transitions (222)
# # The next best model is 05C RelGenSpec + DW_roots 
# 
# tidy (PCA_axes_fit06, effects = c("fixed"))
# 
# 
# 
# 
# 
# ### get posteriors from the model####
# 
# # 
# dataNL_sample %>%
#   #group_by(PlaSpe) %>%
#   #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
#   add_epred_draws(PCA_axes_fit06) %>%
#   ggplot(aes(x = Dim.3, y = AMF, color = ordered(PlaSpe))) +
#   stat_lineribbon(aes(y = .epred)) +
#   geom_point(data = dataNL_sample) +
#   scale_fill_brewer(palette = "Greys") +
#   scale_color_brewer(palette = "Set2") #+
# facet_wrap(~PlaSpe, scales = "free_x")
# 
# # conditional effects ####
# #conditions = make_conditions (PCA_axes_fit04, "RelGenSpec")
# 
# 
# plot (conditional_effects(PCA_axes_fit06, effects = "DW_roots",  ndraws = 10000, spaghetti = F,  
#                           prob = 0.9), points =T)
# 
# #conditions2 = make_conditions (PCA_axes_fit04, "DW_roots")
# 
# 
# #plot (conditional_effects(PCA_axes_fit04, effects = "RelGenSpec",  ndraws = 10000, spaghetti = F,  prob = 0.9, conditions = conditions2), points = T)
# 
# plot (conditional_effects(PCA_axes_fit01, effects = "Dim.3",  ndraws = 10000, 
#                           spaghetti = F,  prob = 0.5
# ), points = T)
# 
# 
# 
# plot (conditional_effects(PCA_axes_fit06, effects = "Dim.3:DW_roots",  ndraws = 10000, spaghetti = F,  
#                           prob = 0.9), points = T)
# 
# launch_shinystan(PCA_axes_fit06 )
# 



### BEST MODEL and FIGURE ######

## model 
priorAMFNew <- get_prior (   AMF  ~  Dim.3 +  DW_roots ,
                             ### is only needed to account for diff between the species 
                             #+ OTHER than phylogenetic (environmental factors, niches)
                             family = gaussian (),
                             data = dataNL_sample, data2 = list(A = A))



## Rerun with more chains etc##
fitAMF06New <- brm(
  AMF  ~ Dim.3 + DW_roots,
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = priorAMFNew, 
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)


fitAMF06New <- add_criterion(fitAMF06New, "loo")

loo (fitAMF06New)
fitAMF06New %>%  bayes_R2 ()

## plot conditional effects model fit2 ####

plot_1 <- conditional_effects(fitAMF06New, effects = "Dim.3:DW_roots",  ndraws = 50000, spaghetti = F,  prob = 0.5)


##make  ggplot #####

Plot_RG_AMF <- 
  ggplot (plot_1$`Dim.3:DW_roots`, aes (x = Dim.3, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_bw ( ) + 
  theme (axis.title = element_text(size =12), 
         legend.position = "top") +
  xlab ( "PCA axis 3") +
  ylab(expression("NL 16:1"*omega*"5 in nmol/g soil")) +
  scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
  scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 

launch_shinystan(fitAMF06New)




