
# Part 9 - Bayesian linear models for soil bacterial biomass ##

# Packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library(GGally)
library (rstatix)
library (shinystan)


## Bayesian model for total bacterial biomass  ####
# The bacterial biomass was estimated using PLFA as a proxy. 
# Data cleanup
dataNLPL  <- 
  dataNLPL %>% 
  filter (!is.na (totalBact)) %>% 
  mutate (DWRA = DW_roots / DW_above) %>% 
  left_join (PCA_metric_data_sample) 



# #Run some tests in linear mixed models ####
# 
# PLLmer1  <- lmerTest::lmer (totalBact ~ AMF +Dim.1 + Dim.2 + Dim.3 + DW_above + DW_roots + (1|PlantSpeciesfull) , data = dataNLPL               )
# Anova (PLLmer1, type = 3)
# 
# summary (PLLmer1)


# 0A Specify the Nullmodel 0 ####
## get default  priors  #### 
# prior <- get_prior ( totalBact ~ Dim.1 + Dim.2 + Dim.3 + DW_above + DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull)
#                      #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
#                      #+ OTHEr than phylogenetic (environmental factors, niches)
#                      , data = dataNLPL, data2 = list(A = A))

# Fit the model with only the random structure
PLm0 <- brm(
  totalBact ~  1+ (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) ,
  data = dataNLPL,
  family = gaussian(),
  data2 = list(A = A),
    sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

#
PLm0B <- brm(
  totalBact ~  1+  (1|PlantSpeciesfull) ,
  data = dataNLPL,
  family = gaussian(),
  data2 = list(A = A),
  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

# Set up hypothesis #
hyp0 <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"
hyp1 <- "sd_PlantSpeciesfull__Intercept^2 / (sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"

hyp <- hypothesis(PLm0, hyp1, class = NULL)

hyp

PLm0B <-add_criterion(PLm0B, "loo")
loo (PLm0B)



#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(PLm0)


# check convergence #
launch_shinystan(PLm0)

# posterior predictive checks #

pp_check (PLm0, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

bayes_R2 (PLm0)


# 1 Update to model 01 ####

#and then fit it again - adding more variables

PLm01 <- brm(
  totalBact ~  Dim.1  ,
  data = dataNLPL,
  family = gaussian(),
 chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

summary (PLm01)
PLm01 %>%  bayes_R2()

# check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc


# check convergence #
launch_shinystan(PLm01)

# posterior predictive checks #

pp_check (PLm01, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###
PLm0 <- add_criterion(PLm0, "loo")

PLm01 <- add_criterion(PLm01, "loo")
loo_compare (PLm0, PLm01)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

PLm02 <- update(
  PLm01, formula = ~ . +Dim.2  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

# check sampling quality of model 02 ##

# Look for rhat, ESS, Sd etc
summary (PLm02)

# check convergence #
launch_shinystan(PLm02)
bayes_R2(PLm02)
# posterior predictive checks #

pp_check (PLm02, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

# Compare models ###

PLm02 <- add_criterion(PLm02, "loo")
loo_compare (PLm0, PLm02, PLm01)



# 3 Update to model 03 ####

#and then fit it again - adding more variables

PLm03 <- update(
  PLm01, formula = ~ .  +Dim.3 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
bayes_R2(PLm03)

# check sampling quality of model 03 ###

# Look for rhat, ESS, Sd etc
summary (PLm03)

# check convergence #
launch_shinystan(PLm03)

# posterior predictive checks #

pp_check (PLm03, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLm03 <- add_criterion(PLm03, "loo")
loo_compare (PLm02, PLm03, PLm01, PLm0)
### Model 02 is the best ###


# 4 Update to model 04 ####

#and then fit it again - adding more variables

PLm04 <- update(
  PLm01, formula = ~ .  + DW_roots  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm04)

PLm04 <- add_criterion(PLm04, "loo")

loo_compare(PLm04, PLm01)
bayes_R2 (PLm04)


# 5 Update to model 05 ####

#and then fit it again - adding more variables

PLm05 <- update(
  PLm01, formula = ~ .  - Dim.1 +Dim.1:DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

bayes_R2(PLm05)




# check sampling quality of model 05###

# Look for rhat, ESS, Sd etc
summary (PLm05)

# check convergence #
launch_shinystan(PLm05)

# posterior predictive checks #

pp_check (PLm05, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ##

PLm05 <- add_criterion(PLm05, "loo")
loo_compare (PLm02, PLm03, PLm05, PLm04)



# 6 Update to model 06 ####

#and then fit it again - adding more variables

PLm06 <- update(
  PLm05, formula = ~ .   +  DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# check sampling quality of model 06 ###

# Look for rhat, ESS, Sd etc
summary (PLm06)

# check convergence #
launch_shinystan(PLm06)

# posterior predictive checks #

pp_check (PLm06, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLm06 <- add_criterion(PLm06,"loo")
loo_compare (PLm02, PLm03, PLm06, PLm05)

bayes_R2 (PLm06)
#  7 Update to model 7 ####
#and then fit it again - adding more variables

PLm07 <- update(
  PLm01, formula = ~ .   -Dim.1 + Dim.1:DW_roots:DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm07)
bayes_R2(PLm07)

PLm07 <- add_criterion(PLm07, "loo")

loo_compare(PLm05, PLm02, PLm07, 
            PLm06)

#  8 Update to model 8 ####
#and then fit it again - adding more variables

PLm08 <- update(
  PLm01, formula = ~ .   +Dim.2 + Dim.3 + DW_above + DW_roots,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm08)


PLm08 <- add_criterion(PLm08, "loo")

loo_compare(PLm05, PLm02, PLm07, 
            PLm06, PLm08)

bayes_R2(PLm08)

#  9 Update to model 9 ####
#and then fit it again - adding more variables

PLm09 <- update(
  PLm01, formula = ~ .   +DW_roots:DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm09)


PLm09 <- add_criterion(PLm09, "loo")

loo_compare(PLm05, PLm02, PLm07, 
            PLm06, PLm09)
bayes_R2(PLm09)
#  10 Update to model 10 ####
#and then fit it again - adding more variables

PLFm10 <- update(
  PLm04, formula = ~ .   +DW_above  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFm10)


PLFm10 <- add_criterion(PLFm10, "loo")

loo_compare(PLm05, PLm02, PLm07, 
            PLm06, PLm09, PLFm10)

bayes_R2(PLFm10)

#  11 Update to model 11 ####
#and then fit it again - adding more variables

PLm11 <- update(
  PLm01, formula = ~ .   +DW_above  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PLm11 %>% bayes_R2()


PLm11 <- add_criterion(PLm11, "loo")

loo_compare(PLm05, PLFm10, PLm07, 
            PLm06, PLm09, PLm11)

#  12 Update to model 12 ####
#and then fit it again - adding more variables

PLm12 <- update(
  PLm02, formula = ~ .   - Dim.1 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_12)


PLm12 <- add_criterion(PLm12, "loo")

loo_compare(PLm12, PLm02, PLm07, 
            PLm06, PLm09)

bayes_R2 (PLm12)


#  13 Update to model 13 ####
#and then fit it again - adding more variables

PLm13 <- update(
  PLm06, formula = ~ .   - Dim.1 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm13)


PLm13 <- add_criterion(PLm13, "loo")

loo_compare(PLm12, PLm02, PLm05, 
            PLm06, PLm09, PLm13)

bayes_R2 (PLm13)



#  14 Update to model 14 ####
#and then fit it again - adding more variables

PLm14 <- update(
  PLm06, formula = ~ .   -DW_above + DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm14)


PLm14 <- add_criterion(PLm14, "loo")

loo_compare(PLm14, PLm02, PLm05, 
            PLm06, PLm09, PLm13)

bayes_R2 (PLm14)

#  15 Update to model 15 ####
#and then fit it again - adding more variables

PLm15 <- update(
  PLm02, formula = ~ .   - DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm15)


PLm15 <- add_criterion(PLm15, "loo")

loo_compare(PLm14, PLm02, PLm05, 
            PLm06, PLm09, PLm15)

bayes_R2 (PLm15)



#  16 Update to model 16 ####
#and then fit it again - adding more variables

PLm16 <- update(
  PLm15, formula = ~ .   - DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm16)


PLm16 <- add_criterion(PLm16, "loo")

loo_compare(PLm16, PLm02, PLm05, 
            PLm06, PLm09, PLm13)

bayes_R2 (PLm16)



#  17 Update to model 17 ####
#and then fit it again - adding more variables

PLm17 <- update(
  PLm15, formula = ~ .   - (1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm17)


PLm17 <- add_criterion(PLm17, "loo")

loo_compare(PLm12, PLm02, PLm05, 
            PLm06, PLm09, PLm17)

bayes_R2 (PLm17)

#  18 Update to model 18 ####
#and then fit it again - adding more variables

PLm18 <- update(
  PLm16, formula = ~ .   -  (1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLm18)


PLm18 <- add_criterion(PLm18, "loo")

loo_compare(PLm12, PLm02, PLm05, 
            PLm06, PLm09, PLm18)

bayes_R2 (PLm18)

# rerun best model #####
## get default  priors  #### 
prior06 <- get_prior ( totalBact ~ Dim.1 +DW_roots + DW_above 
                       , data = dataNLPL, data2 = list(A = A))

# model 06New
loo_compare (PLm0B, PLm11)

PLm10New <- brm(
  totalBact ~ Dim.1 + DW_roots  + DW_above  ,
  data = dataNLPL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior06,  sample_prior = TRUE, chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)
bayes_R2 (PLm10New)
#  NOTE: Additional models might have been run and discarded from the code if they did not converge #


# Figure 3 for paper  ####
conditions2 = make_conditions (PLm10New, "DW_above")

plot_2 <- conditional_effects(PLm10New, effects = "Dim.1:DW_roots",  ndraws = 50000, 
                              spaghetti = F,  prob = 0.5, conditions = conditions2)


label_DWA <- c("shoot biomass = 0.52 g","shoot biomass = 0.71 g","shoot biomass = 0.90 g" )

names (label_DWA) <- c ("DW_above = 0.52", "DW_above = 0.71", "DW_above = 0.9")
ggplot (plot_2$`Dim.1:DW_roots`, aes (x = Dim.1, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_bw ( ) + 
  geom_point (data = dataNLPL, aes (x= Dim.1, y = totalBact), inherit.aes = F)+
  theme (axis.title = element_text(size =16, color = "black" ), 
         axis.text = element_text(size = 14, color = "black"),
         legend.position = "right", 
         legend.text = element_text(size = 14, color = "black"), 
         legend.title = element_text(size = 16, color = "black"), 
         panel.grid = element_blank()) +
  facet_wrap(~cond__, labeller = labeller (cond__ = label_DWA)) +
  theme ( strip.text.x = element_text(size = 14),
         legend.position = "right") +
  xlab ( "Principal component 1") +
  ylab(expression("Bacterial PLFA in nmol/g soil")) +
  scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
  scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 

ggsave ( "figures/Bact_PC.png", height = 6, width =12, dpi = 300)    






# Appendix S2: Figure S8 pp check ####
pp_check (PLm10New, ndraws = 100) +
  xlab ("Bacterial biomass in soil")

# Appendix S2: Table S9 ####

launch_shinystan(PLm10New)                 
## Appendix S2: Table S6 #####
dataNL %>%  
  filter (sampleID != "R18")  %>%   ## plant lost shoot biomass as it was dead at harvest, removed from analysis
left_join ()
dataNLPL %>%  select (sampleID, PlantSpeciesfull, AMFroot, AMF, totalBact, DW_roots, DW_above) %>% 
  write_csv("results/tableS6.csv")

# Appendix ####
# Hypothesis test ####
h <- "Dim.1 <0"

 hypothesis (PLm10New, h)
 launch_shinystan(PLm10New)
# Full model hypothesis tests ####
# 
 PLmFULL <- brm(
   totalBact ~ Dim.1 + Dim.2 + Dim.3 + DW_roots  + DW_above  ,
   data = dataNLPL,
   family = gaussian(),
   data2 = list(A = A),
   prior = prior06,  sample_prior = TRUE, chains = 8, cores = 8,
   iter = 10000, warmup = 3000
 )
 
 h<- "Dim.3>0"
 hypothesis (PLmFULL, h)
 
# Richness 

