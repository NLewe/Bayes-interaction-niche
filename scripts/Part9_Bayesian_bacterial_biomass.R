
# Part 9 ##

## Bayesian model for total bacterial biomass  ####
# The bacterial biomass was estimated using PLFA as a proxy. 

dataNLPL  <- 
  dataNLPL %>% 
  filter (!is.na (totalPLFA)) %>% 
  mutate (DWRA = DW_roots / DW_above) %>% 
  left_join (PCA_metric_data_sample) 



#Run some tests in linear mixed models ####

PLLmer1  <- lmerTest::lmer (totalBact ~ AMF +Dim.1 + Dim.2 + Dim.3 + DW_above + DW_roots + (1|PlantSpeciesfull) , data = dataNLPL               )
Anova (PLLmer1, type = 3)

summary (PLLmer1)


# 0A Specify the  model 0 ####
## get default  priors  #### 
prior <- get_prior ( totalBact ~ Dim.1 + Dim.2 + Dim.3 + DW_above + DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull)
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     , data = dataNLPL, data2 = list(A = A))

# model 0 


PLFAtot_PCA_0 <- brm(
  totalBact ~  Dim.1 + Dim.2 + Dim.3 + DW_above + DW_roots  + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) ,
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
summary(PLFAtot_PCA_0)


# check convergence #
launch_shinystan(PLFAtot_PCA_0)

# posterior predictive checks #

pp_check (PLFAtot_PCA_0, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

bayes_R2 (PLFAtot_PCA_0)


# 1 Update to model 01 ####

#and then fit it again - adding more variables

PLFAtot_PCA_01 <- update(
  PLFAtot_PCA_0, formula = ~ . - Dim.2 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_01)

PLFAtot_PCA_01 %>%  bayes_R2()

# check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc


# check convergence #
launch_shinystan(PLFAtot_PCA_01)

# posterior predictive checks #

pp_check (PLFAtot_PCA_01, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###
PLFAtot_PCA_0 <- add_criterion(PLFAtot_PCA_0, "loo")

PLFAtot_PCA_01 <- add_criterion(PLFAtot_PCA_01, "loo")
loo_compare (PLFAtot_PCA_0, PLFAtot_PCA_01)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

PLFAtot_PCA_02 <- update(
  PLFAtot_PCA_01, formula = ~ . -Dim.3  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

# check sampling quality of model 02 ##

# Look for rhat, ESS, Sd etc
summary (PLFAtot_PCA_02)

# check convergence #
launch_shinystan(PLFAtot_PCA_02)

# posterior predictive checks #

pp_check (PLFAtot_PCA_02, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

# Compare models ###

PLFAtot_PCA_02 <- add_criterion(PLFAtot_PCA_02, "loo")
loo_compare (PLFAtot_PCA_0, PLFAtot_PCA_02, PLFAtot_PCA_01)



# 3 Update to model 03 ####

#and then fit it again - adding more variables

PLFAtot_PCA_03 <- update(
  PLFAtot_PCA_01, formula = ~ .  - (1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
bayes_R2(PLFAtot_PCA_03)

# check sampling quality of model 03 ###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_PCA_03)

# check convergence #
launch_shinystan(PLFAtot_PCA_03)

# posterior predictive checks #

pp_check (PLFAtot_PCA_03, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLFAtot_PCA_03 <- add_criterion(PLFAtot_PCA_03, "loo")
loo_compare (PLFAtot_PCA_02, PLFAtot_PCA_03, PLFAtot_PCA_01, PLFAtot_PCA_0)
### Model 02 is the best ###


# 4 Update to model 04 ####

#and then fit it again - adding more variables

PLFAtot_PCA_04 <- update(
  PLFAtot_PCA_03, formula = ~ .  - Dim.3 - DW_roots  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_04)

PLFAtot_PCA_04 <- add_criterion(PLFAtot_PCA_04, "loo")

loo_compare(PLFAtot_PCA_04, PLFAtot_PCA_0)
bayes_R2 (PLFAtot_PCA_04)


# 5 Update to model 05 ####

#and then fit it again - adding more variables

PLFAtot_PCA_05 <- update(
  PLFAtot_PCA_04, formula = ~ .  - (1|gr(PlaSpe, cov = A))  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

bayes_R2(PLFAtot_PCA_05)




# check sampling quality of model 05###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_PCA_05)

# check convergence #
launch_shinystan(PLFAtot_PCA_05)

# posterior predictive checks #

pp_check (PLFAtot_PCA_05, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ##

PLFAtot_PCA_05 <- add_criterion(PLFAtot_PCA_05, "loo")
loo_compare (PLFAtot_PCA_02, PLFAtot_PCA_03, PLFAtot_PCA_05, PLFAtot_PCA_0)



# 6 Update to model 06 ####

#and then fit it again - adding more variables

PLFAtot_PCA_06 <- update(
  PLFAtot_PCA_05, formula = ~ .   +  DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# check sampling quality of model 06 ###

# Look for rhat, ESS, Sd etc
summary (PLFAtot_PCA_06)

# check convergence #
launch_shinystan(PLFAtot_PCA_06)

# posterior predictive checks #

pp_check (PLFAtot_PCA_06, ndraws= 100) +
  xlab ("totalPLFA biomass in soil")

#  Compare models ###

PLFAtot_PCA_06 <- add_criterion(PLFAtot_PCA_06,"loo")
loo_compare (PLFAtot_PCA_02, PLFAtot_PCA_03, PLFAtot_PCA_06, PLFAtot_PCA_05)

bayes_R2 (PLFAtot_PCA_06)
#  7 Update to model 7 ####
#and then fit it again - adding more variables

PLFAtot_PCA_07 <- update(
  PLFAtot_PCA_06, formula = ~ .   -Dim.1 - DW_roots + Dim.1:DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_07)


PLFAtot_PCA_07 <- add_criterion(PLFAtot_PCA_07, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06)

#  8 Update to model 8 ####
#and then fit it again - adding more variables

PLFAtot_PCA_08 <- update(
  PLFAtot_PCA_05, formula = ~ .   - DW_above + DW_roots:DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_08)


PLFAtot_PCA_08 <- add_criterion(PLFAtot_PCA_08, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06, PLFAtot_PCA_08)

bayes_R2(PLFAtot_PCA_08)

#  9 Update to model 9 ####
#and then fit it again - adding more variables

PLFAtot_PCA_09 <- update(
  PLFAtot_PCA_03, formula = ~ .   -Dim.3 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_09)


PLFAtot_PCA_09 <- add_criterion(PLFAtot_PCA_09, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06, PLFAtot_PCA_09)
bayes_R2(PLFAtot_PCA_09)
#  9 Update to model 9B ####
#and then fit it again - adding more variables

PLFAtot_PCA_09B <- update(
  PLFAtot_PCA_09, formula = ~ .    + (1|gr(PlaSpe, cov = A)) -(1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_09B)


PLFAtot_PCA_09B <- add_criterion(PLFAtot_PCA_09B, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_09B, 
            PLFAtot_PCA_06, PLFAtot_PCA_09)

bayes_R2(PLFAtot_PCA_09B)

#  10 Update to model 10 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_10 <- update(
  PLFAtot_PCA_02, formula = ~ .   -DW_roots  ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_10)


PLFAtot_relGfit_10 <- add_criterion(PLFAtot_relGfit_10, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_relGfit_10)


#  11 Update to model 11 ####
#and then fit it again - adding more variables

PLFAtot_relGfit_11 <- update(
  PLFAtot_PCA_06, formula = ~ .   - DW_above - Dim.1 + Dim.1:DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

PLFAtot_relGfit_11 %>% bayes_R2()


PLFAtot_relGfit_11 <- add_criterion(PLFAtot_relGfit_11, "loo")

loo_compare(PLFAtot_PCA_05, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_relGfit_11)

#  12 Update to model 12 ####
#and then fit it again - adding more variables

PLFAtot_PCA_12 <- update(
  PLFAtot_PCA_02, formula = ~ .   - Dim.1 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_relGfit_12)


PLFAtot_PCA_12 <- add_criterion(PLFAtot_PCA_12, "loo")

loo_compare(PLFAtot_PCA_12, PLFAtot_PCA_02, PLFAtot_PCA_07, 
            PLFAtot_PCA_06, PLFAtot_PCA_09)

bayes_R2 (PLFAtot_PCA_12)


#  13 Update to model 13 ####
#and then fit it again - adding more variables

PLFAtot_PCA_13 <- update(
  PLFAtot_PCA_06, formula = ~ .   - Dim.1 ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_13)


PLFAtot_PCA_13 <- add_criterion(PLFAtot_PCA_13, "loo")

loo_compare(PLFAtot_PCA_12, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_13)

bayes_R2 (PLFAtot_PCA_13)



#  14 Update to model 14 ####
#and then fit it again - adding more variables

PLFAtot_PCA_14 <- update(
  PLFAtot_PCA_06, formula = ~ .   -DW_above + DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_14)


PLFAtot_PCA_14 <- add_criterion(PLFAtot_PCA_14, "loo")

loo_compare(PLFAtot_PCA_14, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_13)

bayes_R2 (PLFAtot_PCA_14)

#  15 Update to model 15 ####
#and then fit it again - adding more variables

PLFAtot_PCA_15 <- update(
  PLFAtot_PCA_02, formula = ~ .   - DW_above ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_15)


PLFAtot_PCA_15 <- add_criterion(PLFAtot_PCA_15, "loo")

loo_compare(PLFAtot_PCA_14, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_15)

bayes_R2 (PLFAtot_PCA_15)



#  16 Update to model 16 ####
#and then fit it again - adding more variables

PLFAtot_PCA_16 <- update(
  PLFAtot_PCA_15, formula = ~ .   - DW_roots ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_16)


PLFAtot_PCA_16 <- add_criterion(PLFAtot_PCA_16, "loo")

loo_compare(PLFAtot_PCA_16, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_13)

bayes_R2 (PLFAtot_PCA_16)



#  17 Update to model 17 ####
#and then fit it again - adding more variables

PLFAtot_PCA_17 <- update(
  PLFAtot_PCA_15, formula = ~ .   - (1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_17)


PLFAtot_PCA_17 <- add_criterion(PLFAtot_PCA_17, "loo")

loo_compare(PLFAtot_PCA_12, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_17)

bayes_R2 (PLFAtot_PCA_17)

#  18 Update to model 18 ####
#and then fit it again - adding more variables

PLFAtot_PCA_18 <- update(
  PLFAtot_PCA_16, formula = ~ .   -  (1|PlantSpeciesfull) ,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (PLFAtot_PCA_18)


PLFAtot_PCA_18 <- add_criterion(PLFAtot_PCA_18, "loo")

loo_compare(PLFAtot_PCA_12, PLFAtot_PCA_02, PLFAtot_PCA_05, 
            PLFAtot_PCA_06, PLFAtot_PCA_09, PLFAtot_PCA_18)

bayes_R2 (PLFAtot_PCA_18)

# rerun best model #####
## get default  priors  #### 
prior06 <- get_prior ( totalBact ~ Dim.1 +DW_roots + DW_above 
                       , data = dataNLPL, data2 = list(A = A))

# model 06New


PLFAtot_PCA_06New <- brm(
  totalBact ~ Dim.1 + DW_roots  + DW_above  ,
  data = dataNLPL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior06,  sample_prior = TRUE, chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)




# Figure 3 for paper  ####
conditions2 = make_conditions (PLFAtot_PCA_06New, "DW_above")

plot_2 <- conditional_effects(PLFAtot_PCA_06New, effects = "Dim.1:DW_roots",  ndraws = 50000, 
                              spaghetti = F,  prob = 0.5, conditions = conditions2)


label_DWA <- c("shoot biomass = 0.52 g","shoot biomass = 0.71 g","shoot biomass = 0.90 g" )

names (label_DWA) <- c ("DW_above = 0.52", "DW_above = 0.71", "DW_above = 0.9")
ggplot (plot_2$`Dim.1:DW_roots`, aes (x = Dim.1, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_bw ( ) + 
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





# plot (conditional_effects(PLFAtot_PCA_06New, effects = "Dim.1",  ndraws = 10000, 
#                           spaghetti = F,  prob = 0.9), points = T)
# 

# Appendix S2: Figure S8 pp check ####
pp_check (PLFAtot_PCA_06New, ndraws = 100) +
  xlab ("Bacterial biomass in soil")

# Appendix S2: Table S9 ####

launch_shinystan(PLFAtot_PCA_06New
                 )
## Appendix S2: Table S6 #####
dataNL %>%  
  filter (sampleID != "R18")  %>%   ## plant lost shoot biomass as it was almost dead at harvest, removed from analysis
left_join ()
dataNLPL %>%  select (sampleID, PlantSpeciesfull, AMFroot, AMF, totalBact, DW_roots, DW_above) %>% 
  write_csv("results/tableS6.csv")
