### Part 11 - Additional Bayesian models for Appendix S2  ####


# Using diversity metrics, the models for describing changes in the AMF and bacterial biomass were refitted.
# 
# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library (shinystan)


# get data ####
 data_metrics<-
  dataNL_sample %>%
  left_join (All_Metrics_E2_sample) %>%
  filter (!is.na (AMF)) %>% 
   rename ("Richness" = `Richness S`) %>% 
   unique ()
  

# 0A Specify the MPD model 0####
## get default  priors  #### 
priorB <- get_prior ( AMF ~ MPD  + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     data = data_metrics, data2 = list(A = A))


# Null model MPD + roots ###### 


mpd_fit0 <- brm(
  AMF ~ MPD + DW_roots,
  data = data_metrics,
  family = gaussian(),
  data2 = list(A = A),
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)


#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(mpd_fit0)

bayes_R2 (mpd_fit0)
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

bayes_R2 (mpd_fit0)

# 1A Update to model 01 BEST ####

#and then fit it again - adding more variables

mpd_fit01 <- update(
  mpd_fit0, formula = ~ . -MPD + MPD:DW_roots ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (mpd_fit01)
bayes_R2 (mpd_fit01)



# 1B check sampling quality of model 01 ####


# check convergence #
launch_shinystan(mpd_fit01)

# posterior predictive checks #

pp_check (mpd_fit01, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models #####
mpd_fit0 <- add_criterion(mpd_fit0, "loo")

mpd_fit01 <- add_criterion(mpd_fit01, "loo")
loo_compare (mpd_fit0, mpd_fit01)
# best performing model will be named at top
#
bayes_R2 (mpd_fit01)


# 2A Update to model 02 ####

#and then fit it again - adding more variables

mpd_fit02 <- update(
  mpd_fit0, formula = ~ . + DW_above ,
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
loo_compare (mpd_fit01, mpd_fit02, mpd_fit0)
bayes_R2(mpd_fit02)


# 3A Update to model 03 ####

#and then fit it again - adding more variables

mpd_fit03 <- update(
  mpd_fit01, formula = ~ .    ,
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

# 3C Compare models #####

mpd_fit03 <- add_criterion(mpd_fit03, "loo")
loo_compare (mpd_fit02, mpd_fit0, mpd_fit03)
### Model 02 is the best ###




# 4A Update to model 04 ####

#and then fit it again - adding more variables

mpd_fit04 <- update(
  mpd_fit02, formula = ~ .   + Richness ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 4B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit04)

# check convergence #
launch_shinystan(mpd_fit04)

# posterior predictive checks #

pp_check (mpd_fit04, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 4C Compare models #####

mpd_fit04 <- add_criterion(mpd_fit04, "loo")
loo_compare (mpd_fit02, mpd_fit03, mpd_fit04)
bayes_R2 (mpd_fit04)
### Model 02 is the best ###



# Hypothesis testing MPD model #####
# 
h <- "Richness<0"
hypothesis (mpd_fit04,h)

### get posteriors from the model

# 5A Update to model 05 ####

#and then fit it again - adding more variables

mpd_fit05 <- update(
  mpd_fit02, formula = ~ .   + Richness -MPD ,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 5B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit05)

# check convergence #
launch_shinystan(mpd_fit05)

# posterior predictive checks #

pp_check (mpd_fit05, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 5C Compare models #####

mpd_fit05 <- add_criterion(mpd_fit05, "loo")
loo_compare (mpd_fit02, mpd_fit03, mpd_fit04, mpd_fit05)
bayes_R2 (mpd_fit05)


#6A and then fit it again #####

mpd_fit06 <- update(
  mpd_fit05, formula = ~ .   -DW_above,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 6B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit06)

# check convergence #
launch_shinystan(mpd_fit06)

# posterior predictive checks #

pp_check (mpd_fit06, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 6C Compare models #####

mpd_fit06 <- add_criterion(mpd_fit06, "loo")
loo_compare (mpd_fit02, mpd_fit03, mpd_fit06, mpd_fit05, mpd_fit04)
bayes_R2 (mpd_fit06)

#7 A and then fit it again#####
#
mpd_fit07 <- update(
  mpd_fit0, formula = ~ .   + PlantSpeciesfull,
  newdata = data_metrics, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


# 7B check sampling quality of model 03 ####

# Look for rhat, ESS, Sd etc
summary (mpd_fit07)

# check convergence #
launch_shinystan(mpd_fit07)

# posterior predictive checks #

pp_check (mpd_fit06, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 7C Compare models #####

mpd_fit07 <- add_criterion(mpd_fit07, "loo")
loo_compare (mpd_fit02, mpd_fit03, mpd_fit06, mpd_fit05, mpd_fit07)
bayes_R2 (mpd_fit07)

h <- "Richness>0"
hypothesis (mpd_fit06,h)




# Figure AMF MPD Bayes ######

plot3 <- conditional_effects(mpd_fit0, effects = "MPD:DW_roots", ndraws = 1000, spaghetti = F, mean = T, prob = 0.5)
ggplot (plot3$`MPD:DW_roots`, aes (x = MPD, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_classic ( ) + 
  geom_point (data = data_metrics, aes (x= MPD, y = AMF), inherit.aes = F) +
  theme (axis.title = element_text(size =16, color = "black" ), 
         axis.text = element_text(size = 14, color = "black"),
         legend.position = "right", 
         legend.text = element_text(size = 14, color = "black"), 
         legend.title = element_text(size = 16, color = "black")) +
  xlab ( "Mean phylogenetic distance (MPD)") +
  ylab(expression("NLFA 16:1"*omega*"5 in nmol/g soil")) +
  scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
  scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 
ggsave ( "figures/MPD_AMF_Bayes.png", height = 6, width =9, dpi = 300)    

pp_check (mpd_fit02, ndraws= 100) +
  xlab ("AMF biomass in soil")

launch_shinystan (mpd_fit0)



# Richness Bacterial Biomass model #####
# 
data_metricsPL<-
  dataNLPL  %>%
  left_join (All_Metrics_E2_sample) %>%
  filter (!is.na (totalBact)) %>% 
  rename (Richness = 'Richness S')

# Note that more models as shown here were run.

Richness_fit0 <- brm(
  totalBact ~ Richness,
  data = data_metricsPL,
  family = gaussian(),
  data2 = list(A = A),
  chains = 4, cores = 8,
  iter = 2000, warmup = 500
)
bayes_R2 (Richness_fit0)
summary (Richness_fit0)
Richness_fit0 <- add_criterion(Richness_fit0, "loo")

Richness_fit01 <- update(
  Richness_fit0, formula = ~ .  + DW_roots ,
  newdata = data_metricsPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (Richness_fit01)
bayes_R2 (Richness_fit01)
Richness_fit01 <- add_criterion(Richness_fit01, "loo")


loo_compare(Richness_fit0, Richness_fit01)



Richness_fit02 <- update(
  Richness_fit0, formula = ~ .  + DW_roots + DW_above ,
  newdata = data_metricsPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (Richness_fit02)
bayes_R2 (Richness_fit02)
Richness_fit02 <- add_criterion(Richness_fit02, "loo")


loo_compare(Richness_fit0, Richness_fit01, Richness_fit02)

Richness_fit03 <- update(
  Richness_fit0, formula = ~ .  + DW_roots:DW_above ,
  newdata = data_metricsPL, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

summary (Richness_fit03)
bayes_R2 (Richness_fit03)
Richness_fit03 <- add_criterion(Richness_fit03, "loo")


loo_compare( PLm10New, Richness_fit03)
PLm10New <- add_criterion(PLm10New, "loo")

# Hypothesis testing Richnes #####
# 
 h <- "Richness <0"

hypothesis (Richness_fit02, h)

Richness_fit03 <- brm(
  totalBact ~ Richness + DW_roots + DW_above,
  data = data_metricsPL,
  family = gaussian(),
  data2 = list(A = A),
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)

launch_shinystan(Richness_fit03)

pp_check(Richness_fit03, ndraws= 100) +
  xlab ("Bacterial biomass in soil")

# Figure ###
# 
# Figure 3 for paper  ####
conditions3 = make_conditions (Richness_fit03, "DW_above")

plot_4 <- conditional_effects(Richness_fit03, effects = "Richness:DW_roots",  ndraws = 50000, 
                              spaghetti = F,  prob = 0.5, conditions = conditions3)


label_DWA <- c("shoot biomass = 0.52 g","shoot biomass = 0.71 g","shoot biomass = 0.90 g" )

names (label_DWA) <- c ("DW_above = 0.52", "DW_above = 0.71", "DW_above = 0.9")
ggplot (plot_4$`Richness:DW_roots`, aes (x = Richness, y = estimate__, color = effect2__)) +
  geom_smooth  (method = lm) + 
  geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
  theme_bw ( ) + 
  geom_point (data = data_metricsPL, aes (x= Richness, y = totalBact), inherit.aes = F)+
  theme (axis.title = element_text(size =16, color = "black" ), 
         axis.text = element_text(size = 14, color = "black"),
         legend.position = "right", 
         legend.text = element_text(size = 14, color = "black"), 
         legend.title = element_text(size = 16, color = "black"), 
         panel.grid = element_blank()) +
  facet_wrap(~cond__, labeller = labeller (cond__ = label_DWA)) +
  theme ( strip.text.x = element_text(size = 14),
          legend.position = "right") +
  xlab ( "Richness") +
  ylab(expression("Bacterial PLFA in nmol/g soil")) +
  scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
  scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 

ggsave ( "figures/Bact_Richness_PC.png", height = 6, width =12, dpi = 300)    







