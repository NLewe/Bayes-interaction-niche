

# Part 8 Bayesian models for interaction generalism ####
# The following models are examples for variable constellation to find the best model explaining 
#  the AMF biomass in the soil of plant species.
# As variables describing plants' interaction generalism, we used the loadings
# from a principal component analysis.
# 
# 
# Packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library(GGally)
#library (lmerTest)
library (rstatix)
library (shinystan)

# Tutorials can be found here: 

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# https://mc-stan.org/loo/articles/online-only/faq.html


# Data cleanup  ####

dataNL_sample <- 
  dataNL %>%  
  filter (!is.na (AMF)) %>% 
  filter (sampleID != "R18")  %>%   ## plant lost shoot biomass, dead at harvest, removed from analysis
  left_join (PCA_metric_data_sample) %>% 
  left_join (metaM0 %>%  select (PlaSpe, PlantType) %>%  unique ()) %>% 
  # dplyr::select (!starts_with ("Dim")) %>% 
  mutate (DWRA = DW_roots /DW_above ) %>% 
  mutate (totRootAMF = AMFroot * DW_roots)

# AMF is the AMF biomass in the soil, calculated using the fatty acid biomarker for AMF
# AMF root is the AMF biomass in the plant roots


# some data checks #
# e.g. what is the distribution of the data?
(hist <- ggplot(dataNL_sample, aes(x = DWRA )) +
    geom_histogram(bins = 40) +
    theme_classic())


# Note that all combinations of the variables have been run and the results compared.
# Model selection was done by different criteria (see publication)

# 0A Specify  model 0####
## get default  priors  #### 
# prior0 <- get_prior ( AMF ~ Dim.1 + Dim.2 +Dim.3 + DW_roots + DW_above + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
#                       data = dataNL_sample, data2 = list(A = A))
# 
# # model 0 
#control = list(adapt_delta = 0.8) ## then rerun!!


mNull <- brm(
  AMF ~ 1 ,
  data = dataNL_sample,
  data2 = list(A = A),
  sample_prior = "yes",
  chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
summary (mNull)

m0 <- brm(
  AMF ~ 1 +     
    (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) + (1|PlantType),
  data = dataNL_sample,
  data2 = list(A = A),
  sample_prior = "yes",
   chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
summary (m0)
# 
# m0B <- brm(
#   AMF ~    1 + (1|PlantSpeciesfull),
#   data = dataNL_sample,
#   data2 = list(A = A),
#   chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# summary (m0B)
# m0C <- brm(
#   AMF ~    1 +  (1|gr(PlaSpe, cov = A)),
#   data = dataNL_sample,
#   data2 = list(A = A),
#   sample_prior = "yes",
#   chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# summary (m0C)

m0D <- brm(
  AMF ~    1 +  (1|PlantType),
  data = dataNL_sample,
  data2 = list(A = A),
  sample_prior = "yes",
  chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
summary (m0D)

## Add criteria for model fit comparison 
m0 <- add_criterion(m0, c("loo", "waic"))
m0B <- add_criterion(m0B, c("loo", "waic"))
m0C <- add_criterion(m0C, c("loo", "waic"))
m0D <- add_criterion(m0D, c("loo", "waic"))
mNull <- add_criterion(mNull, c("loo", "waic"))
# Compare loo
loo_compare(m0, mNull)

# Set up hypotheses for testing 
hyp0 <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"
hyp1 <- "sd_PlantSpeciesfull__Intercept^2 / (sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
hyp2 <- "sd_PlantType__Intercept^2 / (sd_PlantType__Intercept^2 + sigma^2) = 0"

hyp <- hypothesis(m0D, hyp2, class = NULL)


hyp


##Full model ######
##
mFull <- brm(
AMF ~ Dim.1+ Dim.2 +Dim.3 + DW_roots + DW_above  ,
data = dataNL_sample,
data2 = list(A = A),
sample_prior = "yes",
chains = 4, cores = 8,
iter = 2000, warmup = 500
)

# Read out the priors
get_prior(mFull)


# Hypothesis testing 

h6 <- "Dim.1 +Dim.2 +Dim.3>0"
h5 <- "Dim.3>0"
h <- "Dim.3>0"
h <- "DW_roots>0"

hypothesis (fitAMF06New,h)

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

# Obtain model summaries.
summary(m)
bayes_R2 (m)
 
# check convergenceetc #
launch_shinystan(m)

# posterior predictive checks #

pp_check (m, ndraws= 100) +
  xlab ("AMF biomass in soil")


# 1 Update to model 01 ####

#and then fit it again - adding more variables

m1 <- brm(
  AMF ~    Dim.1 ,
  data = dataNL_sample,
  data2 = list(A = A),
  chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
summary (m1)

hyp4 <- "Dim.3>1"
  hypothesis ( m3, hyp4)

  
##
m0 <- add_criterion(m0, c("loo", "waic"))

m1 <- add_criterion(m1, c("loo", "waic"))

bayes_R2 (m1)


loo_compare (mrandom, m1)
#  check sampling quality of model 01 ##

# Look for rhat, ESS, Sd etc
summary (m1)

# check convergence #
launch_shinystan(m1)

# posterior predictive checks #

pp_check (m1, ndraws= 100) +
  xlab ("AMF biomass in soil")



m2 <- update(
  m1, formula = ~ . -Dim.1 + Dim.2  ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)




# 3 Compare models ###

m2 <- add_criterion(m2, "loo")

loo_compare (m0, m1, m2)
# best performing model will be named at top
#



# 2 Update to model 02 ####

#and then fit it again - adding more variables

m2 <- update(
  m1, formula = ~ .  - Dim.1 + Dim.3,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m2 <- add_criterion(m2, "loo")


#  check sampling quality of model 02 ###

# Look for rhat, ESS, Sd etc
summary (m2)

bayes_R2 (m2)
# check convergence #
launch_shinystan(m2)

# posterior predictive checks #

pp_check (m2, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 3 Compare models ###

m2 <- add_criterion(m2, "loo")
loo_compare (m, m1, m2)




# 3 Update to model 03 ####
m3 <- update(
  m1, formula = ~ .  -Dim.1  + Dim.3  ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m3 <- add_criterion(m3, "loo")


loo_compare (m2,  m0, 
             m1, m3)

   
# Look for rhat, ESS, Sd etc
summary (m3)
bayes_R2(m3)

# 4 Update to model 04 ####

#and then fit it again - adding more variables

m4 <- update(
  m3, formula = ~ .  +  DW_above ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m4 <- add_criterion(m4, "loo")
loo_compare (m0,  m1, 
             m2, m3, 
             m4)


# Look for rhat, ESS, Sd etc
summary (m4)
pp_check (m4, ndraws = 100)


bayes_R2 (m4)

# 5 Update to model 05 ####

#and then fit it again - adding more variables

m5 <- update(
  m1, formula = ~ .  - Dim.1 + Dim.3:DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)


m5 <- add_criterion(m5, "loo")
loo_compare (m5,  m1, m3, 
             m4)
pp_check (m5, ndraws = 50)
bayes_R2 (m5)


# Look for rhat, ESS, Sd etc
summary (m5)


# 6 Update to model 06 ####

#and then fit it again - adding more variables

m6 <- update(
  m1, formula = ~ .   - Dim.1 + Dim.3 + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m6 <- add_criterion(m6, "loo")
loo_compare (m2,  m6, m3)
bayes_R2 (m6)


# Look for rhat, ESS, Sd etc
summary (m6)


# 7 Update to model 07 ####

#and then fit it again - adding more variables

m7 <- update(
  m1, formula = ~ .   - Dim.1 + DW_roots,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m7 <- add_criterion(m7, c("loo", "waic"))
loo_compare (m6, m7, m1, m5)

#  check sampling quality of model  ###

# Look for rhat, ESS, Sd etc
summary (m7)
bayes_R2 (m7)


# 8 Update to model 08 ####

#and then fit it again - adding more variables

m8 <- update(
  m1, formula = ~ . + Dim.3:DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m8 <- add_criterion(m8, "loo")
loo_compare (m6, m2, m3, m5, m4, m8, m7, m1)
launch_shinystan(m8)
summary (m8)

# 9 Update to model 09 ####

#and then fit it again - adding more variables

m9 <- update(
  m8, formula = ~ .  +Dim.2 - Dim.1,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m9 <- add_criterion(m9, "loo")
loo_compare (m6, m9, m5   , m8, m4)


m9 %>%  bayes_R2
m9 %>%  summary ()

# 10 Update to model 10 ####

#and then fit it again - adding more variables
dataNL_sample <- dataNL_sample %>% 
  left_join (metaM0 %>% select (PlaSpe, PlantType ))

m10 <- update(
  m3, formula = ~ .  + PlantType ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m10 <- add_criterion(m10, "loo")
loo_compare ( m10, m2, m6)
bayes_R2 (m10)

# 11 Update to model 11 ####
# 
m11 <- update(
  m10, formula = ~ .  -Dim.3 ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m11 <- add_criterion(m11, "loo")
loo_compare (  m10, m11)


bayes_R2 (m11)
# # 
# 12 Update to model 12 ####
m12 <- update(
  m3, formula = ~ .  + (1|PlantType)  ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m12 <- add_criterion(m12, "loo")
loo_compare ( m6, m12)
#
m13 <- update(
  m12, formula = ~ .  + DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)
#
m13 <- add_criterion(m13, "loo")
loo_compare ( m13, m6, m12)
# # 
bayes_R2 (m13)

m14 <- update(
  m13, formula = ~ .  -Dim.3 ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m14 <- add_criterion(m14, "loo")
loo_compare (m14, m10, m12, m13,m6)

m15 <- update(
  m6, formula = ~ .  + PlantType ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

m15 <- add_criterion(m15, "loo")
loo_compare ( m15, m13)

# 
# 
# PCA_axes_fit16 <- update(
#   PCA_axes_fit12, formula = ~ .  - DW_roots:RelGenSpec + RelGenSpec:DW_roots ,
#   newdata = dataNL_sample, chains = 4, cores = 8,
#   iter = 5000, warmup = 2000
# )
# 
# PCA_axes_fit16 <- add_criterion(PCA_axes_fit16, "loo")
# loo_compare (PCA_axes_fit16, PCA_axes_fit15, m, PCA_axes_fit12)
# 
# 
# 
# 

# ### Model  is the best ###
# 
# loo_compare (m,  m1, 
#              m2, m3, 
#              m4, m8, 
#              m5B, m5C, 
#              m1B, m5, 
#              m7, m9, 
#              PCA_axes_fit10, m6)

# While model 5B only including RelGenSpec is the best based on loo criterion, it is calculated with large amount of divergent transitions (222)
# The next best model is 05C RelGenSpec + DW_roots 

tidy (m6, effects = c("fixed")
### get posteriors from the model##
      
#
      
      
dataNL_sample %>%
        #group_by(PlaSpe) %>%
        #data_grid(PlaSpe = seq_range(PlaSpe, n = 51)) %>%
        add_epred_draws(fitAMF06New) %>%
        ggplot(aes(x = Dim.3, y = AMF, color = ordered(PlaSpe))) +
        stat_lineribbon(aes(y = .epred)) +
        geom_point(data = dataNL_sample) +
        scale_fill_brewer(palette = "Greys") +
        scale_color_brewer(palette = "Set2") #+
        facet_wrap(~PlaSpe, scales = "free_x")
      
      
### BEST MODEL and FIGURE ######
      
## Run model again with more chains and iterations ###
priorAMFNew <- get_prior (   AMF  ~  Dim.3 +  DW_roots ,
                                   ### is only needed to account for diff between the species 
                                   #+ OTHER than phylogenetic (environmental factors, niches)
                                   family = gaussian (),
                                   data = dataNL_sample, data2 = list(A = A))
      
      fitAMF06New <- brm(
        AMF  ~ Dim.3 + DW_roots,
        data = dataNL_sample,
        family = gaussian(),
        data2 = list(A = A),
        chains = 8, cores = 8,
        iter = 10000, warmup = 3000
      )
      launch_shinystan(fitAMF06New)
      
      fitAMF06New <- add_criterion(fitAMF06New, "loo")
      pp_check (fitAMF06New, ndraws = 100)
      loo (fitAMF06New)
      fitAMF06New %>%  bayes_R2 ()
      loo_compare (mpd_fit02, fitAMF06New)

## plot conditional effects model  ####
      
plot_1 <- conditional_effects(fitAMF06New, effects = "DW_roots:Dim.3",  ndraws = 5000, spaghetti = F,  prob = 0.5)
      
      
##Prepare plot Figure 2 #####
      
        ggplot (plot_1$`Dim.3:DW_roots`, aes (x = Dim.3, y = estimate__, color = effect2__)) +
        geom_smooth  (method = lm) + 
        geom_ribbon (aes (ymin = lower__, ymax = upper__, fill = effect2__ ), linewidth = 0, alpha = 0.5) +
        theme_classic ( ) + 
        geom_point (data = dataNL_sample, aes (x= Dim.3, y = AMF), inherit.aes = F) +
        theme (axis.title = element_text(size =16, color = "black" ), 
               axis.text = element_text(size = 14, color = "black"),
               legend.position = "right", 
               legend.text = element_text(size = 14, color = "black"), 
               legend.title = element_text(size = 16, color = "black")) +
        xlab ( "Principal component 3") +
        ylab(expression("NLFA 16:1"*omega*"5 in nmol/g soil")) +
        scale_color_manual(values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") +
        scale_fill_manual (values = c ("darkslategray4", "coral2", "blueviolet"), name = "Root biomass") 
# Save plot
ggsave ( "figures/AMF_PC.png", height = 6, width =9, dpi = 300)    
      
  
  
  plot_1b <- conditional_effects(fitAMF06New, effects = "Dim.3",  ndraws = 5000, spaghetti = F,  prob = 0.5)
  
  
# ##Prepare plot using ggplot #####
#   
# ggplot (plot_1$`Dim.3`, aes (x = Dim.3, y = estimate__)) +
#     geom_smooth  (method = lm) + 
#     geom_ribbon (aes (ymin = lower__, ymax = upper__), linewidth = 0, alpha = 0.5) +
#     theme_classic ( ) + 
#     geom_point (data = dataNL_sample, aes (x= Dim.3, y = AMF), inherit.aes = F) +
#     theme (axis.title = element_text(size =16, color = "black" ), 
#            axis.text = element_text(size = 14, color = "black"),
#            legend.position = "right", 
#            legend.text = element_text(size = 14, color = "black"), 
#            legend.title = element_text(size = 16, color = "black")) +
#     xlab ( "Principal component 3") +
#     ylab(expression("NLFA 16:1"*omega*"5 in nmol/g soil"))
# 
# ggsave ( "figures/AMF_PC_one_line.png", height = 6, width =9, dpi = 300)    
  
launch_shinystan(fitAMF06New)
      
# Appendix S2: Figure S6
      pp_check (fitAMF06New, ndraws = 100) +
        xlab ("AMF biomass in soil")
      
      
      