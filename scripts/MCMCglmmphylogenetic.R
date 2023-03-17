## Phylogenetic informed analysis####


# packages #####
library(tidyverse)
#library (ape)
library (broom.mixed)
library (brms)
library (shinystan)
library (tidybayes)
# data  ####

dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")
metaM0 <- read_xlsx("dataCh2/meta/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")
PL_FA_conc_Soil <- read_rds("data/PL_FA_conc_soil.rds")
All_metrics_M1_8Sp <- read_rds ("data/All_metrics_M1_8Sp.rds")
PC_data_M1roots  <-  read_rds ("data/PC_data_M1roots.rds")


# cleanup data NL ####
Metrics_exp2  <-
  adiv_richness_M1 %>% select (sampleID, Observed, Shannon) %>% 
  left_join(meta_M1 %>%  select (sampleID, PlaSpe))  %>%  
  left_join (stand_pd_Glo_all_M1 %>% select (sampleID, pd.obs)) %>% 
  left_join(ses.MPD_Glo_M1 %>%  select (sampleID, mpd.obs))  %>%  
  drop_na () %>% 
  add_column (exp = "Exp2")

saveRDS(Metrics_exp2, "data/Metrics_exp2.rds")

dataNL <- 
  dataNLrootstosoil %>%  
  select (sampleID, AMF, AMFroot, PlaSpe, 
          PlantSpeciesfull, DW_above, DW_roots, totalAMF, 
          Dim.1, Dim.2, Dim.3) %>% 
  left_join (Metrics_exp2) %>% 
  mutate (AMF = 100 * AMF, AMFroot = 100 * AMFroot, totalAMF = 100 * totalAMF)



# check out data ###


(plot_hist <- ggplot(dataNL, aes(x = AMF)) +
    geom_histogram(binwidth = 0.1) +
    theme_classic())


## slightly right skewed

# Poisson?? 


# BRMS package MCMCglmm #####

## Gutes tutoria: ###
#   https://ourcodingclub.github.io/tutorials/brms/   #

# data ##
head (dataNL)

#To enable parallel computing, 
#you can run this line of code and 
#then later on in the model code, you can specify how many cores you want to use.
options(mc.cores = parallel::detectCores())

# get variance structure for model ###
#A <- ape::vcv.phylo(treeML) 

# test first model #
model_simple <- brm(
  AMF ~ mpd.obs + (1|gr(PlaSpe, cov = A)),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A)
)
# change from 0.8 to higher (0.8 to 1) if divergent transitions ##
control = list(adapt_delta = 0.9)

summary(model_simple)

#On the top of the output, some general information on the model is given, such as family, formula, number of iterations and chains. 

#Next, group-level effects are displayed separately for each grouping factor in terms of standard deviations and 
#(in case of more than one group-level effect per grouping factor; not displayed here) correlations between group-level effects. 

#On the bottom of the output, population-level effects (i.e. regression coefficients) are displayed.
#If incorporated, autocorrelation effects and family specific parameters (e.g., the residual standard deviation ‘sigma’ in normal models) are also given.

# If Rhat is considerably greater than 1 (i.e., > 1.1), the chains have not yet converged and it is
 #necessary to run more iterations and/or set stronger priors
#in general only fully trust the sample if R-hat is less than 1.01. In early workflow, R-hat below 1.1 is often sufficien

# EFF more than 1000 is a good sign


# One way to assess model convergence is by visually examining the trace plots. 
# They should be fuzzy with no big gaps, breaks or gigantic spikes.

plot(model_simple, N = 2, ask = FALSE)

# Effects of population-level predictors can also be visualized with the
 #conditional_effects method
plot(conditional_effects(model_simple), points = TRUE)


hyp <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"
hyp <- hypothesis(model_simple, hyp, class = NULL)

#Phylogenetic Model with Repeated Measurements ####

# add the means to the table
dataNL$mean_DW_roots <-
  with(dataNL, sapply(split(DW_roots, PlaSpe), mean)[PlaSpe])

dataNL$mean_DW_above <-
  with(dataNL, sapply(split(DW_above, PlaSpe), mean)[PlaSpe])
#The variable mean_xx just contains the mean of the cofactor for each species. 
head(dataNL)
#The code for the repeated measurement phylogenetic model looks as follows:


#The most important reason to use control is to decrease (or eliminate at best) the number of divergent transitions
#that cause a bias in the obtained posterior samples. Whenever you see the warning
#"There were x divergent transitions after warmup.", you should really think about
#increasing adapt_delta. To do this, write control = list(adapt_delta = <x>), where
#<x> should usually be a value between 0.8 (current default) and 1. Increasing adapt_delta
#will slow down the sampler but will decrease the number of divergent transitions threatening
#the validity of your posterior samples.

# Model AMF ~  Dim.1 ####

## Set appropriate priors using get_prior function #### 
prior <- get_prior ( AMF ~ mpd.obs + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull), data = dataNL, data2 = list(A = A))

# model 1
 ##
model_repeat1 <- brm(
  AMF ~ mpd.obs + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(model_repeat1)

# estimate of phylogenetic signal!!
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(model_repeat1, hyp, class = NULL)

hyp
### means we ahave a phylogenetc signal 
plot (hyp)


#So far, we have completely ignored the variability of the cofactor within species. 
#To incorporate this into the model, we define

dataNL$within_spec_DW_roots <- dataNL$DW_roots - dataNL$mean_DW_roots

dataNL$within_spec_DW_above <- dataNL$DW_above - dataNL$mean_DW_above

#and then fit it again using within_spec_xx as an additional predictor.

model_repeat2 <- update(
  model_repeat1, formula = ~ . + within_spec_DW_roots + within_spec_DW_above,
  newdata = dataNL, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

summary (model_repeat2)
# compare with first model
summary (model_repeat1)


# estimate phylogenetic signal ###
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(model_repeat1, hyp, class = NULL)

hyp
plot (hyp)




# fit new model, 

model_repeat3 <- update(model_repeat2, formula. = ~ . - mean_DW_above - within_spec_DW_above)

summary (model_repeat3)


pp_check (model_repeat3, ndraws = 100)

launch_shinystan(model_repeat3)


plot (conditional_effects(model_repeat3, effects = "mpd.obs", method = "posterior_predict"), points = T)



