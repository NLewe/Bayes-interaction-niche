


prior0 <- get_prior (   AMF  ~ RelGenSpec : DWRA + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) , 
                      ### is only needed to account for diff between the species 
                      #+ OTHER than phylogenetic (environmental factors, niches)
                         family = gaussian(),
                      data = dataNL_sample, data2 = list(A = A))
# model 0 
control = list(adapt_delta = 0.93) ## then rerun!!


fit1 <- brm(
  AMF  ~ RelGenSpec :DWRA + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull) ,
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior0,
  chains = 4, cores = 8,
  iter = 2000, warmup = 500
)


fit1B <-  update(
  fit1, formula = ~ . - (1|gr(PlaSpe, cov = A)) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)


fit1 <- add_criterion (fit1, c ("waic", "loo"))

fit1B <- add_criterion (fit1B, c ("waic", "loo"))
loo_compare (fit1, fit1B) 


fit1BSIM <-  update(
  fit1B, sample_prior = "only", 
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)



# Does RG explain it alone
fit1C <-  update(
  fit1B, formula = ~ . -   RelGenSpec:DWRA + RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)


## Or is only the biomass of relevance?
fit1D <-  update(
  fit1B, formula = ~ . -   RelGenSpec:DWRA + DWRA ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)

## Is additive formula better as the interaction?
fit1E <-  update(
  fit1D, formula = ~ . +   RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)
fit1C <- add_criterion (fit1C, c ("waic", "loo"))

fit1B <- add_criterion (fit1B, c ("waic", "loo"))

fit1D <- add_criterion (fit1D, c ("waic", "loo"))
fit1E <- add_criterion (fit1E, c ("waic", "loo"))

loo_compare (fit1B, fit1D, fit1C, fit1E, fit1) 


## best model with interaction but without any random effect? 
plot (conditional_effects (fit1B, effects = "RelGenSpec:DWRA", prob =0.5), points = T)
 fit1F <-  update(
  fit1B, formula = ~ . -   (1 | PlantSpeciesfull) ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)

 fit1F <- add_criterion (fit1F, c ("waic", "loo"))

## explais less (R2, but no problems) 


## additive instead of interaction?
fit1G <-  update(
  fit1F, formula = ~ . -   RelGenSpec:DWRA + DWRA + RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)
## maybe slightly better? Similar loo
fit1G <- add_criterion (fit1G, c ("waic", "loo"))
loo_compare (fit1G, fit1F)


## only DWRA??
fit1H <-  update(
  fit1G, formula = ~ . -   RelGenSpec ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)
fit1H <- add_criterion (fit1H, c ("waic", "loo"))
loo_compare (fit1G,  fit1H)


# model without any variable?

fit1J <-  update(
  fit1H, formula = ~ . -   DWRA + 1 ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)

fit1J <- add_criterion (fit1J, c ("waic", "loo"))
loo_compare (fit1G,  fit1H, fit1J, criterion ="waic")


fit1E <- add_criterion (fit1E, c ("waic", "loo"))





pp_check (fit1B, ndraws = 100)

launch_shinystan(fit1B)

bayes_R2(fit1F)


plot (conditional_effects (fit1F, prob = 0.9), points = T)

### FIT2 IS BEST ######

## model 2 additional
prior2 <- get_prior (   AMF  ~ RelGenSpec: DW_roots + (1|PlantSpeciesfull),
                        ### is only needed to account for diff between the species 
                        #+ OTHER than phylogenetic (environmental factors, niches)
                        family = gaussian (),
                        data = dataNL_sample, data2 = list(A = A))

fit2 <- brm(
  AMF  ~ RelGenSpec: DW_roots + (1|PlantSpeciesfull),
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior2, 
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)


fit2 <- add_criterion(fit2, "loo")

loo_compare (fit2, fitBEST)
             
             
             fit3 <-  update(
  fit2, formula = ~ . -  DW_above,
  newdata = dataNLPL, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)

fit3SIM <-  update(
  fit2, formula = ~ . -  DW_above, sample_priors = "only", 
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 3000, warmup = 1000
)


pp_check (fit3, ndraws = 100)

launch_shinystan(fit3)

fit3 <- add_criterion(fit3, "loo")
conditional_effects(fit3, effects = "RelGenSpec:DW_roots",  ndraws = 10000, 
                    spaghetti = F,  prob = 0.5)

# Best model Refit #####
priorBest <- get_prior (   AMF  ~ RelGenSpec : DWRA  + (1|PlantSpeciesfull) , 
                        ### is only needed to account for diff between the species 
                        #+ OTHER than phylogenetic (environmental factors, niches)
                        family = gaussian(),
                        data = dataNL_sample, data2 = list(A = A))
fitBEST <- brm(
  AMF  ~ RelGenSpec :DWRA +  (1|PlantSpeciesfull) ,
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = priorBest,
  chains = 8, cores = 8,
  iter = 10000, warmup = 3000
)

fitBEST <- add_criterion(fitBEST, "loo")

## plot conditional effects model fit2 ####
conditional_effects(fitBEST, effects = "RelGenSpec:DWRA",  ndraws = 10000, spaghetti = F,  prob = 0.5)

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


