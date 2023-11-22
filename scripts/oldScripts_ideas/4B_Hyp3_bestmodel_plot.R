
# packages ####
library (tidyverse)
library (tidybayes)
library (brms)
library(modelr)
library (broom.mixed)
library(GGally)

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# https://mc-stan.org/loo/articles/online-only/faq.html







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

