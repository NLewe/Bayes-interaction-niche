

# Data and start on Bayes modelling fo janelle ### 

library (tidyverse)
library (tidybayes)
library (brms)


library (phangorn)

# This is how the covariance structure was calculated: '

## Build tree ####
dnaTRNL<-read.dna("data/trnlFasta8.fasta",format="fasta")

dnaTRNL
labels(dnaTRNL) #check, change
rownames(dnaTRNL)<-c("AchMil" ,"CicInt","PlaLan", "HolLan" ,"PoaCit", "BroWil","SchAru", "AgrCap") #### YOU CAN NAME YOUR PLANTS HERE

## Build  tree
# I used JalView (software) to generate the tree and the pairwise alignment
tree_jalView<-read.tree("data/NJ_tree_8")
tree_jalView$tip.label # check
plot(tree_jalView)
# 
# # rename 
tree_jalView$tip.label <- c("PlaLan","CicInt","AchMil", "BroWil","AgrCap","SchAru", "PoaCit","HolLan"  ) 

## optimise tree using maximum likelihood
dna2 <- as.phyDat(dnaTRNL) 
class(dna2)

tre.ini<-nj(dist.dna(dnaTRNL, model="TN93"))
tre.ini

#To  initialize  the  optimization  procedure,  we  need  an  initial  fit  for  the model chosen.  
#This is computed using pml
fit.ini<-pml(tre.ini, dna2,k=4)

fit.ini<-pml(tree_jalView, dna2,k=4)

#Optimize tree
fit<-optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE,optGamma = TRUE)
fit
class(fit)
names(fit) ##all the useful parameters of the model and the optimal tree
#optimal tree is stored in fit$tree

##test if optimized tree is better than original
anova(fit.ini,fit)
AIC(fit.ini)
AIC(fit) # YES

treeML<-root(fit$tree,1)
plot(treeML)

# get variance structure for model ###
A <- ape::vcv.phylo(treeML) 
A <- saveRDS ("data/covariance_str_tree_trnl_plants.rds") #covaraince structure based on trnl tree , for glmm models




# Bayes LMM ####
# phylogenetic covariance
# data 
A <- read_rds ("data/covariance_str_tree_trnl_plants.rds") 

RelGen_E1_E2_sample <- readRDS("data/Exp1_Exp2_relGeneralims_NL.rds")

### maybe using the means per plant species??
## Add to your data 

# # 0A Specify the  model 0####

## get default  priors  #### 

prior0 <- get_prior ( AMF ~ RelGenSpec * DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
                      ### is only needed to account for diff between the species 
                      #+ OTHER than phylogenetic (environmental factors, niches)
                      data = YOURDATA, data2 = list(A = A)) # data2 is based on the covariance, you need to have the same names in A and in YOURDATA

# model 0 


Gen_samples_fit0 <- brm(
  AMF ~ RelGenSpec * DW_roots  + (1 |gr(PlaSpe, cov = A))+ (1|PlantSpeciesfull),
  data = dataNL_sample,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior0,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(Gen_samples_fit0)

# check convergence , look at estimates!
launch_shinystan(Gen_samples_fit0)

# posterior predictive checks #

pp_check (Gen_samples_fit0, ndraws= 100) +
  xlab ("AMF biomass in soil")

# 1 Update to model 01 ####

#and then fit it again - adding or removing more variables - there is no specific rule

Gen_samples_fit01 <- update(
  Gen_samples_fit0, formula = ~ . - DW_roots ,
  newdata = dataNL_sample, chains = 4, cores = 8,
  iter = 5000, warmup = 2000
)

## Model selection/ cross validation
#https://mc-stan.org/loo/articles/online-only/faq.html

Gen_samples_fit0 <- add_criterion(Gen_samples_fit0, "loo")

Gen_samples_fit01 <- add_criterion(Gen_samples_fit01, "loo")


loo_compare (Gen_samples_fit0, Gen_samples_fit01)



