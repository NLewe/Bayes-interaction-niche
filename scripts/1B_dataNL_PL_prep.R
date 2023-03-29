## Phylogenetic informed analysis####


# packages #####
library(tidyverse)
library(readxl)
# data  ####

## all  data is prepared in the R_interaction_niches project ##
# only copied into folder data/  #

#meta files ##
metaM0 <- read_xlsx("data/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")

# NL data ##
dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")# NLFA_data for exp 2 #

# tree data ##
A <- read_rds ("data/covariance_str_tree_trnl_plants.rds") #covaraince structure based on trnl tree , for glmm models

# PLFA data ###






# check out data ###


(plot_hist <- ggplot(dataNL, aes(x = AMF)) +
    geom_histogram(binwidth = 0.001) +
    theme_classic())


## slightly right skewed

# Poisson?? 


# summary of data 

dataNL %>%  filter (!is.na (AMF)) %>% 
  group_by (PlantSpeciesfull) %>%  
  summarize (meanAMFsoil = mean (AMF), SDAMFsoil = sd (AMF))


dataNL %>%  filter (!is.na (AMF)) %>% 
  ungroup () %>% 
  summarize (meanMPD = mean (mpd.obs), SDMPDsoil = sd (mpd.obs))
