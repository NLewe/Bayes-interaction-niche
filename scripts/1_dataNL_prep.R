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

dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")# NLFA_data for exp 2 #
All_metrics_M1_8Sp <- read_rds ("data/All_metrics_M1_8Sp.rds")# mean metrics for exp 2 interaction niches #
Metrics_exp2 <- read_rds("data/Metrics_exp2.rds")  # samplewise results of interaction niche properties #

PC_data_M1roots  <-  read_rds ("data/PC_data_M1roots.rds")# results from PCA of interaction niche properties #

A <- read_rds ("data/covariance_str_tree_trnl_plants.rds") #covaraince structure based on trnl tree , for glmm models


# data 
dataNL <- 
  dataNLrootstosoil %>%  
  select (sampleID, AMF, AMFroot, PlaSpe, 
          PlantSpeciesfull, DW_above, DW_roots, totalAMF, # both AMF and AMFroot are in µmol per g DW soil or roots
          Dim.1, Dim.2, Dim.3) %>% 
  left_join (Metrics_exp2) %>% 
  mutate (AMF = AMF*1000)



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
