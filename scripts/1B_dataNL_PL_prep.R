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
meta_M1Wcontrol <- read_xlsx ("data/meta_M1Wcontrol.xlsx")

# NL data ##
dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")# NLFA_data for exp 2 #
dataNL <- dataNLrootstosoil %>%  select (sampleID, AMF, AMFroot, ID, PlaSpe, PlantSpeciesfull, DW_above, DW_roots, totalAMF, PlantFamily)
# tree data ##
A <- read_rds ("data/covariance_str_tree_trnl_plants.rds") #covaraince structure based on trnl tree , for glmm models

# PLFA data ###

PL_FA_conc_Soil <-  readRDS ("data/PL_FA_conc_Soil.rds") 


PLFA_Soil <-
  PL_FA_conc_Soil %>%
  filter (Biom2Soil !="IS") %>%
  filter (Biom2Soil != "NAFA") %>%  
  group_by(sampleID, Biom2Soil) %>%  
  summarise (group = sum (conc_FA)) %>%  
  left_join( meta_M1 %>%  select (sampleID, ID), by='sampleID') %>%  # add meta again
  pivot_wider(names_from = Biom2Soil, values_from = group ) %>% 
  mutate (totalPLFA = Actinobacteria +AMF+ bacteria+Fungi+gram.neg+gram.pos) %>% 
  mutate  (RatioFtoB  = (AMF +Fungi)/(bacteria + gram.neg + gram.pos + Actinobacteria)) %>% 
  mutate (GPtoGN = gram.pos / gram.neg)
  
# Nl and PLFa data in one tibble ###
dataNLPL <- 
  dataNLrootstosoil %>% 
  left_join(PLFA_Soil %>%  
              ungroup () %>%  
              select (!c(AMF, sampleID))) %>% 
  filter (!(is.na (AMF)))  # removed samples without soil measurements



                                