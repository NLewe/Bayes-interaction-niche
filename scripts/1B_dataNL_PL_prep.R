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
dataNL <- dataNLrootstosoil %>%  select (sampleID, AMF, AMFroot, ID, PlaSpe, PlantSpeciesfull, 
                                         DW_above, DW_roots, totalAMF, PlantFamily) %>% 
  mutate (AMF = 1000 * AMF)  # soil AMF is now in nmol

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
  mutate (GPtoGN = gram.pos / gram.neg) %>% 
  mutate (totalBact = Actinobacteria + bacteria+gram.neg+gram.pos)
  
  
# Nl and PLFa data in one tibble ###
dataNLPL <- 
  dataNLrootstosoil %>% # amf soil is in µmol
  mutate (AMF = AMF * 1000) %>%  # AMF is now in nmol per g 
  left_join(PLFA_Soil %>%  
              ungroup () %>%  
              select (!c(AMF, sampleID))) %>% 
  filter (!(is.na (AMF)))  # removed samples without soil measurements

# Tables for paper #### 
# Table S5 - all values metrics per plant species  ###
All_metrics_E1_8Sp %>%  rbind(All_metrics_E2) %>% 
  mutate  (Shannon = round (Shannon,2 ), 
           Richness = round (Richness,2), CU = round (CU, 2), 
           PD = round (PD,2), MPD = round (MPD, 2)) %>%  
  select (!PlantFamily) %>% write_excel_csv("results/All_metrics_E1_8Sp.csv")



# table S 6 alpha diversity per sample  #
adiv_richness %>%  
  left_join(metaM0 %>%  select (sampleID, PlantSpeciesfull)) %>%  
  mutate (Shannon = round (Shannon, 3)) %>% 
  select (-Chao1, -se.chao1) %>% write_excel_csv("results/tableS6_alpha_div_E1.csv")

adiv_richness_M1 %>%  select (-PlaSpe, -PlantFamily) %>%  
  mutate (Shannon = round (Shannon, 3)) %>% 
  write_excel_csv("results/tableS6_alpha_div_E2.csv")

                                