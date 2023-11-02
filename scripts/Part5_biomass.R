## Part 5 Addition of biomass data from PLFA and NLFA  ####


# packages #####
library(tidyverse)
library(readxl)
# data  ####

#meta files ##
metaM0 <- read_xlsx("data/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <- read_xlsx ("data/meta_M1Wcontrol.xlsx")

# NLFA data: neutral lipid analysis of the AMF biomarker as a proxy for the AMF biomass in the soil and roots ##
dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")# NLFA_data for exp 2 #
dataNL <- dataNLrootstosoil %>%  select (sampleID, AMF, AMFroot, ID, PlaSpe, PlantSpeciesfull, 
                                         DW_above, DW_roots, totalAMF, PlantFamily) %>% 
  mutate (AMF = 1000 * AMF)  # soil AMF is now in nmol per g soil 

# tree data ##
A <- read_rds ("data/covariance_str_tree_trnl_plants.rds") #covaraince structure based on trnl tree , for glmm models

# PLFA data upload, phospholid data is used as a proxy for the bacterial biomass in the soil ###

PL_FA_conc_Soil <-  readRDS ("data/PL_FA_conc_Soil.rds") 

# remove IS = internal standard, and NAFA = non-identified fatty acids from the dataset #
PLFA_Soil <-
  PL_FA_conc_Soil %>%
  filter (Biom2Soil !="IS") %>%
  filter (Biom2Soil != "NAFA") %>%  
  group_by(sampleID, Biom2Soil) %>%  
  summarise (group = sum (conc_FA)) %>%  # get the sum per bacterial group
  left_join( meta_M1 %>%  select (sampleID, ID), by='sampleID') %>%  # add meta again
  pivot_wider(names_from = Biom2Soil, values_from = group ) %>% 
  mutate (totalPLFA = Actinobacteria +AMF+ bacteria+Fungi+gram.neg+gram.pos) %>% 
  mutate  (RatioFtoB  = (AMF +Fungi)/(bacteria + gram.neg + gram.pos + Actinobacteria)) %>% # calculate interesting ratios that might be of use later
  mutate (GPtoGN = gram.pos / gram.neg) %>% 
  mutate (totalBact = Actinobacteria + bacteria+gram.neg+gram.pos) # total bacterial biomass


# Put NL and PLFa data in one tibble ###
dataNLPL <- 
  dataNLrootstosoil %>% # amf soil is in µmol
  mutate (AMF = AMF * 1000) %>%  # AMF is now in nmol per g 
  left_join(PLFA_Soil %>%  
              ungroup () %>%  
              select (!c(AMF, sampleID))) %>% 
  filter (!(is.na (AMF)))  # removed samples without soil measurements

                            