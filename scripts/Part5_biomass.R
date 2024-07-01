## Part 5 Addition of biomass data from PLFA and NLFA  ####


# Packages #####
library(tidyverse)
library(readxl)

# Import data ####

## Import meta files ##
metaM0 <- read_xlsx("data/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <- read_xlsx ("data/meta_M1Wcontrol.xlsx")

## NLFA data #### 
## neutral lipid fatty acid analysis data of the AMF biomarker 
# as a proxy for the AMF biomass in the soil and roots ##
dataNLrootstosoil <- 
  read_rds ("data/dataNLrootstosoil.rds")# NLFA_data for exp 2 #

# NLFA data cleanup #
dataNL <- 
  dataNLrootstosoil %>%  
  select (sampleID, AMF, AMFroot, ID, PlaSpe, PlantSpeciesfull, 
          DW_above, DW_roots, totalAMF, PlantFamily) %>% 
  mutate (AMF = 1000 * AMF)  # soil AMF is now in nmol per g soil 

##PLFA data ####
## phospholid fatty acid analysis data ##
# is used as a proxy for the bacterial biomass in the soil ###

PL_FA_conc_Soil <-  readRDS ("data/PL_FA_conc_Soil.rds") 

# PLFA data cleanup #
# remove IS = internal standard, and NAFA = non-identified fatty acids from the dataset #
PLFA_Soil <-
  PL_FA_conc_Soil %>%
  filter (Biom2Soil !="IS") %>%
  filter (Biom2Soil != "NAFA") %>%  
  group_by(sampleID, Biom2Soil) %>%  
  summarise (group = sum (conc_FA)) %>%  # get the sum per bacterial group
  left_join( meta_M1 %>%  select (sampleID, ID), by='sampleID') %>%  # add meta data again
  pivot_wider(names_from = Biom2Soil, values_from = group ) %>% 
  mutate (totalPLFA = Actinobacteria +AMF+ bacteria+Fungi+gram.neg+gram.pos) %>% 
  mutate  (RatioFtoB  = (AMF +Fungi)/(bacteria + gram.neg + gram.pos + Actinobacteria)) %>% # calculate interesting ratios that might be of use later
  mutate (GPtoGN = gram.pos / gram.neg) %>% 
  mutate (totalBact = Actinobacteria + bacteria+gram.neg+gram.pos) # total bacterial biomass


# Put NL and PLFA data together in one tibble ###
dataNLPL <- 
  dataNLrootstosoil %>% # AMF soil is in µmol
  mutate (AMF = AMF * 1000) %>%  # AMF is now in nmol per g 
  left_join(PLFA_Soil %>%  
              ungroup () %>%  
              select (!c(AMF, sampleID))) %>% 
  filter (!(is.na (AMF)))  # removed samples without soil measurements


## NLFA & PLFA control soil ####
## Add information on control soils (without plant) ##
## data ###
meta_control <- read_xlsx ("data/meta_soil_control.xlsx")
Biomarker<- read_xlsx("data/FAME_RRF_biomarker.xlsx", sheet = "BM")

path_controlNL <- "data/SoilControlNL.xlsx"
path_controlPL <- "data/SoilControlPL.xlsx"

# NL Soil control cleanup
NL_SoiCon<- path_controlNL%>%
  excel_sheets()%>%
  set_names() %>% 
  map_df(~read_excel(path=path_controlNL, sheet = .x), .id= "sheet")  %>%  ###here, it would be possible to give a range for the sheets etc
  dplyr::rename("sampleID"="sheet") %>%   #if using dplyr, add dplyr:: because of interference with plyr
  mutate (sampleID = str_replace_all (sampleID, "NL_", "")) %>% 
  add_column (PL_NL= "NLFA") %>% 
  select(sampleID, Name, area, PL_NL)


# add Biomarker info to data
NL_SoiCon<-merge(NL_SoiCon, Biomarker %>%  select(-c(Biom1Soil, Soil_A1_kons, MW, RI)), by="Name")

# add meta data
NL_SoiCon<-left_join(NL_SoiCon, meta_control, by="sampleID") 


# Calculate the loss of sample and correct areas of all FAMEs for loss
correction<- NL_SoiCon %>% select(c(sampleID,area, Name))%>% 
  filter(Name=="19:0")%>%
  mutate(IS_area=1*area) %>%  
  select(c(sampleID,IS_area )) 

# Full table for all NL fatty acids
NL_FA_conc_SoiCon <-  left_join(NL_SoiCon, correction, by ="sampleID")%>% 
    mutate (conc_FA = 5 * (area/IS_area) * (1/RRF) * (1/m_soil))  %>%  ## correction for sample loss via IS, and use of biomarker RRF (=relative retention factor) to calculate conc of each FA in nmol pe  g soil sample
  filter (Name != "19:0")   # remove IS from dataset
#Note RRF are only valid for the analytical instrument they were calculated on ###

dataNL_SoiCon <- 
  NL_FA_conc_SoiCon %>% 
  filter (Biom2Soil == "AMF") %>% 
  select  (sampleID, conc_FA, PlantFamily, PlantSpeciesfull, treatment)

dataNL_SoiCon %>% 
  add_row(sampleID = "Sc1", conc_FA= 0, treatment = "control") %>% 
  group_by (treatment) %>% 
  summarise (mean = mean (conc_FA), Sd = sd (conc_FA))




# PL Soil controls cleanup ###
PL_SoiCon<- path_controlPL%>%
  excel_sheets()%>%
  set_names() %>% 
  map_df(~read_excel(path=path_controlPL, sheet = .x), .id= "sheet")  %>%  
  dplyr::rename("sampleID"="sheet") %>%   
  mutate (sampleID = str_replace_all (sampleID, "PL_", "")) %>% 
  add_column (PL_NL= "PLFA") %>% 
  select(sampleID, Name, area, PL_NL)

# add Biomarker info to data
PL_SoiCon<-merge(PL_SoiCon, Biomarker %>%  select(-c(Biom1Soil, Soil_A1_kons, MW, RI)), by="Name")

# add meta
PL_SoiCon<-left_join(PL_SoiCon, meta_control, by="sampleID")

# Calculate the loss of sample and correct areas of all FAMEs for loss
correction<- PL_SoiCon %>% select(c(sampleID,area, Name))%>% 
  filter(Name=="19:0")%>%
  mutate(IS_area=1*area) %>%  
  select(c(sampleID,IS_area )) 


PL_FA_conc_SoiCon <-  left_join(PL_SoiCon, correction, by ="sampleID")%>% 
  mutate (conc_FA = 5* (area/IS_area) * (1/RRF) * (1/m_soil))  %>%  ## correction for sample loss via IS, and use of biomarker RR factors to calculate conc of each FA in nmol pe  g soil sample
  filter (Name != "19:0")   # remove IS from dataset


dataPL_SoiCon <- 
  PL_FA_conc_SoiCon %>% 
  filter (Biom2Soil !="IS") %>%
  filter (Biom2Soil != "NAFA") %>%  
  group_by(sampleID, Biom2Soil) %>%  
  summarise (group = sum (conc_FA)) %>%  # get the sum per bacterial group
  pivot_wider(names_from = Biom2Soil, values_from = group ) %>% 
  replace (is.na (.),0) %>% 
  mutate (totalPLFA = Actinobacteria +AMF+ bacteria+Fungi+gram.neg+gram.pos) %>% 
  mutate (totalBact = Actinobacteria + bacteria+gram.neg+gram.pos) # total bacterial biomass

dataPL_SoiCon %>% 
  add_column(treatment="PLFA") %>% 
  group_by(treatment) %>% 
  summarise (mean = mean (totalBact), Sd = sd (totalBact))
                            