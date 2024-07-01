# Part 4 Calculate diversity metrics #######

## packages ####
library (tidyverse)
library(readxl)
library(phyloseq)
library (vegan)
library (picante)
library (readxl)
library (rstatix)
library(ggpubr)
library (phylosmith)
library(metagMisc)

# Load data #####
# These objects include only the data for AMF , all other fungal sequences have been removed##
ps_E1 <- read_rds ("data/ps_E1.rds")  %>% prune_taxa( taxa_sums (.)>0, .) %>% prune_samples(sample_sums (.) > 0, .)
ps_E2 <- read_rds ("data/ps_E2.rds") %>% prune_taxa( taxa_sums (.)>0, .) %>% prune_samples(sample_sums (.) > 0, .) ##



# Load meta tables
metaM0 <- read_xlsx("data/M0_meta.xlsx")
meta_M1Wcontrol <- read_xlsx ("data/meta_M1Wcontrol.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1 <- read.csv ("data/meta_M1.csv")



# change meta in phyloseq object #

sample_data (ps_E1) <- sample_data (metaM0 %>% data.frame (row.names = "sampleID"))


# E1 -  All metrics calculation ####
#### Get a ASV table (all AMF ASVs)##
ASV_table_Glo  <- 
  otu_table (ps_E1) %>% 
  data.frame() %>%  
  as_tibble (rownames = "ASV_ID") %>%   
  left_join ((tax_table (ps_E1) %>% data.frame() %>%  as_tibble (rownames = "ASV_ID")), by = "ASV_ID")

## mean Richness per plant species####
meannumberASVsper_species   <-
  ASV_table_Glo  %>% 
  select(-c(Kingdom, Phylum, Class, Family, Order, Genus, Species))  %>% 
  gather (!ASV_ID, key=sampleID, value = ASV_counts)  %>%
  left_join(metaM0, by="sampleID") %>% 
  filter (ASV_counts!=0) %>% 
  group_by (sampleID) %>%  
  tally(name="numberASVs") %>% 
  left_join(metaM0, by="sampleID") %>% 
  group_by (PlaSpe) %>% 
  mutate(meannumberASVs=mean(numberASVs)) %>% 
  arrange(meannumberASVs)   %>% 
  select (-c(sampleID, numberASVs)) %>% 
  group_by (PlantSpeciesfull, PlantFamily, PlantType)  %>% 
  summarize (mean_n_ASV_per_species  = mean(meannumberASVs))

# total ASVs in dataset###
totalASVs  <-ASV_table_Glo %>%  tally()

##y- diversity = unique ASVs per plant species ####
uniqueASVsperPlaSpe   <-  ASV_table_Glo  %>% 
  select(-c(Kingdom, Phylum, Class, Family, Order, Genus, Species))  %>% 
  gather (!ASV_ID, key=sampleID, value = ASV_counts)  %>%
  left_join(metaM0, by="sampleID") %>% 
  filter (ASV_counts!=0)  %>% 
  group_by(PlantSpeciesfull, ASV_ID)  %>% 
  tally ()  %>% 
  group_by(PlantSpeciesfull)  %>% 
  tally (name ="uniqueASVsperPlSpe")




## Richness S, Shannon' H' #####

adiv_richness  <- estimate_richness(ps_E1, measures = c("Observed",  "Shannon")) %>% 
  as_tibble (rownames = "sampleID")



## beta (CU) - compositional units  ####
#Number of replicates remaining after filtering of sequences gotr Glomeromycotoina is calculated ###
replicates_samples_after_Glo  <-
  ASV_table_Glo  %>% 
  select(-c(Kingdom, Phylum, Class, Family, Order, Genus, Species))  %>% 
  gather (!ASV_ID, key=sampleID, value = ASV_counts)  %>%
  left_join(metaM0, by="sampleID") %>% 
  filter (ASV_counts!=0) %>% 
  group_by (sampleID) %>%  
  tally(name="numberASVs") %>% 
  left_join(metaM0, by="sampleID") %>%  group_by(PlantSpeciesfull) %>%  tally() %>% 
  dplyr:: rename ("repl"= "n")

## beta (CU) ##
comp_units   <- 
  meannumberASVsper_species  %>% 
  left_join(uniqueASVsperPlaSpe)  %>% 
  left_join (replicates_samples_after_Glo) %>% 
  group_by (PlantSpeciesfull) %>% 
  mutate (unique = mean(uniqueASVsperPlSpe)) %>% 
  mutate (CUnits = 1* unique/(mean_n_ASV_per_species*repl)) %>% # repl adjusted for , times 5 to get back to number of 5 replicates ## I don't think that should be used ### check and correct ###
  mutate ("1-CU" = 1 - CUnits) 


##ß core ####

data_E1_melted <- melt_phyloseq(ps_E1) # make phyloseq to data table

### Test different cutoffs for beta (core) for Supplementary material ####
ASVs_in_1_samples <-
  data_E1_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  left_join (replicates_samples_after_Glo) %>% # add information on number of replicates
  filter (n/repl>=0.2) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()


ASVs_in_2_samples <-
  data_E1_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  left_join (replicates_samples_after_Glo) %>% 
  filter (n/repl>=0.4) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()


ASVs_in_3_samples <-
data_E1_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>%   
  left_join (replicates_samples_after_Glo) %>% 
  filter (n/repl>=0.6) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

ASVs_in_4_samples <-
  data_E1_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  left_join (replicates_samples_after_Glo) %>% 
  filter (n/repl>=0.8) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

ASVs_in_5_samples <-
  data_E1_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>%   
  left_join (replicates_samples_after_Glo) %>% 
  filter (n/repl>=1) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

beta_core_E1_1 <- 
  uniqueASVsperPlaSpe %>% 
  left_join (ASVs_in_1_samples) %>% 
  mutate (beta_core1 = round (1-(n/uniqueASVsperPlSpe),2))

beta_core_E1_2 <- 
  uniqueASVsperPlaSpe %>% 
  left_join (ASVs_in_2_samples) %>% 
  mutate (beta_core2 =round ( 1-(n/uniqueASVsperPlSpe), 2))

beta_core_E1_3 <- 
  uniqueASVsperPlaSpe %>% 
  left_join (ASVs_in_3_samples) %>% 
  mutate (beta_core3 = round (1-(n/uniqueASVsperPlSpe), 2))

beta_core_E1_4 <- 
  uniqueASVsperPlaSpe %>% 
  left_join (ASVs_in_4_samples) %>% 
  mutate (beta_core4 = round (1-(n/uniqueASVsperPlSpe),2))

beta_core_E1_5 <- 
  uniqueASVsperPlaSpe %>% 
  left_join (ASVs_in_5_samples) %>% 
  mutate (beta_core5 = round (1-(n/uniqueASVsperPlSpe),2))

### Table for supplementary material ####
beta_core_E1_1 %>% 
  left_join(beta_core_E1_2 %>%  select ( PlantSpeciesfull, beta_core2)) %>% 
  left_join(beta_core_E1_3 %>%  select ( PlantSpeciesfull, beta_core3 )) %>% 
  left_join(beta_core_E1_4 %>%  select ( PlantSpeciesfull, beta_core4)) %>% 
  left_join(beta_core_E1_5 %>%  select ( PlantSpeciesfull, beta_core5)) %>% 
  write_csv("results/beta_core_all_cutoffs.csv")


## MPD ####
dist_Glo  <- cophenetic(phy_tree(ps_E1))
comm_df_Glo  <- data.frame (t(otu_table (ps_E1))) # vegan expects samples as rows and ASVs species as columns

ses.MPD_Glo  <-  as_tibble (ses.mpd (comm_df_Glo, dist_Glo, null.model = "independentswap" ), rownames = "sampleID")

## phylogenetic γ-diversity ####
#agglomerate to Plant Species
ps_Glo_perPlaSpe  <- merge_samples(ps_E1, group = "PlantSpeciesfull")

# make df for data 
comm_df_Glo_PlaSpe  <- data.frame (otu_table (ps_Glo_perPlaSpe)) # vegan expects samples as rows and ASVs species as columns
# calculate MPD per plant species (gamma-diversity)
ses.MPD_Glo_PlaSpe  <-  as_tibble (ses.mpd (comm_df_Glo_PlaSpe, dist_Glo, null.model = "independentswap" ), rownames = "sampleID")# tree is the same as before

# Unifrac #####

Unifrac_E1 <- 
  ps_E1 %>%  phyloseq_sep_variable("PlantSpeciesfull") %>% 
  map (~phyloseq::distance(., method = "uunifrac", type = "samples")) %>% 
  map (~mean (.)) %>% 
  as_tibble () %>% 
  add_column (Exp = "E1") %>% 
  pivot_longer(!Exp ,  values_to = "Unifrac", names_to = "PlantSpeciesfull")

# Table of all diversity metrics ####
# All metrics are put into one table and the mean per plant species is calculated ##
All_metrics_E1 <- 
  adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull)  %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID,  Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(beta_core_E1_3 %>% select (!n)) %>% 
  left_join(Unifrac_E1) %>% 
  select (!c(PlantType,  unique))  %>% 
  left_join (ses.MPD_Glo  %>%  
               left_join (metaM0)  %>% 
               drop_na ()  %>% 
               group_by (PlantSpeciesfull)  %>%   
               mutate ( mean_mpd = mean (mpd.obs)) %>%  
               select (PlantSpeciesfull,   mean_mpd, -PlantFamily)   %>% 
               unique ())  %>% 
  left_join (ses.MPD_Glo_PlaSpe %>%  select (sampleID, mpd.obs), by = c("PlantSpeciesfull" = "sampleID")) %>% 
  dplyr::rename("Richness S"= "mean_n_ASV_per_species","Shannon's H'" = "Shannon", "β(CU)" = "CUnits", "β(core)"  = "beta_core3",
                 "MPD" = "mean_mpd", "γ-diversity" = "uniqueASVsperPlSpe", "Phyl. γ-diversity" = "mpd.obs" ) 

write_csv(All_metrics_E1, "results/20240501_All_metrics_E1.csv")

# dataframe of the table #
All_metrics_E1_df <- 
  All_metrics_E1 %>% 
  as.data.frame(row.names = NULL)

rownames(All_metrics_E1_df)  <- All_metrics_E1_df$PlantSpeciesfull
All_metrics_E1_df  <- All_metrics_E1_df[,-1]


# All metrics E1 - without soil  ###
All_metrics_E1  <-
  All_metrics_E1 %>% 
  right_join(meta_plants %>%  filter (PlaSpe != "soil") %>% 
               select (PlantSpeciesfull) , by = "PlantSpeciesfull") 

# dataframe ##
All_metrics_E1_df <- 
  All_metrics_E1 %>% 
  as.data.frame(row.names = NULL)

rownames(All_metrics_E1_df)  <- All_metrics_E1_df$PlantSpeciesfull
All_metrics_E1_df  <- All_metrics_E1_df[,-1]

# All metrics per sample (individual) for E1 ####

All_metrics_E1_samples <- 
  adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull, PlaSpe)  %>% 
  select (sampleID,  Shannon, Observed) %>% 
  left_join (comp_units %>% select (!c(repl,mean_n_ASV_per_species )) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(beta_core_E1_3 %>% select (!n)) %>% 
  left_join(Unifrac_E1) %>% 
  select (!c(PlantType,  unique))  %>% 
  left_join (ses.MPD_Glo %>%  select (sampleID, mpd.obs) )  %>% 
  left_join (ses.MPD_Glo_PlaSpe %>%  select (sampleID, mpd.obs), by = c("PlantSpeciesfull" = "sampleID")) %>% 
  dplyr::rename("Richness S" = "Observed","Shannon's H'" = "Shannon", "β(CU)" = "CUnits",  "β(core)"  = "beta_core3",
                "γ-diversity" = "uniqueASVsperPlSpe", "MPD" = "mpd.obs.x",   "Phyl. γ-diversity" = "mpd.obs.y") 
  
write_csv (All_metrics_E1_samples, "results/20240501_All_metrics_E1_sample.csv")




# E2 - all metrics calculation ######

## Richness S,  Shannon's H####
adiv_richness_M1  <- 
  estimate_richness(ps_E2, measures = c("Observed", "Shannon")) %>% 
  as_tibble (rownames = "sampleID") %>% 
  left_join (meta_M1 %>%  select (sampleID, PlaSpe, PlantSpeciesfull, PlantFamily), by = "sampleID")  


## Mean richness per plant  species ####
meannumberASVsper_species_M1  <-
  adiv_richness_M1 %>% 
  group_by (PlantSpeciesfull) %>% 
  summarize (meanASV = mean (Observed))


ASV_table_Glo_roots_M1 <- 
  otu_table (ps_E2) %>%  
  data.frame() %>%  t () %>%  
  as_tibble (rownames = "sampleID") 

ASV_table_Exp2_Glo  <- 
  otu_table (ps_E2) %>%  
  data.frame() %>%  
  as_tibble (rownames = "ASV_ID") %>%   
  left_join ((tax_table (ps_E2) %>% data.frame() %>%  as_tibble (rownames = "ASV_ID")), by = "ASV_ID")


## y-diversity ####
unique_M1  <- 
  ASV_table_Glo_roots_M1  %>% 
  pivot_longer(!sampleID, names_to = "ASV_ID", values_to = "ASV_count") %>% 
  filter (ASV_count != 0) %>%  
  left_join(meta_M1 %>%  select (sampleID, PlantSpeciesfull)) %>%  
  group_by (ASV_ID,PlantSpeciesfull ) %>%  
  tally () %>% 
  group_by(PlantSpeciesfull)  %>% 
  tally (name ="uniqueASVsperPlSpe")






## CU ####
comp_units_M1   <- meannumberASVsper_species_M1  %>% 
  left_join(unique_M1, by = "PlantSpeciesfull")  %>% 
  mutate (CUnits = uniqueASVsperPlSpe/(meanASV *5))%>% # repl adjustment nor needed here, because all 5 replicates are available
  mutate ("1-CU" = 1 - CUnits)


# meannumberASVsper_species  %>% 
#   left_join(uniqueASVsperPlaSpe)  %>% 
#   left_join (replicates_samples_after_Glo) %>% 
#   group_by (PlantSpeciesfull) %>% 
#   mutate (unique = mean(uniqueASVsperPlSpe)) %>% 
#   mutate (CUnits = 5* unique/(mean_n_ASV_per_species*repl)) %>% # repl adjusted for , times 5 to get back to number of 5 replicates 
#   mutate ("1-CU" = 5 - CUnits)

## Cores####

# ps_test <- ps_E2  # make a new ps to change its sample_data to get taxa_core and conglomerate_samples to work
# sample_data (ps_test)  <- sample_data (ps_E2) %>% data.frame () %>%  select (PlantSpeciesfull) %>%  sample_data()
# ps_focal_roots_M1_core  <- 
#   phylosmith::taxa_core (ps_test, treatment =  "PlantSpeciesfull", 
#                          frequency =0.6) #see above
# ps_focal_roots_M1_core <- 
#   phylosmith::conglomerate_samples(ps_focal_roots_M1_core, treatment = "PlantSpeciesfull")  %>%  
#     prune_taxa( taxa_sums (.)>0, .) %>% prune_samples(sample_sums (.) > 0, .)

data_E2_melted <- melt_phyloseq(ps_E2)

ASVs_in_1_samples <-
  data_E2_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  filter (n>=1) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()


ASVs_in_2_samples <-
  data_E2_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  filter (n>=2) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()


ASVs_in_3_samples <-
  data_E2_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>%   
  filter (n>=3) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

ASVs_in_4_samples <-
  data_E2_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>% 
  left_join (replicates_samples_after_Glo) %>% 
  filter (n>=4) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

ASVs_in_5_samples <-
  data_E2_melted %>% 
  as_tibble () %>% 
  filter (Abundance >0) %>% 
  group_by (PlantSpeciesfull, OTU) %>% 
  count () %>%   
  left_join (replicates_samples_after_Glo) %>% 
  filter (n>=5) %>%  ## Change here 
  select (PlantSpeciesfull, n) %>% 
  group_by (PlantSpeciesfull) %>% 
  count ()

beta_core_E2_1 <- 
  unique_M1 %>% 
  left_join (ASVs_in_1_samples) %>% 
  mutate (beta_core1 = round (1-(n/uniqueASVsperPlSpe),2))

beta_core_E2_2 <- 
  unique_M1 %>% 
  left_join (ASVs_in_2_samples) %>% 
  mutate (beta_core2 =round ( 1-(n/uniqueASVsperPlSpe), 2))

beta_core_E2_3 <- 
  unique_M1 %>% 
  left_join (ASVs_in_3_samples) %>% 
  mutate (beta_core3 = round (1-(n/uniqueASVsperPlSpe), 2))

beta_core_E2_4 <- 
  unique_M1 %>% 
  left_join (ASVs_in_4_samples) %>% 
  mutate (beta_core4 = round (1-(n/uniqueASVsperPlSpe),2))

beta_core_E2_5 <- 
  unique_M1 %>% 
  left_join (ASVs_in_5_samples) %>% 
  mutate (beta_core5 = round (1-(n/uniqueASVsperPlSpe),2))

beta_core_E2_1 %>% 
  left_join(beta_core_E2_2 %>%  select ( PlantSpeciesfull, beta_core2)) %>% 
  left_join(beta_core_E2_3 %>%  select ( PlantSpeciesfull, beta_core3 )) %>% 
  left_join(beta_core_E2_4 %>%  select ( PlantSpeciesfull, beta_core4)) %>% 
  left_join(beta_core_E2_5 %>%  select ( PlantSpeciesfull, beta_core5)) %>% 
  write_csv("results/beta_core_all_cutoffs_E2.csv")


#calculate beta diversity as percentage 
# cores_roots_M1  <- 
#   as.data.frame (otu_table (ps_focal_roots_M1_core )) %>%   as_tibble(rownames = "ASV_ID") %>% 
#   gather (!ASV_ID, key=PlantSpeciesfull, value = ASV_counts)  %>%  
#   group_by (PlantSpeciesfull) %>%  
#   dplyr::filter (ASV_counts != "0")  %>% 
#   dplyr::tally(name="core") %>% 
#   left_join(unique_M1, by = "PlantSpeciesfull")  %>% 
#   mutate (perc_core = 100- (100* core/uniqueASVsperPlSpe)) %>% 
#   mutate (perc_core = round (perc_core,1) ) %>% 
#   select (PlantSpeciesfull, perc_core) 
# 

## MPD ####

dist_Glo_M1  <- cophenetic(phy_tree(ps_E2))
comm_df_Glo_M1  <- data.frame (t(otu_table (ps_E2)) )# vegan expects samples as rows and ASVs species as columns

ses.MPD_Glo_M1  <-  as_tibble (ses.mpd (comm_df_Glo_M1, dist_Glo_M1, null.model = "independentswap" ), rownames = "sampleID")

## phylogenetic γ-diversity ####

#agglomerate to Plant Species to calculate phylogenetic y-diversity #
ps_Glo_perPlaSpe_M1  <- merge_samples(ps_E2, group = "PlantSpeciesfull")

#make df for Glo_PlaSpe data (empty samples  pruned )
comm_df_Glo_PlaSpe_M1  <- data.frame (otu_table (ps_Glo_perPlaSpe_M1)) # vegan expects samples as rows and ASVs species as columns

ses.MPD_Glo_PlaSpe_M1  <-  as_tibble (ses.mpd (comm_df_Glo_PlaSpe_M1, dist_Glo_M1, null.model = "independentswap" ), rownames = "sampleID")# tree is the same as before

## Unifrac #####
## 
Unifrac_E2 <- 
  ps_E2 %>%  phyloseq_sep_variable("PlantSpeciesfull") %>% 
  map (~phyloseq::distance(., method = "uunifrac", type = "samples")) %>% 
  map (~mean (.)) %>% 
  as_tibble () %>% 
  add_column (Exp = "E2") %>% 
  pivot_longer(!Exp ,  values_to = "Unifrac", names_to = "PlantSpeciesfull")






# All metrics E2  ####
#Get values for all metrics per plant species in one tibble

All_metrics_E2 <-
  adiv_richness_M1 %>% 
  group_by(PlantFamily, PlantSpeciesfull)  %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(beta_core_E2_3%>% select (!n)) %>% 
  left_join (Unifrac_E2) %>% 
  left_join (ses.MPD_Glo_M1  %>%  
               left_join (meta_M1) %>% 
               drop_na ()  %>% 
               group_by (PlantSpeciesfull)  %>%   
               mutate ( mean_mpd = mean (mpd.obs)) %>%  
               select (PlantSpeciesfull,   mean_mpd, -PlantFamily)   %>% 
               unique ()) %>% 
  left_join (ses.MPD_Glo_PlaSpe_M1 %>%  select (sampleID, mpd.obs), by = c("PlantSpeciesfull" = "sampleID")) %>% 
  dplyr::rename("Richness S"= "meanASV","Shannon's H'" = "Shannon", "β(CU)" = "CUnits", "β(core)"  = "beta_core3",
                "MPD" = "mean_mpd", "γ-diversity" = "uniqueASVsperPlSpe", "Phyl. γ-diversity" = "mpd.obs" )  


# dataframe of table #
All_metrics_E2_df  <- All_metrics_E2 %>% 
  as.data.frame(row.names = NULL) %>% 
  relocate (where (is.numeric))

rownames(All_metrics_E2_df)  <- All_metrics_E2_df$PlantSpeciesfull

## All metrics per sample  ####

All_Metrics_E2_sample  <- 
  adiv_richness_M1 %>% 
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(beta_core_E2_3 %>% select (!n)) %>% 
  left_join(Unifrac_E2) %>% 
  left_join (ses.MPD_Glo_M1 %>%  select (sampleID, mpd.obs) )  %>% 
  select (!meanASV) %>% 
  left_join (ses.MPD_Glo_PlaSpe_M1 %>%  select (sampleID, mpd.obs), by = c("PlantSpeciesfull" = "sampleID")) %>% 
  dplyr::rename("Richness S" = "Observed","Shannon's H'" = "Shannon", "β(CU)" = "CUnits",  "β(core)"  = "beta_core3",
                "γ-diversity" = "uniqueASVsperPlSpe", "MPD" = "mpd.obs.x",   "Phyl. γ-diversity" = "mpd.obs.y") 

saveRDS(All_Metrics_E2_sample, "results/20240501_All_metrics_E2_sample.rds")


### Tree for the plant species using two genes: matK and trnL #####

# sequences were sourced from NCBI database
# alignment using MAFFT was applied with software JalView
# genes were concatenated
# NJ tree was build in JalView 

# This is how the covariance structure was calculated:
library (phangorn)
## Build tree ####
dnaTRNL<-read.dna("data/Exp2Fasta/matK genes/concat_align.fa",format="fasta") ## these are the aligned sequences

dnaTRNL
labels(dnaTRNL) #check, change
rownames(dnaTRNL)<-c("AchMil" ,"CicInt","PlaLan", "HolLan" ,"PoaCit", "BroWil","SchAru", "AgrCap") #### YOU CAN NAME YOUR PLANTS HERE, i checked those names manually

## Build  tree from plant gene sequences ###
# I used JalView (software) to generate the tree and the pairwise alignment
tree_jalView<-read.tree("data/Exp2Fasta/matK genes/concat_tree_mafft")
tree_jalView$tip.label # check
plot(tree_jalView)
# 
# # rename tree tips with names of plant species (metadata column: PlaSpe)
tree_jalView$tip.label <- c("PlaLan","CicInt","AchMil", "BroWil","PoaCit","SchAru", "AgrCap", "HolLan"  ) # manual chec of names 

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
saveRDS (A, "data/covariance_str_tree_concat_plants.rds") 





## Calculations for the tables in the Supplementary Material ####


## Mean phylogenetic distance for experiment E2 ##
ses.MPD_Glo_M1 %>% 
  mutate (mpd.obs = round (mpd.obs,2) , mpd.rand.mean = round (mpd.rand.mean,2), 
          mpd.rand.sd = round (mpd.rand.sd,2), mpd.obs.z = round (mpd.obs.z,2), 
          mpd.obs.p = round (mpd.obs.p,4)) %>%  
  select (-runs) %>% write_excel_csv("results/mpd_E2.csv")


# All_Metrics_E2_sample %>% ggplot (aes (x= PlantSpeciesfull, y = `Richness S`))  + geom_boxplot () +coord_flip()
# 
# All_Metrics_E2_sample %>% ggplot (aes (x= PlantSpeciesfull, y = MPD))  + geom_boxplot () +coord_flip()
