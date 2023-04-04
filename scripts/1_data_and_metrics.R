## packages ##


library (tidyverse)
library(readxl)
library(phyloseq)
library (vegan)
library (picante)
# data# ####




#meta
metaM0 <- read_xlsx("data/M0_meta.xlsx")

meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1 <- read.csv ("data/meta_M1.csv")
ps_Glo <- read_rds("resultsCh2/ps_Glo.rds")

# change meta

sample_data (ps_Glo) <- sample_data (metaM0 %>% data.frame (row.names = "sampleID"))


# get only 8 species 

ps_Glo_E1_8Sp <- subset_samples (ps_Glo, PlaSpe == "AchMil" | PlaSpe == "CicInt"| PlaSpe == "AgrCap"|
                  PlaSpe == "BroWil"| PlaSpe == "PoaCit"| PlaSpe == "HolLan"| 
                  PlaSpe == "PlaLan"| PlaSpe  == "SchAru")

ps_Glo_E1_8Sp <- prune_taxa (taxa_sums (ps_Glo_E1_8Sp) > 1, ps_Glo_E1_8Sp)


# E1 -  All metrics calculation ####
#### Glo Table ASVs##
ASV_table_Glo  <- 
  otu_table (ps_Glo_E1_8Sp) %>% 
  data.frame() %>%  
  as_tibble (rownames = "ASV_ID") %>%   
  left_join ((tax_table (ps_Glo_E1_8Sp) %>% data.frame() %>%  as_tibble (rownames = "ASV_ID")), by = "ASV_ID")


# means ##
#mean ASVs per sample####

meannumberASVspersample   <-
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
  arrange(meannumberASVs) 

#mean ASVs per plant species####
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

# total ASVs in dataset####
totalASVs  <-ASV_table_Glo %>%  tally()

##Unique --- unique ASVs per plant species ####
uniqueASVsperPlaSpe   <-  ASV_table_Glo  %>% 
  select(-c(Kingdom, Phylum, Class, Family, Order, Genus, Species))  %>% 
  gather (!ASV_ID, key=sampleID, value = ASV_counts)  %>%
  left_join(metaM0, by="sampleID") %>% 
  filter (ASV_counts!=0)  %>% 
  group_by(PlantSpeciesfull, ASV_ID)  %>% 
  tally ()  %>% 
  group_by(PlantSpeciesfull)  %>% 
  tally (name ="uniqueASVsperPlSpe")

## Richness, Shannon, Chao #####

adiv_richness  <- estimate_richness(ps_Glo_E1_8Sp, measures = c("Observed",  "Shannon", "Chao1")) %>% 
  as_tibble (rownames = "sampleID")



## CU - compositional unity  ####
#repl after selecting Glo ###
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

comp_units   <- meannumberASVsper_species  %>% 
  left_join(uniqueASVsperPlaSpe)  %>% 
  left_join (replicates_samples_after_Glo) %>% 
  group_by (PlantSpeciesfull) %>% 
  mutate (unique = mean(uniqueASVsperPlSpe)) %>% 
  mutate (CUnits = 1* unique/(mean_n_ASV_per_species*repl)) %>% # repl adjusted for , times 5 to get back to number of 5 replicates 
  mutate ("1-CU" = 1 - CUnits) 


#ßeta core ####

# taxa core
ps_Glo_core  <- phylosmith::taxa_core (ps_Glo_E1_8Sp, treatment = "PlaSpe", frequency =0.6)#see above
ps_core_PlaSpe  <- merge_samples(ps_Glo_core, group = "PlantSpeciesfull", fun = sum)

#calculate beta diversity as percentage 
cores  <- 
  as.data.frame (otu_table (ps_core_PlaSpe )) %>%   as_tibble(rownames = "PlantSpeciesfull") %>% 
  pivot_longer (!PlantSpeciesfull, names_to="ASV_ID", values_to  = "ASV_counts")  %>%  
  group_by (PlantSpeciesfull) %>%  
  dplyr::filter (ASV_counts != "0")  %>% 
  dplyr::tally(name="core") %>% 
  left_join(uniqueASVsperPlaSpe)  %>% 
  mutate (perc_core = 100 -(100* core/uniqueASVsperPlSpe)) %>% 
  mutate (perc_core = round (perc_core,1) ) %>% 
  select (PlantSpeciesfull, perc_core)




## PD####


#make df for Glo data (empty samples samples pruned )
comm_df_Glo  <- data.frame (t(otu_table (ps_Glo_E1_8Sp))) # vegan expects samples as rows and ASVs species as columns

#hier total randomised community data frame as null model
stand_pd_Glo_all  <- as_tibble (ses.pd(comm_df_Glo, phy_tree(ps_Glo_E1_8Sp), include.root = FALSE, null.model = "independentswap", runs=100, iterations=999), rownames="sampleID")

#In addition, the observed and standardised Faith's PD~PlaSpe~ were calculated for each plant species after agglomeration of the individual samples into one sample per plant species. Comparison of the agglomerated Faith's PD PlaSpe and Faith's PD I reveal ... yeah, the same as in ordination actually.

#agglomerate to Plant Species
ps_Glo_perPlaSpe  <- merge_samples(ps_Glo_E1_8Sp, group = "PlantSpeciesfull")

#make df for Glo_PlaSpe data (empty samples samples pruned )
comm_df_Glo_PlaSpe  <- data.frame (otu_table (ps_Glo_perPlaSpe)) # vegan expects samples as rows and ASVs species as columns

stand_pd_Glo_PlaSpe  <- as_tibble (ses.pd(comm_df_Glo_PlaSpe, phy_tree(ps_Glo_E1_8Sp), 
                                          include.root = FALSE, null.model = "independentswap", 
                                          runs=1000, iterations=999), rownames="sampleID")


## MPD ####
dist_Glo  <- cophenetic(phy_tree(ps_Glo_E1_8Sp))
ses.MPD_Glo  <-  as_tibble (ses.mpd (comm_df_Glo, dist_Glo, null.model = "independentswap" ), rownames = "sampleID")

#agglomerate to Plant Species
ps_Glo_perPlaSpe  <- merge_samples(ps_Glo_E1_8Sp, group = "PlantSpeciesfull")

#make df for Glo_PlaSpe data (empty samples samples pruned )
comm_df_Glo_PlaSpe  <- data.frame (otu_table (ps_Glo_perPlaSpe)) # vegan expects samples as rows and ASVs species as columns

ses.MPD_Glo_PlaSpe  <-  as_tibble (ses.mpd (comm_df_Glo_PlaSpe, dist_Glo, null.model = "independentswap" ), rownames = "sampleID")# tree is the same as before

# All metrics E1 ####


# all 20 species ##
All_metrics_E1 <- 
adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull)  %>% 
  mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Chao1, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(cores) %>% 
  select (!c(PlantType, Chao1, unique))  %>% 
  left_join (stand_pd_Glo_all %>% 
               left_join(ses.MPD_Glo, by ="sampleID")  %>%  
               left_join (metaM0)  %>% 
               drop_na ()  %>% 
               group_by (PlantSpeciesfull)  %>%   
               mutate (mean_pd = mean (pd.obs), mean_mpd = mean (mpd.obs)) %>%  
               select (PlantSpeciesfull,   mean_pd, mean_mpd, -PlantFamily)   %>% 
               unique ())  %>% 
  dplyr::rename("Richness"= "mean_n_ASV_per_species", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd", "uniqueASV" = "uniqueASVsperPlSpe" )  %>% 
  add_column (Exp = "E1")


# df of  
All_metrics_E1_df <- 
  All_metrics_E1 %>% 
  as.data.frame(row.names = NULL)

rownames(All_metrics_E1_df)  <- All_metrics_E1_df$PlantSpeciesfull
All_metrics_E1_df  <- All_metrics_E1_df[,-1]


# All metrics E1 - 8 species #####
All_metrics_E1_8Sp  <-
 All_metrics_E1 %>% 
  right_join(meta_plants %>%  filter (PlaSpe != "soil") %>% 
               select (PlantSpeciesfull) , by = "PlantSpeciesfull") 

# df 
All_metrics_E1_8Sp_df <- 
  All_metrics_E1_8Sp %>% 
  as.data.frame(row.names = NULL)

rownames(All_metrics_E1_8Sp_df)  <- All_metrics_E1_8Sp_df$PlantSpeciesfull
All_metrics_E1_8Sp_df  <- All_metrics_E1_8Sp_df[,-1]

# All metrics E1 per sample ####

All_metrics_E1_samples <- 
  adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull, PlaSpe)  %>% 
  select (sampleID, Chao1, Shannon) %>% 
  left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(cores) %>% 
  select (!c(PlantType, Chao1, unique))  %>% 
  left_join (stand_pd_Glo_all %>% select (sampleID, pd.obs)) %>% 
  left_join (ses.MPD_Glo %>%  select (sampleID, mpd.obs) )  %>% 
  dplyr::rename("Richness"= "mean_n_ASV_per_species", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "pd.obs", "MPD" = "mpd.obs", "uniqueASV" = "uniqueASVsperPlSpe" )  %>% 
  add_column (Exp = "E1") 



# E2 - all metrics calculation ######
# data #
ps_M1_roots_Glo_rar <- readRDS("data/ps_M1_roots_Glo_rar.rds")


## richness, chao, shannon ####
adiv_richness_M1  <- 
  estimate_richness(ps_M1_roots_Glo_rar, measures = c("Observed", "Shannon")) %>% 
  as_tibble (rownames = "sampleID") %>% 
  left_join (meta_M1 %>%  select (sampleID, PlaSpe, PlantSpeciesfull, PlantFamily), by = "sampleID")  


## richness per species -average
meannumberASVsper_species_M1  <-
  adiv_richness_M1 %>% 
  group_by (PlantSpeciesfull) %>% 
  summarize (meanASV = mean (Observed))


ASV_table_Glo_roots_M1 <- 
  otu_table (ps_M1_roots_Glo_rar) %>%  
  data.frame ()  %>%  t () %>%  
  as_tibble (rownames = "sampleID") 

## unique asvs ####
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

ps_test <- ps_M1_roots_Glo_rar  # make a new ps to change its sample_data (core and conglomerte have problems when sample-data can not conglomerated)
sample_data (ps_test)  <- sample_data (ps_M1_roots_Glo_rar) %>% data.frame () %>%  select (PlantSpeciesfull, PlaSpe) %>%  sample_data()
ps_focal_roots_M1_core  <- phylosmith::taxa_core (ps_test, treatment =  "PlantSpeciesfull", frequency =0.6)#see above
ps_focal_roots_M1_core <- phylosmith::conglomerate_samples(ps_focal_roots_M1_core, treatment = "PlantSpeciesfull")

#calculate beta diversity as percentage 
cores_roots_M1  <- 
  as.data.frame (otu_table (ps_focal_roots_M1_core )) %>%   as_tibble(rownames = "ASV_ID") %>% 
  gather (!ASV_ID, key=PlantSpeciesfull, value = ASV_counts)  %>%  
  group_by (PlantSpeciesfull) %>%  
  dplyr::filter (ASV_counts != "0")  %>% 
  dplyr::tally(name="core") %>% 
  left_join(unique_M1, by = "PlantSpeciesfull")  %>% 
  mutate (perc_core = 100- (100* core/uniqueASVsperPlSpe)) %>% 
  mutate (perc_core = round (perc_core,1) ) %>% 
  select (PlantSpeciesfull, perc_core) 

## PD####
comm_df_Glo_M1 <- data.frame (t(otu_table (ps_M1_roots_Glo_rar))) # vegan expects samples as rows and ASVs species as columns

#hier total randomised community data frame as null model
stand_pd_Glo_all_M1  <- as_tibble (ses.pd(comm_df_Glo_M1, phy_tree(ps_M1_roots_Glo_rar), include.root = FALSE, null.model = "independentswap", runs=100, iterations=999), rownames="sampleID")

stand_pd_Glo_all_M1  <- 
  stand_pd_Glo_all_M1 %>%  
  left_join (meta_M1 %>% select (sampleID, PlantSpeciesfull)) 

## MPD ####

dist_Glo_M1  <- cophenetic(phy_tree(ps_M1_roots_Glo_rar))
ses.MPD_Glo_M1  <-  as_tibble (ses.mpd (comm_df_Glo_M1, dist_Glo_M1, null.model = "independentswap" ), rownames = "sampleID")

#agglomerate to Plant Species
ps_Glo_perPlaSpe_M1  <- merge_samples(ps_M1_roots_Glo_rar, group = "PlantSpeciesfull")

#make df for Glo_PlaSpe data (empty samples samples pruned )
comm_df_Glo_PlaSpe_M1  <- data.frame (otu_table (ps_Glo_perPlaSpe_M1)) # vegan expects samples as rows and ASVs species as columns

ses.MPD_Glo_PlaSpe_M1  <-  as_tibble (ses.mpd (comm_df_Glo_PlaSpe_M1, dist_Glo_M1, null.model = "independentswap" ), rownames = "sampleID")# tree is the same as before


# together PD, MPD
PD_MPD_M1 <-  
  stand_pd_Glo_all_M1 %>%  
  select (sampleID, pd.obs) %>% 
  left_join(ses.MPD_Glo_M1, by = "sampleID") %>%  
  left_join (meta_M1) %>% 
  select (sampleID, pd.obs, mpd.obs, PlaSpe, PlantSpeciesfull) %>% 
  drop_na() %>% 
  group_by (PlantSpeciesfull) %>% 
  mutate (mean_pd = mean (pd.obs), mean_mpd = mean (mpd.obs)) %>%  
  select (PlantSpeciesfull,  mean_pd, mean_mpd)   %>% 
  unique ()  





# All metrics E2  ####
#Get values for all metrics in one tibble

All_metrics_E2 <-
  adiv_richness_M1 %>% 
  group_by(PlantFamily, PlantSpeciesfull)  %>% 
  #mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(cores_roots_M1) %>% 
  left_join (PD_MPD_M1) %>% 
  dplyr::rename("Richness"= "meanASV", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd" , "uniqueASV" = "uniqueASVsperPlSpe") %>% 
  add_column (Exp = "E2")


  # df 
All_metrics_E2_df  <- All_metrics_E2 %>% 
    as.data.frame(row.names = NULL) %>% 
  relocate (where (is.numeric))

rownames(All_metrics_E2_df)  <- All_metrics_E2_df$PlantSpeciesfull

## All metrics per sample  ####

All_Metrics_E2_sample  <- 
  adiv_richness_M1 %>% 
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(cores_roots_M1) %>% 
  left_join (stand_pd_Glo_all_M1 %>%  select (sampleID, pd.obs)) %>%
  left_join (ses.MPD_Glo_M1 %>%  select (sampleID, mpd.obs) )  %>% 
  select (!meanASV) %>% 
  dplyr::rename("Richness" = "Observed","CU" = "CUnits", "β(core)"  = "perc_core",
                 "uniqueASV" = "uniqueASVsperPlSpe", "MPD" = "mpd.obs", "PD" = "pd.obs") %>% 
  add_column (Exp = "E2")

saveRDS(All_Metrics_E2_sample, "results/All_metrics_E2_sample.rds")

