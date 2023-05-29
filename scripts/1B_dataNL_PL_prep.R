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

# Table S  dada resuls ####

dada <- read_csv("results/track_reads_M0_PP.csv")


dada_E1_8Sp <- meta_M1 %>%  select (PlaSpe) %>% unique  () %>% left_join (dada, by = "PlaSpe")  %>%
  tibble () %>%  select (-sortedPlantSpecies, - PlantFamily, -PlantType, -value)

# all reads ##

dada_E1_8Sp %>%  summarize (sum (input))

dada_E1_8Sp %>%  summarize (sum (nonchim))



## get number of ASVs,  etc  ##

ps_ALL_E1 <- readRDS("data/ps_ALL_E1.rds") ## needs ps of experiment 1 

ps_ALL_E1_8Sp <- subset_samples (ps_ALL_E1, PlaSpe == "AchMil"| PlaSpe == "CicInt" |
                                   PlaSpe == "SchAru"| PlaSpe == "PoaCit"|
                                   PlaSpe == "PlaLan"| PlaSpe == "HolLan"|
                                   PlaSpe == "BroWil"| PlaSpe == "AgrCap") # filterd for the plants f interest
ps_ALL_E1_8Sp  <-  prune_taxa (taxa_sums(ps_ALL_E1_8Sp)>=1, ps_ALL_E1_8Sp) ## ASVs removed that are not found in those 8 species

## 
estimate_richness(ps_ALL_E1_8Sp) %>% summarize (n = mean (Observed),sd = sd (Observed))  # to get average ASV per smalple
estimate_richness(ps_ALL_E1_8Sp) %>%  arrange( Observed)

# `rarefaction-ALL, fig.cap = "A) Rarefaction (solid line segments) and extrapolation (dotted line segments) 
#sampling curves of fungal ASV richness for each plant species and soil. B) Sample completeness curves by sequence reads. c) Sample coverage curves. Shaded areas in each plot show the 95% confidence interval of the extrapolation. ", fig.height=8, fig.width=12, dpi = 300 }
#agglomerate to Plant Species
ps_ALL_E1_8Sp_perPlaSpe  <- merge_samples(ps_ALL_E1_8Sp, group = "PlantSpeciesfull")

#make df for Glo_PlaSpe data (empty samples samples pruned )
comm_df_ALL_E1_8Sp_PlaSpe  <- data.frame (otu_table (ps_ALL_E1_8Sp_perPlaSpe))### needs ps with all asvs only 8plants)) 

#rarefaction in iNEXT 
iNext_ALL_E1_8_sp_PlaSpe_richness <- iNEXT (t (comm_df_ALL_E1_8Sp_PlaSpe), 
                                    q= 0, # richness
                                    datatype =  "abundance", 
                                    endpoint = NULL, 
                                    knots = 400 )
# plot 
p1 <- ggiNEXT(iNext_ALL_E1_8_sp_PlaSpe_richness, type = 1) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_y_continuous (name = "ASV richness") +
  scale_x_continuous(name = "Number of sequence reads")+
  theme (axis.text = element_text(size =14), 
         axis.title = element_text(size =18),
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))


p2 <- ggiNEXT(iNext_ALL_E1_8_sp_PlaSpe_richness, type = 2) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_x_continuous(name = "Number of sequence reads")+
  theme (axis.text = element_text(size =14), 
         axis.title = element_text(size =18), 
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))

p3 <- ggiNEXT(iNext_ALL_E1_8_sp_PlaSpe_richness, type = 3) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_y_continuous(name = "ASV richness")+
  scale_x_continuous (labels = scales::percent) + 
  theme (axis.text = element_text(size =14), 
         axis.title = element_text(size =18),  
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))

ggarrange (p1,p2,p3, ncol = 3, common.legend = T, labels = "AUTO", legend = "bottom")

### AMF families description ####
Glo_Fam_perc_reads_E1_8Sp  <- 
  ASV_table_Glo  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  group_by (Family)  %>% 
  summarize(Fam_sum=sum(counts)) %>% 
  add_column (Glo_total_counts)  %>% 
  mutate (Glo_Fam_perc = 100  * Fam_sum/total)  %>% 
  arrange(desc(Fam_sum))  %>% 
  as.data.frame()


                                