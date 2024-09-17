# Part 10 Supplementary results ####

library (tidyverse)
library (iNEXT)
library (phyloseq)
library (ggpubr)
library(vegan)


#Appendix S2 ##
ps_ALL_E1_8Sp  <-  readRDS("data/ps_AllASVs_E1.rds")  



#Section S1 #
# Figure S1 rarefaction curves for E1 ####

# A) Rarefaction (solid line segments) and extrapolation (dotted line segments) 
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

write_rds(iNext_ALL_E1_8_sp_PlaSpe_richness, "results/iNext_ALL_E1_8_sp_PlaSpe_richness.rds")


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


ggarrange (p1, p2, p3, common.legend = T, nrow = 1, ncol = 3, legend = "bottom")


# Figure S2 rarefaction curves for E2 ####
ps_roots_allASVs <- readRDS ("data/ps_M1_allASVs.rds")  %>%  
  subset_samples(roots_soil == "roots") %>% 
  prune_taxa( taxa_sums (.)>0, .) %>% 
 prune_samples(sample_sums (.) > 0, .)
# 
ps_roots_M1_allASVs_perPlaSpe  <- merge_samples(ps_roots_allASVs, group = "PlantSpeciesfull")

comm_df_M1_PlaSpe  <- data.frame (otu_table (ps_roots_M1_allASVs_perPlaSpe)) 

# Rarefaction in iNEXT ##
iNext_M1_PlaSpe_richness <- iNEXT (t (comm_df_M1_PlaSpe), 
                                   q= 0, # richness
                                   datatype =  "abundance", 
                                   endpoint = NULL, 
                                   knots = 400 )
# Plot rarefaction curve ###
p1 <- ggiNEXT(iNext_M1_PlaSpe_richness, type = 1) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_y_continuous (name = "ASV richness") +
  scale_x_continuous(name = "Number of sequence reads")+
  theme (axis.text = element_text(size =16), 
         axis.title = element_text(size =20),
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))

p2 <- ggiNEXT(iNext_M1_PlaSpe_richness, type = 2) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_x_continuous(name = "Number of sequence reads")+
  theme (axis.text = element_text(size =16), 
         axis.title = element_text(size =20),
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))

p3 <- ggiNEXT(iNext_M1_PlaSpe_richness, type = 3) +
  theme_classic() + 
  scale_shape_manual(values = c(0,1, 2,3,4,5,6,7)) + 
  scale_y_continuous(name = "ASV richness")+
  scale_x_continuous (labels = scales::percent) + 
  theme (axis.text = element_text(size =16), 
         axis.title = element_text(size =20),
         legend.text = element_text(face = "italic"))+
  guides(     colour=guide_legend(title="Plant species"), 
              fill=guide_legend(title="Plant species"), 
              shape=guide_legend(title="Plant species"))

ggarrange (p1, p2, p3, common.legend = T, nrow = 1, ncol = 3, legend = "bottom")

-------------------------------------------------------------------------------------------------------------


# Appendix S2: Table S2  sequencing results s ####

dada <- read_csv("results/track_reads_E1.csv") %>% as_tibble () %>% 
  select (!...1)#  use the 


dada_E1 <- meta_M1 %>%  select (PlaSpe) %>% unique  () %>% left_join (dada, by = "PlaSpe")  %>%
  tibble () 

# all reads ##

dada_E1 %>%  summarize (sum (input))

dada_E1 %>%  summarize (sum (nonchim))



## get number of ASVs,  etc  ##

ps_ALL_E1 <- readRDS("data/ps_ALL_E1.rds") ## needs ps of experiment 1 or experiment 2, before removal of the non-Glomeromycotan ASVs

## Calculate the average number of ASVs per sample #
estimate_richness(ps_ALL_E1, measures = "Observed") %>% summarize (n = mean (Observed),sd = sd (Observed))  # to get average ASV per smalple
estimate_richness(ps_ALL_E1, measures =  "Observed") %>%  arrange( Observed)



###
# data##
ps_E1 <- read_rds ("data/ps_E1.rds")
ps_E2 <- read_rds ("data/ps_E2.rds")

# Load meta tables
metaM0 <- read_xlsx("data/M0_meta.xlsx")
meta_M1Wcontrol <- read_xlsx ("data/meta_M1Wcontrol.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1 <- read.csv ("data/meta_M1.csv")

sample_data (ps_E1) <- sample_data (metaM0 %>% data.frame (row.names = "sampleID"))


# E1 -  All metrics calculation ####
#### Get a ASV table (all AMF ASVs)##
ASV_table_Glo  <- 
  otu_table (ps_E1) %>% 
  data.frame() %>%  
  as_tibble (rownames = "ASV_ID") %>%   
  left_join ((tax_table (ps_E1) %>% data.frame() %>%  as_tibble (rownames = "ASV_ID")), by = "ASV_ID")


### AMF sequence stats ####

# Total AMF sequence counts in E1##
Glo_total_counts_E1_8Sp  <-
  ASV_table_Glo  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  summarize(total=sum(across(counts)))


#Percentage of counts per AMF family ##
Glo_Fam_perc_reads_E1_8Sp  <- 
  ASV_table_Glo  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  group_by (Family)  %>% 
  summarize(Fam_sum=sum(counts)) %>% 
  add_column (Glo_total_counts_E1_8Sp)  %>% 
  mutate (Glo_Fam_perc = 100  * Fam_sum/total)  %>% 
  arrange(desc(Fam_sum))  %>% 
  as.data.frame()

# Number of different ASVs per AMF family ###
ASV_table_Glo  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  group_by (Family) %>%  select (-sampleID) %>% filter (counts!= 0) %>% 
  select (-counts) %>%  unique () %>%
  group_by(Family ) %>%  tally ()

# Mean number of AMF ASVs over all samples ##
ASV_table_Glo  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts") %>% 
  select (-Family) %>% filter (counts!=0) %>% 
  group_by (sampleID)  %>% 
  tally () %>% 
  summarize (meanASV = mean (n))


## Experiment 2 ##

ps_M1_allASVs_rarefied <- readRDS("data/ps_M1_allASVs_rarefied.rds") 
ps_M1_roots_allASVs_rarefied <- ps_M1_allASVs_rarefied %>%  
  subset_samples(roots_soil =="roots") %>% 
  prune_taxa( taxa_sums (.)>0, .)

ps_M1_roots_aftreDADA <-
ps_M1_allASVs %>%  
  subset_samples(roots_soil =="roots") %>% 
  prune_taxa( taxa_sums (.)>0, .)

ave_allASV <- 
  sample_sums (ps_M1_roots_aftreDADA) %>%  
  as_tibble (rownames = "sampleID") %>% 
  summarize (mean= mean (value))
sum_allASV <- 
  sample_sums (ps_M1_roots_aftreDADA) %>%  
  as_tibble (rownames = "sampleID") %>%  
  summarize (sum= sum (value))


ave_allASV_rar <- 
  sample_sums (ps_M1_roots_allASVs_rarefied) %>%  
  as_tibble (rownames = "sampleID") %>% 
  summarize (mean= mean (value))
sum_allASV_rar <- 
  sample_sums (ps_M1_roots_allASVs_rarefied) %>%  
  as_tibble (rownames = "sampleID") %>%  
  summarize (sum= sum (value))

# E2 -  All metrics calculation ####
#### Get a ASV table (all AMF ASVs)##
ASV_table_Glo_E2  <- 
  otu_table (ps_E2) %>% 
  data.frame() %>%  
  as_tibble (rownames = "ASV_ID") %>%   
  left_join ((tax_table (ps_E2) %>% data.frame() %>%  as_tibble (rownames = "ASV_ID")), by = "ASV_ID")


### AMF sequence stats ####

# Total AMF sequence counts in E2##
Glo_total_counts_E2  <-
  ASV_table_Glo_E2  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  summarize(total=sum(across(counts)))


#Percentage of counts per AMF family ##
Glo_Fam_perc_reads_E2  <- 
  ASV_table_Glo_E2  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  group_by (Family)  %>% 
  summarize(Fam_sum=sum(counts)) %>% 
  add_column (Glo_total_counts_E2)  %>% 
  mutate (Glo_Fam_perc = 100  * Fam_sum/total)  %>% 
  arrange(desc(Fam_sum))  %>% 
  as.data.frame()

# Number of different ASVs per AMF family ###
ASV_table_Glo_E2  %>%  
  select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts")  %>%
  group_by (Family) %>%  select (-sampleID) %>% filter (counts!= 0) %>% 
  select (-counts) %>%  unique () %>%
  group_by(Family ) %>%  tally ()

# Mean number of AMF ASVs over all samples ##
ASV_table_Glo_E2 %>% 
   select (-c(Kingdom, Phylum, Class, Order, Genus, Species))  %>%  
  pivot_longer (!c(Family, ASV_ID), names_to = "sampleID", values_to="counts") %>% 
  select (-Family) %>% filter (counts!=0) %>% 
  group_by (sampleID)  %>% 
  tally () %>% 
  summarize (meanASV = mean (n))

#
## Figure S3: y-diversity (unique ASVs) ####
#Plot unique ASVs per plant##
# and disstribution of AMF genera per plant species 
p1  <-  
  ASV_table_Glo  %>%  
  select(-c(Kingdom, Phylum, Class, Order, Species))  %>% 
  gather (-c(ASV_ID, Genus, Family), key=sampleID, value = ASV_counts)  %>%
  left_join(metaM0, by="sampleID") %>% 
  filter (ASV_counts!=0)  %>% 
  group_by (ASV_ID, PlantSpeciesfull, Genus) %>% 
  tally () %>%
  group_by(Genus, PlantSpeciesfull )  %>% 
  tally () %>% 
  replace_na(list (Genus="unidentified"))   %>% 
  left_join (meannumberASVsper_species, by = "PlantSpeciesfull")  %>% 
  add_column (Exp = "E1") %>% 
  rbind(ASV_table_Glo_E2  %>%  
          select(-c(Kingdom, Phylum, Class, Order, Species))  %>% 
          gather (-c(ASV_ID, Genus, Family), key=sampleID, value = ASV_counts)  %>%
          left_join(meta_M1Wcontrol, by="sampleID") %>% 
          filter (ASV_counts!=0)  %>% 
          group_by (ASV_ID, PlantSpeciesfull, Genus) %>% 
          tally () %>%
          group_by(Genus, PlantSpeciesfull )  %>% 
          tally () %>% 
          replace_na(list (Genus="unidentified")) %>% 
          left_join (meannumberASVsper_species, by = "PlantSpeciesfull")  %>% 
          add_column (Exp = "E2")) %>% 
  ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species), y= n, fill=Genus))  +
  geom_bar(position="stack", stat ="identity") +
  facet_wrap(~Exp, scales = "free_x")


  p1+ theme_light()+ 
  coord_flip () +
  ylab("Unique AMF ASVs per plant species") +
  xlab ("Plant species") +
  theme_light () +
  theme(axis.text.y= element_text(family= "sans", face= "italic", size = 12)) +
  theme (axis.text.x=element_text(family= "sans", size = 12))  +
  theme (axis.title = element_text(family = "sans", size = 14 ),  
         legend.title=element_text(size=13), 
         legend.text=element_text(size=11, face ="italic"),
         legend.position = "bottom",
         strip.background = element_rect(fill = "white", colour = "grey"), 
         strip.text = element_text(colour = "black")) +
    scale_fill_manual(values = c( "unidentified" ="#C7EAE5", "Archaeospora"  = "#287D8EFF" , 
                                "Acaulospora" =  "#C1A363" , "Diversispora" =  "darksalmon",
                                "Funneliformis" = "#F6E8C3", "Cetraspora"  = "#DFC27D", 
                                "Claroideoglomus" =  "#20A386FF", "Rhizophagus" = "darkgray",
                                "Glomus" = "#80CDC1" , "Paraglomus" = "#35978F", 
                                "Scutellospora"= "#01665E" ), name = "AMF genus"    )  

  



#PermANOVA to test if soils have significantly different AMF communities #####
#
  
ps_SoilE1 <- subset_samples(ps_ALL, PlaSpe == "Soil")  ## Exp1
  
ps_SoilE2 <- subset_samples(ps_ALL_rar_E2, PlaSpe == "SoiCon")    
  
  
  
ps_Soil_combined <- merge_phyloseq(
    phyloseq::otu_table(ps_SoilE1), phyloseq::otu_table(ps_SoilE2),
    phyloseq::sample_data(ps_SoilE1), phyloseq::sample_data(ps_SoilE2),
    phyloseq::tax_table(ps_SoilE1), phyloseq::tax_table(ps_SoilE2)
  )
ps_Soil  <-  subset_taxa(ps_Soil_combined, Phylum == "Glomeromycota")
  
  
  
# Extract the distance matrix (e.g., Bray-Curtis)
dist_matrix <- distance(ps_Soil, method = "bray")
  
# Extract the sample data 
sample_data_df <- as(sample_data(ps_Soil), "data.frame")
  
# Run PERMANOVA (the variable PlaSpe contains the different infomration o the different soil origin)
permanova_result <-  adonis2(dist_matrix ~ PlaSpe, data = sample_data_df, permutations = 9999)
  
# View the PERMANOVA result
print(permanova_result)
  

  

# Tables for paper - supporting information #### 
# Appendix S2: Table S3  all values metrics per plant species  ###
All_metrics_E1 %>%  rbind(All_metrics_E2) %>% 
    mutate  (`Shannon's H'` = round (`Shannon's H'`,2 ), 
             `Richness S` = round ( `Richness S`,2), CU = round (`β(CU)`, 2), 
             MPD = round (MPD, 2),
             `Phyl. γ-diversity` = (round(`Phyl. γ-diversity`,2))) %>%  
    select (!PlantFamily) %>% write_excel_csv("results/All_metrics_E1_E2.csv")
  
  
  
 # Appendix S2: Table S4  alpha diversity per sample  #
adiv_richness %>%  
    left_join(metaM0 %>%  select (sampleID, PlantSpeciesfull)) %>%  
    mutate (`Shannon's H'` = round (Shannon, 3)) %>% 
    write_excel_csv("results/tableS6_alpha_div_E1.csv")
  
adiv_richness_M1 %>%  select (-PlaSpe, -PlantFamily) %>%  
    mutate (`Shannon's H'` = round (Shannon, 3)) %>% 
    write_excel_csv("results/tableS6_alpha_div_E2.csv")
  
  
# Table S6: AMF per root and soil sample, total bacterial biomass (PLFA) , root and shoot biomass ############
  
dataNLrootstosoil %>% # amf soil is in µmol here
    mutate (AMF = AMF * 1000) %>%  # AMF is now in nmol per g now
    left_join(PLFA_Soil %>%  
                ungroup () %>%  
                select (!c(AMF, sampleID)))  %>%  
    select (sampleID, PlantSpeciesfull, AMFroot, AMF, totalBact, DW_above, DW_roots) %>%  #PLFA_soil is in nmol per g soil
    write_csv("results/TableS6.csv")

                          