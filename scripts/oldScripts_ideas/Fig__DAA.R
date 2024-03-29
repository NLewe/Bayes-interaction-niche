# packages s####
library (edgeR)
library (phyloseq)
library (tidyverse)
library(RVAideMemoire)
library (ggpubr)

##data 
ps_M1_allASVs <- readRDS ("data/ps_M1_allASVs.rds")

### Add a new meta column to the data
# combination of PlaSpe and root/soil, e.g. AchMil-soil, AchMil-roots
# for use as grouping factor in edgeR!
ps_edgeR <- ps_M1_allASVs %>% subset_samples(roots_soil =="soil") 



# include better identification of taxa##
# check ASVs not identified at family level on UNITE database#

# use tsv  file, copy sequence, blast on UNITE 

# tax_table (ps_edgeR) %>%  
#   data.frame () %>%  
#   as_tibble (rownames = "ASV_ID") %>% 
#   write_csv ("results/testtableTaxonomy.csv")


# used blast to identify the following  AVs
#386, 1320, 1484, 1541, 1106 , added the result to csv table


corrected_taxTable <- 
  read_csv ("results/testtableTaxonomy.csv") %>% 
  data.frame (row.names = "ASV_ID") %>%  as.matrix() %>%   tax_table ()


# add corrected taxa data tp ps object ####

tax_table(ps_edgeR) <- corrected_taxTable


## add sample data
new_sampledata <- 
  sample_data (ps_edgeR) %>% 
  data.frame () %>% 
  unite (col =PlaSpe_RS, c(PlaSpe, roots_soil), sep= "_" ) %>% 
  sample_data()
sample_data(ps_edgeR)  <- new_sampledata

## different version of medelling - generalism as variable 



# at genus level - needed for Tree building and taxon names
ps_edgeR_genus <- tax_glom(ps_edgeR, taxrank = "Genus", NArm = F)

#Object for edgeR ##
# rows OTUs, columns = samples
# get the count table as  df
df <- otu_table(ps_edgeR) %>%  
  data.frame ()

#Get taxon table for the ASVs - needed as meta table for the DAA results
taxa_edgeR <- 
  tax_table (ps_edgeR) %>% 
  data.frame () %>%  
  mutate(GenusLabel = ifelse(!is.na(Genus), paste(Genus), 
                             ifelse(!is.na(Family), paste('unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('unid. ', Class, sep = ""), paste("Unid. ", Phylum, sep = "")))))) %>% 
  as_tibble (rownames = "ASV_ID" ) %>%  
  select (ASV_ID, GenusLabel,Phylum, Class, Family, Genus)


## Grouping : PlaSpe in roots and PlaSpe in soil = 16 treatments
group = factor (sample_data (ps_edgeR)$PlaSpe_RS)

## prepare edgeR object y
y <- DGEList(counts = df, group = group)
DAA_dim_before_filter <- dim(y)
# Filtering ####
#filter small OTUs- that removes about 2000 OTUs
keep <- filterByExpr(y)
# filter the dataset for the ASVs to keep
y <- y[keep, , keep.lib.sizes=FALSE]
DAA_dim_after_filter <-  dim (y)

#Normalisation####
#The calcNormFactors function normalises the library sizes by finding a set of scaling factors
#for the library sizes that minimizes the log-fold changes between the samples for most genes.
y<- calcNormFactors(y)
#y$samples
#A normalization factor below one indicates that a small number of high count genes
#are monopolizing the sequencing, causing the counts for other genes to be lower than would
#be usual given the library size.
y$samples$group  <- relevel (y$samples$group, ref = "SoiCon_soil")

design <- model.matrix (~ group , data = y$samples)

#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design)

#add sample information
y$samples$PlaSpe <- factor (sample_data (ps_edgeR)$PlaSpe)
y$samples$PlantSpeciesfull <- factor (sample_data (ps_edgeR)$PlantSpeciesfull)
y$samples$PlantFamily <- factor (sample_data (ps_edgeR)$PlantFamily)


# relevel reference level

# DESIGN: Decide on grouping ####
# see Law 2020, A guide to creating design matrices ...
## Model design 

colnames(design) <- levels(y$samples$group)
#check the names of the columns with 
# colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)

# testing for differentially abundant (DA) OTUs
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().


#GLM ####
fit <- glmQLFit(y, design)






#GLM F-test #####
#compare the coefficients of the glm 
# i.e. compare treatments
AchMil_soil <- glmQLFTest(fit, coef=2) 
#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_AchMil_soil <-
  AchMil_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AchMil", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

# Agr Cap soil ##
# i.e. compare treatments
AgrCap_soil <- glmQLFTest(fit, coef=3) 

# Ach Mil data soil
DAA_AgrCap_soil <-
  AgrCap_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AgrCap", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

## BroWIl soil
BroWil_soil <- glmQLFTest(fit, coef=4) 

DAA_BroWil_soil <-
  BroWil_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "BroWil", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns"))

# CicInt data soil


CicInt_soil <- glmQLFTest(fit, coef=5) 

DAA_CicInt_soil <-
  CicInt_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "CicInt", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

# HolLan data soil


HolLan_soil <- glmQLFTest(fit, coef=6) 

DAA_HolLan_soil <-
  HolLan_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "HolLan", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

# PlaLan data soil

PlaLan_soil <- glmQLFTest(fit, coef=7) 

DAA_PlaLan_soil <-
  PlaLan_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PlaLan", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

# PoaCit data soil


PoaCit_soil <- glmQLFTest(fit, coef=8) 

DAA_PoaCit_soil <-
  PoaCit_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PoaCit", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 


# SchAru data soil

SchAru_soil <- glmQLFTest(fit, coef=9) 

DAA_SchAru_soil <-
  SchAru_soil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "SchAru", roots_soil = "soil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) 

## Subset the genus level data  for Glo
ps_edgeR_glo_genus <- subset_taxa(ps_edgeR_genus, Phylum == "Glomeromycota")
ps_edgeR_glo_genus <- prune_taxa (taxa_sums(ps_edgeR_glo_genus)>=1, ps_edgeR_glo_genus)  %>%  
  taxa_prune("ASV_3007")

# get the tree data
MyTree <-    ps_edgeR_glo_genus %>% phy_tree

# Save names of taxa in tree
TreeTax <-  taxa_names(ps_edgeR_glo_genus)


# create new taxonomy labels
df.tax <-  ps_edgeR_glo_genus %>% tax_table %>% as.data.frame
df.tax$ASV_ID   <-  df.tax %>% row.names # add column with the ASV iD 

df.tax <-  df.tax  %>% mutate( TaxLabel = paste(Family, Genus, sep = "_")) %>%
  select(ASV_ID, TaxLabel, Phylum, Class, Order, Family, Genus)

# Change the NA in the taxon table to the nearest identified taxon
df.tax = df.tax %>%
  mutate(GenusLabel = ifelse(!is.na(Genus), paste(Genus), 
                             ifelse(!is.na(Family), paste( 'unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('unid. ', Class, sep = ""), paste("Unid. ", Phylum, sep = "")))))) 

# get a tibble of the whole taxa table incl new GenusLabel
taxa_names <- df.tax %>% as_tibble ()

library (ggtree)
library (ape)
## check tree ##

p  = ggtree(MyTree, ladderize = F) %<+% df.tax

p_tree <- 
  p +  geom_tiplab(aes(label = GenusLabel, color = Order), align = TRUE, size = 4) +
  # geom_tippoint(aes(color= Order), size=3, show.legend = F) +
  theme_tree2() +
  theme (plot.margin = unit (c(8,5,6.5,5), "mm")) + 
  xlim(NA, 0.7) +
  theme (legend.position = "none")

# Tree plot ####
test_tree <- rotateConstr(MyTree, constraint = c("ASV_3904" , "ASV_3003", "ASV_3844" ,
                                                 "ASV_14" , "ASV_713"  ,"ASV_357",
                                                 "ASV_988", "ASV_247","ASV_386", "ASV_1106" , 
                                                 "ASV_270", "ASV_550", "ASV_1334", 
                                                 "ASV_978",    "ASV_2852", 
                                                 "ASV_1320","ASV_1540"))

p  = ggtree(test_tree, ladderize = F) %<+% df.tax

p_tree <- 
  p +  geom_tiplab(aes(label = GenusLabel, color = Order), align = TRUE, size = 4) +
  # geom_tippoint(aes(color= Order), size=3, show.legend = F) +
  theme_tree2() +
  theme (plot.margin = unit (c(8,5,6.5,5), "mm")) + 
  xlim(NA, 0.7) +
  theme (legend.position = "none")

#  get the order of the tree tips for use in the DAA plot
order_taxa <- get_taxa_name(p) %>%  
  as_tibble () %>% 
  dplyr::rename ("ASV_ID" = "value") %>%  
  left_join (taxa_names) %>% 
  add_column (order = c(1:17))

# get relative abundances for the ASVs for later use in the DAA plot
ps_edgeR_glo_ASV <- subset_taxa(ps_edgeR, Phylum == "Glomeromycota")
ps_edgeR_glo_ASV_rel <- relative_abundance(ps_edgeR_glo_ASV)

rel_ASV_abundance_soil <- 
  data.frame (otu_table(ps_edgeR_glo_ASV_rel))  %>%  
  as_tibble(rownames ="ASV_ID") %>%
  pivot_longer(!ASV_ID, names_to = "sampleID", values_to = "rel_ASV_per_sample") %>% 
  left_join((tax_table(ps_edgeR_glo_ASV_rel) %>% data.frame () %>% as_tibble (rownames = "ASV_ID")), by= "ASV_ID")  %>% 
  #filter (rel_ASV_per_sample!=0)  %>% 
  left_join (meta_M1Wcontrol) %>% ## these are the relative abundances per sample
  group_by(roots_soil,PlaSpe, ASV_ID) %>%  # group to seperate soil from roots for each  PlaSpe , to get rel abu of ASVs per PlaSpe
  summarise (mean_rel_ASV =mean(rel_ASV_per_sample)) %>% 
  filter (roots_soil =="soil") # for Soil data only


#DAA ALL ####
#get data  for all DAA results into one tibble
# add information on relative abundance of the ASVs per Glo community!
# add genus label
DAA_PlaSpe <- bind_rows(DAA_AchMil_soil, DAA_AgrCap_soil, DAA_BroWil_soil, DAA_CicInt_soil, 
                 DAA_HolLan_soil, DAA_PlaLan_soil, DAA_PoaCit_soil,  DAA_SchAru_soil) %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil)  %>% 
  filter (!is.na (mean_rel_ASV))  # remove all rel abundances and logFC that have 0 abundance in the respective PlaSpe

## make a dummy df to add a blank geom to the DAA plot (without it, some of the 
# AMF genera are mssing in the DAA plot
# and it cannot be aligned with tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 




# Plot DAA ####
plotDAA <-
  DAA_PlaSpe %>%   #mutate (across (PlaSpe, factor, levels = c("AchMil", "CicInt", "PlaLan","PoaCit", "SchAru", "BroWil", "AgrCap", "HolLan" ))) %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV_per_PlaSpe)) + 
  geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order), size =mean_rel_ASV_per_PlaSpe)) +  ##this adds blank data - to include All AMF that are in the tree!!
  geom_jitter(width =0.4) +  # geom_jitter
  facet_wrap(~PlaSpe, nrow =1) +
  theme_classic() + 
  theme (axis.title.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.line.y = element_blank(),
         axis.text.y = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
  scale_alpha_discrete(range = c(0.3, 1)) +
  theme (legend.position = "right" ) +
  xlim(-9,13) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance")


# Plot tree and DAA together #####
ggarrange (p_tree, plotDAA, nrow = 1, widths = c(0.4,1))

DAA_Dim1 <- DAA
### compare DAA generlaism and plaSpe #### 
 #Boxplot ####
DAA_Dim1 %>%  
  add_column(model = "Dim1") %>%
  bind_rows (DAA_Dim3 %>%  add_column (model = "Dim3")) %>% 
  bind_rows (DAA_Dim2 %>%  add_column (model = "Dim2")) %>% 
  bind_rows(DAA_PlaSpe2 %>%  add_column (model ="Plant") %>% filter (PlaSpe != "RelGen")) %>% 
  filter (!is.na (Family)) %>% 
  filter (PlaSpe != "Dim1") %>% 
  left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (Family, Phylum, GenusLabel, Order, order)) %>% 
  #filter (Sign == "sign.") %>% 
  group_by (Family, model, PlantSpeciesfull) %>% 
  summarize (meanLogFC = mean (logFC)) %>%  #### how to tae a mean of log fold values??
  ggplot (aes (x = PlantSpeciesfull, y = meanLogFC, color = model)) + 
  geom_col(aes (fill = model),position = position_dodge ()) +
  #geom_point () +
  facet_wrap(~  (Family), nrow = 3 )  +
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 
  
#eom_point(position =  "dodge2" )


### Plot comparison log FC ####
DAA_dim123 %>% filter (!is.na (Family)) %>% 
  #filter (PlaSpe != "RelGen") %>% 
  #left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (GenusLabel, order)) %>% 
  #filter (Sign == "sign.") %>% 
  ggplot (aes (x = PlaSpe, y = logFC)) + 
  #geom_violin () +
  geom_boxplot() +
  geom_point () +
  #facet_wrap(~ reorder (Family, order), nrow = 3 )+
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 


# t tests for the boxplots ###  
# t-test on model 

# Perm_DAA_boxplot_tTest <- 
# DAA %>%  add_column(model = "RelGen") %>% 
#   bind_rows(DAA_PlaSpe %>%  add_column (model ="Plant")) %>% 
#   filter (!is.na (GenusLabel)) %>% 
#   filter (Family!= "Diversisporaceae") %>% 
#   filter (PlaSpe != "RelGen") %>% 
#   left_join (meta_plants) %>%  
#   group_split(PlantSpeciesfull, GenusLabel) %>% 
#   map (~unite (.,"PlaSpe_Genus", PlaSpe, GenusLabel, sep= "_", remove = F  )) %>% 
#   purrr::set_names(purrr::map_chr(., ~.x$PlaSpe_Genus[1])) %>% 
#   map (~perm.t.test (.$logFC ~ .$model, nperm= 9999, paired = T, .id = .$PlaSpe_Genus)) 

#%>% 
Perm_DAA_boxplot_tTest %>% map_dbl ("p.value")


# ## Permutational t test ~ PlaSpe, paired  ####
# #to compare if the log FC differ between the models (for each plant species) ##
# Perm_DAA_boxplot_tTest_plaspe <- 
#   DAA %>%  add_column(model = "RelGen") %>% 
#   bind_rows(DAA_PlaSpe %>%  add_column (model ="Plant")) %>% 
#   filter (!is.na (GenusLabel)) %>% 
#   #filter (Family!= "Diversisporaceae") %>% 
#   filter (PlaSpe != "RelGen") %>% 
#   left_join (meta_plants) %>%  
#   group_split(PlantSpeciesfull) %>% 
#   purrr::set_names(purrr::map_chr(., ~.x$PlaSpe[1])) %>% 
#   map (~perm.t.test (.$logFC ~ .$model, nperm= 999, paired = T, .id = .$PlaSpe)) 
# 
# # Permutational t-test ~ GenusLabel ####
# 
# Perm_DAA_boxplot_tTest_GenusLabel <- 
#   DAA %>%  add_column(model = "RelGen") %>% 
#   bind_rows(DAA_PlaSpe %>%  add_column (model ="Plant")) %>% 
#   filter (!is.na (GenusLabel)) %>% 
#   filter (Family!= "Diversisporaceae") %>% 
#   filter (PlaSpe != "RelGen") %>% 
#   left_join (meta_plants) %>%  
#   group_split(GenusLabel) %>% 
#   purrr::set_names(purrr::map_chr(., ~.x$GenusLabel[1])) %>% 
#   map (~perm.t.test (.$logFC ~ .$model, nperm= 9999, paired = T, .id = .$GenusLabel)) 

## add model name to DAA_empty so that plot identifies it ##

DAA_empty_fields <-DAA_empty_fields %>%  add_column(model = "Plant")

### both models in one DAA plot #####
plotDAA_bothM <- DAA %>%  add_column(model = "RelGen") %>%  bind_rows(DAA_PlaSpe %>%  add_column (model ="Plant")) %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color = Order , 
               shape =model, alpha = Sign, size = mean_rel_ASV)) + 
  geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order), color = Order , 
                                                     shape =model,
                                                   size =mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
  geom_jitter(width =0.4) +  # geom_jitter
  #geom_point()+
  facet_wrap(~PlaSpe , nrow =1) +
  theme_classic() + 
  theme (axis.title.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.line.y = element_blank(),
         axis.text.y = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
  scale_alpha_discrete(range = c(0.3, 1)) +
  theme (legend.position = "right" ) +
 # xlim(-9,13) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", shape = "Model")


# Plot tree and DAA together #####
ggarrange (p_tree, plotDAA_bothM, nrow = 1, widths = c(0.4,1))

