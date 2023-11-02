# DAA boxplot comprison dimenion 3 vs plant species  #### Figure ###

library (tidyverse)

library (edgeR)
library (phyloseq)
library (ggtree)
library (phylosmith)
library (ggpubr)


##data 
ps_M1_allASVs <- readRDS ("data/ps_M1_allASVs.rds")


test2 <- 
  PCA_metric_data_sample %>%  
  left_join(meta_M1Wcontrol %>% select (sampleID, ID, pot)) %>% 
  select (-PlantSpeciesfull, -PlantFamily) %>%  
  # left_join (meta_M1Wcontrol %>%  select (sampleID,  ID, DW_roots, DW_above), by = "sampleID") %>% 
  select (-sampleID)


### Add a new meta column to the data
# combination of PlaSpe and root/soil, e.g. AchMil-soil, AchMil-roots
# for use as grouping factor in edgeR!
ps_edgeR_relGenSpec <- ps_M1_allASVs %>% subset_samples(roots_soil =="soil") #%>%  subset_samples (AC.BC !="control")

ps_edgeR_relGenSpec<- prune_taxa(taxa_sums (ps_edgeR_relGenSpec)>1,ps_edgeR_relGenSpec)


new_sampledata_gen <- 
  sample_data (ps_edgeR_relGenSpec) %>% 
  data.frame () %>% 
  select (-DW_roots, -DW_above) %>% 
  as_tibble (rownames= "sampleID") %>% 
  left_join (test2, by = "ID") %>% 
  mutate (Dim.1 = replace_na ( Dim.1,0), Dim.2 = replace_na ( Dim.2,0), 
          Dim.1 = replace_na ( Dim.1,0), Dim.4 = replace_na ( Dim.4,0)) %>% 
  data.frame (row.names = "sampleID") %>% 
  sample_data()
sample_data(ps_edgeR_relGenSpec)  <- new_sampledata_gen

## different version of medelling - Dim.1 as variable 

Gen_table  <-  (new_sampledata_gen$Dim.1)

# at genus level - needed for Tree building and taxon names
ps_edgeR_genus <- tax_glom(ps_edgeR_relGenSpec, taxrank = "Genus", NArm = F)

#Object for edgeR ##
# rows OTUs, columns = samples
# get the count table as  df
df <- otu_table(ps_edgeR_relGenSpec) %>%  
  data.frame ()

#Get taxon table for the ASVs - needed as meta table for the DAA results
taxa_edgeR_gen <- 
  tax_table (ps_edgeR_relGenSpec) %>% 
  data.frame () %>%  
  mutate(GenusLabel = ifelse(!is.na(Genus), paste(Genus), 
                             ifelse(!is.na(Family), paste('Unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('Unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('Unid. ', Class, sep = ""), paste("Unid. ", Phylum, sep = "")))))) %>% 
  as_tibble (rownames = "ASV_ID" ) %>%  
  select (ASV_ID, GenusLabel,Phylum, Class, Family, Genus)



## Grouping : PlaSpe in roots and PlaSpe in soil = 16 treatments
group = factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe) ## The plant species identity is used to group the samples.

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
#add sample information
y$samples$PlaSpe <- factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe)
y$samples$PlantSpeciesfull <- factor (sample_data (ps_edgeR_relGenSpec)$PlantSpeciesfull)
y$samples$PlantFamily <- factor (sample_data (ps_edgeR_relGenSpec)$PlantFamily)
y$samples$Dim.1 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.1)
y$samples$Dim.1 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.1)
y$samples$Dim.2 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.2)

# relevel reference level
y$samples$group  <- relevel (y$samples$group, ref = "SoiCon")

# DESIGN: Decide on grouping ####
# see Law 2020, A guide to creating design matrices ...
## Model design 
design <- model.matrix (~  as.numeric (Dim.1) + group  , data = y$samples)


#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design)






colnames(design)
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
DAA_dim1_PS_dim1 <- glmQLFTest(fit, coef = 2) 

summary(decideTests(DAA_dim1_PS_dim1))



#Soil PlaSpe results####
# dim1 data soil #####
DAA_dim1_PS_dim1_table <-
  DAA_dim1_PS_dim1$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "Dim.1") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

## boxplot - changes in AMF taxa foe change of dim1 = 1
DAA_dim1_PS_dim1_table %>% ggplot (aes (x = Family , y = logFC)) + geom_violin() + geom_point()



#GLM F-test Achillea millefolium - fit test  #####
#compare the coefficients of the glm 
# i.e. compare treatments
# coefficent 2 represents the differnces in ASV abundances between the ASVs of the rhizosphere from A. millefolium 
# and the control soil (plant-free) etc.

DAA_dim1_PS_AchMil <- glmQLFTest(fit, coef = 3) 

DAA_dim1_PS1_table <-
  DAA_dim1_PS_AchMil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AchMil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Agrostis capillaris - fit 3 #####
DAA_dim1_PS_AgrCap <- glmQLFTest(fit, coef = 4) 

DAA_dim1_PS3_table <-
  DAA_dim1_PS_AgrCap$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AgrCap") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Bromus willdenowii - fit test 4 #####
DAA_dim1_PS4_vsC <- glmQLFTest(fit, coef = 5) 

DAA_dim1_PS4_table <-
  DAA_dim1_PS4_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "BroWil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Cichorium intybus - fit test 5#####
DAA_dim1_PS5_vsC <- glmQLFTest(fit, coef = 6) 

DAA_dim1_PS5_table <-
  DAA_dim1_PS5_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "CicInt") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test Holcus lanatus - fit test 6#####
DAA_dim1_PS6_vsC <- glmQLFTest(fit, coef = 7) 

DAA_dim1_PS6_table <-
  DAA_dim1_PS6_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "HolLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Platago lanceolata - fit test 7 #####
DAA_dim1_PS7_vsC <- glmQLFTest(fit, coef = 8) 

DAA_dim1_PS7_table <-
  DAA_dim1_PS7_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PlaLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Poa cita - fit test 8 #####

DAA_dim1_PS8_vsC <- glmQLFTest(fit, coef = 9) 

DAA_dim1_PS8_table <-
  DAA_dim1_PS8_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PoaCit") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Schedonorus arundinaceum - fit test 9 #####
DAA_dim1_PS9_vsC <- glmQLFTest(fit, coef = 10) 

DAA_dim1_PS9_table <-
  DAA_dim1_PS9_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "SchAru") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


# Prepare the data for tree building ##
## Subset the genus level data  for Glo
ps_edgeR_glo_genus <- subset_taxa(ps_edgeR_genus, Phylum == "Glomeromycota")
ps_edgeR_glo_genus <- prune_taxa (taxa_sums(ps_edgeR_glo_genus)>=1, 
                                  ps_edgeR_glo_genus) %>%  
  taxa_prune("ASV_3007") # this taxon is an unidentified Glomeromycota, therefore non- informative 

# get the tree data ####
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
                             ifelse(!is.na(Family), paste('Unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('Unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('Unid. ', Class, sep = ""), paste("Unid. ", Phylum, sep = "")))))) 

# get a tibble of the whole taxa table incl new GenusLabel
taxa_names <- df.tax %>% as_tibble ()

library (ape)

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
  theme_tree2() +
  theme (plot.margin = unit (c(8,5,6.5,5), "mm")) + 
  xlim(NA, 0.7) +
  theme (legend.position = "none")

#  get the order of the tree tips for use in the DAA plot
order_taxa <- get_taxa_name(p) %>%  
  as_tibble () %>% 
  dplyr::rename ("ASV_ID" = "value") %>%  
  left_join (taxa_names) %>% 
  add_column (order = c(1:18))

# get relative abundances for the ASVs for later use in the DAA plot
ps_edgeR_glo_ASV <- subset_taxa(ps_edgeR_relGenSpec, Phylum == "Glomeromycota")
ps_edgeR_glo_ASV_rel <- relative_abundance(ps_edgeR_glo_ASV)

rel_ASV_abundance_soil_PlaSpe <- 
  data.frame (otu_table(ps_edgeR_glo_ASV_rel))  %>%  
  as_tibble(rownames ="ASV_ID") %>%
  pivot_longer(!ASV_ID, names_to = "sampleID", values_to = "rel_ASV_per_sample") %>% 
  left_join((tax_table(ps_edgeR_glo_ASV_rel) %>% data.frame () %>% as_tibble (rownames = "ASV_ID")), by= "ASV_ID")  %>% 
  #filter (rel_ASV_per_sample!=0)  %>% 
  left_join (meta_M1Wcontrol) %>% ## these are the relative abundances per sample
  group_by(roots_soil, PlaSpe,  ASV_ID) %>%  # group to seperate soil from roots for each  PlaSpe , to get rel abu of ASVs per PlaSpe
  summarise (mean_rel_ASV =mean(rel_ASV_per_sample)) %>% 
  filter (roots_soil =="soil") # for Soil data only

rel_ASV_abundance_soil_all  <- 
  data.frame (otu_table(ps_edgeR_glo_ASV_rel))  %>%  
  as_tibble(rownames ="ASV_ID") %>%
  pivot_longer(!ASV_ID, names_to = "sampleID", values_to = "rel_ASV_per_sample") %>% 
  left_join((tax_table(ps_edgeR_glo_ASV_rel) %>% data.frame () %>% as_tibble (rownames = "ASV_ID")), by= "ASV_ID")  %>% 
  #filter (rel_ASV_per_sample!=0)  %>% 
  left_join (meta_M1Wcontrol) %>% ## these are the relative abundances per sample
  group_by(roots_soil,  ASV_ID) %>%  # group to seperate soil from roots for each  PlaSpe , to get rel abu of ASVs per PlaSpe
  summarise (mean_rel_ASV =mean(rel_ASV_per_sample)) %>% 
  filter (roots_soil =="soil") %>% 
  add_column ("PlaSpe" = "Dim.1")


rel_ASV_abundance_soil <- 
  rel_ASV_abundance_soil_PlaSpe %>% rbind (rel_ASV_abundance_soil_all) 

#DAA ALL ####
#get data  for all DAA results into one tibble
# add information on relative abundance of the ASVs per AMF community!
# add genus labels
DAA_dim1_PS_all <- 
  bind_rows ( DAA_dim1_PS1_table, DAA_dim1_PS_dim1_table, DAA_dim1_PS3_table, DAA_dim1_PS4_table,
             DAA_dim1_PS5_table, DAA_dim1_PS6_table, DAA_dim1_PS7_table, DAA_dim1_PS8_table, DAA_dim1_PS9_table)  %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil) 


## make a dummy data frame to add a blank geom to the DAA plot (without it, some of the 
# AMF genera would be removed from each DAA plot
# and it would not be possible to align them with the phylogenetic tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 

# plot only the panel showing the relative interaction generalism ###
plotDAA_gen <-
  DAA_dim1_PS_all %>% filter (PlaSpe == "Dim.1")  %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV)) + 
  # geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order) , size = mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
  geom_jitter(width =0.4) +  
  theme_classic() + 
  theme (axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
  scale_alpha_discrete(range = c(0.3, 1)) +
  theme (legend.position = "right" ) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
        title = "Relative interaction generalism") 


# Plot DAA for all plant species####
# This plot shows the logFC of the rhizosphere soils from each plant species in comparison to the control soil #
plotDAA <-
  DAA_dim1_PS_all %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV)) + 
  geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order) , size = mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
  geom_jitter(width =0.4) +  # geom_jitter
  facet_wrap(~PlaSpe, nrow = 1) +#, scales = "free_x") +
  theme_classic() + 
  theme (axis.title.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.line.y = element_blank(),
         axis.text.y = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
  scale_alpha_discrete(range = c(0.3, 1)) +
  theme (legend.position = "right" ) +
  #xlim(-20,24) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
        title = "Dim1 + PlaSpe") 



# Arrange the tree and DAA plots into one figure #####
ggarrange (p_tree, plotDAA_gen, nrow = 1, widths = c(0.6,1)) ## depending on the plots, this might need to be adjusted
# manual adjustments of the figures can also be done in inkscape if needed.



# COMPARISON DAA Plant species with / without dim 3 ######

# DESIGN####
# see Law 2020, A guide to creating design matrices ...
## Model design 
design <- model.matrix (~   group  , data = y$samples)


#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design)






colnames(design)
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




#GLM F-test Achillea millefolium - fit test  #####
#compare the coefficients of the glm 
# i.e. compare treatments
# coefficent 2 represents the differnces in ASV abundances between the ASVs of the rhizosphere from A. millefolium 
# and the control soil (plant-free) etc.

DAA_PS_AchMil <- glmQLFTest(fit, coef = 2) 

DAA_PS1_table <-
  DAA_PS_AchMil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AchMil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Agrostis capillaris - fit 3 #####
DAA_PS_AgrCap <- glmQLFTest(fit, coef =3) 

DAA_PS3_table <-
  DAA_PS_AgrCap$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AgrCap") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Bromus willdenowii - fit test 4 #####
DAA_PS4_vsC <- glmQLFTest(fit, coef = 4) 

DAA_PS4_table <-
  DAA_PS4_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "BroWil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Cichorium intybus - fit test 5#####
DAA_PS5_vsC <- glmQLFTest(fit, coef = 5) 

DAA_PS5_table <-
  DAA_PS5_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "CicInt") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test Holcus lanatus - fit test 6#####
DAA_PS6_vsC <- glmQLFTest(fit, coef = 6) 

DAA_PS6_table <-
  DAA_PS6_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "HolLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Platago lanceolata - fit test 7 #####
DAA_PS7_vsC <- glmQLFTest(fit, coef = 7) 

DAA_PS7_table <-
  DAA_PS7_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PlaLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Poa cita - fit test 8 #####

DAA_PS8_vsC <- glmQLFTest(fit, coef = 8) 

DAA_PS8_table <-
  DAA_PS8_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PoaCit") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Schedonorus arundinaceum - fit test 9 #####
DAA_PS9_vsC <- glmQLFTest(fit, coef = 9) 

DAA_PS9_table <-
  DAA_PS9_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "SchAru") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#DAA ALL plant species ####
#get data  for all DAA results into one tibble
# add information on relative abundance of the ASVs per AMF community!
# add genus labels
DAA_PS <- 
  bind_rows ( DAA_PS1_table, DAA_PS3_table, DAA_PS4_table,
              DAA_PS5_table, DAA_PS6_table, DAA_PS7_table, DAA_PS8_table, DAA_PS9_table)  %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil) 


## make a dummy data frame to add a blank geom to the DAA plot (without it, some of the 
# AMF genera would be removed from each DAA plot
# and it would not be possible to align them with the phylogenetic tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 

### compare DAA generlaism and plaSpe #### 
#Barplot Family * Plant species ####
DAA_dim1_PS_all %>%  
  add_column(model = "PS_removed_Dim.1") %>%
  bind_rows(DAA_PS %>%  add_column (model ="Plant")) %>% 
  filter (!is.na (Family)) %>% 
  filter (PlaSpe != "Dim.1") %>% 
  left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (Family, Phylum, GenusLabel, Order, order)) %>% 
  #filter (Sign == "sign.") %>% 
  group_by (Family, model, PlantSpeciesfull) %>% 
  summarize (meanLogFC = mean (logFC)) %>%  #### how to take a mean of log fold values?? That is actually okay to show trends like that
  ggplot (aes (x = PlantSpeciesfull, y = meanLogFC, color = model)) + 
  geom_col(aes (fill = model),position = position_dodge ()) +
  #geom_point () +
  facet_wrap(~  (Family), nrow = 3 )  +
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 



#Barplot AMF family####
DAA_dim1_PS_all %>%  
  add_column(model = "PS_removed_Dim.1") %>%
  bind_rows(DAA_PS %>%  add_column (model ="Plant")) %>% 
  filter (!is.na (Family)) %>% 
  filter (PlaSpe != "Dim.1") %>% 
  left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (Family, Phylum, GenusLabel, Order, order)) %>% 
  #filter (Sign == "sign.") %>% 
  group_by (Family, model) %>% 
  summarize (meanLogFC = mean (logFC)) %>%  #### how to take a mean of log fold values?? That is actually okay to show trends like that
  ggplot (aes (x = model, y = meanLogFC, color = model)) + 
  geom_col(aes (fill = model),position = position_dodge ()) +
  #geom_point () +
  #facet_wrap(~  (Family), nrow = 3 )  +
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 


#Boxplot AMF family####
DAA_dim1_PS_all %>%  
  add_column(model = "PS_removed_Dim.1") %>%
  bind_rows(DAA_PS %>%  add_column (model ="Plant")) %>% 
  filter (!is.na (Family)) %>% 
  filter (PlaSpe != "Dim.1") %>% 
  left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (Family, Phylum, GenusLabel, Order, order)) %>% 
  #filter (Sign == "sign.") %>% 
  group_by (Family, model) %>% 
  ggplot (aes (x = PlantSpeciesfull, y = logFC, color = model)) + 
  #geom_col(aes (fill = model),position = position_dodge ()) +
  geom_boxplot () +
  facet_wrap(~  (Family), nrow = 3 )  +
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 
#eom_point(position =  "dodge2" )


#Boxplot compare effect dim 1 and dim 3 family####
DAA_dim1_PS_all %>%  
  add_column(model = "PS_removed_Dim.1") %>%
  bind_rows(DAA_Dim3_PS_all %>%  add_column (model ="PS_removed_dim3")) %>% 
  filter (!is.na (Family)) %>% 
  filter (PlaSpe != "Dim.1", PlaSpe != "Dim.3") %>% 
  left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (Family, Phylum, GenusLabel, Order, order)) %>% 
  filter (Sign == "sign.") %>% 
  group_by (Family, model) %>% 
  ggplot (aes (x = model, y = logFC, color = model)) + 
  #geom_col(aes (fill = model),position = position_dodge ()) +
  geom_boxplot () +
  #facet_wrap(~  (Family), nrow = 3 )  +
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



