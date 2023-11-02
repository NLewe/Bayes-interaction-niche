



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
          Dim.3 = replace_na ( Dim.3,0), Dim.4 = replace_na ( Dim.4,0)) %>% 
  data.frame (row.names = "sampleID") %>% 
  sample_data()
sample_data(ps_edgeR_relGenSpec)  <- new_sampledata_gen

## different version of medelling - Dim.3 as variable 

#Gen_table  <-  (new_sampledata_gen$Dim.3)

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
y$samples$Dim.3 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.3)
y$samples$Dim.1 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.1)
y$samples$Dim.2 <- as.numeric (sample_data (ps_edgeR_relGenSpec)$Dim.2)

# relevel reference level
y$samples$group  <- relevel (y$samples$group, ref = "SoiCon")

# Design dimension 1 ####
design_Dim1 <- model.matrix (~ as.numeric (Dim.1)  + group  , data = y$samples)

# see Law 2020, A guide to creating design matrices ...
## Model design 

colnames(design_Dim1)
#check the names of the columns with 
# colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)

#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design_Dim1)

# testing for differentially abundant (DA) OTUs
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().


#GLM ####
fit <- glmQLFit(y, design_Dim1)

#GLM F-test #####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fit_Dim1 <- glmQLFTest(fit, coef = 2) 

summary(decideTests(DAA_fit_Dim1))

# Dim1 data table soil #####
DAA_fit_Dim1_table <-
  DAA_fit_Dim1$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "Dim1") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

## boxplot - changes in AMF taxa for change of Dim1 = 1
DAA_fit_Dim1_table %>% ggplot (aes (x = Family , y = logFC)) + geom_violin() + geom_point()



# Design dimension 2 ####
design_Dim2 <- model.matrix (~ as.numeric (Dim.2)  + group  , data = y$samples)
# DESIGN: Decide on grouping ###
# see Law 2020, A guide to creating design matrices ...
## Model design 

colnames(design_Dim2)
#check the names of the columns with 
# colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)


#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design_Dim2)



# testing for differentially abundant (DA) OTUs
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().


#GLM ####
fit <- glmQLFit(y, design_Dim2)

#GLM F-test #####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fit_Dim2 <- glmQLFTest(fit, coef = 2) 

summary(decideTests(DAA_fit_Dim2))

# Dim2 data table soil #####
DAA_fit_Dim2_table <-
  DAA_fit_Dim2$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "Dim2") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

## boxplot - changes in AMF taxa for change of Dim1 = 1
DAA_fit_Dim2_table %>% ggplot (aes (x = Family , y = logFC)) + geom_violin() + geom_point()


# Design dimension 3 ####
design_Dim3 <- model.matrix (~ as.numeric (Dim.3)  + group  , data = y$samples)
# DESIGN: Decide on grouping ###
# see Law 2020, A guide to creating design matrices ...
## Model design 

colnames(design_Dim3)
#check the names of the columns with 
# colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)


#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design_Dim3)



# testing for differentially abundant (DA) OTUs
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().


#GLM ####
fit <- glmQLFit(y, design_Dim3)

#GLM F-test #####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fit_Dim3 <- glmQLFTest(fit, coef = 2) 

summary(decideTests(DAA_fit_Dim3))

# Dim3 data table soil #####
DAA_fit_Dim3_table <-
  DAA_fit_Dim3$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "Dim3") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

## boxplot - changes in AMF taxa for change of Dim1 = 1
DAA_fit_Dim3_table %>% ggplot (aes (x = Family , y = logFC)) + geom_violin() + geom_point()


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
  filter (roots_soil =="soil") 
#%>% 
 # add_column ("PlaSpe" = "Dim1")


rel_ASV_abundance_soil <- 
  rel_ASV_abundance_soil_PlaSpe %>% rbind (rel_ASV_abundance_soil_all) 

#DAA ALL ####
#get data  for all DAA results into one tibble
# add information on relative abundance of the ASVs per AMF community!
# add genus labels
DAA_dim123 <- 
  bind_rows (DAA_fit_Dim1_table, DAA_fit_Dim2_table, DAA_fit_Dim3_table)  %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil_all) 


## make a dummy data frame to add a blank geom to the DAA plot (without it, some of the 
# AMF genera would be removed from each DAA plot
# and it would not be possible to align them with the phylogenetic tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 

# plot only the panel showing the log FC correlated to the three PCA dimensions  ###
plotDAA_dim123 <-
  DAA_dim123  %>% #filter (PlaSpe == "Dim3")  %>% 
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
  facet_wrap ( ~ PlaSpe) +
  theme (legend.position = "right" ) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
        title = "Relative interaction generalism") 


# Plot DAA for all plant species####
# This plot shows the logFC of the rhizosphere soils from each plant species in comparison to the control soil #
# plotDAA <-
#   DAA_Dim_all %>% 
#   ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV)) + 
#   geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order) , size = mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
#   geom_jitter(width =0.4) +  # geom_jitter
#   facet_wrap(~PlaSpe, nrow = 1) +#, scales = "free_x") +
#   theme_classic() + 
#   theme (axis.title.y = element_blank(),
#          axis.ticks.y = element_blank(), 
#          axis.line.y = element_blank(),
#          axis.text.y = element_blank()) +
#   geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
#   scale_alpha_discrete(range = c(0.3, 1)) +
#   theme (legend.position = "right" ) +
#   #xlim(-20,24) +
#   labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
#         title = "Dim1 + PlaSpe") 



### Plot comparison log FC ####
DAA_dim123 %>% filter (!is.na (Family)) %>% 
  #filter (PlaSpe != "RelGen") %>% 
  #left_join (meta_plants) %>% 
  select (!order) %>% 
  left_join (order_taxa %>%  select (GenusLabel, order)) %>% 
  #filter (Sign == "sign.") %>% 
  ggplot (aes (x = PlaSpe, y = logFC)) + 
  #geom_violin () +
  geom_boxplot() +   #### boxplot do not make sense with log fold change !!!1
  geom_point () +
  facet_wrap(~ reorder (Family, order), nrow = 3 )+
  theme_classic() +
  theme (axis.text.x = element_text (face ="italic", angle = 45, hjust =1 )) 




# Arrange the tree and DAA plots into one figure #####
ggarrange (p_tree, plotDAA_gen, nrow = 1, widths = c(0.6,1)) ## depending on the plots, this might need to be adjusted
# manual adjustments of the figures can also be done in inkscape if needed.





# 
# tab <- as.data.frame (topTags(DAA_fit_Dim3, n = 310))
# 
# tab2 <-  
#   tab %>% as_tibble (rownames = "ASV_ID") %>%  left_join (taxa_edgeR_gen) %>%  
#   filter (Phylum == "Glomeromycota") %>% data.frame (row.names = "ASV_ID")
# 
# 
# logCPM.obs <- cpm(y, log=TRUE, prior.count=fit$prior.count)
# 
# logCPM.fit <- cpm(DAA_fit_Dim3, log=TRUE)
# 
# 
# 
# 


# ### None of the ASVs show any relevant changes ##
# par(mfrow=c(5,5))
# 
# for(i in 1:26) {
#   ASV_ID <- row.names(tab2)[i]
#   logCPM.obs.i <- logCPM.obs[ASV_ID,]
#   logCPM.fit.i <- logCPM.fit[ASV_ID,]
#   plot(new_sampledata_gen$Dim.3, logCPM.obs.i, ylab="log-CPM", main=ASV_ID, pch=16)
#  lines(new_sampledata_gen$Dim.3, logCPM.fit.i, col="red", lwd=2)
#   }



# PlaLan data soil

PlaLan_soil <- glmQLFTest(fit, coef=13) 




## Subset the genus level data  for Glo
ps_edgeR_glo_genus <- subset_taxa(ps_edgeR_genus, Phylum == "Glomeromycota")
ps_edgeR_glo_genus <- prune_taxa (taxa_sums(ps_edgeR_glo_genus)>=1, ps_edgeR_glo_genus)

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
                             ifelse(!is.na(Family), paste('Unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('Unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('Unid. ', Class, sep = ""), paste("Unid. ", Phylum, sep = "")))))) 

# get a tibble of the whole taxa table incl new GenusLabel
taxa_names <- df.tax %>% as_tibble ()

library (ape)
# Tree plot ####
test_tree <- rotateConstr(MyTree, constraint = c("ASV_3904" ,"ASV_1540", "ASV_1320",   "ASV_1856", 
                                                 "ASV_2852","ASV_1705","ASV_1334", "ASV_550" ,
                                                 "ASV_270","ASV_386",  "ASV_662" , "ASV_988" , 
                                                 "ASV_247" , "ASV_713"  , "ASV_357", "ASV_253" ,  "ASV_14" ,   "ASV_3003"))

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
  add_column (order = c(1:18))

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
  group_by(roots_soil, ASV_ID) %>%  # group to seperate soil from roots for each  PlaSpe , to get rel abu of ASVs per PlaSpe
  summarise (mean_rel_ASV =mean(rel_ASV_per_sample)) %>% 
  filter (roots_soil =="soil") # for Soil data only


#DAA ALL ####
#get data  for all DAA results into one tibble
# # add information on relative abundance of the ASVs per Glo community!
# # add genus label
# DAA <- (DAA_fit066_table) %>%  
#   left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
#   left_join (rel_ASV_abundance_soil)  %>% 
#   filter (!is.na (mean_rel_ASV))  # remove all rel abundances and logFC that have 0 abundance in the respective PlaSpe

## make a dummy df to add a blank geom to the DAA plot (without it, some of the 
# AMF genera are mssing in the DAA plot
# and it cannot be aligned with tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 




# Plot DAA ####
plotDAA <-
  DAA %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV)) + 
  geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order) , size = mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
  geom_jitter(width =0.4) +  # geom_jitter
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
# ggarrange (p_tree, plotDAA, nrow = 1, widths = c(0.4,1))

