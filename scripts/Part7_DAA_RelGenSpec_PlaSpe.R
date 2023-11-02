# Part 7  Differential abundance analysis  DAA ####

# packages ####
library (tidyverse)
library (edgeR)
library (phyloseq)
library (ggtree)
library (phylosmith)
library (ggpubr)


##data 
ps_allASVs_E2 <- readRDS ("data/ps_allASVs_E2.rds")

# Subset the ps object for the soil samples (rhizosphere)
ps_edgeR_relGenSpec <- ps_allASVs_E2 %>% subset_samples(roots_soil =="soil") 
# prune taxa that are not present in those subsetted samples
ps_edgeR_relGenSpec<- prune_taxa(taxa_sums (ps_edgeR_relGenSpec)>1,ps_edgeR_relGenSpec)
# add taxon information to ps object
corrected_taxTable <- 
  read_csv ("results/testtableTaxonomy.csv") %>% 
  data.frame (row.names = "ASV_ID") %>% 
  as.matrix() %>%   
  tax_table ()
tax_table(ps_edgeR_relGenSpec) <- corrected_taxTable


# Prepare the sample data - information of ID needed instead of sampleID, 
# and the values of the relative interaction generalism need to be added
test2 <- 
  RelGen_E1_E2_sample %>%  
  filter (Exp =="E2") %>%  
  arrange (RelGenSpec) %>% 
  select (-PlantSpeciesfull, - PlaSpe, -Exp, -PlantFamily) %>%    # these are already present in table
  left_join (meta_M1Wcontrol %>%  
               select (sampleID,  ID, DW_roots, DW_above), by = "sampleID") %>% 
  select (-sampleID) %>% 
  right_join(meta_M1Wcontrol %>%  filter (roots_soil == "soil") %>% select (ID, sampleID)) %>% 
  mutate (DW_roots = as.numeric (DW_roots), DW_above = as.numeric (DW_above))


# Change the sample data of the ps object to sample data that includes the relative interaction generalism
new_sampledata_gen <- 
  sample_data (ps_edgeR_relGenSpec) %>% 
  data.frame () %>% 
  select (-DW_roots, -DW_above) %>% 
  as_tibble (rownames= "sampleID") %>% 
  left_join (test2) %>% 
  mutate (RelGenSpec = replace_na ( RelGenSpec,0)) %>% 
data.frame (row.names = "sampleID") %>% 
  sample_data()
sample_data(ps_edgeR_relGenSpec)  <- new_sampledata_gen

# Prepare a ps object for building a ree for a figure
# for that, the taxa are agglomerated at the genus level.
ps_edgeR_genus <- tax_glom(ps_edgeR_relGenSpec, taxrank = "Genus", NArm = F)

# Prepare a dataframe from the ASV table for the use in package edgeR ##
# rows are OTUs, columns = samples
df <- otu_table(ps_edgeR_relGenSpec) %>%  
  data.frame ()

# Get a taxon table for the ASVs - needed as meta table for the DAA results
# NAs are changed to "unidentified" etc
taxa_edgeR_gen <- 
  tax_table (ps_edgeR_relGenSpec) %>% 
  data.frame () %>%  
  mutate(GenusLabel = ifelse(!is.na(Genus), paste(Genus), 
                             ifelse(!is.na(Family), paste('unid. ', Family,  sep = ""), 
                                    ifelse(!is.na(Order), paste('unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('unid. ', Class, sep = ""), paste("unid. ", Phylum, sep = "")))))) %>% 
  as_tibble (rownames = "ASV_ID" ) %>%  
  select (ASV_ID, GenusLabel,Phylum, Class, Family, Genus)

# Differential abundance analysis -  DAA#####
## Grouping : 
group = factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe) ## The plant species identity is used to group the samples.

## prepare edgeR object y
y <- DGEList(counts = df, group = group)

# Add sample information to the edgeR object
y$samples$PlaSpe <- factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe) # This is relevantt information 
y$samples$PlantSpeciesfull <- factor (sample_data (ps_edgeR_relGenSpec)$PlantSpeciesfull) # additional sample information
y$samples$PlantFamily <- factor (sample_data (ps_edgeR_relGenSpec)$PlantFamily)
y$samples$RelGenSpec <-  sample_data (ps_edgeR_relGenSpec)$RelGenSpec  ## this is relevant information

# Calculate the number of ASVs before filtering 
DAA_dim_before_filter <- dim(y)
# Filtering ####
# Remove ASVs of insufficient read counts from the dataset 
keep <- filterByExpr(y)
# Subset the dataset for the ASVs to keep
y <- y[keep, , keep.lib.sizes=FALSE]
DAA_dim_after_filter <-  dim (y)

#Normalisation####
#The calcNormFactors function normalises the library sizes by finding a set of scaling factors
#for the library sizes that minimizes the log-fold changes between the samples for most genes.
y<- calcNormFactors(y)

#A normalization factor below one indicates that a small number of high count genes
#are monopolizing the sequencing, causing the counts for other genes to be lower than would
#be usual given the library size.

# The data is releveled so that the soil control samples are the level of reference. 
# I.e. The comparison of all other levels  obtained by the grouping of the samples happens against those reference samples.
y$samples$group  <- relevel (y$samples$group, ref = "SoiCon")

# MODEL ## 
# The model design includes the grouping (plant species identity) and the relative interaction generalism as variables.
design <- model.matrix (~  group + RelGenSpec  , data = y$samples)  
design <- model.matrix (~ group, data = y$samples)

# see Law 2020, A guide to creating design matrices ... for how to design models and more for DAA. ##

#check the names of the columns with 
 colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)

# Testing for differentially abundant (DA) ASVs ######
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().

#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design)

#GLM ####
fit <- glmQLFit(y, design)

#GLM F-test Achillea millefolium - fit test 1 #####
#compare the coefficients of the glm 
# i.e. compare treatments
 # coefficent 2 represents the differnces in ASV abundances between the ASVs of the rhizosphere from A. millefolium 
 # and the control soil (plant-free) etc.
 
DAA_fit_test_AchMil <- glmQLFTest(fit, coef = 2) 

DAA_fit_test1_table <-
  DAA_fit_test_AchMil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AchMil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test  Relative interaction generalism - fit_test 2#####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fit_test2_vsC <- glmQLFTest(fit, coef = 10) # The coefficient 10 is the coefficient 
#representing the relative interaction generalism
# The resulting values are based on a comparison with the reference level, which we defined using the soil control samples.

#DA ASVs due the change in the relative interaction generalism
DAA_fit_test2_table <-
  DAA_fit_test2_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "RelGen") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test Agrostis capillaris - fit_test 3 #####
DAA_fit_test3_vsC <- glmQLFTest(fit, coef = 3) 

DAA_fit_test3_table <-
  DAA_fit_test3_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "AgrCap") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Bromus willdenowii - fit test 4 #####
DAA_fit_test4_vsC <- glmQLFTest(fit, coef = 4) 

DAA_fit_test4_table <-
  DAA_fit_test4_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "BroWil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Cichorium intybus - fit test 5#####
DAA_fit_test5_vsC <- glmQLFTest(fit, coef = 5) 

DAA_fit_test5_table <-
  DAA_fit_test5_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "CicInt") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test Holcus lanatus - fit test 6#####
DAA_fit_test6_vsC <- glmQLFTest(fit, coef = 6) 

DAA_fit_test6_table <-
  DAA_fit_test6_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "HolLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Platago lanceolata - fit test 7 #####
DAA_fit_test7_vsC <- glmQLFTest(fit, coef = 7) 

DAA_fit_test7_table <-
  DAA_fit_test7_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PlaLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Poa cita - fit test 8 #####

DAA_fit_test8_vsC <- glmQLFTest(fit, coef = 8) 

DAA_fit_test8_table <-
  DAA_fit_test8_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (PlaSpe = "PoaCit") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test Schedonorus arundinaceum - fit test 9 #####
DAA_fit_test9_vsC <- glmQLFTest(fit, coef = 9) 

DAA_fit_test9_table <-
  DAA_fit_test9_vsC$table %>%  
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
                             ifelse(!is.na(Family), paste('unid. ', Family, sep = ""), 
                                    ifelse(!is.na(Order), paste('unid. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('unid. ', Class, sep = ""), paste("unid. ", Phylum, sep = "")))))) 

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
  add_column (order = c(1:17))

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
  add_column ("PlaSpe" = "RelGen")


rel_ASV_abundance_soil <- 
  rel_ASV_abundance_soil_PlaSpe %>% rbind (rel_ASV_abundance_soil_all) 

#DAA ALL ####
#get data  for all DAA results into one tibble
# add information on relative abundance of the ASVs per AMF community!
# add genus labels
DAA <- 
  bind_rows (DAA_fit_test1_table, DAA_fit_test2_table, DAA_fit_test3_table, DAA_fit_test4_table,
                 DAA_fit_test5_table, DAA_fit_test6_table, DAA_fit_test7_table, DAA_fit_test8_table, DAA_fit_test9_table)  %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil) 


DAA_PlaSpe2 <- 
  bind_rows (DAA_fit_test1_table, DAA_fit_test2_table, DAA_fit_test3_table, DAA_fit_test4_table,
             DAA_fit_test5_table, DAA_fit_test6_table, DAA_fit_test7_table, DAA_fit_test8_table, DAA_fit_test9_table)  %>%  
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil) 

## make a dummy data frame to add a blank geom to the DAA plot (without it, some of the 
# AMF genera would be removed from each DAA plot
# and it would not be possible to align them with the phylogenetic tree)
DAA_empty_fields <- order_taxa %>%  add_column(logFC = 0, Sign = "ns", mean_rel_ASV = 0) 

# plot only the panel showing the relative interaction generalism ###
plotDAA_gen <-
  DAA %>% filter (PlaSpe == "RelGen")  %>% 
  ggplot (aes (x= logFC,  y= reorder (GenusLabel, -order), color =Order, alpha = Sign, size = mean_rel_ASV)) + 
  geom_blank (data =DAA_empty_fields, mapping = aes (x= logFC, y = reorder (GenusLabel, -order) , size = mean_rel_ASV)) +  ##this adds blank data - to include All AMF that are in the tree!!
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
  DAA %>% 
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
 xlim(-20,24) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
        title = "RelGenSpec + PlaSpe") 
   


# Arrange the tree and DAA plots into one figure #####
ggarrange (p_tree, plotDAA_gen, nrow = 1, widths = c(0.6,1)) ## depending on the plots, this might need to be adjusted
# manual adjustments of the figures can also be done in inkscape if needed.


