
# packages ####
library (tidyverse)
library (edgeR)
library (phyloseq)
library (ggtree)
library (phylosmith)
library (ggpubr)


##data 
ps_M1_allASVs <- readRDS ("data/ps_M1_allASVs.rds")

#repair meta
# repairSampledata <- 
#   sample_data (ps_M1_allASVs) %>%  
#   data.frame () %>%  
#   as_tibble (rownames ="sampleID") %>%  
#   select (!ID) %>% 
#   left_join (meta_M1Wcontrol %>%  select (sampleID, ID)) %>% 
#   data.frame (row.names = "sampleID") %>% 
#   sample_data()
# 
# sample_data (ps_M1_allASVs)  <- repairSampledata
# 
# saveRDS (ps_M1_allASVs, "data/ps_M1_allASVs.rds")

### Add a new meta column to the data
# combination of PlaSpe and root/soil, e.g. AchMil-soil, AchMil-roots
# for use as grouping factor in edgeR!
ps_edgeR_relGenSpec <- ps_M1_allASVs %>% subset_samples(roots_soil =="soil") #%>%  subset_samples (PlaSpe != "SoiCon")

ps_edgeR_relGenSpec<- prune_taxa(taxa_sums (ps_edgeR_relGenSpec)>1,ps_edgeR_relGenSpec)



# prep sample data - information of ID needed
test2 <- 
  RelGen_E1_E2_sample %>%  
  filter (Exp =="E2") %>%  
  arrange (RelGenSpec) %>% 
 # mutate (GenBins = paste0 ("Bin", rep(1:8, each = 5))) %>% 
  select (-PlantSpeciesfull, - PlaSpe, -Exp, -PlantFamily) %>%    # these are alread n other tbl
  left_join (meta_M1Wcontrol %>%  
               select (sampleID,  ID, DW_roots, DW_above), by = "sampleID") %>% 
  select (-sampleID) %>% 
  right_join(meta_M1Wcontrol %>%  filter (roots_soil == "soil") %>% select (ID, sampleID)) %>% 
  mutate (DW_roots = as.numeric (DW_roots), DW_above = as.numeric (DW_above))



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

## different version of medelling - generalism as variable 

#Gen_table  <-  factor(new_sampledata_gen$GenBins)


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

# DAA#####
## Grouping : 
group = factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe)

## prepare edgeR object y
y <- DGEList(counts = df, group = group)

#add sample information
y$samples$PlaSpe <- factor (sample_data (ps_edgeR_relGenSpec)$PlaSpe)
y$samples$PlantSpeciesfull <- factor (sample_data (ps_edgeR_relGenSpec)$PlantSpeciesfull)
y$samples$PlantFamily <- factor (sample_data (ps_edgeR_relGenSpec)$PlantFamily)
#y$samples$Gen <- factor (sample_data (ps_edgeR_relGenSpec)$Gen)
y$samples$RelGenSpec <-  sample_data (ps_edgeR_relGenSpec)$RelGenSpec


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
y$samples$group  <- relevel (y$samples$group, ref = "SoiCon")
design <- model.matrix (~  group + RelGenSpec  , data = y$samples)

#Dispersions - GLM####
#For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25].
y <- estimateDisp(y, design)




# relevel reference level

# DESIGN: Decide on grouping ####
# see Law 2020, A guide to creating design matrices ...
## Model design 

#colnames(design) <- levels(y$samples$group)
#check the names of the columns with 
 colnames(design)
# to get the needed comparisons (either by writing contrasts or by direct use of the coefficient's results)

# testing for differentially abundant (DA) OTUs
#While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small.
#The QL dispersion estimation and hypothesis testing can be done by using the functions
#glmQLFit() and glmQLFTest().

 # my_contrasts <- makeContrasts (
 #  AchMil = groupAchMil - groupSoiCon - RelGenSpec, AgrCap = groupAgrCap - groupSoiCon,
 #  BroWil = groupBroWil - groupSoiCon, CicInt = groupCicInt - groupSoiCon,
 #  HolLan = groupHolLan - groupSoiCon, PlaLan = groupPlaLan - groupSoiCon,
 #  PoaCit = groupPoaCit - groupSoiCon,SchAru = groupSchAru - groupSoiCon,  
 #  levels = colnames (design)
 # )

#GLM ####
fit <- glmQLFit(y, design)

#GLM F-test Bin 1 #####
#compare the coefficients of the glm 
# i.e. compare treatments
#DAA_fitBin1_vsC <- glmQLFTest(fit, contrast = my_contrasts [,"B1"]) 
DAA_fittest_RelGenSpec <- glmQLFTest(fit, coef = 10) 


summary(decideTests(DAA_fittest1_vsC))

#Soil test Rel Gen Spec results####
# RelGen data soil ##
DAA_fittest_RelGenSpec_table <-
  DAA_fittest_RelGenSpec$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "RelGen") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test AchMil #####
#compare the coefficients of the glm 
# i.e. compare treatments
#DAA_fitBin1_vsC <- glmQLFTest(fit, contrast = my_contrasts [,"AchMil"]) 
DAA_fittest_AchMil <- glmQLFTest(fit, coef = 2) 


#Soil test 1 results####
# Ach Mil data soil ##
DAA_fittest1_table <-
  DAA_fittest_AchMil$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "AchMil") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

## boxplot - changes in AMF taxa foe change of RelGnSpec = 1
#DAA_fitGen03_table %>% ggplot (aes (x = Genus , y = logFC)) + geom_violin() + geom_point()


#GLM F-test  test 2#####
#compare the coefficients of the glm 
# i.e. compare treatments
#DAA_fittest2_vsC <- glmQLFTest(fit, contrast = my_contrasts [,"B2"]) ## for design matrix ~ 0+ group
DAA_fittest2_vsC <- glmQLFTest(fit, coef = 10) 

summary(decideTests(DAA_fittest2_vsC))

#Soil test 2 results####
# Ach Mil data soil ##
DAA_fittest2_table <-
  DAA_fittest2_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "RelGen") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test test 3 #####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest3_vsC <- glmQLFTest(fit, coef = 3) 

summary(decideTests(DAA_fittest3_vsC))

#Soil test 3 table ####
# Ach Mil data soil ##
DAA_fittest3_table <-
  DAA_fittest3_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "Agr") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test test 4 #####
#compare the coefficients of the glm 
# i.e. compare treatments
#DAA_fitgen06_vsC <- glmQLFTest(fit, contrast = my_contrasts [,"G6"]) 
DAA_fittest4_vsC <- glmQLFTest(fit, coef = 4) 

summary(decideTests(DAA_fittest4_vsC))

#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_fittest4_table <-
  DAA_fittest4_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "Bro") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test Bi 5#####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest5_vsC <- glmQLFTest(fit, coef = 5) 

summary(decideTests(DAA_fittest5_vsC))

#Soil test 5 results####
# Ach Mil data soil ##
DAA_fittest5_table <-
  DAA_fittest5_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "Cic") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)


#GLM F-test test 6#####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest6_vsC <- glmQLFTest(fit, coef = 6) 

summary(decideTests(DAA_fittest6_vsC))

#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_fittest6_table <-
  DAA_fittest6_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "Hol") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)



#GLM F-test test7#####
#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest7_vsC <- glmQLFTest(fit, coef = 7) 

summary(decideTests(DAA_fittest7_vsC))

#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_fittest7_table <-
  DAA_fittest7_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "PlaLan") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test test8#####

#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest8_vsC <- glmQLFTest(fit, coef = 8) 

summary(decideTests(DAA_fittest8_vsC))

#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_fittest8_table <-
  DAA_fittest8_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "PoaCit") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)

#GLM F-test test9#####

#compare the coefficients of the glm 
# i.e. compare treatments
DAA_fittest9_vsC <- glmQLFTest(fit, coef = 9) 

summary(decideTests(DAA_fittest9_vsC))

#Soil PlaSpe results####
# Ach Mil data soil ##
DAA_fittest9_table <-
  DAA_fittest9_vsC$table %>%  
  as_tibble (rownames = "ASV_ID") %>% 
  left_join (taxa_edgeR_gen) %>% 
  filter (Phylum == "Glomeromycota" ) %>% 
  add_column (Gen = "Sch") %>% ### change here 
  mutate (Sign = case_when(PValue<=0.05 ~ "sign.", PValue > 0.05 ~ "ns")) %>%  
  arrange (PValue)





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
ps_edgeR_glo_ASV <- subset_taxa(ps_edgeR_relGenSpec, Phylum == "Glomeromycota")
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
# add information on relative abundance of the ASVs per Glo community!
# add genus label
DAA <- 
  bind_rows (DAA_fittest_RelGenSpec_table, DAA_fittest1_table, DAA_fittest2_table, DAA_fittest3_table, DAA_fittest4_table,
                 DAA_fittest5_table, DAA_fittest6_table, DAA_fittest7_table, DAA_fittest8_table, DAA_fittest9_table) %>%  
 # DAA_fitBin1_table %>% 
  left_join (order_taxa %>% select (Order, GenusLabel, order), by = "GenusLabel") %>% 
  left_join (rel_ASV_abundance_soil)  %>% 
  filter (!is.na (mean_rel_ASV))  # remove all rel abundances and logFC that have 0 abundance in the respective PlaSpe

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
  facet_wrap(~Gen, nrow = 1, scales = "free_x") +
  theme_classic() + 
  theme (axis.title.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.line.y = element_blank(),
         axis.text.y = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color ="darkgrey") + 
  scale_alpha_discrete(range = c(0.3, 1)) +
  theme (legend.position = "right" ) +
 # xlim(-9,13) +
  labs (size = "Mean relative abundance of ASV ", color = "AMF order", alpha = "Significance", 
        title = "RelGenSpec + PlaSpe") 
   


# Plot tree and DAA together #####
# ggarrange (p_tree, plotDAA, nrow = 1, widths = c(0.4,1))

