#### Part3 Assign fungal taxonomy 

## Note that this script needs to be run twice, once for each experiment. 
# For that, use the files seqtab.nochim.E1.rds and seqtab.nochim.E2.rds instead of the placeholder "seqtab.nochim" #

# Packages ####
library(dada2)
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(tidyverse)
library(readxl)
library(phyloseq); packageVersion("phyloseq")

library (metagMisc)

seqtab.nochim  <-readRDS ("seqtab.nochim.rds")

###Assign taxonomy####

#DADA2 supports fungal taxonmic assignment using the UNITE database! 
#The DADA2 package provides a native implementation of the naive Bayesian classifier method 
#for taxonomic assignment. The assignTaxonomy function takes as input a set of
#sequences to be classified, and a training set of reference sequences with known taxonomy,
#and outputs taxonomic assignments with at least minBoot bootstrap confidence. For fungal taxonomy, the General Fasta release files
#from the UNITE ITS database can be downloaded and used as the reference.
#Reference: When using this resource, please cite it as follows:
# Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2020): UNITE general FASTA release for Fungi. Version 04.02.2020. UNITE Community. https://doi.org/10.15156/BIO/786368

unite.ref <- "data/sh_general_release_dynamic_10.05.2021_dev.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = FALSE, tryRC = TRUE) #multithread switched off bec Windows


##Save seqtab and taxa####
saveRDS(seqtab.nochim, "data/seqtab_final.rds") 
saveRDS(taxa, "data/taxa.rds") 


# The final result, the count matrix of samples (rows) by non-chimeric sequence variants (columns), 
# is stored as as serialized R object. 


# save  ASV tables ####
# Table using the full sequences as rownames:

asv_table_full  <- t(seqtab.nochim)
write.table(asv_table_full, "data/ASV_table_fullASVs.tsv", sep="\t", quote=F, col.names=NA)

# Save taxa table -full_ASVs (full names of sequence )
write.table(taxa, "data/ASVs_taxonomy_fullASVs.tsv", sep = "\t", quote=F, col.names=NA)


# Save mit simplified names(ASV_1 etc#)##

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# name ASV automatically
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# Writing out a fasta file of the final ASV seqs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "data/ASV_newnames.fa") #Fasta file of the new (short) names and the full sequences

# Count table with the new fasta names
asv_table <- t(seqtab.nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "data/ASV_counts_simple.tsv", sep="\t", quote=F, col.names=NA)

# Taxa table with simple names ASV_1 etc.#
taxa_simple  <- taxa
rownames(taxa_simple) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(taxa_simple, "data/ASVs_taxonomy_simple.tsv", sep ="\t", quote = F, col.names = NA)


# Same taxa table after removal of the prefixes of all taxa for easier reading
taxa_2 <- as.matrix(read.table("data/ASVs_taxonomy_simple.tsv", header=T,
                               row.names=1, check.names=F, sep="\t")) %>%  
  as_tibble(rownames = "ASV_ID") %>% 
  mutate (Phylum=str_remove(Phylum,"p__"))  %>% 
  mutate (Class=str_remove(Class,"c__")) %>% 
  mutate (Order=str_remove(Order,"o__"))  %>% 
  mutate (Family=str_remove(Family,"f__"))  %>% 
  mutate (Genus=str_remove(Genus,"g__")) %>% 
  mutate (Species=str_remove(Species,"s__")) %>%  
  as.data.frame ( ) 
rownames (taxa_2) <- taxa_2$ASV_ID
taxa_2  <- taxa_2[,-1]
taxa_2  <- as.matrix(taxa_2)



# Build phylogenetic tree ######
#https://rstudio-pubs-static.s3.#amazonaws.com/345955_fba1ccbdcd8f424aa5505c15bfd75bf7.html
## Packages ####
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(seqinr)
library(msa)

# Note: The building of the phyloseq objects and subsequently of the phylogenetic trees was done separately
# for the two experiments of the mentioned publication. 
# The reason for that was that these experiments were run in different years and the bioinformatics were calculated after sequencing.
# 


meta_table <- read_xlsx("data/meta_table.xlsx")

# A phyloseq object is build 
# Get a sample table for the phyloseq oblect from meta table
# get meta data
meta_phyloseq.df   <-   meta_table %>%  
  column_to_rownames(var="sampleID")  %>% 
  as.data.frame()


dna<- readDNAStringSet ("data/ASV_newnames.fa")


ASVs_align  <-msa(dna, "ClustalW")
ASVs_align_forseqinr  <- msaConvert(ASVs_align, type= "seqinr::alignment")

# Ctation of the package
# U. Bodenhofer, E. Bonatesta, C. Horejš-Kainrath, and S. Hochreiter (2015). 
# msa: an R package for multiple sequence alignment. Bioinformatics 31(24):3997-3999. 
# DOI: 10.1093/bioinformatics/btv494.
#   ClustalW:
# J. D. Thompson, D. G. Higgins, and T. J. Gibson (1994).CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap 
#penalties and weight matrix choice. Nucleic Acids Res., 22(22):4673 4680. DOI: 10.1093/nar/22.22.4673. 
# MUSCLE:
#  R. C. Edgar (2004) MUSCLE: a multiple sequence alignment method with reduced 
#time and space complexity. BMC Bioinformatics, 5(5):113. DOI: 10.1186/1471-2105-5-113.


# Build distance matrix as preparation for tree building

d <- dist.alignment(ASVs_align_forseqinr, "identity")  #seqinr package
# Build a phylogenetic tree
# by neighbor joining
ASV_tree  <- nj(d)  #ape package

# To check if neighbor joining is the appropriate approach, run:
# x  <- as.vector(d)
# y  <- as.vector(as.dist(cophenetic(ASV_tree)))
# 
# plot (x, y, xlab = "original distance", ylab = "distance in the tree",
#       main = "Is NJ appropriate?", pch = 20, col = transp("black", 0.1), cex = 3)
# abline(lm((y ~ x), col = "red"))
# cor(x, y)^2
# If points fall approximately on the red line, then the approach was appropriate.





# Build PHYLOSEQ ####
# Use the tables made before and combine them into a phyloseq objct #
ps_ALL <- phyloseq(otu_table (asv_table, taxa_are_rows = T), 
                   tax_table (taxa_2), 
                   sample_data( meta_phyloseq.df)                   , 
                   ASV_tree)

saveRDS (ps_ALL, "data/ps_ALL.rds")


# Rarefy to coverage ####
# Run 99 iterations of the rarefaction function, which produces 99 phyloseq objects #
# This is a very time consuming step!!

 ps_ALL_iter <-  phyloseq_coverage_raref(ps_ALL, iter =99) #produces 99 ps objects in list
 write_rds (ps_ALL_iter, "data/ps_ALL_iter99.rds") ## save this object 

 
# ps_ALL_iter <- read_rds("data/ps_ALL_iter99.rds")
 # From those 99 phyloseq objects, we need to get one mean result
 #get mean result of 99 iterations of resampling
 ALL_coverage_tables <-
   map(ps_ALL_iter, otu_table) %>% # all ps objects are changed into tables
   map (data.frame) %>%   map (t) %>% # are tables are changes into data frames and transposed
   map (function (df) as_tibble (df,rownames = "sampleID" )) %>%
   map_dfr (bind_rows, .id = "TableNumber") %>%
   pivot_longer(!c(TableNumber, sampleID), names_to =  "ASV_ID", values_to = "counts") %>%
   group_by (sampleID, ASV_ID) %>% select (!TableNumber) %>% 
   summarize (mean_counts = ceiling (mean(counts))) %>% # we take the mean
   unique () %>%
   pivot_wider(id_cols = , names_from = "sampleID", values_from = "mean_counts") %>% # we rebuild the tabel again 
   data.frame()
 rownames (ALL_coverage_tables)<- ALL_coverage_tables$ASV_ID
 ALL_coverage_tables <- ALL_coverage_tables[,-1]
# 
# 
# Now our original ps_ALL object will be changed to the resampled version we just produced

 otu_table (ps_ALL) <- otu_table(ALL_coverage_tables, taxa_are_rows = T)



# Subset the data to include sequences from the Glomeromycota (AMF) only 
ps_glo  <- subset_taxa(ps_ALL, Phylum == "Glomeromycota" )


# Remove samples that do not include any AMF sequences 
# This is the example for experiment 1, change E1 to E2 for the data of experiment 2
# The following analysis does not make use of the soil samples from experiment 1, they can be removed with

ps_E1 <- subset_samples (ps_glo, roots_soil == "roots" ) %>% prune_taxa( taxa_sums (.)>0, .)
ps_E1 <- prune_samples(sample_sums (ps_glo) > 0, ps_glo)
# 



# Subset the phyloseq objects
# At each subsetting step, taxa need to be pruned from the dataset, as not all data subsets include 
# the same sequences. Sequences with less than 1 count over all samples are removed.
# For experiment 2.
ps_E2 <- subset_samples (ps_E2, roots_soil == "roots") %>% prune_taxa( taxa_sums (.)>0, .)
#ps_E2_soil <- subset_samples (ps_E2, roots_soil == "soil") %>%  prune_taxa (taxa_sums (.)>0, .)  # if needed


saveRDS(ps_E1, "data/ps_E1.rds")

saveRDS(ps_E2, "data/ps_E2.rds")










