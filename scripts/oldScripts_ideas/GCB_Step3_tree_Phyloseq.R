# Prepare Tree, add funGuild, and save a phyloseq object####


#Here we build a tree - which is also a very time consuming step ###


library (msa)
library(seqinr)
library (ape)
library(adegenet)
library (phyloseq)
library(tidyverse)
library(readxl)

#Tree building ####
#Import fasta file as dna sequences for later alignment
dna_All  <- readDNAStringSet("results/Pipits_2021/3_GCB_process/repseqs.fasta")
#dna_All  <- fasta2DNAbin("results/ASV_sequences.fa")


#Alignment with ClustalW###
ASVs_align  <-msa(dna_All, "ClustalW")
#Convert for different package, here seqinr
ASVs_align_forseqinr  <- msaConvert(ASVs_align, type= "seqinr::alignment")
#
#Build distance matrix as prep for tree###
d_mat <- seqinr::dist.alignment(ASVs_align_forseqinr, "identity")  #seqinr package

#d_mat <- dist.dna (dna_All, model = "TN93")
#phylogenetic tree with neighbor joining
ASV_tree  <- njs(d_mat)
write.tree (ASV_tree, "results/OTU_tree.tre")

# validation of NJ tree
#check if NJ was appropriate:
#Careful, this takes a long time!
# test parts of data instead, e.g. plot (x[1:5000],y[1:5000], ......)
x  <- as.vector(d_mat)
y  <- as.vector(as.dist(cophenetic(ASV_tree)))

plot (x, y, xlab = "original distance", ylab = "distance in the tree",
      main = "Is NJ appropriate?", pch = 20,col = transp("black", 0.1), cex = 3)
abline(lm(y ~ x))
cor(x, y)^2
# points along line? ->okay

# Phyloseq####
#Import otu table and representative seq into phyloseq
ps_ALL <- import_biom(BIOMfilename = "results/Pipits_2021/3_GCB_process/otu_table.biom" , 
                      refseqfilename = "results/Pipits_2021/3_GCB_process/repseqs.fasta",
                      treefilename = "results/OTU_tree.tre")

#add meta table to phyloseq (change to data frame first)
meta.df <- meta %>% data.frame(row.names ="sampleID")
ps_ALL  <- merge_phyloseq(ps_ALL, sample_data(meta.df))


saveRDS (ps_ALL, "results/Pipits_2021/GCB_phyloseq_ALL.rds")


## All preparations are done

# Go to LW.Rmd



