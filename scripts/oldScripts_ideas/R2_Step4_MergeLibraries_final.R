##

##Libraries####
library(dada2)
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(tidyverse)
library(readxl)
library(phyloseq); packageVersion("phyloseq")

meta <- read_xlsx("meta_data/meta_table.xlsx")


### merge runs/libraries####
# Merge multiple runs (if necessary)
st03 <- readRDS("data/seqtab_Lib03.rds")
st04 <- readRDS("data/seqtab_Lib04.rds")

st05 <- readRDS("data/seqtab_Lib05.rds")
st06 <- readRDS("data/seqtab_Lib06.rds")

st07 <- readRDS("data/seqtab_Lib07.rds")
st08 <- readRDS("data/seqtab_Lib08.rds")

st09 <- readRDS("data/seqtab_Lib09.rds")
st10 <- readRDS("data/seqtab_Lib10.rds")

st11 <- readRDS("data/seqtab_Lib11.rds")
st12 <- readRDS("data/seqtab_Lib12.rds")

st13 <- readRDS("data/seqtab_Lib13.rds")
st14 <- readRDS("data/seqtab_Lib14.rds")

st15 <- readRDS("data/seqtab_Lib15.rds")
st16 <- readRDS("data/seqtab_Lib16.rds")

st17 <- readRDS("data/seqtab_Lib17.rds")
st18 <- readRDS("data/seqtab_Lib18.rds")

st19 <- readRDS("data/seqtab_Lib19.rds")

st01 <- readRDS("data/seqtab_Lib01.rds")
st02 <- readRDS("data/seqtab_Lib02.rds")



st.all <- mergeSequenceTables(st01, st02, st03, st04, st05, st06, st07, st08, st09, st10, st11, st12, st13, st14, st15, st16, st17, st18,st19)  #etc ergänzen!



###Remove chimeras####
#he core dada method corrects substitution and indel errors, but chimeras remain.
#Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler 
#than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed 
#by combining a left-segment and a right-segment 
#from two more abundant “parent” sequences.

seqtab.nochim<- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)   

old_seqtab_names <- rownames(seqtab.nochim)  #  to get names 
new_seqtab_names <- sapply(strsplit(basename(old_seqtab_names), "_1.fastq.gz"), `[`, 1)

rownames(seqtab.nochim) <- new_seqtab_names

###Inspect distribution of sequence lengths:####
  
table(nchar(getSequences(seqtab.nochim)))



###Track reads through the pipeline#####

#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:
  
track2 <-  as_tibble(as.data.frame( rowSums(seqtab.nochim)), rownames="sampleID") #CHANGE

TR04  <- read_csv("results/track_reads_Lib04.csv",  col_types = cols(sampleID = col_character()))  
TR03  <- read_csv("results/track_reads_Lib03.csv",  col_types = cols(sampleID = col_character()))
TR05  <- read_csv("results/track_reads_Lib05.csv",  col_types = cols(sampleID = col_character()))
TR06  <- read_csv("results/track_reads_Lib06.csv",  col_types = cols(sampleID = col_character()))
TR07  <- read_csv("results/track_reads_Lib07.csv",  col_types = cols(sampleID = col_character()))
TR08  <- read_csv("results/track_reads_Lib08.csv",  col_types = cols(sampleID = col_character()))  
TR09  <- read_csv("results/track_reads_Lib09.csv",  col_types = cols(sampleID = col_character()))  
TR10  <- read_csv("results/track_reads_Lib10.csv",  col_types = cols(sampleID = col_character()))  
TR11  <- read_csv("results/track_reads_Lib11.csv",  col_types = cols(sampleID = col_character()))  
TR12  <- read_csv("results/track_reads_Lib12.csv",  col_types = cols(sampleID = col_character()))  
TR13  <- read_csv("results/track_reads_Lib13.csv",  col_types = cols(sampleID = col_character()))
TR14  <- read_csv("results/track_reads_Lib14.csv",  col_types = cols(sampleID = col_character()))
TR15  <- read_csv("results/track_reads_Lib15.csv",  col_types = cols(sampleID = col_character()))
TR16  <- read_csv("results/track_reads_Lib16.csv",  col_types = cols(sampleID = col_character()))
TR17  <- read_csv("results/track_reads_Lib17.csv",  col_types = cols(sampleID = col_character()))
TR18  <- read_csv("results/track_reads_Lib18.csv",  col_types = cols(sampleID = col_character()))
TR19  <- read_csv("results/track_reads_Lib19.csv",  col_types = cols(sampleID = col_character()))

TR01  <- read_csv("results/track_reads_Lib01.csv",  col_types = cols(sampleID = col_character()))
TR02  <- read_csv("results/track_reads_Lib02.csv",  col_types = cols(sampleID = col_character()))





TR.all <- bind_rows(TR01, TR02,TR03, TR04, TR05, TR06, TR07, TR08, TR09, TR10, TR11, TR12, TR13, TR14, TR15, TR16, TR17, TR18, TR19)  %>% 
     left_join(track2)   %>% filter (!is.na(sampleID)) %>% 
       as.data.frame ()  %>%
  write_csv  ("results/track_reads_ALL.csv")


  
#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
 
  
   
###Assign taxonomy####
  
#DADA2 supports fungal taxonmic assignment using the UNITE database! 
#The DADA2 package provides a native implementation of the naive Bayesian classifier method 
#for taxonomic assignment. The assignTaxonomy function takes as input a set of
#sequences to be classified, and a training set of reference sequences with known taxonomy,
#and outputs taxonomic assignments with at least minBoot bootstrap confidence. For fungal taxonomy, the General Fasta release files
#from the UNITE ITS database can be downloaded and used as the reference.
#Reference: When using this resource, please cite it as follows:
# Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2020): UNITE general FASTA release for Fungi. Version 04.02.2020. UNITE Community. https://doi.org/10.15156/BIO/786368

unite.ref <- "meta_data/sh_general_release_dynamic_10.05.2021_dev.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = FALSE, tryRC = TRUE) #multithread switched off bec Windows
saveRDS(taxa, "data/taxa.rds") 


##Save seqtab and taxa####
saveRDS(seqtab.nochim, "data/seqtab_final.rds") # CHANGE ME to where you want sequence table saved



#The final result, the count matrix of samples (rows) by non-chimeric sequence variants (columns), 
#is stored as as serialized R object. 
#Read it back into R with foo <- readRDS("path/to/run1/output/seqtab.rds").

# save results ####
#mit den ASVs als Namen - brauche ich ertstmal so nicht mehr

#save otu table -full sequence names
asv_table_full  <- t(seqtab.nochim)
write.table(asv_table_full, "results/ASV_table_fullASVs_M2.tsv", sep="\t", quote=F, col.names=NA)

#save taxa table -full_ASVs (full names of sequence )
write.table(taxa, "results/ASVs_taxonomy_fullASVs.tsv", sep = "\t", quote=F, col.names=NA)

#sample table is from meta table
#get meta data
meta_phyloseq.df   <-   meta %>%  
  select (c(sampleID, GSPlaSpe, AC.BC, focal.sp, Lib, m.aboveground, treatment))  %>%  
  filter (Lib != "NA") %>% 
  column_to_rownames(var="sampleID")  %>% 
  as.data.frame()


#save mit veienfachten Namen (ASV_1 etc#)##

#get ASV seq, prep for fasta  table
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# name ASV automatically
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of  final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "results/ASV_newnames.fa") #Fasta file mit ASV_1 und der sequences dazu

# count table with new fasta names
asv_table <- t(seqtab.nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "results/ASV_counts_simple.tsv", sep="\t", quote=F, col.names=NA)

#taxa table with simple names ASV_1 etc.#
taxa_simple  <- taxa
rownames(taxa_simple) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(taxa_simple, "results/ASVs_taxonomy_simple.tsv", sep ="\t", quote = F, col.names = NA)


##same taxa table with removed letters for better figures etc
taxa_2 <- as.matrix(read.table("results/ASVs_taxonomy_simple.tsv", header=T,
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



##Build phylogenetic TREE######
#https://rstudio-pubs-static.s3.#amazonaws.com/345955_fba1ccbdcd8f424aa5505c15bfd75bf7.html
#Lib
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(seqinr)
library(msa)

#Import data - -das sind hier ALLE ASVs! nicht nur Glo.

dna<- readDNAStringSet ("results/ASV_newnames.fa")


ASVs_align  <-msa(dna, "ClustalW")
ASVs_align_forseqinr  <- msaConvert(ASVs_align, type= "seqinr::alignment")

# Citing this package
# If you use this package for research that is published later,
#you are kindly asked to cite it as follows:
#   U. Bodenhofer, E. Bonatesta, C. Horejš-Kainrath, and S. Hochreiter (2015). 
#msa: an R package for multiple sequence alignment. Bioinformatics 31(24):3997-3999. DOI: 10.1093/bioinformatics/btv494.

# Moreover, we insist that, any time you use/cite the package, 
#you also cite the original paper in which the algorithm/method/package
#that you have been using has been introduced:
#   ClustalW:
#   J. D. Thompson, D. G. Higgins, and T. J. Gibson (1994).CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap 
#penalties and weight matrix choice. Nucleic Acids Res., 22(22):4673 4680. DOI: 10.1093/nar/22.22.4673. 
# MUSCLE:
#   R. C. Edgar (2004) MUSCLE: a multiple sequence alignment method with reduced 
#time and space complexity. BMC Bioinformatics, 5(5):113. DOI: 10.1186/1471-2105-5-113.


#Build distance matrix as prep for tree

d <- dist.alignment(ASVs_align_forseqinr, "identity")  #seqinr package
#tree
#neighbor joining
ASV_tree  <- nj(d)  #ape package

#check if nJ was appropriate
# x  <- as.vector(d)
# y  <- as.vector(as.dist(cophenetic(ASV_tree)))
# 
# plot (x, y, xlab = "original distance", ylab = "distance in the tree",
#       main = "Is NJ appropriate?", pch = 20, col = transp("black", 0.1), cex = 3)
# abline(lm((y ~ x), col = "red"))
# cor(x, y)^2
# if points on red line, then fine. looked good enough






#PHYLOSEQ####
ps_ALL <- phyloseq(otu_table (asv_table, taxa_are_rows = T), 
                   tax_table (taxa_2), 
                   sample_data( meta_phyloseq.df)                   , 
                   ASV_tree)
#only roots 
ps_roots_ALL <- subset_samples(ps_ALL, Lib != "02Lib")
saveRDS(ps_ALL, "results/ps_ALL.rds")

# Get phyloseq for Glo only ####

#subset for Glo
psALL_glo  <- subset_taxa(ps_ALL, Phylum == "Glomeromycota" )


ps_M1  <- metagMisc::phyloseq_sep_variable(psALL_glo, "AC.BC")$M1 ## seperates into treatment 
ps_control  <- metagMisc::phyloseq_sep_variable(psALL_glo, "AC.BC")$control ## seperates into treatment 

saveRDS(ps_M1, "results/ps_M1.rds")

saveRDS(ps_control, "results/ps_control.rds")






  



