#Renaming samples and dada2 - for all Libraries#####
# Note: The description and steps shown here are almost entirely copied from the tutorial 
# 
#####Library####
#see installation of dada2 on webpage
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(tidyverse)
library(readxl)


## LIB01####
#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib01_merged/CA_primer_01"
list.files(path)
meta <- read_xlsx("meta_data/meta_table.xlsx")

###get meta of Lib####
metaLib01<-  meta %>%  
  select(sampleID, filename, Lib) %>%  
  filter (Lib== "Lib01") %>% 
  select (sampleID, filename)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib01, by=c("value"="filename"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib01, by=c("value"="filename"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)

### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


###Check primer removal####
#man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
#define primer GCATCGATGAAGAACGCAGC
FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
REV <- "TCCTSCGCTTATTGATATGC"
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

#Show the orientations
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients



###Check number of primer in one forward and one reverse read sample
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

##Check for two samples (the ones in brackets)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))


#Note: Orientation mixups are a common trip-up. 
#If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
#then replace REV with its reverse-complement orientation 
#(REV <- REV.orient[["RevComp"]]) before proceeding

###


###Inspect read profiles ####

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
#quality profile plot is a gray-scale heatmap of the frequency 
#each quality score at each base position. 
#The median quality score at each position is shown by the green line, 
#the quartiles of the quality score distribution by the orange 
#The red line shows the scaled proportion of 
#that extend to at least that position. (i.e. length of read)



###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)



###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib01.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib01.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.



####LIB02####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib02_merged/CA_primer_02_test"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib02<-  meta %>%  
  select(sampleID, filename, Lib) %>%  
  filter (Lib== "Lib02") %>% 
  select (sampleID, filename)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib02, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib02, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
REV <- "TCCTSCGCTTATTGATATGC"
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

#Show the orientations
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

##Check for two samples (the ones in brackets)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib02.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib02.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.


####LIB03####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib03_merged/CA_primer_03"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib03<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "03Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib03, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib03, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib03.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib03.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB04####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib04_merged/CA_primer_04"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib04<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "04Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib04, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib04, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib04.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib04.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB05####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib05_merged/CA_primer_05"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib05<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "05Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib05, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib05, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib05.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib05.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB06####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib06_merged/CA_primer_06"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib06<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "06Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib06, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib06, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib06.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib06.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB07####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib07_merged/CA_primer_07"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib07<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "07Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib07, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib07, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib07.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib07.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB08####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib08_merged/CA_primer_08"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib08<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "08Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib08, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib08, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib08.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib08.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.



####LIB09####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib09_merged/CA_primer_09"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib09<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "09Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib09, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib09, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib09.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib09.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.



####LIB10####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib10_merged/CA_primer_10"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib10<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "10Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib10, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib10, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib10.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib10.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB11####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib11_merged/CA_primer_11"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib11<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "11Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib11, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib11, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib11.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib11.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB12####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib12_merged/CA_primer_12"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib12<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "12Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib12, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib12, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib12.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib12.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB13####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib13_merged/CA_primer_13"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib13<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "13Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib13, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib13, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib13.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib13.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB14####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib14_merged/CA_primer_14"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib14<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "14Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib14, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib14, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib14.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib14.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB15####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib15_merged/CA_primer_15"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib15<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "15Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib15, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib15, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib15.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib15.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

####LIB16####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib16_merged/CA_primer_16"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib16<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "16Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib16, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib16, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib16.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib16.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB17####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib17_merged/CA_primer_17"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib17<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "17Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib17, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib17, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib17.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib17.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB18####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib18_merged/CA_primer_18"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib18<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "18Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib18, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib18, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib18.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib18.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.
####LIB19####

#Files anschauen'
path = "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data_demux_merged/Lib19_merged/CA_primer_19"
list.files(path)
#meta <- read_xlsx("meta_data/meta_table.xlsx")

## meta of Lib####
metaLib19<-  meta %>%  
  select(sampleID, Primer.comb, Lib) %>%  
  filter (Lib== "19Lib") %>% 
  select (sampleID, Primer.comb)

###generate name list####
setwd(path)
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files(path, pattern = ".R1.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)#anpassen an filesa
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib19, by=c("value"="Primer.comb"))

new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =old.names, to=new.names)

## For R2 files#
File_names <- sort(list.files(path, pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)#anpassen an files
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, metaLib19, by=c("value"="Primer.comb"))
new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz")
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =old.names, to=new.names)


### DADA2 ###--------------------------------------
###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)
head(sample.names)


# ###Check primer removal####
# #man kann so auch checken, wo sich die Primer im read befinden: in reverse hinten? etc
# #define primer GCATCGATGAAGAACGCAGC
# FWD <- "CAHCGATGAAGAACGYRG" #okay, das ist nur der primer ohne barcodes
# REV <- "TCCTSCGCTTATTGATATGC"
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# 
# #Show the orientations
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# REV.orients
# 
# 
# 
# ###Check number of primer in one forward and one reverse read sample
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# ##Check for two samples (the ones in brackets)
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[2]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[2]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[2]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[2]]))
# 
# 
# #Note: Orientation mixups are a common trip-up. 
# #If, for example, the REV primer is matching the Reverse reads in its RevComp orientation, 
# #then replace REV with its reverse-complement orientation 
# #(REV <- REV.orient[["RevComp"]]) before proceeding
# 
# ###Inspect read profiles ####
# 
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# #quality profile plot is a gray-scale heatmap of the frequency 
# #each quality score at each base position. 
# #The median quality score at each position is shown by the green line, 
# #the quartiles of the quality score distribution by the orange 
# #The red line shows the scaled proportion of 
# #that extend to at least that position. (i.e. length of read)


###Prepare filepath for filtering####
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# on windows, set multithread = FALSE
View(out)

#standard filtering paraments: maxN=0 
#(DADA2 requires sequences contain no Ns), truncQ = 2,
#rm.phix = TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here,
#to get rid of spurious very low-length sequences. 
#f too few reads are passing the filter, consider relaxing maxEE, 
#especially on the reverse reads (eg. maxEE=c(2,5)

###Error rates####
errF <- learnErrors(filtFs,nbases = 1e8, multithread=FALSE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=FALSE)

#plotErrors(errF, nominalQ=TRUE)
#
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#Are the estimated error rates (black line) a good fit to the observed rates (points)? 
#Do the error rates drop with increased quality as expected?
#If everything looks reasonable, we proceed with confidence.


###!!Dereplicate identical reads####

#Dereplication combines all identical sequencing reads into into “unique sequences” with a 
#corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#DADA2 retains a summary of the quality information associated with each unique sequence. 
#The consensus quality profile of a unique sequence is the average of the positional 
#qualities from the dereplicated reads. These quality profiles inform the error model 
#of the subsequent sample inference step, 
#significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

###Sample Inference####
#DADA2 infers sample sequences exactly and resolves differences of as little as 1 nucleotide.
#At this step, the core sample inference algorithm is applied to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool ="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")
names (derepFs[-31])
#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant

###Merge paired reads####
#merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs[-31], derepFs[-31], dadaRs[-31], derepRs[-31], verbose=TRUE)

###Construct Sequence Table####

#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#equence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

saveRDS(seqtab, "C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/data/seqtab_Lib19.rds") # CHANGE ME to where you want sequence table saved

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))


# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track.reads <- as_tibble (track, rownames="sampleID")  
View(track.reads)

## Save
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/PhD/Results/Marsden2/M2_seq/R_M2_seq/R_M2/results/track_reads_Lib19.csv")

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

## From here, we go into the script mergeLibraries.R 


