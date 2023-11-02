#Step 2: Changing the file names to show the sample IDs instead of primer combination.



# working in R again.
#Rename files #####

library(tidyverse)
library(readxl)

getwd()


# eg. C:/Users/RProject_InteractionNiches"
path = "path/to/directory/Lib_primer_removed"  # e.g. C:/Users/RProject_InteractionNiches/Lib_primer_removed
#have a look at the file names
list.files(path)  ## add the appropriate path to the directory containing the resulting files from the last step (see script Part1_primer removal)


# we read the meta file that contains the information about the primer combination used into R for experiment 1 & 2#
meta_names <- read_csv("data/meta_filenames.csv")   


### in some case, we want to use only the information for a specific library for naming the files### For example:
# meta01<-  meta %>% 
#   select (sampleID, filename) ### Do these step for each folder containing files to be named if that is the case.


###generate name list##
# setwd to directory containing the fastq files of one of the experiments (i.e. do this separately for each experiment)
setwd("path/to/directory/Lib_primer_removed")


##R1
#read files are in the format SAMPLENAME.R2.fastq.gz and SAMPLENAME.R2.fastq.gz
File_names <- sort(list.files( pattern = ".R1.fastq.gz", full.names = TRUE))# CHANGE
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)## CHANGE
#als tibble für left join
Primer.name <- as_tibble(Primer.name)


Pr_sample <- left_join(Primer.name, meta_names, by=c("value"="filename")) 


new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz") 
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")

file.rename(from =File_names, to=new.names) ## this renames all files that end with ".R1.fastq.gz" with the sample names followed by the new ending "_R1.fastq.gz2

## R2
# now change all instances of "R1" to "R2" and do the same steps again for the R2 files
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files( pattern = ".R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)
# chnage to a tibble for next step (a left join)
Primer.name <- as_tibble(Primer.name)

#Find correct sample names for each primer combination
Pr_sample <- left_join(Primer.name, meta_names, by=c("value"="filename"))# 


new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz") 
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")

file.rename(from =File_names, to=new.names)

# set working directory 
setwd("your/workingdirectory")

# check that the renaming worked
list.files("path/to/directory/Lib_primer_removed")

## If the files were in separate folder Lib01_primer_removed etc, these steps need to be repeated for each folder and the code changed accordingly.
## IMPORTANT! th raw data file from Lib19 contained sequnces from a different project - these have to be removed for the following scripts to work - all fastq.gz files that 
# were NOT renamed but start with Lib19_ can be deleted
# From now, we assume, we have all fastq files in one folder and the files are named with their sample ID, R01, S01 etc
# 

##


# Note: The description and steps shown here are almost entirely copied from the tutorial 
# https://benjjneb.github.io/dada2/ITS_workflow.html by Benjamin Callahan
# The authors thank Benjamin Callahan for the work as developer and for maintaining excellent tutorials.

# 
#####packages ####
#see installation of dada2 on webpage   
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(tidyverse)
library(readxl)
library(phyloseq); packageVersion("phyloseq")




### DADA2 ###--------------------------------------


###generate name list####
#generate matched lists of the forward and reverse read files, 
#as well as parsing out the sample name. Here we assume forward and reverse 
#read files are in the format SAMPLENAME.R1.fastq.gz and SAMPLENAME.R2.fastq.gz
fnFs <- sort(list.files(path, pattern = ".1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), ".1.fastq.gz"), `[`, 1)
head(sample.names)


###Check primer removal####
#define primer 
FWD <- "CAHCGATGAAGAACGYRG" # ITS3 and ITS4 primers without barcodes
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

## Apply this check to two samples (for that change the values in brackets [[]], now the first file in the list will be used)
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



###Prepare a new filepath for filtering####
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
#the resulting sequence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants.

# For later use, we can save this file in R.
saveRDS(seqtab, "your/workingdirectory/seqtab_Lib.rds") # CHANGE ME to where you want sequence table saved



##Libraries####


meta <- read_xlsx("data/meta_table.xlsx")


###Remove chimeras####
#he core dada method corrects substitution and indel errors, but chimeras remain.
#Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler 
#than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed 
#by combining a left-segment and a right-segment 
#from two more abundant “parent” sequences.

seqtab.nochim<- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)   

old_seqtab_names <- rownames(seqtab.nochim)  #  to get names 
new_seqtab_names <- sapply(strsplit(basename(old_seqtab_names), ".1.fastq.gz"), `[`, 1)

rownames(seqtab.nochim) <- new_seqtab_names

###Inspect distribution of sequence lengths:####

table(nchar(getSequences(seqtab.nochim)))



###Track reads through the pipeline#####

#As a  check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,  getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

track2 <-  as_tibble(as.data.frame( rowSums(seqtab.nochim)), rownames="sampleID") 
colnames (track2) <- c ("sampleID", "no chimeras")
#rownames(track2) <- sample.names

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

track.reads <- as_tibble (track, rownames="sampleID")  %>% 
  left_join (track2)

View(track.reads)

## Save the information on read numbers per step into a file.
track.reads  %>% 
  as.data.frame ()  %>% 
  write_csv  ("your/workingdirectory/track_reads.csv")
# If you run this script for both experiments separately, save the track reads as track_read_E1 and track_reads_E2

#Outside of filtering (depending on how stringent you want to be) there should no step in which a majority 
#of reads are lost.

saveRDS(seqtab.nochim, "your/workingdirectory/seqtab.nochim.rds") # CHANGE ME to where you want sequence table saved

## This script needs to be run for both experiments separately, therefore, use distinguishing names for the
#resulting files to save, for example
# seqtab.nochim.E1.rds
# seqtab.nochim.E2.rds
