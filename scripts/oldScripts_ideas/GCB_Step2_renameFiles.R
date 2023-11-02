#Step 2: Changing the file names toshow the sample IDs instead of primer combination.

# working in R again.
#Rename files

library(tidyverse)
library(readxl)

#Lib01 #####
#have a look at the file names
list.files("docs/SubmissionData/fasrqDaaExp2/Lib02merged")

# we need the meta file
meta <- read_excel("data/meta_table_seq.xlsx") %>% 
  mutate (Libr = Lib) %>% 
  unite ("filename", c(Libr, Primer.comb), sep = "_" )

###get meta of appropriate Lib###
meta02<-  meta %>%  #CHANGE Lib'
  select(sampleID, filename, Lib) %>%  
  filter (Lib== "Lib02") %>% ##CHANGE
  select (sampleID, filename)

###generate name list##
# setwd to directory containing the fastq files
setwd("C:/Users/Nat/OneDrive - Victoria University of Wellington - 
      STAFF/Publications/Paper_Ch2_interaction_niches/Interaction_niche_AMFbiomass/docs/SubmissionData/fasrqDaaExp2/Lib02merged")
##DO FOR BOTH R1 and R2 files!! ###

#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files( pattern = ".R1.fastq.gz", full.names = TRUE))# CHANGE
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R1.fastq.gz"), `[`, 1)## CHANGE
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Zuordnung der Primer combination zu sample names
Pr_sample <- left_join(Primer.name, meta02, by=c("value"="filename"))# CHANGE meta01, 02 etc


new.names <-  paste0(Pr_sample$sampleID, "_1.fastq.gz") # CHANGE
old.names <-  paste0(Pr_sample$value, ".R1.fastq.gz")# CHANGE

file.rename(from =File_names, to=new.names)


# now change all R1 to R2 and do again
#read files are in the format SAMPLENAME_R2.fastq.gz and SAMPLENAME_R2.fastq.gz
File_names <- sort(list.files( pattern = ".R2.fastq.gz", full.names = TRUE))# CHANGE
# Extract sample names, assuming filenames have format: 
Primer.name <- sapply(strsplit(basename(File_names), ".R2.fastq.gz"), `[`, 1)## CHANGE
#als tibble für left join
Primer.name <- as_tibble(Primer.name)

#Finf correct sample names for each primer combination
Pr_sample <- left_join(Primer.name, meta02, by=c("value"="filename"))# CHANGE meta01, 02 etc


new.names <-  paste0(Pr_sample$sampleID, "_2.fastq.gz") # CHANGE
old.names <-  paste0(Pr_sample$value, ".R2.fastq.gz")# CHANGE

file.rename(from =File_names, to=new.names)

# set working directory 
setwd("C:/Users/Nat/OneDrive - Victoria University of Wellington - STAFF/Publications/Paper_Ch2_interaction_niches/Interaction_niche_AMFbiomass")

# check that the renaming worked
list.files("path/to/file/directory")

## These renaming steps can be either done library wise or for all samples together,
#if the fastq files are in the same folder - and the meta file includes the samples from all libraries.





### send data to Rapoi:
# open Git bash Locally!
#scp /C/Users/Nat/'OneDrive - Victoria University of Wellington - STAFF'/PhD/'Lab groups'/R_Lesleys_GCB/data/GCB_All_demultiplexed/*.gz lewena@raapoi.vuw.ac.nz:/nfs/home/lewena/GCB/GCB_seq/

#scp lewena@raapoi.vuw.ac.nz:/nfs/home/lewena/*sh /C/Users/Nat/'OneDrive - Victoria University of Wellington - STAFF'/PhD/'Lab groups'/R_Lesleys_GCB/scripts/