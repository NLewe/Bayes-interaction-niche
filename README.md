# Plant interaction traits determine the biomass of arbuscular mycorrhizal fungi and bacteria in soil - R-code 

code written by: Natascha Lewe

The recommended start point to use the data is from Part4, i.e. downloading the phyloseq files directly. 


This repository contains the meta data, the R scripts and intermediate results/data files for the accepted manuscript:

Lewe, Natascha; Keyzers, Robert A.; Tylianakis, Jason M.; Deslippe, Julie R. (2024). Plant interaction traits determine the biomass of arbuscular mycorrhizal fungi and bacteria in soil. Ecology



Our goal was to compare the plant-AMF interactions among plants and experimental conditions using sequencing of the AMF community in roots and soil.
We ran two glasshouse experiments, each with eight grassland/pasture plant species in replicates of five. We sourced the soil for the two experiments from a former pasture and mixed it in different ratios with sand and other growth media, depending on the experiment. This was to ensure different growth conditions and AMF pools.
We extracted DNA from roots and soil and sequenced the ITS2 region to identify the fungal communities. We filtered for Glomeromycotina (arbuscular mycorrhizal fungi, AMF) and we determined the biomass of the AMF and bacteria in the soil by phospholipid fatty acid analysis and neutral lipid fatty acid analysis.
By calculating measures of AMF diversity, we were able to estimate the plants' interaction generalism with AMF. We then used a metric of generalism to examine the relationship between the interaction trait and the biomass of the soil microbial community by Bayesian analysis.


The repository contains the following parts:

## Part 1 Primer removal
- This file contains code needed to remove primer sequences from raw fastq.gz files. 
- Linux environment needed
- output: fastq files

#### Files needed:
- Part1_primer_removal.txt (script to remove primers from sequence data)
- all fastq.gz datafiles from https://www.ncbi.nlm.nih.gov/home/download/ , using the BioProject ID:      PRJNA997080 



## Part 2 Renaming and dada2 pipeline
- This file contains the R-code to first name the files based on their sample name, followed by inferring amplicon sequence variants (ASVs) using the DADA2 pipeline. Sequences are merged and the steps of the dada2 are tracked by counting the number of sequence reads after each step.
- output: a sequence table for each experiment, saved as "seqtab.rds" , a table ("track_reads.csv") tracking the read numbers per sample for each step of the pipeline.
The code is given in a form that allows the user to change it to its own directory structures: change the mention of "path/to/directory" or similar to their directory structure.

#### Files needed
- Part2_RenameSamples_Inference_ASVs.R
- fastq.gz files resulting from Step1
- meta_filenames.csv



## Part 3 Identification and tree building
- The R scripts contains the identification of the fungal sequences using the UNITE database.
- A phylogenetic tree for the plant species is built. 
- A phyloseq object is build as output, including the phylogenetic tree. Phyloseq is a well known and often used package for analysis of phylogenetic datasets. 
- A coverage based rarefaction is applied to each dataset.


#### Files needed
- Part3_Identification_Tree
- seqtab.nochim.E1.rds (resulting from former step: Part2)
- seqtab.nochim.E2.rds (resulting from former step: Part2)
- meta_tables: 
    M0_meta.xlsx (meta for E1)
    meta_M1wcontrol.xlsx and meta_M1.csv (meta for E2)
- unite database file (downloadable at UNITE website)



## Part 4 - Diversity metrics & interaction generalism ###
All diversity metrics are calculated both at the level of an individual plant and at the plant species level.


#### Files needed 
- Part4_diversity_metrics_generalism.R
- ps_E1.rds
- ps_E2.rds
- meta_plants
- meta_M1wcontrol.xlsx
- meta_M1.csv
- M0_meta.xlsx (M0 and M1 are historically used names for experiment E1 and E2)
- fasta files for plant genes: concat_align.fa 


## Part 5 - Biomass of soil bacteria and AMF ####
The biomass data for the AMF in the soil and for the bacteria in the soil is calculated.

#### Files needed
- Part5_biomass.R
- dataNLrootstosoil.rds
- meta files:
  - meta_plants
  - meta_M1wcontrol.xlsx
  - meta_M1.csv
  - meta_soil_control.xlsx
- PL_FA_conc_Soil.rds
- SoilControlPL.xlsx
- SoilControlNL.xlsx
- FAME_RRF_biomarker.xlsx


## Part 6  -  Radarplots and Procrustes analysis ####
Figure 1 of the publication is calculated as well as the Procrustes analysis comparing the similarity of interaction niches as defined by diversity metrics. For each plant , the absolute values of the diversity metrics are used and the values compared between the two experiments. This is done by permutational Procrustes analysis.


#### Files needed 
- Part6_Figure1_Radarplots.R
- output from part 4 (All_metrics_E1, All_metrics_E2)



## Part 7 Principal component analysis (PCA) ####
based on the per sample values for the diversity metrics, principal components are calculated and extracted for use in bayesian models.

#### Files needed
- Part7_PCA.R
- output from parts 4


## Part 8 Bayesian model and figure ###
The effect of the interaction generalism on the soil AMF biomass is explored using Bayesian models. 
A plot showing the conditional effect of the plant biomass is calculated.

#### Files needed
- Part8_AMFbiomass_Bayes.R
- outputs from parts 5 and 7


## Part 9 ####
The impact of the interaction generalism on the soil bacterial biomass is explored by Bayesian models. 
A plot describing the best resulting model is shown.
Note that not all possible variations of the models are shown.

#### Files needed
- Part9_bacterial_biomass_Bayes.R
- output from parts 5 and 6


## Part 10 ####
- R code for tables etc. of the Appendix S2.

#### Files needed 
- Part10_AppendixS2.R
- ps_E1.rds
- ps_E2.rds
- track_reads_E1.csv (output from part 1&2)
- ps_AllASVs_E1.rds (output from part 3)
- iNext_ALL_E1_8_sp_PlaSpe_richness.rds (output from part 3)


## Part 11 Bayesian models using diversity metrics - for Appendix S2###

#### Files needed
- output from parts 5 and 7
- Part11_Biomass_diversity_Bayes.R


