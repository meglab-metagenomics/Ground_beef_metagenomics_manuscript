## Start with this staging file to set up your analysis. 
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
#setwd("")

# Set the output directory for graphs:
graph_output_dir = 'graphs'
# Set the output directory for statistics:
stats_output_dir = 'stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

####################
## File locations ##
####################
## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine. 

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

# Where is the metadata file stored?
amr_metadata_filepath = 'data/resistome_amr++/AMR_metadata.csv'
amr_count_matrix_filepath = 'data/resistome_amr++/strict_SNP_confirmed_AMR_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/resistome_amr++/megares_annotations_v1.03.csv'

#################################
## Microbiome - 16S or kraken? ##
#################################

# Where is the metadata file for the microbiome samples stored?
microbiome_temp_metadata_file = "data/microbiome-qiime2/Microbiome_metadata.csv"

# Location of results extracted from qiime2
biom_file <- "data/microbiome-qiime2/otu_table_json.biom"
tre_file <- "data/microbiome-qiime2/tree.nwk"
tax_fasta <- "data/microbiome-qiime2/aligned-dna-sequences.fasta"
taxa_file <- "data/microbiome-qiime2/taxonomy.tsv"


###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# The following is a list of analyses based on variables in 
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
# NOTE: Exploratory variables cannot be numeric. 

AMR_exploratory_analyses = list(
  # Analysis Store
  # Description: 
  list(
    name = 'Store',
    subsets = list(),
    exploratory_var = 'Blinded_Store',
    order = ''
  ),  
  # Analysis Dilution
  # Description: 
  list(
    name = 'Dilution',
    subsets = list(),
    exploratory_var = 'Dilution',
    order = ''
  ),  
  # Analysis ID
  # Description: 
  list(
    name = 'ID',
    subsets = list(),
    exploratory_var = 'ID',
    order = ''
  ),
  # Analysis 2
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    exploratory_var = 'Treatment',
    order = ''
  ),
  # Analysis 3
  # Description:
  list(
    name = 'Packaging',
    subsets = list(),
    exploratory_var = 'Packaging',
    order = ''
  ),
  # Analysis 3
  # Description:
  list(
    name = 'sample',
    subsets = list(),
    exploratory_var = 'sample',
    order = ''
  )
)



microbiome_exploratory_analyses = list(
  # Analysis Store
  # Description: 
  list(
    name = 'Store',
    subsets = list(),
    exploratory_var = 'Blinded_Store',
    order = ''
  ), 
  # Analysis ID
  # Description: 
  list(
    name = 'ID',
    subsets = list(),
    exploratory_var = 'ID',
    order = ''
  ),
  # Analysis Treatment
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    exploratory_var = 'Treatment',
    order = ''
  ),
  # Analysis Packaging
  # Description:
  list(
    name = 'Packaging',
    subsets = list(),
    exploratory_var = 'Packaging',
    order = ''
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
AMR_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Dilution',
    subsets = list(),
    model_matrix = '~ 0 + Dilution',
    contrasts = list('DilutionNone - DilutionHalf'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  )
)


## Run the analysis
#
## Pick the correct script that handles resistome data and/or microbiome data. 
source('scripts/metagenomeSeq_qiime2.R')
source('scripts/metagenomeSeq_megares.R')

## Run code to make some exploratory figures, output count matrices.
source('scripts/print_AMR_figures.R')
source('scripts/print_microbiome_figures.R')

## Run code to output zero inflated gaussian model results
source('scripts/print_AMR_ZIG_results.R')
source('scripts/print_microbiome_ZIG_results.R')


# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"


##############################################################
# Code for manuscript figures and other statistical analysis #
##############################################################

# Combine metadata with "melted counts"
# set keys first
setkey(amr_melted_raw_analytic,ID)
setkey(amr_melted_analytic,ID)
setkey(microbiome_melted_analytic,ID)
setkey(microbiome_melted_raw_analytic,ID)
setkey(metadata,ID)
setkey(microbiome_metadata,ID)
# Merge objects
microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
microbiome_melted_raw_analytic <- microbiome_melted_raw_analytic[microbiome_metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]

# Explore the 'scripts/ANOSIM.R' script for ordination and ANOSIM tests
source("scripts/ANOSIM.R")
