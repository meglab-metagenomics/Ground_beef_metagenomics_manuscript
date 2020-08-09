## Script to summarize counts
### Code for dataset analysis 
library(data.table)
library(tidyr)
#setwd("~/Dropbox/WRITING/NCBA2_JAN2018/NCBA2_analysis/") 
#source('scripts/Frontiers_NCBA2_analysis.R')
options(scipen = 999) # to decrease the use of scientific notation


full_microbiome_counts <- as.data.table(MRcounts(microbiome), keep.rownames = TRUE)
names(full_microbiome_counts)[names(full_microbiome_counts) == "rn"] <- "feature"
setkey(full_microbiome_counts, feature)

#write.csv(full_microbiome_counts[microbiome_taxonomy],"full_microbiome_counts.csv",row.names = FALSE)

##
###
#### Comparison of Sequencing read counts
###
##

Treatment_read_count_comparison <- metadata[,.(total_raw_reads = sum(Raw_paired_reads),mean_raw_reads = mean(Raw_paired_reads), sd_raw_reads = sd(Raw_paired_reads)), by=Treatment]

Packaging_read_count_comparison <- metadata[,.(total_raw_reads = sum(Raw_paired_reads),mean_raw_reads = mean(Raw_paired_reads), sd_raw_reads = sd(Raw_paired_reads), min_raw_reads = min(Raw_paired_reads),max_raw_reads = max(Raw_paired_reads)), by=Packaging]

ggplot(microbiome_metadata, aes(x=Packaging, y=Raw_paired_reads, color=Packaging)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1) , dotsize=1)
ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supp_Microbiome_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")

ggplot(metadata, aes(x=Packaging, y=Raw_paired_reads, color=Packaging)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1) , dotsize=1)
ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supp_AMR_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")


##
### AMR majority of counts
##

gene_counts = amr_melted_analytic[Level_ID=='Gene',]
setkey(gene_counts,Name)
setkey(annotations,header)
gene_counts <- gene_counts[annotations][ID != "NA"]
setkey(gene_counts,ID)

## Total AMR percentage for study
total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class) ][,total:=sum(class_sum)][,class_percentage:= class_sum/total * 100]

total_amr_percentage[, sum(class_percentage)]

## Mechanism proportion by class
mech_percentage_byclass <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(mechanism,class) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(class)]
mech_percentage_byclass

## Total AMR abundance by treatment group
AMR_counts_by_Treatment <- gene_counts[,.(class_sum = sum(Normalized_Count), mean_class = mean(Normalized_Count), sd_class = sd(Normalized_Count)), by= .(Treatment) ]
AMR_counts_by_Treatment

## Total AMR abundance by Packaging
AMR_counts_by_Packaging <- gene_counts[,.(class_sum = sum(Normalized_Count), mean_class = mean(Normalized_Count), sd_class = sd(Normalized_Count)), by= .(Packaging) ]
AMR_counts_by_Packaging

# AMR percentage by Treatment
treatment_class_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(Treatment, class) ][,class_percentage:= class_sum/sum(class_sum) * 100, by = Treatment]
treatment_class_amr_percentage

# Total AMR percentage
total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class,mechanism) ][,class_percentage:= class_sum/sum(class_sum) * 100]
total_amr_percentage
total_amr_percentage[, .(sum = sum(class_percentage))]

## AMR percentage of each mechanism per class, in complete dataset
mech_percentage_byclass <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(mechanism,class) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(class)]
mech_percentage_byclass[,.(test=sum(class_sum)), by=.(class, mechanism)]
## AMR classes by sample group
AMR_classes_percent_bygroup <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class,Group) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(Group)]






########################
##############
########
###
# MICROBIOME
###
########
##############
########################

## Microbiome and resistome Total, mean, std, and median counts per sample

# Microbiome mapped reads
micro_count_summary_per_sample = microbiome_melted_raw_analytic[Level_ID=="Domain", .(total_counts = sum(Normalized_Count), mean_counts = mean(Normalized_Count), std_dev_counts = sd(Normalized_Count)), by = ID]
micro_count_summary_project = micro_count_summary_per_sample[,.(total_reads = sum(total_counts),mean_per_sample=mean(total_counts), std_dev_counts= sd(total_counts),min_per_sample=min(total_counts),max_per_sample=max(total_counts))]
micro_count_summary_project

# Number of features
microbiome
microbiome_analytic_data

# Calculate taxa classification proportions
microbiome_melted_raw_analytic <- as.data.table(microbiome_melted_raw_analytic)


taxa_level_counts <- microbiome_melted_raw_analytic[Normalized_Count > 0 ,.(sum_total_level = sum(Normalized_Count), count_uniq = uniqueN(Name)), by = .(Level_ID)]
taxa_level_counts[,proportion_of_total := sum_total_level/2496913]

taxa_level_counts

## Taxa proportion of counts at each level by Treatment
micro_counts <- microbiome_melted_analytic[,.(sum_taxa = sum(Normalized_Count) ), by = .(Level_ID, Treatment)]
micro_counts[Treatment =="RWA",percentage := (sum_taxa/1166527.90) * 100, by = .(Level_ID)]
micro_counts[Treatment =="CONV",percentage := (sum_taxa/1190749.69) * 100, by = .(Level_ID)]
micro_counts




# Group by Order
order_counts = microbiome_melted_analytic[Level_ID=='Order',]
setkey(order_counts,Name)
order_annotations = microbiome_taxonomy[,id := NULL]
setkey(order_annotations, Order)  # Data tables are SQL objects with optional primary keys
order_annotations <- unique(order_annotations, by = key(order_annotations))
setkey(order_annotations,Order)
order_counts <- order_counts[order_annotations][ID != "NA"]
setkey(order_counts,ID)

## Total AMR percentage for study
microbiome_percentage_perPhylum <- order_counts[,.(sum_hits = sum(Normalized_Count)), by= .(Phylum, Class) ][,Class_percentage:= sum_hits/sum(sum_hits) * 100, by = .(Phylum)]
microbiome_percentage_perPhylum

total_phylum <- microbiome_melted_analytic[Level_ID=="Phylum"][,.(sum_phylum = sum(Normalized_Count)), by = Name][,total := sum(sum_phylum)][,percentage := sum_phylum/total * 100 , by = Name]
total_phylum 

total_class <- microbiome_melted_analytic[Level_ID=="Class"][,.(sum_class = sum(Normalized_Count)), by = Name][,total := sum(sum_class)][,percentage := sum_class/total * 100, by = Name]
total_order <- microbiome_melted_analytic[Level_ID=="Order"][,.(sum_order = sum(Normalized_Count)), by = Name][,total := sum(sum_order)][,percentage := sum_order/total * 100, by = Name]



########################
##############
########
###
# METADATA
###
########
##############
########################




