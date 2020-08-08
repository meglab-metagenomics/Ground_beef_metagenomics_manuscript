## Script to summarize counts
### Code for dataset analysis 
library(data.table)
library(tidyr)
#setwd("~/Dropbox/WRITING/NCBA2_JAN2018/NCBA2_analysis/") 
#source('scripts/Frontiers_NCBA2_analysis.R')
options(scipen = 999) # to decrease the use of scientific notation


### Looking for clinically important AMR genes

amr_group_check <- amr_group_raw[!group %in% snp_regex, ]
important_AMR_regex = c('OXA',
                        'SME',
                        'sme',
                        'IMI',
                        'NDM',
                        'GES',
                        'KPC',
                        'CPHA',
                        'TEM',
                        'SHV',
                        'CTX',
                        'CMY',
                        'VGA',
                        'VGAB',
                        'VGAD',
                        'VATA',
                        'VATB',
                        'VATC',
                        'VATD',
                        'VATE',
                        'CFRA')

amr_raw_important_AMR <- amr_group_check[group %in% important_AMR_regex, ]

View(amr_raw_important_AMR[group %in% important_AMR_regex,  ])

melted_important_AMR <- amr_melted_raw_analytic[ Level_ID =='Group' & Name %in% important_AMR_regex,  ][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Packaging_samples)]#[order(-sum_class )]
melted_important_AMR[,.(sum_important_genes = sum(Normalized_Count))]

ggplot(data = melted_important_AMR , aes(x = Packaging_samples, y = Name)) +
  geom_tile(aes(fill = log(Normalized_Count))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=5) 

#ggsave("~/Dropbox/WRITING/FC_meat_2019/FC_meat_manuscript/FCmeat_figures/Supp_AMR_important_genes.jpeg", width = 30, height = 20, units = "cm")



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
total_amr_percentage

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

# Microbiome
micro_count_summary_per_sample = microbiome_melted_analytic[Level_ID=="Phylum", .(total_counts = sum(Normalized_Count), mean_counts = mean(Normalized_Count), std_dev_counts = sd(Normalized_Count)), by = ID]
micro_count_summary_project = micro_count_summary_per_sample[,.(total_reads = sum(total_counts),mean_per_sample=mean(total_counts), std_dev_counts= sd(total_counts))]
micro_count_summary_project


order_counts = microbiome_melted_analytic[Level_ID=='Order',]
setkey(order_counts,Name)


order_annotations = microbiome_taxonomy[,id := NULL]
setkey(order_annotations, Order)  # Data tables are SQL objects with optional primary keys
order_annotations <- unique(order_annotations, by = key(order_annotations))

setkey(order_annotations,Order)
order_counts <- order_counts[order_annotations][ID != "NA"]
setkey(order_counts,ID)

## Total AMR percentage for study
total_microbiome_percentage <- order_counts[,.(Order_sum = sum(Normalized_Count)), by= .(Name) ][,total:=sum(Order_sum)][,Order_percentage:= signif(Order_sum/total * 100,digits=3)]

microbiome_percentage_perPhylum <- order_counts[,.(sum_hits = sum(Normalized_Count)), by= .(Phylum, Class) ][,Class_percentage:= sum_hits/sum(sum_hits) * 100, by = .(Phylum)]

total_phylum <- microbiome_melted_analytic[Level_ID=="Phylum"][,.(sum_phylum = sum(Normalized_Count)), by = Name][,total := sum(sum_phylum)][,proportion := signif(sum_phylum/total * 100,digits=3) , by = Name]
total_class <- microbiome_melted_analytic[Level_ID=="Class"][,.(sum_class = sum(Normalized_Count)), by = Name][,total := sum(sum_class)][,proportion := signif(sum_class/total * 100,digits=3), by = Name]
total_order <- microbiome_melted_analytic[Level_ID=="Order"][,.(sum_order = sum(Normalized_Count)), by = Name][,total := sum(sum_order)][,proportion := signif(sum_order/total * 100,digits=3), by = Name]



########################
##############
########
###
# METADATA
###
########
##############
########################




