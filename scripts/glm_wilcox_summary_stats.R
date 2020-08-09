#
##
### GLM for read counts by Treatment, Dilution, Packaging
##
#

summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Packaging)))

summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Packaging)))

summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Packaging)))

summary(glm(metadata$resistome_raw_mapped_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$resistome_raw_mapped_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$resistome_raw_mapped_reads ~ as.factor(metadata$Packaging)))

## Resistome sequencing results
summary(glm(microbiome_metadata$Raw_paired_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$Raw_paired_reads ~ as.factor(microbiome_metadata$Packaging)))
summary(glm(microbiome_metadata$Raw_paired_reads ~ as.factor(microbiome_metadata$Dilution)))


summary(glm(microbiome_metadata$QC_filtered_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$QC_filtered_reads ~ as.factor(microbiome_metadata$Packaging)))
summary(glm(microbiome_metadata$QC_filtered_reads ~ as.factor(microbiome_metadata$Dilution)))

summary(glm(microbiome_metadata$nonhost_filtered_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$nonhost_filtered_reads ~ as.factor(microbiome_metadata$Dilution)))
summary(glm(microbiome_metadata$nonhost_filtered_reads ~ as.factor(microbiome_metadata$Packaging)))

summary(glm(microbiome_metadata$microbiome_raw_mapped_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$microbiome_raw_mapped_reads ~ as.factor(microbiome_metadata$Packaging)))


#
##
### Wilcox testing of diversity
##
#

#### AMR diversity comparisons

## Test by treatment, Class and Mech
wilcox.test(metadata$AMR_class_Richness ~ metadata$Treatment)
wilcox.test(metadata$AMR_class_Shannon ~ metadata$Treatment)
wilcox.test(metadata$AMR_mech_Richness ~ metadata$Treatment)
wilcox.test(metadata$AMR_mech_Shannon ~ metadata$Treatment)

## Test by Dilution, Class and Mech
wilcox.test(metadata$AMR_class_Richness ~ metadata$Dilution)
wilcox.test(metadata$AMR_class_Shannon ~ metadata$Dilution)
wilcox.test(metadata$AMR_mech_Richness ~ metadata$Dilution)
wilcox.test(metadata$AMR_mech_Shannon ~ metadata$Dilution)



#### Microbiome diversity comparisons
###
##
#
## Test by treatment, Class and Mech
wilcox.test(microbiome_metadata$microbiome_phylum_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_phylum_Shannon ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_class_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_class_Shannon ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_order_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_order_Shannon ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_family_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_family_Shannon ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_genus_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_genus_Shannon ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_species_Richness ~ microbiome_metadata$Treatment)
wilcox.test(microbiome_metadata$microbiome_species_Shannon ~ microbiome_metadata$Treatment)

