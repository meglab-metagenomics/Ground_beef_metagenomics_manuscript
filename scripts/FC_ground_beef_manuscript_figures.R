### final figures
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
library(ggpubr)

#
##
###
#### Figure 1 - NMDS comparing AMR samples by unique sample (dilution and normal)
###
##
#



#
##
###
#### Figure 2 showing resistome composition and microbiome composition
###
##
#

AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID,sample, Name, Packaging, Treatment)][order(-Packaging )]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=c("Sulfonamides","Rifampin","Trimethoprim","Triclosan","Fluoroquinolones","Aminocoumarins","Nitrofuran","Fosfomycin" , "Elfamycins" ,"Phenicol","Bacitracin","Cationic antimicrobial peptides","Aminoglycosides", "MLS" ,"betalactams" , "Multi-drug resistance" , "Tetracyclines"))

AMR_class_sum$sample = factor(AMR_class_sum$sample ,levels=c("Sample 1.1","Sample 1.2","Sample 2.1","Sample 2.2","Sample 3.1","Sample 3.2","Sample 4.1",
                                                             "Sample 4.2","Sample 5.1","Sample 5.2","Sample 6.1","Sample 6.2","Sample 7.1","Sample 7.2","Sample 8.1","Sample 8.2",
                                                             "Sample 9.1","Sample 9.2","Sample 10.1","Sample 10.2","Sample 11.1","Sample 11.2","Sample 12.1","Sample 12.2","Sample 13.1","Sample 13.2","Sample 14.1",
                                                             "Sample 14.2","Sample 15.1","Sample 15.2","Sample 16.1","Sample 16.2"))

AMR_class_sum$Treatment = factor(AMR_class_sum$Treatment ,levels=c("RWA","CONV"))

AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
fig1_A <- ggplot(AMR_class_sum, aes(x = sample, y = percentage, fill = Class)) +
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Treatment, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=24, angle=30, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("Resistome composition by sample") +
  xlab('') +
  ylab('Relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1)
fig1_A

#
## Microbiome
#

microbiome_phylum_sum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID,sample, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_phylum_sum <- microbiome_phylum_sum[percentage < .01, Name := 'Low Abundance Phyla (< 1%)' ]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]

microbiome_phylum_sum$Name = droplevels(microbiome_phylum_sum$Name)
microbiome_phylum_sum$Name = factor(microbiome_phylum_sum$Name ,levels=c("Low Abundance Phyla (< 1%)","Planctomycetes","Gemmatimonadetes","Chloroflexi","Verrucomicrobia", "Acidobacteria","Actinobacteria" ,"Bacteroidetes" , "Proteobacteria" , "Firmicutes"))

microbiome_phylum_sum$Sample = factor(microbiome_phylum_sum$sample ,levels=c("Sample 1","Sample 2","Sample 3","Sample 4","Sample 5","Sample 6","Sample 7","Sample 8","Sample 9","Sample 10","Sample 11",
                                                                             "Sample 12","Sample 13","Sample 14","Sample 15","Sample 16"))

microbiome_phylum_sum$Treatment = factor(microbiome_phylum_sum$Treatment ,levels=c("RWA","CONV"))



microbiome_phylum_sum$Phylum <- microbiome_phylum_sum$Name
#microbiome_phylum_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
fig1_B <- ggplot(microbiome_phylum_sum, aes(x = Sample, y = percentage, fill = Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Treatment, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=24, angle=30, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20", direction = "1") +
  #ggtitle("Microbiome composition in by treatment (only taxa > 1% per sample)") +
  xlab('') +
  ylab('Relative abundance')
fig1_B

## Combine figures
figure <- ggarrange(fig1_A, fig1_B,
                    labels = c("A)", "B)"),
                    ncol = 1, nrow = 2, vjust = 1.4, font.label = list(size = 30, color = "black", face =
                                                                         "bold", family = NULL))


# Output jpeg figure
jpeg("FC_ground_beef_manuscript_figures/Figure1-FC_meat_resistome_and_microbiome_composition_by_Treatment.jpeg", width =1850, height = 1250)
figure
dev.off()


#
##
###
#### Figure 3 - Resistome ordination by treatment
###
##
#

# Use figure "graphs/AMR/Treatment/NMDS_Treatment_Class.png"


#
##
###
#### Figure 4 - Resistome Diversity figures
###
##
#

par(mfrow = c(2, 2))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

boxplot(metadata$AMR_class_Richness ~ metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 3, "AMR Class", line = 1, cex=2)
mtext(side = 2, "Richness", line = 4, cex=2)

boxplot(metadata$AMR_mech_Richness ~ metadata$Treatment,col=c("Red","Grey"),names=FALSE)
#title("AMR mechanism richness by group")
mtext(side = 3, "AMR Mechanism", line = 1, cex=2)

boxplot(metadata$AMR_class_Shannon ~ metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 2, "Shannon\'s Diversity", line = 4, cex=2)


boxplot(metadata$AMR_mech_Shannon ~ metadata$Treatment, col=c("Red","Grey"),names=FALSE)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", c("CONV","RWA"),x.intersp = .3, xpd = TRUE, horiz = TRUE, inset = c(9, 0), bty = "n", pch = 15, fill = c("Red","Grey"), cex = 1.5)

# Export as png 1200 x 800


########################
##############
########
###
# MICROBIOME
###
########
##############
########################


#
##
###
#### Figure 5 - Microbiome ordination by treatment at the phylum level
###
##
#

# Use figure "graphs/Microbiome/Treatment/NMDS_Treatment_Phylum.png"


#
##
###
#### Figure 6 - Resistome Diversity figures
###
##
#

par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

boxplot(microbiome_metadata$microbiome_phylum_Richness ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 3, "Phylum", line = 1, cex=2)
mtext(side = 2, "Richness", line = 4, cex=2)
boxplot(microbiome_metadata$microbiome_class_Richness ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 3, "Class", line = 1, cex=2)
boxplot(microbiome_metadata$microbiome_order_Richness ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 3, "Order", line = 1, cex=2)
boxplot(microbiome_metadata$microbiome_phylum_Shannon ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)
mtext(side = 2, "Shannon's Diversity", line = 4, cex=2)
boxplot(microbiome_metadata$microbiome_class_Shannon ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)
boxplot(microbiome_metadata$microbiome_order_Shannon ~ microbiome_metadata$Treatment, col=c("Red","Grey"),names=FALSE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("CONV","RWA"),x.intersp = .5, xpd = TRUE, horiz = TRUE, inset = c(9, 0), bty = "n", pch = 15, fill = c("Red","Grey"), cex = 1.5)

# Export as png 1200 x 800

#
##
###
#### Figure 7 - Microbiome ordination by store at the phylum level
###
##
#

# Use figure "graphs/Microbiome/Store/NMDS_Blinded_Store_Phylum.png"

###################################################
##                                               ##
##             Supplemental figures              ##
##                                               ##
###################################################

#
##
###
#### Supplemental Figure 1 - Resistome ordination at the phylum level by dilution
###
##
#


#
##
###
#### Supplemental Figure 1 - Heatmap
###
##
#

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


melted_important_AMR <- amr_melted_raw_analytic[ Level_ID =='Group' & Name %in% important_AMR_regex,  ][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Packaging_samples)]#[order(-sum_class )]
melted_important_AMR[,.(sum_important_genes = sum(Normalized_Count))]

ggplot(data = melted_important_AMR , aes(x = Packaging_samples, y = Name)) +
  geom_tile(aes(fill = log(Normalized_Count))) +
  scale_fill_gradient(low = "lightgrey", high = "steelblue") +
  geom_text(aes(label = num_samples), size=5) +
  labs(fill = "Number of samples", y = "Gene name", x = '') + 
  theme(panel.background = element_blank())

#ggsave("~/Dropbox/WRITING/FC_meat_2019/FC_meat_manuscript/FCmeat_figures/Supp_AMR_important_genes.jpeg", width = 30, height = 20, units = "cm")

ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supplemental_Figure_1_AMR_important_genes.png", width = 30, height = 20, units = "cm")


#
##
###
#### Supplemental figure 2 - Microbiome at the Genus level by Packaging
###
##
#

microbiome_genus_sum <- microbiome_melted_analytic[Level_ID=="Genus", .(sum_genus= sum(Normalized_Count)),by=.(ID,sample, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_genus_sum[,total:= sum(sum_genus), by=.(ID)]
microbiome_genus_sum[,percentage:= sum_genus/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_genus_sum <- microbiome_genus_sum[percentage < .01, Name := 'Low Abundance Genus (< 1%)' ]
microbiome_genus_sum[,total:= sum(sum_genus), by=.(ID)]
microbiome_genus_sum[,percentage:= sum_genus/total ,by=.(ID, Name) ]

microbiome_genus_sum$Name = droplevels(microbiome_genus_sum$Name)
microbiome_genus_sum$Sample = factor(microbiome_genus_sum$sample ,levels=c("Sample 1","Sample 2","Sample 3","Sample 4","Sample 5","Sample 6","Sample 7","Sample 8","Sample 9","Sample 10","Sample 11",
                                                                           "Sample 12","Sample 13","Sample 14","Sample 15","Sample 16"))

microbiome_genus_sum$Genus <- microbiome_genus_sum$Name
ggplot(microbiome_genus_sum, aes(x = Sample, y = percentage, fill = Genus)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=24, angle=30, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20", direction = "1") +
  xlab('') +
  ylab('Relative abundance')

ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supplemental_Figure_2_Microbiome_genus_by_packaging.png", width = 60, height = 40, units = "cm")


#
##
###
#### Supplemental Figure 3 - Microbiome at the phylum level by Packaging
###
##
#

ggplot(microbiome_phylum_sum, aes(x = Sample, y = percentage, fill = Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=24, angle=30, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20", direction = "1") +
  xlab('') +
  ylab('Relative abundance')

ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supplemental_Figure_3_Microbiome_phylum_by_packaging.png", width = 60, height = 40, units = "cm")


#
##
###
#### Supplemental Figure 4 - Resistome composition at the Class level by Packaging
###
##
#


ggplot(AMR_class_sum, aes(x = sample, y = percentage, fill = Class)) +
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=24, angle=30, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('') +
  ylab('Relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1)

ggsave("FC_ground_beef_manuscript_figures/Supplemental_figures/Supplemental_Figure_4_Resistome_class_by_packaging.png", width = 60, height = 40, units = "cm")


