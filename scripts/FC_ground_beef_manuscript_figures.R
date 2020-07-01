### NCBA2 final figures
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
#source('scripts/NCBA2_anosim.r')

### Start of code for figures

setkey(amr_melted_raw_analytic,ID)
setkey(amr_melted_analytic,ID)

setkey(microbiome_melted_analytic,ID)
# Set keys for both metadata files
setkey(metadata,ID)
setkey(microbiome_metadata,ID)
microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]
## Write out counts for figures with tableou



## Figure 1 showing resistome composition and microbiome composition
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=c("Sulfonamides","Rifampin","Trimethoprim","Triclosan","Fluoroquinolones","Aminocoumarins","Nitrofuran","Fosfomycin" , "Elfamycins" ,"Phenicol","Bacitracin","Cationic antimicrobial peptides","Aminoglycosides", "MLS" ,"betalactams" , "Multi-drug resistance" , "Tetracyclines"))

AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
fig1_A <- ggplot(AMR_class_sum, aes(x = ID, y = percentage, fill = Class)) +
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Treatment, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
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



microbiome_phylum_sum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_phylum_sum <- microbiome_phylum_sum[percentage > .01]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]

microbiome_phylum_sum$Name = droplevels(microbiome_phylum_sum$Name)
microbiome_phylum_sum$Name = factor(microbiome_phylum_sum$Name ,levels=c("Planctomycetes","Gemmatimonadetes","Chloroflexi","Verrucomicrobia", "Acidobacteria","Actinobacteria" ,"Bacteroidetes" , "Proteobacteria" , "Firmicutes"))


microbiome_phylum_sum$Phylum <- microbiome_phylum_sum$Name
#microbiome_phylum_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
fig1_B <- ggplot(microbiome_phylum_sum, aes(x = ID, y = percentage, fill = Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20") +
  #ggtitle("Microbiome composition in by treatment (only taxa > 1% per sample)") +
  xlab('Sample ID') +
  ylab('Relative abundance')
fig1_B

jpeg("FC_ground_beef_manuscript_figures/Figure1-FC_meat_resistome_and_microbiome_composition_by_Treatment.jpeg", width =1850, height = 1250)
fig1_plot <- grid.arrange(fig1_A,fig1_B, nrow= 2, widths = 2)
dev.off()

ml <- marrangeGrob(fig1_plot)

#ggsave("FC_ground_beef_manuscript_figures/Figure1-FC_meat_resistome_and_microbiome_composition_by_Treatment.jpeg", width = 30, height = 20, units = "cm")


AMR_mech_sum <- amr_melted_analytic[Level_ID=="Mechanism", .(sum_mech= sum(Normalized_Count),Store),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
AMR_mech_sum[,total:= sum(sum_mech), by=.(ID)]
AMR_mech_sum[,percentage:= sum_mech/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
AMR_mech_sum <- AMR_mech_sum[percentage > .05]
AMR_mech_sum[,total:= sum(sum_mech), by=.(ID)]
AMR_mech_sum[,percentage:= sum_mech/total ,by=.(ID, Name) ]

AMR_mech_sum$Name = droplevels(AMR_mech_sum$Name)
#AMR_mech_sum$Name = factor(AMR_mech_sum$Name ,levels=c("Sulfonamides","Rifampin","Trimethoprim","Triclosan","Fluoroquinolones","Aminocoumarins","Nitrofuran","Fosfomycin" , "Elfamycins" ,"Phenicol","Bacitracin","Cationic antimicrobial peptides","Aminoglycosides", "MLS" ,"betalactams" , "Multi-drug resistance" , "Tetracyclines"))

AMR_mech_sum$mech <- AMR_mech_sum$Name
#AMR_mech_sum[,percentage:= round(sum_mech/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_mech_sum, aes(x = ID, y = percentage, fill = mech)) +
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Store, scales='free') +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("Resistome composition by sample") +
  xlab('Sample ID') +
  ylab('Relative abundance') +
  scale_color_tableau("Tableau 20")

#
microbiome_Genus_sum <- microbiome_melted_analytic[Level_ID=="Genus", .(sum_Genus= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_Genus_sum[,total:= sum(sum_Genus), by=.(ID)]
microbiome_Genus_sum[,percentage:= sum_Genus/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_Genus_sum <- microbiome_Genus_sum[percentage > .05]
microbiome_Genus_sum[,total:= sum(sum_Genus), by=.(ID)]
microbiome_Genus_sum[,percentage:= sum_Genus/total ,by=.(ID, Name) ]

microbiome_Genus_sum$Name = droplevels(microbiome_Genus_sum$Name)
microbiome_Genus_sum$Name = factor(microbiome_Genus_sum$Name ,levels=c("Lactococcus","Leuconostoc","Lactobacillus","Photobacterium","Pseudomonas" ,"DA101","Clostridium" ,"Carnobacterium" , "Streptococcus" , "Bifidobacterium","Psychrobacter"))


microbiome_Genus_sum$Genus <- microbiome_Genus_sum$Name
#microbiome_Genus_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
fig1_B <- ggplot(microbiome_Genus_sum, aes(x = ID, y = percentage, fill = Genus)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20") +
  #ggtitle("Microbiome composition in by treatment (only taxa > 1% per sample)") +
  xlab('Sample ID') +
  ylab('Relative abundance')
fig1_B

jpeg("FC_ground_beef_manuscript_figures/Figure-FC_meat_microbiome_Genus_by_Packaging_core.jpeg", width =1850, height = 1250)
fig1_B
dev.off()





##
###
#### Figure 2 - Ordination comparing dilution
###
##
#




### Figure 3 - Resistome ordination by packaging
jpeg("FC_ground_beef_manuscript_figures/Figure3_AMR_Packaging.jpeg")

par(mfrow = c(2,1))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

##
plot(metaMDS_AMR_class, type="none", display=c("sites"))
points(metaMDS_AMR_class, display="sites", pch=20, col=as.numeric(pData(AMR_analytic_data[[1]])$Packaging))
groupz <- sort(unique(pData(AMR_analytic_data[[1]])$Packaging))
for(i in seq(groupz)) {ordispider(metaMDS_AMR_class, pData(AMR_analytic_data[[1]])$Packaging,font=2, cex=1.5, col=i, show.groups=groupz[i])}

#mtext(side = 3, "AMR Class", line = 1, cex=2)
#mtext(side = 3, "Resistome Ordination", line = 1, cex=2)
mtext(side = 2, "AMR Class", line = 4, cex=2)

plot(metaMDS_AMR_mech, type="none", display=c("sites"))
points(metaMDS_AMR_mech, display="sites", pch=20, col=as.numeric(pData(AMR_analytic_data[[2]])$Packaging))
groupz <- sort(unique(pData(AMR_analytic_data[[2]])$Packaging))
for(i in seq(groupz)) {ordispider(metaMDS_AMR_mech, pData(AMR_analytic_data[[2]])$Packaging,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "AMR Mechanism", line = 4, cex=2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("CONV-Chub","CONV-StoreGrind","CONV-TrayOverwrap","RWA-Vacuum"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(19), col = 1:4, cex = 1.5)
#identify(fig,"sites")
dev.off()



### Figure 4 - Resistome Diversity figures
#jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure4_AMR_richness_diversity_boxplots.jpeg")

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

#title("AMR Class diversity by group")
boxplot(metadata$AMR_mech_Shannon ~ metadata$Treatment, col=c("Red","Grey"),names=FALSE)
#title("AMR mechanism diversity by group")


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", c("CONV","RWA"),x.intersp = .3, xpd = TRUE, horiz = TRUE, inset = c(9, 0), bty = "n", pch = 15, fill = c("Red","Grey"), cex = 1.5)

#legend("bottom", c("CONV-Chub","CONV-StoreGrind","CONV-TrayOverwrap","RWA-Vacuum"),x.intersp = .3, xpd = TRUE, horiz = TRUE, inset = c(9, 0), bty = "n", pch = 15, fill = c("Red","White","deepskyblue3","Grey"), cex = 1.5)


dev.off()



########################
##############
########
###
# MICROBIOME
###
########
##############
########################




## Figure 5

microbiome_phylum_sum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_phylum_sum <- microbiome_phylum_sum[percentage > .01]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]

microbiome_phylum_sum$Name = droplevels(microbiome_phylum_sum$Name)
microbiome_phylum_sum$Name = factor(microbiome_phylum_sum$Name ,levels=c("Planctomycetes","Gemmatimonadetes","Chloroflexi","Verrucomicrobia", "Acidobacteria","Actinobacteria" ,"Bacteroidetes" , "Proteobacteria" , "Firmicutes"))


microbiome_phylum_sum$Class <- microbiome_phylum_sum$Name
#microbiome_phylum_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(microbiome_phylum_sum, aes(x = ID, y = percentage, fill = Class)) +
  geom_bar(stat = "identity")+
  facet_wrap( ~ Treatment, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20") +
  ggtitle("Microbiome composition in by treatment (only taxa > 1% per sample)") +
  xlab('Sample ID') +
  ylab('Relative abundance')

ggsave("FC_ground_beef_manuscript_figures/Figure5-microbiome_phylum_barplot.jpeg", width = 30, height = 20, units = "cm")




# Figure 6 microbiome NMDS


### Figure 6 - Microbiome group ordination ####

#pdf("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure6_Microbiome_Group.pdf")
jpeg("FC_ground_beef_manuscript_figures/Figure6_Microbiome_ordination.jpeg")

par(mfrow = c(3,1))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots
par(tcl = -0.5)
par(mgp = c(2, 0.6, 0))


plot(metaMDS_microbiome_phylum, type="none", display=c("sites"))
points(metaMDS_microbiome_phylum, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[2]])$Treatment))
groupz <- sort(unique(pData(microbiome_analytic_data[[2]])$Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_phylum, pData(microbiome_analytic_data[[2]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

#mtext(side = 3, "Microbiome Ordination", line = 1, cex=2)
mtext(side = 2, "Phylum", line = 4, cex=2)

plot(metaMDS_microbiome_class, type="none", display=c("sites"))
points(metaMDS_microbiome_class, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[3]])$Treatment))
groupz <- sort(unique(pData(microbiome_analytic_data[[3]])$Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_class, pData(microbiome_analytic_data[[3]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "Class", line = 4, cex=2)

plot(metaMDS_microbiome_order, type="none", display=c("sites"))
points(metaMDS_microbiome_order, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[4]])$Treatment))
groupz <- sort(unique(pData(microbiome_analytic_data[[4]])$Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_order, pData(microbiome_analytic_data[[4]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "Order", line = 4, cex=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Day 1-Treated","Day 1-Untreated","Day11-Treated","Day11-Untreated"),xpd = TRUE, horiz = TRUE, inset = c(0,
                                                                                                                            0), bty = "n", pch = c(19), col = 1:4, cex = 1.5)
#identify(fig,"sites")
dev.off()


## Figure 7
jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure7_Micro_richness_diversity.jpeg")

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
dev.off()




































































## Heatmap at group level - by sample
AMR_group_samples <- amr_melted_raw_analytic[Level_ID=="Group"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Group)]#[order(-sum_class )]

group_annotations <- annotations[,-1]
group_annotations <- unique(group_annotations,by = c("group"))

setkey(group_annotations, group)
setkey(AMR_group_samples , Name)

AMR_group_samples  <- group_annotations[AMR_group_samples ]
setkey(AMR_group_samples , ID)
AMR_group_samples  <- metadata[AMR_group_samples]

ggplot(data = AMR_group_samples , aes(x = Group, y = group)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5)



#write.csv(AMR_group_samples, "16S_norm_group_melted_sample_counts_.csv")


## Heatmap at Mechanism level - by sample
AMR_mech_samples <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, ID)]#[order(-sum_class )]

mech_annotations <- annotations[,-1]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))

setkey(mech_annotations, mechanism)
setkey(AMR_mech_samples , Name)

AMR_mech_samples  <- mech_annotations[AMR_mech_samples ]
setkey(AMR_mech_samples , ID)
AMR_mech_samples  <- metadata[AMR_mech_samples]
#write.csv(AMR_mech_samples, "16S_norm_melted_sample_counts.csv")



## Heatmap at Group level
AMR_mech_sample_groups <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Group)]#[order(-sum_class )]
AMR_mech_sample_groups$Group <- factor(AMR_mech_sample_groups$Group,levels = c("Day1-Treated", "Day1-Untreated", "Day11-Treated", "Day11-Untreated"))

mech_annotations <- annotations[,-1]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))

setkey(mech_annotations, mechanism)
setkey(AMR_mech_sample_groups, Name)

AMR_mech_sample_groups <- mech_annotations[AMR_mech_sample_groups]

write.csv(AMR_mech_sample_groups, "16S_norm_melted_Group_counts.csv")


ggplot(data = AMR_mech_sample_groups , aes(x = Group, y = mechanism)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5)



setkey(AMR_mech_samples, ID)
AMR_mech_samples <- AMR_mech_samples[metadata]
