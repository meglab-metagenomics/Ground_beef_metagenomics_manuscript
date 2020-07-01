#FC_meat_anosim code
library(vegan)


#
##
###
####
##### AMR ####
####
###
## AMR by treatment group
#

dist_AMR_class <- vegdist(decostand(t(MRcounts(AMR_analytic_data[[1]], norm=FALSE)), "hell"), "euclidean")
metaMDS_AMR_class <- metaMDS(dist_AMR_class, distance="none",symmetric=TRUE, trymax=100)
#fig <- plot(metaMDS_AMR_class, type="none", display=c("sites"))
#points(metaMDS_AMR_class, dsplay="sites", pch=20, col=as.numeric(pData(AMR_raw_analytic_data[[1]])$Treatment))
#groupz <- sort(unique(pData(AMR_raw_analytic_data[[1]])$Treatment))
#for(i in seq(groupz)) {ordispider(metaMDS_AMR_class, pData(AMR_raw_analytic_data[[1]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}
#identify(fig,"sites")

dist_AMR_mech <- vegdist(decostand(t(MRcounts(AMR_analytic_data[[2]], norm=FALSE)), "hell"), "euclidean")
metaMDS_AMR_mech <- metaMDS(dist_AMR_mech, distance="none",symmetric=TRUE, trymax=100)
#plot(metaMDS_AMR_mech, type="none", display=c("sites"))
#points(metaMDS_AMR_mech, display="sites", pch=20, col=as.numeric(pData(AMR_raw_analytic_data[[1]])$Treatment))
#groupz <- sort(unique(pData(AMR_raw_analytic_data[[2]])$Treatment))
#for(i in seq(groupz)) {ordispider(metaMDS_AMR_mech, pData(AMR_analytic_data[[2]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}


# Anosim testing
#anosim(metaMDS_AMR_class$points,pData(AMR_analytic_data[[1]])$Treatment, distance='euclidean')
#anosim(metaMDS_AMR_mech$points,pData(AMR_analytic_data[[2]])$Treatment, distance='euclidean')



#
##
###
####
##### microbiome ####
####
###
##
#

dist_microbiome_phylum <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[2]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_phylum <- metaMDS(dist_microbiome_phylum, distance="none",symmetric=TRUE, trymax=100)
#plot(metaMDS_microbiome_phylum, type="none", display=c("sites"))
#points(metaMDS_microbiome_phylum, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[2]])$Treatment))
#groupz <- sort(unique(pData(microbiome_analytic_data[[2]])$Treatment))
#for(i in seq(groupz)) {ordispider(metaMDS_microbiome_phylum, pData(microbiome_analytic_data[[2]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

dist_microbiome_class <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[3]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_class <- metaMDS(dist_microbiome_class, distance="none",symmetric=TRUE, trymax=100)
#plot(metaMDS_microbiome_class, type="none", display=c("sites"))
#points(metaMDS_microbiome_class, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[2]])$Treatment))
#groupz <- sort(unique(pData(microbiome_analytic_data[[3]])$Treatment))
#for(i in seq(groupz)) {ordispider(metaMDS_microbiome_class, pData(microbiome_analytic_data[[3]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

dist_microbiome_order <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[4]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_order <- metaMDS(dist_microbiome_order, distance="none",symmetric=TRUE, trymax=100)
#plot(metaMDS_microbiome_order, type="none", display=c("sites"))
#points(metaMDS_microbiome_order, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[2]])$Treatment))
#groupz <- sort(unique(pData(microbiome_analytic_data[[2]])$Treatment))
#for(i in seq(groupz)) {ordispider(metaMDS_microbiome_order, pData(microbiome_analytic_data[[2]])$Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}

dist_microbiome_family <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[5]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_family <- metaMDS(dist_microbiome_family, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_genus <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[6]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_genus <- metaMDS(dist_microbiome_genus, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_species <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[7]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_species <- metaMDS(dist_microbiome_species, distance="none",symmetric=TRUE, trymax=100)

# Anosim testing
#anosim(metaMDS_microbiome_phylum$points,pData(microbiome_analytic_data[[2]])$Treatment, distance='euclidean')
#anosim(metaMDS_microbiome_class$points,pData(microbiome_analytic_data[[3]])$Treatment, distance='euclidean')
#anosim(metaMDS_microbiome_order$points,pData(microbiome_analytic_data[[4]])$Treatment, distance='euclidean')
#anosim(metaMDS_microbiome_family$points,pData(microbiome_analytic_data[[5]])$Treatment, distance='euclidean')
#anosim(metaMDS_microbiome_genus$points,pData(microbiome_analytic_data[[6]])$Treatment, distance='euclidean')
#anosim(metaMDS_microbiome_species$points,pData(microbiome_analytic_data[[7]])$Treatment, distance='euclidean')
