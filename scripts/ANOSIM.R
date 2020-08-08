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


dist_AMR_mech <- vegdist(decostand(t(MRcounts(AMR_analytic_data[[2]], norm=FALSE)), "hell"), "euclidean")
metaMDS_AMR_mech <- metaMDS(dist_AMR_mech, distance="none",symmetric=TRUE, trymax=100)


# Anosim testing
anosim(metaMDS_AMR_class$points,pData(AMR_analytic_data[[1]])$Treatment, distance='euclidean')
anosim(metaMDS_AMR_mech$points,pData(AMR_analytic_data[[2]])$Treatment, distance='euclidean')



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

dist_microbiome_class <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[3]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_class <- metaMDS(dist_microbiome_class, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_order <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[4]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_order <- metaMDS(dist_microbiome_order, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_family <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[5]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_family <- metaMDS(dist_microbiome_family, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_genus <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[6]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_genus <- metaMDS(dist_microbiome_genus, distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_species <- vegdist(decostand(t(MRcounts(microbiome_analytic_data[[7]], norm=FALSE)), "hell"), "euclidean")
metaMDS_microbiome_species <- metaMDS(dist_microbiome_species, distance="none",symmetric=TRUE, trymax=100)

# Anosim testing
anosim(metaMDS_microbiome_phylum$points,pData(microbiome_analytic_data[[2]])$Treatment, distance='euclidean')
anosim(metaMDS_microbiome_class$points,pData(microbiome_analytic_data[[3]])$Treatment, distance='euclidean')
anosim(metaMDS_microbiome_order$points,pData(microbiome_analytic_data[[4]])$Treatment, distance='euclidean')
anosim(metaMDS_microbiome_family$points,pData(microbiome_analytic_data[[5]])$Treatment, distance='euclidean')
anosim(metaMDS_microbiome_genus$points,pData(microbiome_analytic_data[[6]])$Treatment, distance='euclidean')
anosim(metaMDS_microbiome_species$points,pData(microbiome_analytic_data[[7]])$Treatment, distance='euclidean')
