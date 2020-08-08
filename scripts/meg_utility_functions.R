###############
## Functions ##
###############
## Utility functions that can be optionally used later for producing
## graphs and performing reshaping operations
## With some edits for a higher trymax and change to heatmap code

# Call variables from parent scope
`..` <- function (..., .env = sys.parent(2)) {
  get(deparse(substitute(...)), env = .env)
}

# Misc reshape function for data table
melt_dt <- function(D, level_id) {
  temp <- reshape2::melt(D, variable.name='Sample', value.name='Normalized_Count')
  names(temp) <- c('Name', 'ID', 'Normalized_Count')
  temp <- data.table(cbind(rep(level_id, nrow(temp)), temp))
  names(temp)[1] <- 'Level_ID'
  return(temp)
}

# Filter data by quantile
meg_filter_data <- function(data_list,
                            filter_min_threshold) {
  local_obj <- data_list
  for( l in 1:length(local_obj) ) {
    filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
    cumNorm(local_obj[[l]])
  }
  return(local_obj)
}




make_sparse <- function(df, rownames, excludes, filter_min_threshold=0.15) {
  local_df <- df[, .SD, .SDcols=!excludes]
  filter_threshold <- quantile(rowSums(local_df), 0.15)
  if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
  chosen <- which(rowSums(local_df) >= filter_threshold )
  ret <- as.data.frame(local_df[chosen, ])
  rownames(ret) <- df[[rownames]][chosen]
  return(ret)
}

data_subset <- function(data_obj, subsets) {
  local_meta <- data.table(pData(data_obj))
  local_subset <- c()
  for( c in 1:length(subsets) ) {
    
    conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
    conditional_string <- paste('local_meta[[\'', conditional_terms[1],
                                '\']] ', conditional_terms[2],
                                ' \'', conditional_terms[3], '\'',
                                sep='', collapse='')
    if(length(local_subset) > 0) {
      local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
    }
    else {
      local_subset <- which(eval(parse(text=conditional_string)))
    }
  }
  return(data_obj[, local_subset])
}


data_subset_long <- function(data_obj, subsets) {
  local_subset <- c()
  for( c in 1:length(subsets) ) {
    
    conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
    conditional_string <- paste('data_obj[[\'', conditional_terms[1],
                                '\']] ', conditional_terms[2],
                                ' \'', conditional_terms[3], '\'',
                                sep='', collapse='')
    if(length(local_subset) > 0) {
      local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
    }
    else {
      local_subset <- which(eval(parse(text=conditional_string)))
    }
  }
  return(data_obj[local_subset, ])
}


##############################
###                        ###
####     Ordination       ####
###                        ###
##############################
# Calculate the points for convex hulls around data for ordination plots
meg_find_hulls <- function(x) x[chull(x$Ord1, x$Ord2),]

hullPlot <- function(df, grouping, legend = TRUE){
  # requires that your columns are called NMDS1, NMDS2, NMDS3
  
  # requires that the "grouping" column is printed as per the name of the column
  # you want separate hulls for
  
  
  allSites <- sort(as.vector(unique(df[[grouping]])))
  matList <- list()
  hullList <- list()
  cols <- as.character(brewer_pal(type = "qual", palette = 1, direction = 1)(length(allSites)))
  
  
  # this loop creates the points for each site
  # it also calculates a separate hull for each site
  # "site" or whatever grouping variable you are using
  for(thisSite in 1:length(allSites)){
    
    tmp <- df[df[grouping] == allSites[thisSite], ]
    
    plot3d(tmp$NMDS1, tmp$NMDS2, tmp$NMDS3, col = cols[thisSite], box = FALSE,
           type = "s", radius = 0.01, add = ifelse(thisSite > 1, TRUE, FALSE),
           xlab = "", ylab = "", zlab = "")
    
    
    matList[[thisSite]] <- matrix(
      c(tmp[[grep("NMDS1", names(tmp))]],
        tmp[[grep("NMDS2", names(tmp))]],
        tmp[[grep("NMDS3", names(tmp))]]), ncol = 3)
    
    hullList[[thisSite]] <- t(convhulln(matList[[thisSite]]))
    
  }
  
  
  # this will run if you have legend = TRUE (the default)
  if(legend){
    # this slows down the plotting which is necessary otherwise
    # the printing can lag and the legend goes to a funny size
    Sys.sleep(0.5)
    
    # you can change your legend as per 'legend' commands
    legend3d("bottom", legend = allSites,
             # uses the same cols as for plotting
             col = cols,
             # symbol size
             pch = 16,
             inset=c(0.02),
             horiz = TRUE)
    
  }
  
  
  # this loop plots the hulls
  for(hull in seq_along(matList)){
    rgl.triangles(matList[[hull]][hullList[[hull]],1],matList[[hull]][hullList[[hull]],2],matList[[hull]][hullList[[hull]],3],
                  col=cols[hull],
                  
                  # change the alpha to change how see through they are
                  alpha=.6)
    
  }
  
}
# Function for computing ordination plots with convex hulls
# using ggplot2.  You will have to specify the facet variable,
# i.e. which experimental design variable determines the grouping of points.
#
# Function arguments:
#   data_list: a list containing the MRexperiment objects for all levels
#   data_names: a character vector of level names for each MRexperiment in data_list
#   metadata: a data.table of metadata for each sample
#   sample_var: the column name in metadata that specifies the sample IDs (must match MRexperiment columns)
#   hull_var: the metadata column name by which to group data points
#   outdir: the file path of the directory for output of files
#   data_type: the data type being computed on, e.g. AMR or Microbiome
#   method: choice of 'NMDS' or 'PCA' for method of ordination
#
#   outputs:
meg_ordination <- function(data_list,
                           data_names,
                           metadata,
                           sample_var,
                           hull_var,
                           analysis_subset,
                           outdir,
                           data_type,
                           method='NMDS',
                           factor_order) {
  all_ord <- data.table(ID=character(),
                        Level_ID=character(),
                        Ord1=numeric(),
                        Ord2=numeric(),
                        Group_Var=character())
  all_hulls <- data.table(Group_Var=character(),
                          Ord1=numeric(),
                          Ord2=numeric(),
                          Level_ID=character())
  setkey(all_ord, ID)
  
  local_obj <- data_list
  for( l in 1:length(local_obj) ) {
    
    if(length(analysis_subset) > 0) {
      local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
    }
    local_meta <- metadata
    local_meta <- local_meta[local_meta[[sample_var]] %in% colnames(MRcounts(local_obj[[l]])), ]
    
    # Transpose the matrix for NMDS (groups are now in rows and features in columns)
    t_data <- t(MRcounts(local_obj[[l]]))
    local_meta <- local_meta[rowSums(t_data) > 0, ]
    t_data <- t_data[rowSums(t_data) > 0, ]
    t_data <- t_data[which(!is.na(local_meta[[sample_var]]) & local_meta[[sample_var]] != 'NA'), colSums(t_data) > 0]
    
    if( method == 'NMDS' ) {
      # Set parallel to whatever your computer can support in terms of CPU count
      #ord.res <- metaMDS(t_data, autotransform=F, parallel=7, trymax=1000)  #### adding way more try max, from 49 to 1000
      ord.res <- metaMDS(t_data, autotransform=F, trymax=100, stratmax = 0.999999)  #### adding way more try max, from 49 to 1000, removed parallel=7,
      ord_points <- data.table(ord.res$points)
      names(ord_points) <- c('Ord1', 'Ord2')
      ord_points[, ID :=( rownames(ord.res$points) )]
    }
    else if( method == 'PCA' ) {
      ord.res <- prcomp(t_data, center=T, scale=T)
      
      # Format to include metadata for ggplot2
      ord_points <- data.table(ord.res$x[, 1:2])
      names(ord_points) <- c('Ord1', 'Ord2')
      ord_points[, ID :=( rownames(ord.res$x) )]
    }
    else {
      stop('method must be either NMDS or PCA')
    }
    
    ord_points[, Level_ID :=( rep(data_names[l], nrow(ord_points)) )]
    setkey(ord_points, ID)
    ord_points <- metadata[ord_points]
    ord_points <- ord_points[, .SD, .SDcols=c(sample_var, hull_var, 'Ord1', 'Ord2', 'Level_ID')]
    ## Check if needing to change factor order. Need to convert to data frame first, then back to data table
    ord_points <- as.data.frame(ord_points)
    ifelse(factor_order != '', ord_points[[hull_var]] <- factor(ord_points[[hull_var]], levels = eval(factor_order)),1)
    ord_points <- as.data.table(ord_points)
    names(ord_points)[2] <- 'Group_Var'
    all_ord <- rbind(all_ord, ord_points)
    
    hulls <- ord_points[, meg_find_hulls(.SD), .SDcols=c('Ord1', 'Ord2'), by=Group_Var]
    hulls[, Level_ID :=( rep(data_names[l], nrow(hulls)) )]
    hulls[, .SD, .SDcols=!'Level_ID']
    all_hulls <- rbind(all_hulls, hulls)
    
    # Plot graphs with convex hulls
    g_ord <- ggplot(data=ord_points, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
      geom_point(size=2.5) + geom_polygon(data=hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                                          alpha=0.2, show.legend=F)
    #g_ord <- g_ord + scale_fill_tableau("Classic 20", direction = -1)
    #g_ord <- g_ord + scale_color_tableau("Classic 20", direction = -1)
    g_ord <- g_ord +
      ggtitle(paste(method, ' for ', data_type, ' ', 'by ', hull_var, '\nAnnotation Level: ',
                    data_names[l],
                    sep='',
                    collapse='')) +
      labs(color=hull_var) +
      guides(fill=F) +
      theme(strip.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_text(size=24),
            axis.title.y=element_text(size=24, hjust=0.5),
            #legend.position="right",
            legend.title=element_text(size=20),
            legend.text=element_text(size=18),
            plot.title=element_text(size=24, hjust=0.5),
            panel.background = element_rect(fill = "white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.spacing=unit(0.1, "lines"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    if( method == 'NMDS' ) {
      g_ord <- g_ord + xlab('MDS1') + ylab('MDS2')
    }
    else if( method == 'PCA' ) {
      g_ord <- g_ord + xlab('PC1') + ylab('PC2')
    }
    
    # Open the graphics device at the specified location and figure size
    png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                       data_names[l], '.png', sep='', collapse=''),
        width=1024, height=768)
    
    print(g_ord)
    
    # Turn off graphics device to save the graphic
    dev.off()
    
  }
  
  
  all_ord <- within(all_ord, Level_ID
                    <- factor(Level_ID, levels=data_names,
                              ordered=T))
  all_hulls <- within(all_hulls, Level_ID
                      <- factor(Level_ID, levels=data_names,
                                ordered=T))
  
  g_all_ord <- ggplot(data=all_ord, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
    geom_point(size=3) +
    geom_polygon(data=all_hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                 alpha=0.2, show.legend=F) +
    facet_wrap(~Level_ID, nrow=2)
  g_all_ord <- g_all_ord +
    ggtitle(paste(method, ' for ', data_type, ' by ', hull_var, sep='', collapse='')) +
    labs(color=hull_var) +
    guides(fill=F) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, hjust=0.5),
          #legend.position="right",
          legend.title=element_text(size=24, hjust=0.5),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0, "mm"),
          #strip.background = element_rect(size = 2, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  #g_all_ord <- g_all_ord + scale_color_tableau("Classic 20", direction = -1)
  #g_all_ord <- g_all_ord + scale_fill_tableau("Classic 20", direction = -1)
  
  png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                     'AllLevels.png', sep='', collapse=''),
      width=1024, height=768)
  print(g_all_ord)
  dev.off()
  return(all_ord)
}

##############################
###                        ###
####      Heatmap         ####
###                        ###
##############################
heatmap_select_top_counts <- function(X, group_var, sample_var, n) {
  return(X[, tail(.SD, n), by=c(group_var, sample_var)])
}

heatmap_col <- c("#849bed", "#9ce2de","#f7f6ed","#efdfa7","#f77a74")

#heatmap_col <- c("grey","#f7f6ed","#efdfa7","#f77a74")

meg_heatmap <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type,
                        factor_order='') {
  tile_subset <- melted_data[Level_ID == level_var, ]
  colnames(tile_subset)[colnames(tile_subset) == 'ID'] <- sample_var
  setkeyv(tile_subset, sample_var)
  setkeyv(metadata, sample_var)
  
  tile_subset <- metadata[tile_subset]
  tile_subset <- tile_subset[!is.na(tile_subset[[group_var]]), ]
  
  if(length(analysis_subset) > 0) {
    tile_subset <- data_subset_long(tile_subset, analysis_subset)
  }
  
  
  sample_order <- unique(tile_subset[order(group_var), sample_var,with=FALSE]) #Error in unique(tile_subset[order(group_var), sample_var]) : error in evaluating the argument 'x' in selecting a method for function 'unique': Error in `[.data.table`(tile_subset, order(group_var), sample_var) : j (the 2nd argument inside [...]) is a single symbol but column name 'sample_var' is not found. Perhaps you intended DT[,..sample_var] or DT[,sample_var,with=FALSE]. This difference to data.frame is deliberate and explained in FAQ 1.1.
  tile_subset <- within(tile_subset, sample_var
                        <- factor(sample_var,
                                  levels=sample_order,
                                  ordered=T))
  
  setkey(tile_subset, Normalized_Count)
  tile_subset <- tile_subset[, sum(Normalized_Count),
                             by=c(group_var, sample_var, 'Name')]
  names(tile_subset)[length(names(tile_subset))] <- 'Normalized_Count'
  
  numselect <- 20
  tile_names <- heatmap_select_top_counts(tile_subset, group_var,
                                          sample_var, numselect)
  name_count <- length(unique(tile_names$Name))
  
  while(name_count > 25) {
    numselect <- numselect - 1
    tile_names <- heatmap_select_top_counts(tile_subset, group_var,
                                            sample_var, numselect)
    name_count <- length(unique(tile_names$Name))
  }
  
  
  tile_subset <- tile_subset[Name %in% tile_names$Name, ]
  
  ## Check if needing to change factor order. Need to convert to data frame first, then back to data table
  tile_subset <- as.data.frame(tile_subset)
  ifelse(factor_order != '', tile_subset[[group_var]] <- factor(tile_subset[[group_var]], levels = eval(factor_order)),1)
  tile_subset <- as.data.table(tile_subset)
  
  # Plot object
  tile <- ggplot(tile_subset, aes_string(x=sample_var, y='Name')) +
    geom_tile(aes(fill=log2(Normalized_Count+1))) +
    facet_wrap(as.formula(paste('~', group_var)), strip.position ='bottom', scales = 'free_x', nrow = 1) +
    theme(strip.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=20),
          axis.title.y=element_blank(),
          legend.position="bottom",
          legend.title=element_text(size=20),
          plot.title=element_text(size=24, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(paste('Samples by ', group_var, sep='', collapse='')) +
    scale_fill_gradientn(colors=heatmap_col) +
    labs(fill= 'Log2 Normalized Count') +
    ggtitle(paste(data_type, ' top ', name_count, " ", level_var, ' Normalized Counts by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/','Heatmap_', level_var, '_by', group_var,
                     '.png', sep='', collapse=''), width=1400, height=700)
  print(tile)
  dev.off()
}

##############################
###                        ###
####   Alpha diversity    ####
###                        ###
##############################

# Function that returns species count, rarefied species count, and alpha diversity measures
# for each sample in the m x n matrix, m = features, n = samples
alpha_CSS_diversity <- function(X, minlevel, method='invsimpson') {
  S <- specnumber(X, MARGIN=2)
  alphadiv <- vegan::diversity(X, index=method, MARGIN=2)
  return(list(CSS_species_abundance=S,
              CSS_data=X,
              alphadiv=alphadiv))
}


meg_alpha_diversity <- function(data_list,
                                data_names,
                                metadata,
                                sample_var,
                                group_var,
                                analysis_subset,
                                outdir,
                                data_type,
                                factor_order) {
  #Create tables that hold the data
  print(group_var)
  all_alphadiv <- data.table(ID=character(),
                             Level=character(),
                             Value=numeric())
  names(all_alphadiv)[1] <- sample_var
  
  # CSS_species_abundance
  all_species_CSS <- data.table(ID=character(),
                                Level=character(),
                                Value=numeric())
  names(all_species_CSS)[1] <- sample_var
  
  local_data <- data_list
  
  
  all_sample_data <- alpha_CSS_diversity(MRcounts(local_data[[1]]))
  
  
  # I don't currently use the non_zero_sample variable for CSS counts, but it could come in handy for removing those samples or including a report for which samples had 0 counts
  for( l in 1:length(local_data) ) {
    
    if(length(analysis_subset) > 0) {
      local_data[[l]] <- data_subset(local_data[[l]], analysis_subset)
    }
    
    sample_counts <- colSums(MRcounts(local_data[[l]]))
    if(min(sample_counts) == 0) {
      non_zero_sample <- min(sample_counts[sample_counts > 0])
    }
    else {
      non_zero_sample <- 0
    }
    
    # Create main object with alpha rarefaction values.
    # This object includes CSS_species_abundance, CSS_data, and alphadiv
    local_obj <- alpha_CSS_diversity(MRcounts(local_data[[l]]))
    
    # Alpha rarefaction
    temp <- data.table(ID=names(local_obj$alphadiv),
                       Level=rep(data_names[l],
                                 length(local_obj$alphadiv)),
                       Value=as.numeric(local_obj$alphadiv))
    names(temp)[1] <- sample_var
    all_alphadiv <- rbind(all_alphadiv, temp, fill=TRUE)
    
    # Raw species
    temp <- data.table(ID=names(local_obj$CSS_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$CSS_species_abundance)),
                       Value=as.numeric(local_obj$CSS_species_abundance))
    names(temp)[1] <- sample_var
    all_species_CSS <- rbind(all_species_CSS, temp, fill = TRUE)
    
  }
  
  #
  ## Alpha rarefaction
  #
  all_alphadiv <- within(all_alphadiv, Level
                         <- factor(Level, levels=data_names,
                                   ordered=T))
  
  setkeyv(all_alphadiv, sample_var)
  setkeyv(metadata, sample_var)
  all_alphadiv <- metadata[all_alphadiv]
  all_alphadiv <- all_alphadiv[!is.na(all_alphadiv[[group_var]]), ]
  
  alphadiv_type_sums <- all_alphadiv[Level==data_names[2], median(round(Value, digits=0)), by=group_var]
  alphadiv_value_labels <- as.character(alphadiv_type_sums[[group_var]][order(alphadiv_type_sums$V1, decreasing=T)])
  all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]],
                                      levels=alphadiv_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  #print(all_alphadiv)
  #all_alphadiv[[group_var]] <- droplevels( all_alphadiv[[group_var]])
  #all_alphadiv$Level <- droplevels( all_alphadiv$Level)
  # Plot figure
  png(filename=paste(outdir, '/', 'Diversity_', data_type, '_CSS_InvSimpson_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_alphadiv <- ggplot(data=all_alphadiv, aes_string(x=group_var,
                                                     y='Value',
                                                     fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_alphadiv <- g_alphadiv + scale_fill_tableau("Classic 20", direction = -1)
  g_alphadiv <- g_alphadiv +
    ggtitle(paste('Alpha Diversity by ', group_var, ' for CSS data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Inverse Simpson\'s Index\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_alphadiv)
  dev.off()
  
  # Output Diversity counts
  write.table(all_alphadiv, paste(outdir, '/', data_type, '_CSS_InvSimpson_values.csv',
                                  sep='', collapse=''), sep=",")
  
  #
  ## Raw species abundance
  #
  all_species_CSS <- within(all_species_CSS, Level
                            <- factor(Level, levels=data_names,
                                      ordered=T))
  setkeyv(all_species_CSS, sample_var)
  all_species_CSS <- metadata[all_species_CSS]
  all_species_CSS <- all_species_CSS[!is.na(all_species_CSS[[group_var]]), ]
  species_raw_type_sums <- all_species_CSS[Level==data_names[2], median(Value), by=group_var]
  species_raw_value_labels <- as.character(species_raw_type_sums[[group_var]][order(species_raw_type_sums$V1,
                                                                                    decreasing=T)])
  all_species_CSS[[group_var]] <- factor(all_species_CSS[[group_var]],
                                         levels=species_raw_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_species_CSS[[group_var]] <- factor(all_species_CSS[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  
  # Plot figure
  png(filename=paste(outdir, '/', 'Diversity_', data_type, '_CSS_Richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_sraw <- ggplot(data=all_species_CSS, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_sraw <- g_sraw + scale_fill_tableau("Classic 20", direction = -1)
  g_sraw <- g_sraw +
    ggtitle(paste('Species Richness by ', group_var, ' for CSS data\nRichness',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_sraw)
  dev.off()
  
  # Output Diversity counts
  write.table(all_species_CSS, paste(outdir, '/', data_type, '_CSS_richness_values.csv',
                                     sep='', collapse=''), sep=",")
  
  
}



##############################
###                        ###
####   Alpha rarefaction  ####
###                        ###
##############################




alpha_rarefaction <- function(X, minlevel, method='invsimpson') {
  S <- specnumber(X, MARGIN=2)
  raremax <- min(colSums(X))
  if( raremax < minlevel ) raremax <- minlevel
  Srare <- rarefy(X, raremax, MARGIN=2)
  Xrare <- t(rrarefy(t(X), raremax))
  alphadiv <- vegan::diversity(Xrare, index=method, MARGIN=2)
  return(list(raw_species_abundance=S,
              rarefied_species_abundance=Srare,
              rarefied_data=Xrare,
              alphadiv=alphadiv))
}

meg_alpha_rarefaction <- function(data_list,
                                  data_names,
                                  metadata,
                                  sample_var,
                                  group_var,
                                  analysis_subset,
                                  outdir,
                                  data_type,
                                  factor_order) {
  all_alphadiv <- data.table(ID=character(),
                             Level=character(),
                             Value=numeric())
  names(all_alphadiv)[1] <- sample_var
  setkeyv(metadata, sample_var)
  all_species_raw <- data.table(ID=character(),
                                Level=character(),
                                Value=numeric())
  names(all_species_raw)[1] <- sample_var
  
  all_species_rare <- data.table(ID=character(),
                                 Level=character(),
                                 Value=numeric())
  names(all_species_rare)[1] <- sample_var
  
  local_data <- data_list
  
  for( l in 1:length(local_data) ) {
    
    if(length(analysis_subset) > 0) {
      local_data[[l]] <- data_subset(local_data[[l]], analysis_subset)
    }
    
    sample_counts <- colSums(MRcounts(local_data[[l]]))
    if(min(sample_counts) == 0) {
      non_zero_sample <- min(sample_counts[sample_counts > 0])
    }
    else {
      non_zero_sample <- 0
    }
    
    # Create main object with alpha rarefaction values.
    # This object includes raw_species_abundance, rarefied_species_abundance, rarefied_data, and alphadiv
    local_obj <- alpha_rarefaction(MRcounts(local_data[[l]]), minlevel = non_zero_sample)
    
    # Alpha rarefaction
    temp <- data.table(ID=names(local_obj$alphadiv),
                       Level=rep(data_names[l],
                                 length(local_obj$alphadiv)),
                       Value=as.numeric(local_obj$alphadiv))
    names(temp)[1] <- sample_var
    all_alphadiv <- rbind(all_alphadiv, temp)
    
    # Raw species
    temp <- data.table(ID=names(local_obj$raw_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$raw_species_abundance)),
                       Value=as.numeric(local_obj$raw_species_abundance))
    names(temp)[1] <- sample_var
    all_species_raw <- rbind(all_species_raw, temp)
    
    # Rarefied species abundance
    temp <- data.table(ID=names(local_obj$rarefied_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$rarefied_species_abundance)),
                       Value=as.numeric(local_obj$rarefied_species_abundance))
    names(temp)[1] <- sample_var
    all_species_rare <- rbind(all_species_rare, temp)
  }
  
  #
  ## Alpha rarefaction
  #
  all_alphadiv <- within(all_alphadiv, Level
                         <- factor(Level, levels=data_names,
                                   ordered=T))
  setkeyv(all_alphadiv, sample_var)
  all_alphadiv <- metadata[all_alphadiv]
  all_alphadiv <- all_alphadiv[!is.na(all_alphadiv[[group_var]]), ]
  
  alphadiv_type_sums <- all_alphadiv[Level==data_names[2], median(Value), by=group_var]
  alphadiv_value_labels <- as.character(alphadiv_type_sums[[group_var]][order(alphadiv_type_sums$V1, decreasing=T)])
  all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]],
                                      levels=alphadiv_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  
  # Plot figure
  png(filename=paste(outdir, '/', 'Diversity_', data_type, '_rarefied_InvSimpson_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_alphadiv <- ggplot(data=all_alphadiv, aes_string(x=group_var,
                                                     y='Value',
                                                     fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_alphadiv <- g_alphadiv + scale_fill_tableau("Classic 20", direction = -1)
  g_alphadiv <- g_alphadiv +
    ggtitle(paste('Alpha Diversity by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Inverse Simpson\'s Index\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_alphadiv)
  dev.off()
  
  
  # Output Diversity counts
  write.table(all_alphadiv, paste(outdir, '/', data_type, '_rarefied_InvSimpson_values.csv',
                                  sep='', collapse=''), sep=",")
  
  
  #
  ## Raw species abundance
  #
  all_species_raw <- within(all_species_raw, Level
                            <- factor(Level, levels=data_names,
                                      ordered=T))
  setkeyv(all_species_raw, sample_var)
  all_species_raw <- metadata[all_species_raw]
  all_species_raw <- all_species_raw[!is.na(all_species_raw[[group_var]]), ]
  species_raw_type_sums <- all_species_raw[Level==data_names[2], median(Value), by=group_var]
  species_raw_value_labels <- as.character(species_raw_type_sums[[group_var]][order(species_raw_type_sums$V1,
                                                                                    decreasing=T)])
  all_species_raw[[group_var]] <- factor(all_species_raw[[group_var]],
                                         levels=species_raw_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_species_raw[[group_var]] <- factor(all_species_raw[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # Plot figure
  png(filename=paste(outdir, '/', 'Diversity_', data_type, '_raw_richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_sraw <- ggplot(data=all_species_raw, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_sraw <- g_sraw + scale_fill_tableau("Classic 20", direction = -1)
  g_sraw <- g_sraw +
    ggtitle(paste('Species Richness by ', group_var, ' for Raw data\nRichness',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_sraw)
  dev.off()
  
  
  # Output Diversity counts
  write.table(all_species_raw, paste(outdir, '/', data_type, '_raw_data_richness_values.csv',
                                     sep='', collapse=''), sep=",")
  
  
  #
  ## Rarefied species abundance
  #
  all_species_rare <- within(all_species_rare, Level
                             <- factor(Level, levels=data_names,
                                       ordered=T))
  setkeyv(all_species_rare, sample_var)
  all_species_rare <- metadata[all_species_rare]
  all_species_rare <- all_species_rare[!is.na(all_species_rare[[group_var]]), ]
  species_rare_type_sums <- all_species_rare[Level==data_names[2], median(Value), by=group_var]
  species_rare_value_labels <- as.character(species_rare_type_sums[[group_var]][order(species_rare_type_sums$V1,
                                                                                      decreasing=T)])
  
  all_species_rare[[group_var]] <- factor(all_species_rare[[group_var]],
                                          levels=species_rare_value_labels, ordered=T)
  ## Check if needing to change factor order.
  ifelse(factor_order != '', all_species_rare[[group_var]] <- factor(all_species_rare[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  png(filename=paste(outdir, '/', 'Diversity_',data_type, '_rarefied_richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_srare <- ggplot(data=all_species_rare, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_srare <- g_srare + scale_fill_tableau("Tableau 20", direction = -1)
  g_srare <- g_srare +
    ggtitle(paste('Species Richness by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_srare)
  dev.off()
  
  # Output Diversity counts
  write.table(all_species_rare, paste(outdir, '/', data_type, '_rarefied_data_richness_values.csv',
                                      sep='', collapse=''), sep=",")
}






##############################
###                        ###
####     Bar plots        ####
###                        ###
##############################
bar_select_top_counts <- function(X, group_var, n) {
  return(X[, tail(.SD, n), by=group_var])
}



meg_relative_barplot <- function(melted_data,
                                 metadata,
                                 sample_var,
                                 group_var,
                                 level_var,
                                 analysis_subset,
                                 outdir,
                                 data_type,
                                 factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, sum(Normalized_Count) / sample_number,
                                  by=c(group_var, 'Name')])
  ## We should improve how we select the number of taxa to be plotted, for now it's arbitrarily 10 for easy color selection
  numselect <- 21
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 21) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # It's not my favorite thing to call the "by" flag with a column, but I couldn't get the group_var variable to work
  bar_subset[, total := sum(Normalized_Count), by = .(bar_subset[[group_var]])][, percentage := Normalized_Count/total , by = .(Name,bar_subset[[group_var]])]
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='percentage', fill='Name')) +
    geom_bar(stat='identity') +
    #scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Relative abundance\n') +
    #scale_fill_tableau("Classic 20", direction = -1) +
    ggtitle(paste(data_type, ' ', level_var, ' Relative Abundance by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/','RelativeAbundance_BarPlot_', data_type, '_', level_var, '_by_',group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}






meg_barplot <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type,
                        factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out, sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, sum(Normalized_Count) / sample_number,
                                  by=c(group_var, 'Name')])
  
  numselect <- 11
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 11) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  bar_subset$Name <-  droplevels(bar_subset$Name)
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
    geom_bar(stat='identity') +
    scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Mean of Normalized Count\n') +
    ggtitle(paste('Mean ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/','MeanAbundance_BarPlot_', data_type, '_', level_var, '_by_',group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}





meg_median_barplot <- function(melted_data,
                               metadata,
                               sample_var,
                               group_var,
                               level_var,
                               analysis_subset,
                               outdir,
                               data_type,
                               factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out, sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, median(round(Normalized_Count)),
                                  by=c(group_var, 'Name')])
  
  numselect <- 11
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 11) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  bar_subset$Name <-  droplevels(bar_subset$Name)
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
    geom_bar(stat='identity') +
    scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Normalized Count\n') +
    ggtitle(paste('Median ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/','MedianAbundance_BarPlot_', data_type, '_', level_var, '_by_',group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}

##############################
###                        ###
####     ZIG model        ####
###                        ###
##############################

meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval=0.1,
                       top_hits=200) {
  settings <- zigControl(maxit=20000, verbose=F)
  
  local_obj <- data_list
  res <- list()
  for( l in 1:length(local_obj) ) {
    filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
    if(length(analysis_subset) > 0) {
      local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
      amr_sparseFeatures = which(rowSums(MRcounts(local_obj[[l]]) > 0) < 1) ## Check to see if any features in the subset == 0
      if (length(amr_sparseFeatures) != 0) {
        local_obj[[l]] = local_obj[[l]][-amr_sparseFeatures,]
      }
    }
    
    
    col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) > 1))
    local_obj[[l]] <- local_obj[[l]][, col_selection]
    
    mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_obj[[l]]))
    zero_mod_select <- zero_mod[col_selection, ]
    
    cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
    
    tryCatch(
      {
        if( is.na(random_effect_var) ) {
          res[[l]] <- fitZig(obj=local_obj[[l]],
                             mod=mod_select,
                             zeroMod=zero_mod_select,
                             control=settings,
                             useCSSoffset=F)
        }
        else {
          res[[l]] <- fitZig(obj=local_obj[[l]],
                             mod=mod_select,
                             zeroMod=zero_mod_select,
                             control=settings,
                             useCSSoffset=F,
                             useMixedModel=T,
                             block=pData(local_obj[[l]])[, random_effect_var])
        }
      },
      error=function(e) {
        print(paste('Model failed to converge for ', data_type, ' ', data_names[l], ' ', analysis_name,
                    sep='', collapse=''))
      },
      finally={
        if( length(res) != l ) {
          next
        }
      }
    )
    
    local_contrasts <- contrast_list
    local_contrasts[[length(local_contrasts)+1]] <- res[[l]]$fit$design
    names(local_contrasts)[length(local_contrasts)] <- 'levels'
    
    contrast_matrix <- do.call(makeContrasts, local_contrasts)
    colnames(contrast_matrix) <- make.names(contrast_list)
    
    contrast_fit <- contrasts.fit(res[[l]]$fit, contrast_matrix)
    contrast_fit <- eBayes(contrast_fit)
    
    stats_results <- data.table(
      Node.Name=character(),
      Contrast=character(),
      logFC=numeric(),
      CI.L=numeric(),
      CI.R=numeric(),
      AveExpr=numeric(),
      t=numeric(),
      P.Value=numeric(),
      adj.P.Val=numeric(),
      B=numeric()
    )
    
    for( c in 1:ncol(contrast_fit$contrasts) ) {
      tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                          number=top_hits, sort.by='AveExpr', coef=c)
      
      if( nrow(tophits) > 0) {
        temp_res <- data.table(
          Node.Name=rownames(tophits),
          Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
        )
        temp_res <- cbind(temp_res, tophits)
        stats_results <- rbind(stats_results, temp_res)
      }
      else {
        print(paste('No significant results for', data_type,
                    data_names[l], analysis_name,
                    colnames(contrast_fit$contrasts)[c],
                    sep=' ', collapse=''))
      }
    }
    
    if( nrow(stats_results) > 0 ) {
      write.csv(stats_results,
                file=paste(outdir, '/', analysis_name, '_', data_type, '_',
                           data_names[l], '_',
                           contrast_list[1], '_Model_Contrasts.csv',
                           sep='', collapse=''),
                quote=F, row.names=F)
    }
  }
}



myheatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                         distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                            "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                                          w), symm = FALSE, scale = c("none", "row", "column"), 
                         na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                         symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
                         col = "heat.colors", colsep, rowsep, sepcolor = "white", 
                         sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
                         na.color = par("bg"), trace = c("column", "row", "both", 
                                                         "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
                         linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
                         cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                         labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
                                                                                 NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
                         colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5, 
                         density.info = c("histogram", "density", "none"), denscol = tracecol, 
                         symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
                         key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, 
                         key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
                         ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, 
                         ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && any(duplicated(breaks))) 
    stop("breaks may not contian duplicate values")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (all(Colv == "Rowv")) 
    Colv <- Rowv
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
        (dendrogram %in% c("both", "row"))) {
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
      if (dendrogram == "both") 
        dendrogram <- "column"
      else dendrogram <- "none"
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
        (dendrogram %in% c("both", "column"))) {
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
      if (dendrogram == "both") 
        dendrogram <- "row"
      else dendrogram <- "none"
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
                                   nr)) 
      stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr) 
      nr <- length(rowInd)
  }
  else if (is.integer(Rowv)) {
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (!isTRUE(Rowv)) {
    rowInd <- nr:1
    ddr <- as.dendrogram(hclust(dist(diag(nr))))
  }
  else {
    rowInd <- nr:1
    ddr <- as.dendrogram(Rowv)
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd > 
                                   nc)) 
      stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc) 
      nc <- length(colInd)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    distc <- distfun(if (symm) 
      x
      else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    distc <- distfun(if (symm) 
      x
      else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (!isTRUE(Colv)) {
    colInd <- 1:nc
    ddc <- as.dendrogram(hclust(dist(diag(nc))))
  }
  else {
    colInd <- 1:nc
    ddc <- as.dendrogram(Colv)
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (!is.null(colRow)) 
    colRow <- colRow[rowInd]
  if (!is.null(colCol)) 
    colCol <- colCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
      1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) != 
          nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                      1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) != 
          nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                           1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  plot.index <- 1
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!gtools::invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol) && is.null(colCol)) 
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
           offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
         padj = adjCol[2])
  else {
    if (is.null(srtCol) || is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
        adjCol = c(1, NA)
      if (is.null(srtCol)) 
        srtCol <- 90
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
             strheight("M"), labels = labCol, adj = adjCol, 
           cex = cexCol, srt = srtCol, col = colCol)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow) && is.null(colRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
  }
  else {
    if (is.null(srtRow) || is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
           srt = srtRow, col = colRow)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in 1:length(colInd)) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in 1:length(rowInd)) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, 
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  plot.index <- plot.index + 1
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
                                yaxs = "i", leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
        stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
                                leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
        stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }
  else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    mar <- c(5, 4, 2, 1)
    if (!is.null(key.xlab) && is.na(key.xlab)) 
      mar[1] <- 2
    if (!is.null(key.ylab) && is.na(key.ylab)) 
      mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title)) 
      mar[3] <- 1
    par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0) 
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min.breaks
      max.raw <- max.breaks
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
    image(z = matrix(z, nrow = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at = xv, labels = lv, las=2)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 2
    do.call(axis, xargs)
    if (is.null(key.xlab)) {
      if (scale == "row") 
        key.xlab <- "Row Z-Score"
      else if (scale == "column") 
        key.xlab <- "Column Z-Score"
      else key.xlab <- "Value"
    }
    if (!is.na(key.xlab)) {
      mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
            cex = par("cex") * par("cex.lab"))
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE, 
                      from = min.scale, to = max.scale)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                        0.95, labels = pretty(dens$y))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title)) 
        key.title <- "Color Key\nand Density Plot"
      if (!is.na(key.title)) 
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab)) 
        key.ylab <- "Density"
      if (!is.na(key.ylab)) 
        mtext(side = 2, key.ylab, line = par("mgp")[1], 
              padj = 0.5, cex = par("cex") * par("cex.lab"))
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                      labels = pretty(hy))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title)) 
        key.title <- "Color Key\nand Histogram"
      if (!is.na(key.title)) 
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab)) 
        key.ylab <- "Count"
      if (!is.na(key.ylab)) 
        mtext(side = 2, key.ylab, line = par("mgp")[1], 
              padj = 0.5, cex = par("cex") * par("cex.lab"))
    }
    else if (is.null(key.title)) 
      title("Color Key")
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }
  else {
    par(mar = c(0, 0, 0, 0))
    plot.new()
  }
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
  if (!is.null(extrafun)) 
    extrafun()
  invisible(retval)
}
