#----------------------------------------------------------------------------
# Script:  Broughton2_functions.R
# Created: May 2024 from the original Broughton script. EJG
#
# Purpose: Support building and evaluation of Seaweed Clusters for Coastal BC.
# Illustrated at 3 spatial extents using the Broughton Region.
#
# Notes:
#  - 2024/06/05: Updated from the almost final DFO version.


#================================== Load require packages =================================

# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "tidyr","dplyr", "raster", "stringr", "rasterVis",
                        "RColorBrewer", "factoextra", "ggpubr", "cluster", "rmarkdown","lubridate",
                        "knitr", "tinytex", "kableExtra", "e1071")

# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "tinytex", "rmarkdown", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)


#=========================== Data sources and constants =====================================

# proj4 string for albers projection with NAD83 datum
spat_ref <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'


#===================================== Functions =========================================

#---- Loads predictors from specified subdirectory -----
LoadPredictors <- function( pred_dir ) {
  
  print( list.files(path = pred_dir, pattern = '\\.tif$', full.names = FALSE) )
  
  # make a list of predictor rasters
  raster.list <- list.files(path = pred_dir, pattern = '\\.tif$', full.names = TRUE)
  
  
  # make list of raster names (without extension)
  raster.names <- lapply(raster.list, FUN = function(raster.layer){
    substr(basename(raster.layer), 1, nchar(basename(raster.layer)) - 4)
  } )
  
  # create a raster stack from raster_list
  # Go thru one at a time so that we can set NA values ... needed for Romina's rasters
  raster.stack <- stack()
  for (i in raster.list) {
    new.layer <- raster( i )
    if (cellStats(new.layer, min) < -5) {
      values(new.layer)[ values(new.layer) < 0 ] <- NA }
    raster.stack <- addLayer( raster.stack, new.layer)
  }
  return(raster.stack)
}


#---- ClipPredictors: masks away deeper water and limits extents based on a polygon mask
ClipPredictors <- function( stack_in, the_mask){
  z <- stack()
  for (i in 1:dim( stack_in)[[3]] ) {
    x <- raster::crop(stack_in[[i]], the_mask )
    y <- raster::mask(x, the_mask)
    z <- stack( z, y )
    print( paste0('layer ', i, ' clipped.'))
  }
  return( z )
}


#---- Returns a stack of integerized rasters from a raster stack ----
Integerize <- function( in_layers, sig = 1000 ) {
  int_layers <- brick()
  for (i in 1:nlayers(in_layers)) {
    raster_to_scale <- in_layers[[i]]
    int_raster <- as.integer( raster_to_scale * sig )
    int_layers <- stack( int_layers, int_raster )
  }
  
  names( int_layers ) <- names( in_layers )
  return( int_layers )
}


PlotHists <- function( stack_dat ) {
  sk_list <- NULL
  w <- colnames(stack_dat)
  
  for (i in 1:dim( stack_dat )[2]){
    x <- stack_dat[,i]
    y <- skewness(x, na.rm = TRUE)
    sk_list <- c(sk_list, y)
    y <- paste0( "Skew=", round(y, 3))
    hist( x, main = w[i] )
    xy <- par("usr")
    xyrange <- c( xy[2]-xy[1], xy[4]-xy[3] )
    text( xy[1]+xyrange[1]/4, xy[4]-xyrange[2]/10, y )
  }
}

#---- Four layers that needed some work. ----
#Includes 4: julBS_ave, julSS_min, julSSpd_ave, tidal_cur
MakeMoreNormal <- function( the_stack, t_list ){

  x <- the_stack[, "julBS_ave"]
  y <- asin(x/max(x))
  the_stack[, "julBS_ave"] <- y
  
  x <- the_stack[, "julSS_min"]
  y <- x^3
  the_stack[, "julSS_min"] <- y
  
  x <- the_stack[, "julSSpd_ave"]
  y <- x^(1/2)
  the_stack[, "julSSpd_ave"] <- y
  
  x <- the_stack[, "tidal_cur"]
  y <- x^(1/2)
  the_stack[, "tidal_cur"] <- y
  
  return( the_stack) 
}

#---- MakeScreePlot: returns a ggplot. ----
# samp is optional, uses all dat if omitted.
MakeScreePlot <- function( indat, nclust, nrand, maxi, sampsize = 0 ){
  #initialize list for results
  wss <- numeric(nclust) 
  
  #subsample as requested
  if (sampsize > 0) {
    samp <- sample( 1:length( indat[ , 1] ), sampsize )
    dat <- indat[ samp, ]
  } else dat <- indat
  
  for (i in 1:nclust) {
    # Fit the model: km.out
    print( paste0( "centers ",i))
    km.out <- kmeans(dat, centers = i, nstart = nrand, iter.max = maxi)
    # Save the within cluster sum of squares
    wss[i] <- km.out$tot.withinss
    
    # calculate the silhouete width
  }
  
  # Produce the scree plot ... using a tibble, I guess. 
  wss_df <- tibble(clusters = 1:nclust, wss = wss)
  scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    xlab('Number of clusters') +
    ylab('Total within-cluster sum of squares') +
    geom_hline(
      yintercept = wss, 
      linetype = 'dashed')
  
  return( scree_plot)
}


#---- ClusterPCA: Returns a pair of PCA plots showing the separation of the clusters.
# NOTE: Uses a relatively small subset of the overall data so these will change with a different sample.
ClusterPCA <- function( n_samp, clustnum ) {
  
  ssidx <- sample( 1:length( stack_data_clean[ , 1] ), n_samp )
  ssamp <- stack_data_clean[ ssidx, ]
  
  # re-run cluster for smaller sample.
  cluster_result <- kmeans(ssamp, centers = clustnum, nstart = randomz) # less than 10 seconds
  csamp <- cluster_result$cluster
  
  # Create the data structure for PCA profiling (combining predictor data with the clusters). 
  # Put cluster # first so easy to find for PCA
  p_data <- as.data.frame( cbind(cluster = csamp, ssamp ) )
  
  res_pca <- prcomp(p_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
  ind_coord <- as.data.frame(get_pca_ind(res_pca)$coord)
  # Add clusters from the classification
  ind_coord$cluster <- factor( p_data$cluster )
  # Data inspection
  #head(ind_coord)
  
  # Percentage of variance explained by dimensions
  eigenvalue <- round(get_eigenvalue(res_pca), 1)
  var_percent <- eigenvalue$variance.percent
  
  # Look at the clusters for the first 4 PCs
  a <- ggscatter(
    ind_coord, x = "Dim.1", y = "Dim.2", 
    color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 1 (", var_percent[1], "% )" ),
    ylab = paste0("Dim 2 (", var_percent[2], "% )" )
  ) +
    stat_mean(aes(color = cluster), size = 4)

  pca_loadings <- data.frame(res_pca$rotation)
  pca_scores   <- data.frame(res_pca$x)
  pca_loadings <- pca_loadings * max(abs(pca_scores$PC1), abs(pca_scores$PC2))

  arrow_plot <- a +
    geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = pca_loadings, aes(x = PC1, y = PC2, label = rownames(pca_loadings)), 
              hjust = 0, vjust = 1, color = "red")

  b <- ggscatter(
    ind_coord, x = "Dim.3", y = "Dim.4", 
    color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 3 (", var_percent[3], "% )" ),
    ylab = paste0("Dim 4 (", var_percent[4], "% )" )
  ) +
    stat_mean(aes(color = cluster), size = 4)
 
  print( "Done clusters ...")
  return( list(res_pca, arrow_plot, b))
}



#---- PredictClusters: returns cluster assignments for un-clustered pictures. ----
# Assembly required with the classified pixels before plotting. Function by ChatGPT.
PredictClusters <- function(newdata, kmeans_model) {
  centers <- kmeans_model$centers
  dists <- as.matrix(dist(rbind(centers, newdata)))
  dists <- dists[-(1:nrow(centers)), 1:nrow(centers)]
  cluster_assignments <- apply(dists, 1, which.min)
  return(cluster_assignments)
}


#---- TransferClusters: returns new values for the prediction rasters containing predictions. ----
# transferred from the random sample used in the classification ----
# Assembly required with the classified pixels before plotting. Function by ChatGPT.
transferCluster <- function(values_target, c_result){
# uses GLOBALS sidx, stack_data_clean   
  #--- Two steps here: First assign clusters to the rest of the clean stack data, 
  #     THEN put the clean data back in the raster.
  # Uses samp and sidx from lines ~150 above.
  # Find the things not in the sample (using sidx from line ~150 above).
  
  not_sidx <- setdiff( 1:dim(stack_data_clean)[[1]], sidx )
  
  # the data to predict clusters for
  new_dat <- stack_data_clean[ not_sidx, ]
  # Process the new data in chunks of 10k to ensure performance
  chunk_size <- 10000
  n_chunks <- ceiling(nrow(new_dat) / chunk_size)
  
  print( paste0( "predicting clusters for ", dim(new_dat)[[1]], " pixels using ", 
                 n_chunks, " chunks of ", chunk_size, ". Stand by ... ") )
  
  # Where we are putting our new chunks
  predicted_clusters <- vector("integer", nrow(new_dat))
  
  for ( i in seq_len(n_chunks) ) {
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, nrow(new_dat))
    
    chunk <- new_dat[start_index:end_index, ]
    predicted_clusters[start_index:end_index] <- PredictClusters(chunk, c_result)
    cat(i, " ")
  }
  
  #length(predicted_clusters)
  #length(c_result$cluster)
  # Combine the results for the cleaned data ... 
  clean_clusts <- array( 1:dim(stack_data_clean)[1] )
  clean_clusts[ sidx ] <- c_result$cluster
  clean_clusts[ not_sidx ] <- predicted_clusters
  
  # replace clustered values with the cluster results
  values_target[ clean_idx ] <- clean_clusts
  # reassign to the plotting raster
  return( values_target )
}




####------------Depreciated or replaced functions----------------------------

#---- TrimStack: limit extents based on a polygon mask (vestigial)
# NOTE: This function now superseded by the crop/mask approach in ClipPredictors.
TrimStack <- function( stack_in, padsize ) {
  # Assigning extents here worked but results odd - plots appeared but looked empty.
  y <- stack()
  # estimate extents from first raster ... 
  x <- trim(stack_in[[1]], padsize )
  x_ext <- extent(x)
  x_ext <- ceiling( x_ext )
  extent( x ) <- x_ext
  y <- stack(y, x)
  print( paste0('layer 1 trimmed.'))
  
  # use these extents for the remaining layers ...
  for (i in 2:dim( stack_in)[[3]] ) {
    x <- setExtent( stack_in[[i]], x_ext, keepres=T, snap=T )
    y <- stack( y, x )
    print( paste0('layer ', i, ' trimmed.'))
  }
  return( y )  
}


# if (spacesub) {
#   
#   sMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
#   smMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
#   
#   # Data extents No mask and polygons
#   x  <- prepped_layers
#   
#   # Use only one of the loaded shapefiles
#   x  <- raster::mask(x, sMask )
#   # x  <- raster::mask(x, smMask )
#   
#   # Plot first layer in the resulting stack for inspection
#   plot(trim(x[[1]]), main = 'Local')
#   
#   # Plot masks over raster. NB: will only see if raster extends beyond the masks. 
#   sp::plot( sMask, border = "darkgreen", lwd = 2, add=TRUE)
#   sp::plot( smMask, border = "red", lwd = 2, add=TRUE)
#   
#   prepped_layers <- x
#   rm('x')
# }



ScalePredictors <- function( the_stack ){
  # A little shimmy to avoid scaling substrate and name it consistently
  #  scaled_stack <- scale( dropLayer( the_stack, "SUBSTRATE") )
  #  scaled_stack <- stack( scaled_stack, trim_layers$SUBSTRATE )
  # safe to assume its the last layer in the stack
  #  names( scaled_stack[[ dim(scaled_stack)[[3]] ]] ) <- "substrate"
  
} 


# Apply scaling to standardize rasters. Does either min/max or z-score.
# Returns a scaled raster.
Scale.Raster <- function( in.raster, scale.how = "mm", intshift = 'T' ) {
  
  if (scale.how == "z") {
    mean.value <- cellStats( in.raster, stat='mean', na.rm=TRUE )
    sd.value   <- cellStats( in.raster, stat='sd', na.rm=TRUE )  
  } else {
    old.min    <- cellStats( in.raster, stat='min', na.rm=TRUE )  
    old.max    <- cellStats( in.raster, stat='max', na.rm=TRUE )  
  }
  
  if (scale.how == "z") {
    # Perform z-score normalization
    scaled.raster <- (in.raster - mean.value) / sd.value
  } else {
    # Perform min/max scaling
    scaled.raster <- (in.raster - old.min) / (old.max - old.min)
  }
  
  if (intshift) {
    scaled.raster <- as.integer( scaled.raster * 1000 )
  }
  
  return( scaled.raster )
}


# fin
