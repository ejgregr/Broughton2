#################################################################################
# Script:  Broughton_main_LSSM.R - Broughton LSSM version
# Created: February 2024. EJG
# 
# This script sources the necessary libraries and functions, coordinates the 
# analysis, and creates data structures that are then 'knitted' together (I think).
# So the idea is to run bits in here, and then 'Render' the RMD script. 
# Seems straightforward. :)
#
## Updates: 
# 2024/04/29: Steady development past few weeks; alll necessary pieces now developed. 
# 2024/04/29: Git repo created and working version pushed prior to RMarkdown development.
# 2024/05/02: Completed smoothing pass thru code; sorted raster plotting. Pushed.
# 2024/05/07: Another pass thru, adding some controls. Ready for RMD work.Pushed.
# 2024/05/29: After some exploratory work with Romina's tifs from SPECTRAL lab, forked the 
#   code so this version focuses on the LSSM needs.
# 2024/06/05: Updated with relevant code from latest DFO version.
# 2024/08/12: Pre-holidays, worked thru code and began work on rMarkdown script. This script
#   is now ahead of the DFO version. Refamiliarization in progress. 
# 2024/08/30: Back to work on classification. Began by renaming the 2 projects for clarity, and 
#   moved all progress over to DFO version -> finish that first.
# 2024/09/22: DFO work close to done. Nicely developing RMD. Back to here to create 
#  initial clusters for the LSSM. 
#  - 2 key bits to integrate: distributions and scaling, and RMD report.

# TO DO: 
#  1. Reproduce (ish) Romina - use as prototype
#  2. Combine attributes and develop update plan
#  3. Spend a TINY bit of time to see if MSEA data can "improve" fit.
#  Add config section to allow repeated clusterings to be attached to same RMD report
######################################################################################

print('Starting Broughton ... LSSM Version')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "broughton_functions.R" )
# source( "Plot_Functions.R" )

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
raster_dir <- 'C:/Data/SpaceData/Classification/Romina'
data_dir   <- 'C:/Data/Git/Broughton_LSSM/Data'
results_dir<- 'C:/Data/Git/Broughton_LSSM/Results' 

# Processing FLAGS...
loadtifs <- F # If true the data will be re-loaded from TIFs, else it will be loaded from rData.
clipdata <- F # If true a spatial subset of the data will be taken based on a polygon shape file. 
scaledat <- F # If true, imported data will be scaled using raster::scale().
reclust  <- T # If true, re-cluster full data set prior to mapping, else predict to unclassified pixels.

#---- Part 1 of 3: Load, clean, and prepare predictor data.  ----
# If loadtifs == TRUE then run all this else load the processed data.

tif_stack <- stack()
today <- format(Sys.Date(), "%Y-%m-%d")

if (loadtifs) {
  print( "Loading predictors ... ")
  src_stack <- LoadPredictors( raster_dir )
  print( "Data loaded.")
  
  tif_stack <- src_stack
  
  if (clipdata) {
    print( "clipping TIFs ... ")
    #amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
    amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
    tif_stack <- ClipPredictors( tif_stack, amask )
    print('Rasters clipped.')
  }
 
# Romina's Ocean TIFs already trimmed to 40 m and did not contain any land.  
# IF we want to include any DFO data it will need to match Romina's spatial reference.

  if (scaledat) {
    print( "Scaling TIFs ... ")
    tmp_stack <- scale( tif_stack )
    tif_stack <- tmp_stack
    print('Rasters Scaled.')
  }
  
  save( src_stack, file = paste0( data_dir, '/src_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tif_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tifs_SPECTRAL_scaled_', today, '.rData' ))
  
} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
  load( paste0( data_dir, '/tifs_SPECTRAL_scaled_2024-08-14.rData' ))
}

# Standardize rasters spatially (i.e., trim northness and wind to ocean layers)

# Move to matrix space from raster space 
stack_data <- getValues( tif_stack )

#-- Quick correlation across data layers
x <- stack_data
x_clean <- x[ complete.cases(x), ]
cor_table <- cor( x_clean )
cor_table[lower.tri(cor_table)] <- NA
(cor_table >= 0.6) & (cor_table != 1)

# Drop correlated layers layers from the SPECRAL data set
selected_stack <- dropLayer( tif_stack, c("julST_ave", "julSS_ave", "julBSpd_ave", "northness") )
names( selected_stack )

stack_data <- getValues( selected_stack )

# remove any rows with an NA
# NB: This decouples the data from the RasterStack and requires re-assembly for mapping
# THESE are the two key data structures used in subsequent steps
clean_idx <- complete.cases(stack_data)
stack_data_clean <- stack_data[ clean_idx, ]

dim( stack_data )
dim( stack_data_clean )

# Transform and scale the data for k-means classification.

print( "Transforming data  ... ")
  # Transform those predictors with abs(skew) > 1. 
  #Includes 4: julBS_ave, julSS_min, julSSpd_ave, tidal_cur
t_stack_data <- MakeMoreNormal(stack_data_clean)

print("Centering and scaling  ... ")
tmp_stack <- scale(t_stack_data, center = T,  scale = T)
t_stack_data <- tmp_stack
print('Data prepped.')
#save(t_stack_data,
#     file = paste0(data_dir, '/t_stack_data', today, '.rData'))
#print('Scaled data saved.')

#---- Part 2 of 3: Cluster number selection ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

#---- Part 2a: Explore number of clusters using Within-sum-of-squares scree plot ----
# Runs kmeans with increasing number of clusters

nclust   <- 18 # number of clusters for scree plot
nsample  <- 25000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( t_stack_data, nclust, randomz, imax, nsample )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 8 # the number of clusters based on scree plot, above.
nsample <- 500000 # a larger sample for more robust classification

sidx <- sample( 1:length( t_stack_data[ , 1] ), nsample )
samp <- stack_data_clean[ sidx, ]
cluster_result <- kmeans(samp, centers=nclust, nstart=randomz, iter.max=imax) 

#---- Part 2b: Create heat map of within-cluster standard deviations ----

# Define color palette
pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette

profile_data <- as.data.frame( cbind(cluster = cluster_result$cluster, samp ) )

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
head(x)
xm <- melt( x, id.var = "cluster" )

z_heat <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "", x = "Clusters", y = "Predictors", fill = "Value")
z_heat

#---- Part 2c: Examine silhouette plot of the WORKING clusters  ----
# Uses the predictor values and the corresponding assigned cluster
# Need to subsample from the cluster result above as distance matrix take long time.

# Take a subsample of the clusters and the predictors for the silhouette plot. 
sil_n <- 20000
silx <- sample( 1:length( samp[ , 1] ), sil_n )

cs <- cluster_result$cluster[ silx ]
ss <- samp[ silx, ]

# Calculate a distance matrix between the predictors and plot assigned to each cluster.
# both these steps are time consuming, hence a smaller sample.
c_dist <- dist(ss)
sk <- silhouette(cs, c_dist)
#mean( sk[,"sil_width"] )

plot(sk, col = 1:nclust, border=NA, main = "" )


#---- Part 3: Detailed examination of N clusters  ----
#---- Part 3a: Show cluster groupings using PCA ----

#-- Can take some time so it makes its own cluster 
pca_n <- 25000

pca_results <- ClusterPCA( pca_n, nclust ) # uses global variable stack_data_clean
names(pca_results) <- c("loadings", "plot1","plot2")

#Percentage of variance explained by dimensions
#eigenvalue <- round(get_eigenvalue(res_pca), 1)
#var_percent <- eigenvalue$variance.percent

#---- Part 3b: Violins of predictor contributions to WORKING clusters ----

x <- as.data.frame( samp )
x$cluster <- as.factor( cluster_result$cluster )

y <- x %>%
  pivot_longer(cols = -cluster, names_to = "predictor", values_to = "value")

# Create violin plot
vplots <- 
  ggplot(y, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ predictor, scales = "free_y") +
  theme_minimal() +
  labs(title = "Violin Plots of Predictors Across k-means Clusters",
       x = "Cluster",
       y = "Value")
vplots

#---- Part 3c: Spatialize and display the WORKING clusters ----
# NB: To show a comprehensive map, can either:
#     a) re-cluster the entire data set (using imax and randomz from above) or
#     b) Predict to the unsampled portion of the raster. 

# initialize target data structure 
cluster_raster <- selected_stack[[1]]
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

if (reclust == T) {
  # Re-cluster using all the clean data  ... 
  # less than 1 min with iter.max = 20, nstart = 20 for smallest region
  cluster_result <- kmeans(t_stack_data, centers = nclust, nstart = randomz, iter.max = imax)
  
  # Assign the clustered values ... 
  # extract values from the target cluster
  new_values <- values( cluster_raster )
  # replace non-NA values with the cluster results
  new_values[ clean_idx ] <- cluster_result$cluster
  # put the updated values back on the target cluster
  values( cluster_raster ) <- new_values  
} else {
  # Predict values for unclustered cells. Can be more time-consuming than re-classifying everything. 
  values( cluster_raster ) <- transferCluster( values(cluster_raster), cluster_result )
}

#--- Display the results, first as histogram then as map.  
raster::hist( values(cluster_raster ))

# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent") # Max for Accent is 8

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
z_map <- levelplot( cluster_raster, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )
z_map
writeRaster( cluster_raster, paste0( results_dir, "/SPEC_8clusterb.tif"), overwrite=TRUE)


#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
rmarkdown::render( "Broughton_LSSM_PDF.Rmd",   
#rmarkdown::render( "Broughton_test.Rmd",   
                                      output_format = 'pdf_document',
                   output_dir ="C:/Data/Git/Broughton_LSSM/Results",
                   output_file = paste0( "LSSM_gdata_8cluster_", today ))


#---- Some details on correlation analysis ... ----
#-- Correlation across UN-scaled data layers ... 
# foo <- getValues( scaled_layers )
# foo_clean <- na.omit(stack_data)
# pick <- sample( 1:length( foo_clean[ , 1] ), 10000 )
# plot( foo_clean[ pick,'rugosity'] ~ foo_clean[ pick,'standard_deviation_slope'] )

#---- ----

# FIN.






