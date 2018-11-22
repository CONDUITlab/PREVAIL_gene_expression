
#############################
# baseline_clusters.R       #
# Baseline subgroups        #
# Does clustering reveal    #
# any cohesive subtypes?	#	
# Uses COMMUNAL clustering	#
# package (CRAN archive)	#
#############################

source('setup.R')

# source('lactoMainFunctions.R')
# source('plot_range_3d.R')       # For COMMUNAL functions

# Reduce filtered.eset to baseline samples
baseline.eset <- filtered.eset[, filtered.eset$`day:ch1`=='1']
# Reduce to the 10% most variable genes
baseline.99 <- varFilter(baseline.eset, var.cutoff = 0.899) # Chose 89% to get 2500 features
# Setup parameters for COMMUNAL clustering
varRange <- seq(250, 2500, 250)
ks <- 2:8 # look at 2 to 8 clusters
data <- exprs(baseline.99)
# Choice of cluster cohesiveness measures
measures <- c("average.between", "dunn", "widestgap", "dunn2",
              "pearsongamma", "g3", "max.diameter", "avg.silwidth")
comm.results <- clusterRange(dataMtx=data, ks = ks,
                             varRange=varRange,
                             validation=measures, 
                             clus.methods = c('hierarchical', 'sota', 'pam', 'clara', 'agnes'),
                             verbose = T)
algs <- getGoodAlgs(comm.results, algs="all")
monotoneClusterRange(comm.results)
measuresCorr(comm.results)
measures <- getNonCorrNonMonoMeasures(comm.results, goodAlgs=algs, numMeasures = 4)
plot.data <- plotRange3D(comm.results, ks, algs, measures, plot3D=T)
#### plots show that data best support a 3 cluster solution ###

result <- comm.results$all.results$vars_250
clusters <- result$getClustering(k=3)
apply(clusters, 2, table)
mat.key <- clusterKeys(clusters)
examineCounts(mat.key)
# find 'core' clusters
core <- returnCore(mat.key, agreement.thresh=50) 
# Visualize with PCA
pca.coords <- prcomp(Biobase::exprs(baseline.99))$rotation[,1:2]
# cluster labels are in the same order so can simply be added
pca.data <- data.frame(pca.coords, cluster=core)
p <- ggplot(pca.data, aes(x=PC1, y=PC2, colour=cluster)) +
  geom_point(size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(face='bold', size='12')) + 
  theme(axis.title = element_text(face='bold', size=16), axis.text = element_text(size=14)) + 
  theme(strip.text.x = element_text(face='bold', size=12)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# Get details on clustering measures
result$measures
# MEAN SILHOUETTE WIDTH =0.12


