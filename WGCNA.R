# WGCNA
# build gene network, extract MEs, and test MEs for DE 

# last modified: Nov 13rd, 2018

source('setup.R')

options(stringsAsFactors = FALSE)

collapsedExprs <- Biobase::exprs(corrected.eset)
datExprs <- t(collapsedExprs)
samples <- rownames(datExprs)

### determine soft-thresholding power to build network 
# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(datExprs,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

# plot results
sizeGrWindow(10,7)
par(mfrow = c(1,2))
cex1 = 0.9

# scale-free topology fit index as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     labels=candidatePowers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5], labels=candidatePowers, cex=cex1,col="red")

### find modules 
# getModules takes an expression matrix and soft-thresholding power and returns a network
# note: reduce maxBlockSize if running on a computer with <16GB memory
getModules <- function(datExprs,sfPower) {
  bwnet = blockwiseModules(datExprs, maxBlockSize = 15000,
                           power = sfPower, networkType= "signed", TOMType = "signed", minModuleSize = 50,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           corType = "bicor", maxPOutliers = 0.05,
                           verbose = 3)
  return (bwnet)
}

### comparing expression of module eigengenes
# use default soft-thresholding power of 12
net <- getModules(datExprs,20) 
# view number of modules and size of modules 
table(net$colors)
moduleLabelsAutomatic=net$colors
# convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# get data frame with module eigengenes 
MEsAutomatic=net$MEs

allDays <- filtered.eset$`day:ch1`
allCodes <- filtered.eset$`treatment_arm:ch1`

getESet <- function(filteredESet,earlyDays,lateDays) {
  # getESet subsets an expression set based on samples from specific days
  # earlyDays and lateDays are arrays containing the desired timepoints, 
  # ex. c("day 1","day 3")
  earlyLateDays = c(earlyDays,lateDays)
  earlyLateIndices <- logical(ncol(filteredESet))
  for (i in 1:length(earlyLateDays)) {
    earlyLateIndices <- earlyLateIndices | filteredESet$`day:ch1`==earlyLateDays[i]
  }
  earlyLate.eset <- filteredESet[,earlyLateIndices]
  return (earlyLate.eset)
}

# create expression set with only early and late samples
early <- c("1")
late <- c("14")
earlyLate.eset <- getESet(filtered.eset,early,late)
groupA.eset = earlyLate.eset[,earlyLate.eset$`treatment_arm:ch1`=="Placebo"]
groupB.eset = earlyLate.eset[,earlyLate.eset$`treatment_arm:ch1`=="Lactoferrin"]

###########
# group A #
###########
title = "Expression of MEs in control group"
A <- which(allCodes=="Placebo")
early <- which(allDays=="1")
late <- which(allDays=="14")
inds <- intersect(A,c(early,late))
timepointLabels <- sub("14","Late",groupA.eset$`day:ch1`)
timepointLabels <- sub("1","Early",timepointLabels)

# select the MEs from group A or group B 
MEsAutomatic=net$MEs
MEsAutomatic <- MEsAutomatic[inds,] # use inds from groupA OR groupB
# optionally, only look at MEs of interest
# MEsAutomatic <- MEsAutomatic[,c("ME7","ME10","ME12")]
MEsAutomatic$Timepoint <- factor(timepointLabels) # use timepoint labels from groupA OR groupB
MEsAutomatic$Timepoint <- factor(groupA.eset$`day:ch1`) # use groupA.eset OR groupB.eset

# create boxplots
df <- melt(MEsAutomatic)
ggplot(data=df) + geom_boxplot(aes(x=Timepoint,y=value)) + 
  facet_wrap(~variable,scales = "free") + 
  ggtitle(title) + scale_fill_brewer(palette = "Accent") + 
  ylab("Value") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

# identify DE in MEs using t-tests
numModules <- length(MEsAutomatic)-1
ctt <- colttests(data.matrix(MEsAutomatic),MEsAutomatic$Timepoint,tstatOnly = FALSE)
ttests_A <- as_tibble(ctt) %>%
  mutate(M= rownames(ctt)) %>%
  mutate(p.adj = p.adjust(p.value, method='BH'))
filter(ttests_A, p.adj<0.05)

# plot p-values
title = "DE of MEs between Day 1 and Day 14 - Placebo"
barplot(ctt$p.value[1:numModules],names.arg = rownames(ctt)[1:numModules],ylim = c(0,1),main = title,xlab = "Module Eigengenes",ylab = "P-values")
abline(h=0.05,col="red")
legend(x = "topright",legend = "P = 0.05",lty = 1,col = "red")


###########
# group B #
###########

title = "Expression of MEs in lactoferrin group"
B <- which(allCodes=="Lactoferrin")
early <- which(allDays=="1")
late <- which(allDays=="14")
inds <- intersect(B,c(early,late))
timepointLabels <- sub("14","Late",groupB.eset$`day:ch1`)
timepointLabels <- sub("1","Early",timepointLabels)

# select the MEs from group A or group B 
MEsAutomatic=net$MEs
MEsAutomatic <- MEsAutomatic[inds,] # use inds from groupA OR groupB
# optionally, only look at MEs of interest
# MEsAutomatic <- MEsAutomatic[,c("ME7","ME10","ME12")]
MEsAutomatic$Timepoint <- factor(timepointLabels) # use timepoint labels from groupA OR groupB
MEsAutomatic$Timepoint <- factor(groupB.eset$`day:ch1`) # use groupA.eset OR groupB.eset

# create boxplots
df <- melt(MEsAutomatic)
ggplot(data=df) + geom_boxplot(aes(x=Timepoint,y=value)) + 
  facet_wrap(~variable,scales = "free") + 
  ggtitle(title) + scale_fill_brewer(palette = "Accent") + 
  ylab("Value") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

# identify DE in MEs using t-tests
numModules <- length(MEsAutomatic)-1
ctt <- colttests(data.matrix(MEsAutomatic),MEsAutomatic$Timepoint,tstatOnly = FALSE)
ttests_B <- as_tibble(ctt) %>%
  mutate(M= rownames(ctt)) %>%
  mutate(p.adj = p.adjust(p.value, method='BH'))
filter(ttests_B, p.adj<0.05)

# plot p-values
title = "DE of MEs between Day 1 and Day 14 - Lactoferrin"
barplot(ctt$p.value[1:numModules],names.arg = rownames(ctt)[1:numModules],ylim = c(0,1),main = title,xlab = "Module Eigengenes",ylab = "P-values")
abline(h=0.05,col="red")
legend(x = "topright",legend = "P = 0.05",lty = 1,col = "red")



####################
# Get module genes #
####################
# GET MODULE GENES
# extract genes from modules and write them to files for enrichment analysis
# get list of all gene IDs
allGenes <- colnames(datExprs)

# getModuleGenes returns a list of gene IDs from a given module 
# geneList contains all genes in the expression matrix
# geneColours contains an array of module colour assignments for each gene in the expression matrix
# moduleColour is the module of interest's colour 
# note: geneColours may contain numbers instead of colours, in which case moduleColour must also be a number
getModuleGenes <- function(geneList,geneColours,moduleColour) {
  moduleGenes <- (geneColours==moduleColour)
  moduleIDs <- geneList[moduleGenes]
  return (moduleIDs)
}


getProbeBPs <- function(module) {
  probes <- getModuleGenes(geneList=allGenes,geneColours = moduleLabelsAutomatic,moduleColour = module)
  genes <- dplyr::select(filter(annoData, ID %in% probes), Gene.Symbol)
  genes <- unlist(str_split(genes$Gene.Symbol, '///'))
  all_genes <- unlist(str_split(annoData$Gene.Symbol, '///'))
  to_enrich <- data.frame(gene=all_genes, is_candidate=as.integer(all_genes %in% genes))  
  enrich <- go_enrich(to_enrich)
  BPs <- filter(enrich$results, ontology=='biological_process' & FWER_overrep < 0.05)$node_name
  return(BPs)
}

modules_of_interest <- c(11, 23, 16, 3, 5, 6, 7, 15, 10, 8)

module_bps <- list()
for (i in modules_of_interest) {
  module_bps[[i]] <- getProbeBPs(i)
}

names(module_bps) <- as.character(modules_of_interest)










