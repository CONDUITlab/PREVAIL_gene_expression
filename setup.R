
#################################################
# R code for analysis of gene expression        #
# data derived from the PREVAIL study           #
# D. Maslove Nov 2018                           #
#                                               #
# PREVAIL study (10.1097/CCM.0000000000003294)  #
# ###############################################

###############################
# File sequence:              #
#   1. setup.R                #
#   2. baseline_clusters.R    #
#   3. DEG_by_time.R          #
#   4. DEG_by_day.R           #
#   5. Day14.R                #
#   6. Goenrichment.R         #
#   7. WGCNA.R                #
###############################


###############################
# setup.R is used to load     #
# and clean expression data   #
# and generate ExpressionSet  #
# objects.                    #
###############################

library(GEOquery)
library(tidyverse)
library(sva)
library(COMMUNAL)
library(limma)
library(GOfuncR)
library(WGCNA)
library(cluster)

# Download dataset from GEO
# gse <- getGEO("GSE118657", GSEMatrix = TRUE)
# eset.1 <- gse[[1]]
# Uncomment below if working from saved data
load('eset1.rds')

# Exclude the healthy controls
eset.2 <- eset.1[ , eset.1$`treatment_arm:ch1`!='NA']
# correct for batch effects (sva package)
phenoData <- pData(eset.2)
batch <- phenoData$`batch:ch1`
modcombat <- model.matrix(~1, data=phenoData) # no adjustment variables, just fit an intercept term
exprsData <- Biobase::exprs(eset.2)
correctedExprs <- ComBat(dat=exprsData, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# create new corrected expression set
corrected.eset <- eset.2
exprs(corrected.eset) <- correctedExprs
# take the top 50% most variable probes (genefilter package)
# use IQR as variance func - robust to outliers 
filtered.eset <- varFilter(corrected.eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
# Recover unique patient ID
ids <- str_split(filtered.eset$title, ' ')

makeID <- function(x) {
  paste(x[1], x[2], sep='_')
}
patient_id <- unlist(map(ids, makeID))
filtered.eset$patient_id <- patient_id

# Generate Esets for each day
# Look only at downstream samples that also have a Day 1 sample
matchedEset <- function(Eset, second_day) {
  # A function to generate an Eset with matched pairs only
  # i.e. only take a Day "X" sample if that patient also
  # has a Day 1 sample
  newEset <- Eset[ , Eset$`day:ch1` == 1 | 
                       Eset$`day:ch1` == second_day]
  # Use only matched pairs
  patients <- newEset$patient_id
  matched <- table(patients)==2
  return(newEset[ , newEset$patient_id %in% names(matched)[matched]])
}

# The resulting esets below have samples for "day X" as well as "day 1"
D3_eset <- matchedEset(filtered.eset, 3)    # 102 samples
D7_eset <- matchedEset(filtered.eset, 7)    # 62 samples
D14_eset <- matchedEset(filtered.eset, 14)  # 32 samples
D21_eset <- matchedEset(filtered.eset, 21)  # 9 samples

# UNIQUE PATIENTS
all_patients <- tibble(GEO_id=filtered.eset$patient_id, arm=filtered.eset$`treatment_arm:ch1`,
                       title=filtered.eset$title)
included_pats_title <- unique(c(D3_eset$patient_id, D7_eset$patient_id, D14_eset$patient_id, D21_eset$patient_id))
included_pats <- filter(all_patients, GEO_id %in% included_pats_title) %>%
  mutate(title=as.character(title))

# get sample hash table
sample_hash <- read_csv('sample_converter.csv')
patients_used <- left_join(included_pats, sample_hash)

# Get annotation data
annoData <- fData(eset.1)
# Rename the Gene Symbol column for later reference
annoData$Gene.Symbol <- annoData$`Gene Symbol`

#########################
# Functions for single  #
# gene differential     #
# expression analysis   #
#########################

createDesign <- function(timeEset,lateDay) {  
  # createDesign forms a design matrix using the smaller expression set
  index <- as.character(timeEset$`day:ch1`==1) %>%
    str_replace('TRUE', 'early') %>%
    str_replace('FALSE', 'late')
  # create design matrix with the four groups we are interested in (A.early, A.late, etc) 
  designVec <- factor(paste(timeEset$`treatment_arm:ch1`,index,sep="."))
  designMatrix <- model.matrix(~0+designVec)
  colnames(designMatrix) <- levels(designVec)
  return (designMatrix)
}

getDiffs <- function(exprsSet,groupCode,design, corrtype) {
  # get differences in expression between groups  
  # in the Eset passed in as argument, using the design matrix
  # and contrasts passed in as arguments.
  contrastsString <- paste(groupCode,".early-",groupCode,".late",sep="")
  diffs <- runLimma(exprsSet,design,contrastsString, corrtype)
  return(diffs)
}

runLimma <- function(exprsSet,designMatrix,contrastsString, corrtype) {
  # Returns all the genes (regardless of p-value)
  # And allows you to specifiy method for P-value adjustment
  # account for correlations between samples of the same patient
  patientIDs <- exprsSet$patient_id
  corfit <- duplicateCorrelation(exprsSet,designMatrix,block=patientIDs)
  # consensus represents the correlation between measurements made on the same patient 
  print(corfit$consensus)
  # fit a linear model using the design
  fit <- lmFit(exprsSet, designMatrix, block=patientIDs,correlation=corfit$consensus)
  # create a contrast matrix 
  # a positive logFC for a gene means that it's upregulated in the early group
  cont.matrix <- makeContrasts(contrasts=contrastsString, levels=designMatrix)
  # fit a linear model with specified contrasts
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  # adjust p-values based on specified method, which controls family-wide error rate
  # return all genes
  diffs <- topTable(fit2, adjust=corrtype, number=100000)
  return(diffs)
}

firstOnly <- function(gene) {
  # Where a probe represents many genes
  # take the first gene only
  unlist(str_split(gene, ' /// '))[1]
}

timepointComparison <- function(fullEset, lateDay, corrtype) {
  
  # Select the samples with Day 1 and comparator
  time.eset <- fullEset[ , fullEset$`day:ch1`==1 | 
                           fullEset$`day:ch1`==lateDay]
  # Keep only matched samples
  # (i.e. the patient has samples at both time points)
  patients <- time.eset$patient_id
  matched <- table(patients)==2
  time.eset.matched <- time.eset[ , time.eset$patient_id %in% names(matched)[matched]]
  # create design matrix
  design <- createDesign(time.eset.matched, lateDay)
  # Run the limma comparison
  diffs_lacto <- getDiffs(time.eset.matched,groupCode="Lactoferrin",design, corrtype)
  diffs_placebo <- getDiffs(time.eset.matched,groupCode="Placebo",design, corrtype)
  # For multiply annotated genes,
  # take the first one 
  diffs_lacto$gene <- sapply(diffs_lacto$Gene.Symbol, firstOnly)
  diffs_placebo$gene <- sapply(diffs_placebo$Gene.Symbol, firstOnly)
  # collapse values to mins and means
  lacto_tibble <- as.tibble(diffs_lacto) %>%
    dplyr::select(gene, logFC, t, P.Value, adj.P.Val) %>%
    group_by(gene) %>%
    summarise(mean_fc = mean(logFC), mean_t=mean(t), mean_p=mean(P.Value), mean.adj.p=mean(adj.P.Val)) %>%
    dplyr::filter(gene != "") %>%
    dplyr::filter(gene != "---") %>%
    dplyr::filter(is.na(as.numeric(gene)))
  lacto_tibble$arm <- rep("Lactoferrin", nrow(lacto_tibble))
  placebo_tibble <- as.tibble(diffs_placebo) %>%
    dplyr::select(gene, logFC, t, P.Value, adj.P.Val) %>%
    group_by(gene) %>%
    summarise(mean_fc = mean(logFC), mean_t=mean(t), mean_p=mean(P.Value), mean.adj.p=mean(adj.P.Val)) %>%
    dplyr::filter(gene != "") %>%
    dplyr::filter(gene != "---") %>%
    dplyr::filter(is.na(as.numeric(gene)))
  placebo_tibble$arm <- rep("Placebo", nrow(placebo_tibble))
  # Join into master table
  full_T <- rbind(lacto_tibble, placebo_tibble)
  full_T$is_candidate <- as.numeric(full_T$mean.adj.p<0.05)
  return(full_T)
}


analyzeByDay <- function(filteredESet,dayString, corrtype) {
  # analyzeByDay creates an expression set with only samples from that day
  # and uses limma to look for differences in expression between 
  # group A and group B on that day.
  subsetted.eset <- filteredESet[ , filteredESet$`characteristics_ch1.3`==dayString]
  # create vectors indicating whether a sample belongs to group A or B
  groupA <- as.integer(subsetted.eset$`treatment_arm:ch1`=='Placebo')
  groupB <- as.integer(subsetted.eset$`treatment_arm:ch1`=='Lactoferrin') 
  # create design matrix with the two groups we are interested in 
  design <- cbind(groupA, groupB)
  colnames(design) <- c('groupA', 'groupB')
  AvsB <- "groupA-groupB"
  # fit a linear model using the design
  fit <- lmFit(subsetted.eset, design)
  # create a contrast matrix 
  cont.matrix <- makeContrasts(contrasts=AvsB, levels=design)
  # fit a linear model with specified contrasts
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  # adjust p-values based on Holm's method, which controls family-wide error rate
  diffs <- topTable(fit2, adjust=corrtype, number=100000)
  return(diffs)
}

GOwrapper <- function(df) {
  # wrapper function for GO enrichment with GOfuncR
  res <- go_enrich(data.frame(df), n_randsets = 100)
  stats <- res[[1]]
  return(filter(stats, ontology == 'biological_process', FWER_overrep < 0.05))
}

