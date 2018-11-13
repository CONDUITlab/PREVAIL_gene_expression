# PREVAIL_gene_expression
R code for the analysis of gene expression data from the PREVAIL study (GSE118657)

Requires specific R packages including:
-GEOquery
-tidyverse
-sva
-COMMUNAL (archived in CRAN - may need to install from source)
-limma
-GOfuncR
-WGCNA

STEPS TO FOLLOW:
1. Run setup.R to download data from GEO and generate expressionSet
2. baseline_clusters.R does the clustering on baseline samples and returns clustering metrics/plots
3. DEG_by_time.R generates comparisons with baseline at each time point, for each study arm
4. DEG_by_day.R generates comparisons between groups at each time point
5. Goenrichment.R looks at the number of differentially expressed genes (DEGs) between groups, and does functional enrichment for GO biological pathways
6. WGCNA.R does the network-based analysis
