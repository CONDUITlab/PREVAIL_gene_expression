
###############################
# Goenrichment.R              # 
# Gene set enrichment using   #
# the GOfuncR package         #
# to identify Gene Ontology   #
# biological process (BP)     #
# concepts from DEGs          #
###############################

source('setup.R')

# Enrichment of gene sets (WITH BH CORRECTION)
P_diff <- D1vD14_BH %>%
  filter(arm=='Placebo' & is_candidate==1)
L_diff <- D1vD14_BH %>%
  filter(arm=='Lactoferrin' & is_candidate==1)
# how many in each group, and intersect
# Number of DEG in placebo
length(unique(P_diff$gene))       # 2483
# Number of DEG in lactoferrin
length(unique(L_diff$gene))       # 3720
# DE in both
both <- intersect(P_diff$gene, L_diff$gene) # 1916 genes
# DE in placebo only
P_only <- setdiff(P_diff$gene, L_diff$gene) # 567 genes
# DE in lactoferrin only
L_only <- setdiff(L_diff$gene, P_diff$gene) # 1804 genes
length(both)
length(P_only)
length(L_only)
# How many up vs down?
P_up <- pull(filter(P_diff, mean_fc >= 0), gene)
L_up <- pull(filter(L_diff, mean_fc >= 0), gene)
P_down <- pull(filter(P_diff, mean_fc < 0), gene)
L_down <- pull(filter(L_diff, mean_fc < 0), gene)
both_up <- intersect(P_up, L_up) # 1354
P_only_up <- setdiff(P_up, L_up) # 243
L_only_up <- setdiff(L_up, P_up) # 875
both_down <- intersect(P_down, L_down) # 560
P_only_down <- setdiff(P_down, L_down) # 326
L_only_down <- setdiff(L_down, P_down) # 931
length(both_up)
length(both_down)
length(P_only_up)
length(P_only_down)
length(L_only_up)
length(L_only_down)

P_df <- D1vD14_BH %>% 
  filter(arm == 'Placebo', is_candidate==1) %>%
  dplyr::select(gene, is_candidate)
P_bp <- GOwrapper(P_df)

L_df <- D1vD14_BH %>%
  filter(arm == 'Lactoferrin', is_candidate==1) %>%
  dplyr::select(gene, is_candidate)
L_bp <- GOwrapper(L_df)

bp_both <- intersect(P_bp$node_name, L_bp$node_name)
P_only_bp <- setdiff(P_bp$node_name, L_bp$node_name)
L_only_bp <- setdiff(L_bp$node_name, P_bp$node_name)


