
#################################
# ANALYZE BY DAY                #
# Looks at differences between  #
# lactoferrin and placebo       #
# on a particular study day     #
#                               #
# 10/24/18                      #
#################################

# Load the expression sets (from GEO)
source('setup.R')

####################
# Day comparisons  #
####################
diffs_D1 <- analyzeByDay(filtered.eset, 'day: 1', corrtype = 'BH')
diffs_D3 <- analyzeByDay(D3_eset, 'day: 3', corrtype = 'BH')
diffs_D7 <- analyzeByDay(D7_eset, 'day: 7', corrtype = 'BH')
diffs_D14 <- analyzeByDay(D14_eset, 'day: 14', corrtype = 'BH')
diffs_D21 <- analyzeByDay(D21_eset, 'day: 21', corrtype = 'BH')

sample_sizes <- c(length(unique(filtered.eset[, filtered.eset$`day:ch1`==1]$patient_id)),
                  length(unique(D3_eset$patient_id)),
                  length(unique(D7_eset$patient_id)),
                  length(unique(D14_eset$patient_id)),
                  length(unique(D21_eset$patient_id))
                  )

DEGs <- c(
  sum(diffs_D1$adj.P.Val<0.05),
  sum(diffs_D3$adj.P.Val<0.05),
  sum(diffs_D7$adj.P.Val<0.05),
  sum(diffs_D14$adj.P.Val<0.05),
  sum(diffs_D21$adj.P.Val<0.05)
)

diffs_by_day <- tibble(Day=c(1,3,7,14,21), sample_sizes, DEGs)

write_csv(diffs_by_day, 'diffs_by_day.csv')
