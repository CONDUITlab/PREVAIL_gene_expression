
#################################
# DEG_by_time.R                 #
#                               #
# Looks at differences between  #
# Day 1 and Day X               # 
# in lactoferrin vs placebo     #
#                               #
# 10/24/18                      #
#################################

# Load the expression sets (from GEO)
source('setup.R')

# Comparisons for each time span
# (eg. D1 vs D3, D1 vs D7, etc.)
D1vD3 <- timepointComparison(filtered.eset, 3, corrtype = 'BH')
D1vD3$timepoint <- rep(3, nrow(D1vD3))
D1vD7 <- timepointComparison(filtered.eset, 7, corrtype = 'BH')
D1vD7$timepoint <- rep(7, nrow(D1vD3))
D1vD14 <- timepointComparison(filtered.eset, 14, corrtype = 'BH')
D1vD14$timepoint <- rep(14, nrow(D1vD3))
D1vD21 <- timepointComparison(filtered.eset, 21, corrtype = 'BH')
D1vD21$timepoint <- rep(21, nrow(D1vD3))
# combine these
DEG_by_time <- do.call(rbind, list(D1vD3, D1vD7, D1vD14, D1vD21))

# Number of DEGs per time point, by arm
timepoint_summary <- DEG_by_time %>%
  group_by(arm, timepoint) %>%
  summarise(n_DEG=sum(is_candidate==1))

# create figure
p <- ggplot(timepoint_summary, aes(x=timepoint, y=n_DEG, group=arm, colour=arm)) +
  geom_line(size=2) + scale_colour_manual(values = c("navyblue", "darkorange2")) +
  theme_bw() + ylab('Number of differentially expressed genes') +
  theme(axis.title.x = element_text(face='bold', size=16), axis.title.y = element_text(face='bold', size=12)) +
  theme(axis.text = element_text(size=12, face = 'bold')) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.85,0.85)) +
  theme(legend.text = element_text(size=12, face='bold')) +
  scale_x_continuous(breaks = c(3,7,14,21)) + xlab('Day')
