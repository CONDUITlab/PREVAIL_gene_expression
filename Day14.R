
##############################
# Day14.R                    #
# Looks at differences in    #
# expression between groups  #
# on day 14.                 # 
#  11/14/18                  #
##############################

# Load the expression sets (from GEO)
source('setup.R')

# Generate comparisons with p-value adjustments using
# two methods (holm method to get threshold for plot)
D1vD14_BH <- timepointComparison(filtered.eset, 14, corrtype='BH')
D1vD14_holm <- timepointComparison(filtered.eset, 14, corrtype = 'holm')

# Cutoffs
# Uncomment below to re-calculate
BH_cutoff <- filter(D1vD14_BH, mean.adj.p < 0.05)
BH_thresh <- max(BH_cutoff$mean_p)
holm_cutoff <- filter(D1vD14_holm, mean.adj.p < 0.05)
holm_thresh <- max(holm_cutoff$mean_p)

# Volcano plot
p <- ggplot(D1vD14_BH, aes(x=mean_fc, y=-log10(mean.adj.p), colour=arm)) +
  geom_point(size=2, alpha=0.7) +
  scale_colour_manual(values = c('navyblue', 'darkorange2')) +
  facet_wrap(~arm) +
  theme_bw() + xlab('log fold-change') + ylab('-log(adjusted P-value)') +
  theme(axis.title = element_text(face='bold', size=16), axis.text = element_text(size=12, face = 'bold')) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_hline(yintercept = -log10(BH_thresh), linetype='twodash', colour='grey50') +
  geom_hline(yintercept = -log10(holm_thresh), linetype='twodash', colour='grey50') +
  annotate("text", label = "paste(BH, ' ', italic(P), \"= 0.05\")", x= -2.65, y=1.7, size=4, parse=TRUE) +
  annotate("text", label = "paste(Holm, ' ', italic(P), \"= 0.05\")", x= -2.55, y=5.8, size=4, parse=TRUE) +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-3, 3)) +
  theme(strip.text = element_text(face='bold', size=14))

