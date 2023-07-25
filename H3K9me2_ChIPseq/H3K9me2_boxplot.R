
##quantitation data from Seqmonk

l1 = read.delim('young_5kb_L1s.txt')
ran = read.delim('random_5kb_regions.txt')
data = rbind(l1, ran)


##normalisation as per Ebata et al.

data$H3K9me2_VitC.unique = data$H3K9me2_VitC.unique * 0.57


##boxplot ratio (Figure 3B)

ratio = data$H3K9me2_VitC.unique / data$H3K9me2_Untreated.unique

boxplot(ratio ~ data$Probe, outline=F, lty=1, las=1,
	xlab='', ylab='Fold change', col='orange')
abline(h=1, lty=2)
