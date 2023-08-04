library(dplyr)


##RPM counts from Seqmonk

l1 = read.delim('young_5kb_L1s.txt') %>%
	mutate(ctrl=H3K9me2_Untreated.unique/input_DNA_Untreated.unique,
			vitc=H3K9me2_VitC.unique/input_DNA_VitC.unique)
down = read.delim('downstream_5kb.txt') %>%
	mutate(ctrl=H3K9me2_Untreated.unique/input_DNA_Untreated.unique,
			vitc=H3K9me2_VitC.unique/input_DNA_VitC.unique)


##L1/downstream ratio

minr = l1$input_DNA_Untreated.unique>0.1 &
		l1$input_DNA_VitC.unique>0.1 &
		down$input_DNA_Untreated.unique>0.1 &
		down$input_DNA_VitC.unique>0.1

ratio = c(l1$ctrl[minr]/down$ctrl[minr],
		l1$vitc[minr]/down$vitc[minr])
subf = rep(l1$Probe[minr], 2)
cond = rep(c('ctrl','vitC'), each=sum(minr))


##boxplot ratio (Figure 3B)

boxplot(ratio ~ cond + subf, outline=F, lty=1, las=1,
	xlab='', ylab='L1 to downstream ratio', col=c('grey','orange'))

