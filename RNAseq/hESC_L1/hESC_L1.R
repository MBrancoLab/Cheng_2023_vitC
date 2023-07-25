library(tidyverse)


##TeXP data (Supplementary Figure 4C)

texp.files = list.files(pattern='TExP')

texp = tibble()
for (f in texp.files) {
	sub = read_delim(f, delim=' ', skip=1, col_names=c('L1','rpkm')) %>%
		add_column(group = strsplit(f, split='[_.]')[[1]][2])
	texp = rbind(texp, sub)
}

filter(texp, grepl('L1HS', L1) | grepl('L1PA2', L1) | grepl('L1PA3', L1)) %>%
ggplot(aes(x=L1, y=rpkm, fill=group)) +
	geom_bar(stat='identity', position=position_dodge()) +
	theme_classic()



##TEtranscripts data (Supplementary Figure 4B)

tetr.files = list.files(pattern='TEtranscripts')

tetr = tibble()
for (f in tetr.files) {
	sub = read_tsv(f, skip = 1, col_names=c('feature','count')) %>%
		add_column(group = strsplit(f, split='[_.]')[[1]][2]) %>%
		mutate(rpm = count/sum(count)*1e6)
	tetr = rbind(tetr, sub)
}

filter(tetr, grepl('L1HS', feature) | grepl('L1PA2', feature) | grepl('L1PA3', feature)) %>%
ggplot(aes(x=feature, y=rpm, fill=group)) +
	geom_bar(stat='identity', position=position_dodge()) +
	theme_classic()



##SQuIRE subfamily data (Figure 6D)

subf.files = list.files(pattern='subFcounts')

subf = tibble()
for (f in subf.files) {
	sub = read_tsv(f) %>%
		mutate(group = strsplit(f, split='[_]')[[1]][2])
	subf = rbind(subf, sub)
}
colnames(subf)[3] = 'te'

filter(subf, grepl('L1HS', te) | grepl('L1PA2', te) | grepl('L1PA3', te)) %>%
ggplot(aes(x=te, y=fpkm, fill=group)) +
	geom_bar(stat='identity', position=position_dodge()) +
	theme_classic()


##SQuIRE TE data (Figure 6E)

te.files = list.files(pattern='TEcounts')

te = tibble()
for (f in te.files) {
	sub = read_tsv(f) %>%
		mutate(group = strsplit(f, split='[_]')[[1]][2])
	te = rbind(te, sub)
}

l1 = filter(te, grepl('L1HS',TE_name) | grepl('L1PA2',TE_name) | grepl('L1PA3',TE_name),
			TE_stop-TE_start>5000, tx_strand==TE_strand) %>%
	select(TE_ID, TE_chr, TE_start, TE_stop, alignedsize, tot_counts, score, group) %>%
	mutate(rpm = tot_counts/alignedsize*1e6)

av.score = group_by(l1, TE_ID) %>% summarise(score = mean(score))
	
l1.wide = select(l1, TE_ID:TE_stop, group, rpm) %>%
	pivot_wider(names_from=group, values_from=rpm) %>%
	filter(!is.na(ctrl), !is.na(vitC)) %>%
	inner_join(av.score)

#annotate with KDM4A peaks
#kdm = read_tsv('ENCFF021QGZ.bed', col_names=NA)
#has.kdm = logical(nrow(l1.wide))
#for (i in 1:nrow(l1.wide)) {
#	sub = filter(kdm, X1==l1.wide$TE_chr[i],
#		(X3>l1.wide$TE_start[i] & X3<l1.wide$TE_stop[i]) | (X2>l1.wide$TE_start[i] & X2<l1.wide$TE_stop[i]))
#	if (nrow(sub)>0) has.kdm[i] = TRUE
#} #only 3 L1s with KDM4 peaks

#ggplot(l1.wide, aes(x=ctrl, y=vitC)) +
#	geom_point(col='grey') +
#	geom_point(data=filter(l1.wide,score>80), aes(x=ctrl,y=vitC), col='red') +
#	geom_abline(slope=1, intercept=0, linetype='dashed') +
#	theme_classic()

