library(tidyverse)


##read data

d13 = read_tsv('L1seq_13days_TEbreak.txt.gz')
d41.r4 = read_tsv('L1seq_41days_TEbreak_4reads.txt.gz')
d41.r1 = read_tsv('L1seq_41days_TEbreak_1read.txt.gz')


##find insertions unique to one clone

d13.support = gsub('NA,|,NA', '', paste(d13$Sample_support_5p, d13$Sample_support_3p, sep=','))
d13.unique = unlist(lapply(strsplit(d13.support, ','), function(x) length(unique(substr(x,1,1))))) == 1
write_tsv(d13[d13.unique,], 'clonal_insertions_d13.txt')

d41.support = gsub('NA,|,NA', '', paste(d41.r4$Sample_support_5p, d41.r4$Sample_support_3p, sep=','))
d41.unique = unlist(lapply(strsplit(d41.support, ','), function(x) length(unique(substr(x,1,7))))) == 1
write_tsv(d41[d41.unique,], 'clonal_insertions_d41.txt')


##subset vitC- or ctrl-only insertions

l1 = filter(d41.r1, Sample_count==1, !grepl('Day_0', Sample_support_5p), !grepl('Day_0', Sample_support_3p)) %>%
		mutate(group = case_when(grepl('No_VitC',Sample_support_5p) | grepl('No_VitC',Sample_support_3p) ~ 'Ctrl',
								grepl('YesVitC',Sample_support_5p) | grepl('YesVitC',Sample_support_3p) ~ 'VitC'))

l4 = filter(d41.r4, Sample_count==1, !grepl('Day_0', Sample_support_5p), !grepl('Day_0', Sample_support_3p)) %>%
		mutate(group = case_when(grepl('No_VitC',Sample_support_5p) | grepl('No_VitC',Sample_support_3p) ~ 'Ctrl',
								grepl('YesVitC',Sample_support_5p) | grepl('YesVitC',Sample_support_3p) ~ 'VitC'))


##summarise

group_by(l1, Superfamily, group) %>% summarise(count = length(group))
group_by(l4, Superfamily, group) %>% summarise(count = length(group))




