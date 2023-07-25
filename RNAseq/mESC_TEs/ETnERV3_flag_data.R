
##Data

flag.files = list.files(path='flagSubset', full.names=T)

flag = list()
for (i in 1:length(flag.files)) flag[[i]] = read.delim(flag.files[i],as.is=T)
names(flag) = gsub('_flagSubset.txt.gz','',basename(flag.files))


##Keep only sense transcripts

sense = lapply(flag, function(x) x[x$tx_strand==x$TE_strand,])


##Select LTR families

te = lapply(sense, function(x) x[x$TE_name=='ETnERV3-int:ERVK:LTR' |
 x$TE_name=='RLTR13B2:ERVK:LTR' |
 x$TE_name=='RLTR13A3:ERVK:LTR' ,])


##Match elements across samples

all.te = unique(unlist(lapply(te,function(x) x$TE_ID)))
te.match = lapply(te, function(x) x[match(all.te,x$TE_ID),])

for (i in 1:length(te.match)) {
	te.match[[i]]$TE_ID = all.te
	na.lines = is.na(te.match[[i]]$fpkm)
	te.match[[i]]$fpkm[na.lines] = te.match[[i]]$uniq_counts[na.lines] = te.match[[i]]$tot_counts[na.lines] = te.match[[i]]$tot_reads[na.lines] = 0
	te.match[[i]]$alignedsize[na.lines] = te.match[[i]]$alignedsize[!na.lines][1]
}


##Get normalised counts

counts = lapply(te.match, function(x) log2(x$tot_counts/x$alignedsize*1e6 + 0.05))


##Calculate mean fold change

fc1 = counts$'VC10-VC-54801786' - counts$'VC10-Ctrl-54801785'
fc2 = counts$'VC11-VC-54807786' - counts$'VC12-Ctrl-54809768'
fc3 = counts$'VC35B-VC-54809771' - counts$'VC35B-Ctrl-54808775'
fc4 = counts$'VC35C-VC-54809772' - counts$'VC35C-Ctrl-54799820'

mean.counts = rowMeans(matrix(unlist(counts),ncol=length(counts)))
mean.fc = rowMeans(cbind(fc1,fc2,fc3,fc4))


##Subgroups

etn = grep('ETnERV3',all.te)
b2 = grep('RLTR13B2',all.te)
a3 = grep('RLTR13A3',all.te)


##Scatter plot

plot(mean.counts[etn],mean.fc[etn],pch=19,cex=0.5,las=1,
 col='grey',xlab='log2 RPM',ylab='log2 FC',main='ETnERV3-int',
 ylim=c(-2,4.5))
points(mean.counts[b2],mean.fc[b2],pch=19,cex=0.5,col='red')
points(mean.counts[a3],mean.fc[a3],pch=19,cex=0.5,col='blue')
abline(h=0,lty=2)


##Select elements

up = mean.counts>-2 & mean.fc>2
up.id = strsplit(all.te[up], split='\\|')

up.te = data.frame(chr=unlist(lapply(up.id, function(x) x[1])),
			start=unlist(lapply(up.id, function(x) x[2])),
			end=unlist(lapply(up.id, function(x) x[3])),
			strand=unlist(lapply(up.id, function(x) x[6])),
			name=unlist(lapply(up.id, function(x) x[4])),
			log2rpm = mean.counts[up],
			foldChange = mean.fc[up])

write.table(up.te, 'up_ETnERV3.txt', sep='\t', quote=F, row.names=F)
#These 22 elements merge into 10 loci, most of which have a proviral arrangement

