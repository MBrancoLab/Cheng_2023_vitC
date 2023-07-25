
##Data

flag.files = list.files(path='flagSubset', full.names=T)

flag = list()
for (i in 1:length(flag.files)) flag[[i]] = read.delim(flag.files[i],as.is=T)
names(flag) = gsub('_flagSubset.txt.gz','',basename(flag.files))


##Keep only sense transcripts

sense = lapply(flag, function(x) x[x$tx_strand==x$TE_strand,])


##Select L1 families

te = lapply(sense, function(x) x[x$TE_name=='L1Md_T:L1:LINE' | x$TE_name=='L1Md_A:L1:LINE' | x$TE_name=='L1Md_Gf:L1:LINE',])


##Merge elements across samples

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


##Identify ITLs

tx.type = matrix(unlist(lapply(te.match,function(x) x$transcript_type)),ncol=length(te.match))
itl = rowSums(tx.type=='pot_unit_ITL',na.rm=T)>=1


##Full-length elements

te.split = strsplit(all.te,split='\\|')
te.length = unlist(lapply(te.split,function(x) as.numeric(x[3])-as.numeric(x[2])))
fl = te.length>5000


##Confidence values

conf = matrix(unlist(lapply(te.match,function(x) x$avg_conf)),ncol=length(te.match))
mean.conf = rowMeans(conf,na.rm=T)

#hist(mean.conf)
high.conf = mean.conf>50


##Scatter plots

#tf = grepl('L1Md_T',all.te)
#a = grepl('L1Md_A',all.te)
#gf = grepl('L1Md_Gf',all.te)

#plot(mean.counts[itl&fl&tf],mean.fc[itl&fl&tf],pch=19,cex=0.5,las=1,
# col='grey',xlab='log2 RPM',ylab='log2 FC',main='L1Md_T',ylim=c(-4,4))
#points(mean.counts[itl&fl&high.conf&tf],mean.fc[itl&fl&high.conf&tf],pch=19,cex=0.5,col='black')
#abline(h=0,lty=2)

#plot(mean.counts[itl&fl&a],mean.fc[itl&fl&a],pch=19,cex=0.5,las=1,
# col='grey',xlab='log2 RPM',ylab='log2 FC',main='L1Md_A',ylim=c(-4,4))
#points(mean.counts[itl&fl&high.conf&a],mean.fc[itl&fl&high.conf&a],pch=19,cex=0.5,col='black')
#abline(h=0,lty=2)

#plot(mean.counts[itl&fl&gf],mean.fc[itl&fl&gf],pch=19,cex=0.5,las=1,
# col='grey',xlab='log2 RPM',ylab='log2 FC',main='L1Md_Gf',ylim=c(-4,4))
#points(mean.counts[itl&fl&high.conf&gf],mean.fc[itl&fl&high.conf&gf],pch=19,cex=0.5,col='black')
#abline(h=0,lty=2)


##boxplot (Figure 1E)

quartz(w=4.5, h=3)
par(mar=c(4,4,2,5))
boxplot(mean.fc[itl&fl&a & mean.counts>-2],
	mean.fc[itl&fl&tf & mean.counts>-2],
	lty=1, outline=FALSE,
	col='orange', las=1,
	xlab='log2 fold change',
	boxwex=0.6,
	horizontal=TRUE)

y = rep(1:2, c(sum(itl&fl&high.conf&a & mean.counts>-2),
	sum(itl&fl&high.conf&tf & mean.counts>-2)))
points(c(mean.fc[itl&fl&a&high.conf & mean.counts>-2],
	mean.fc[itl&fl&tf&high.conf & mean.counts>-2]),
	y + rnorm(length(y),sd=0.05),
	pch=17, cex=0.7)

abline(v=0,lty=2)

p1 = t.test(mean.fc[itl&fl&a & mean.counts>-2])
p2 = t.test(mean.fc[itl&fl&tf & mean.counts>-2])
p.adjust(c(p1$p.value, p2$p.value), method='BH')
