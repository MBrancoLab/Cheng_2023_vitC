library('DESeq2')


##read subfamily counts

subf.files = list.files(path='subFcounts', full.names=T)

subf = list()
for (i in 1:length(subf.files)) {
	data = read.delim(subf.files[i],as.is=T)
	if (i==1) {
		order = data$Subfamily.Family.Class
		subf[[i]] = data
	} else {
		subf[[i]] = data[match(order,data$Subfamily.Family.Class),]
	}
}
names(subf) = gsub('_subFcounts.txt','',basename(subf.files))


##make total count matrix (unique + squire-assigned fraction of multimappers)

reads = matrix(unlist(lapply(subf, function(x) x$tot_counts)),ncol=length(subf))
rownames(reads) = order


##colData

colData = DataFrame(sample=factor(rep(c('Ctrl','VitC'),4)),replicate=factor(rep(c('10','11','35B','35C'),each=2)))


##Get DE repeats

dds = DESeqDataSetFromMatrix(round(reads),colData,design=~sample+replicate)
sizeFactors(dds) = unlist(lapply(subf,function(x) x$aligned_libsize[1]))/1e6

dds = DESeq(dds)
res = results(dds,contrast=c('sample','VitC','Ctrl'))
sig = as.data.frame(subset(res,padj<0.05))
sig = sig[order(sig$log2FoldChange,decreasing=T),]

write.table(sig,'de_squire_subfamily.txt', sep='\t',quote=F,col.names=NA)


##Volcano plot

plot(res$log2FoldChange,-log10(res$padj),pch=19,cex=0.3,col='grey',
 xlab='log2 FC',ylab='-log10 p-value',las=1)
points(sig$log2FoldChange,-log10(sig$padj),pch=19,cex=0.5,col='red')
lines(c(0,0),c(-10,200),lty=2)


##MA plot (Figure 1C)

quartz(w=4,h=4)
plot(log2(res$baseMean),res$log2FoldChange,pch=19,cex=0.3,col='grey',
 xlab='log2 baseMean',ylab='log2 Fold change',las=1)
points(log2(sig$baseMean),sig$log2FoldChange,pch=19,cex=0.5,col='red')
abline(0,0,lty=2)


