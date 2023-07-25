library('DESeq2')

##Data

reads = read.table('gene_rawCounts.txt', header = T, row.names = 1)
colData = DataFrame(sample=factor(rep(c('Ctrl','VitC'),each=4)),replicate=factor(rep(c('10','11','35B','35C'),2)))


##Expression values (VST)

##dds = DESeqDataSetFromMatrix(reads,colData,design=~sample+replicate)
##vsd = varianceStabilizingTransformation(dds)
##expr = assay(vsd)
expr = read.table('gene_expression_vsd.txt',as.is=T,row.names=1)


##Get DE genes

dds = DESeqDataSetFromMatrix(reads,colData,design=~sample+replicate)
dds = DESeq(dds)
res = results(dds,contrast=c('sample','VitC','Ctrl'))
sig = as.data.frame(subset(res,padj<0.05 & abs(log2FoldChange)>log2(1.5)))
sig = sig[order(sig$log2FoldChange,decreasing=T),]

write.table(sig, 'de_genes_vitC.txt', sep='\t',quote=F,col.names=NA)


##Volcano plot

plot(res$log2FoldChange,-log10(res$padj),pch=19,cex=0.3,col='grey',
 xlab='log2 FC',ylab='-log10 p-value',las=1)
is.sig = res$padj<0.05 & abs(res$log2FoldChange)>log2(1.5)
points(res$log2FoldChange[is.sig],-log10(res$padj[is.sig]),pch=19,cex=0.5,col='red')
lines(c(0,0),c(-10,200),lty=2)


##MA plot (Supplementary Figure 1C)

plot(log2(res$baseMean),res$log2FoldChange,pch=19,cex=0.3,col='grey',
 xlab='log2 expression',ylab='log2 Fold change',las=1)
is.sig = res$padj<0.05 & abs(res$log2FoldChange)>log2(1.5)
points(log2(res$baseMean[is.sig]),res$log2FoldChange[is.sig],pch=19,cex=0.5,col='red')
abline(0,0,lty=2)


