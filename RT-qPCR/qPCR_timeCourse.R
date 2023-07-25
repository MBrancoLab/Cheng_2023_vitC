library(gplots)

x2i = read.csv('2i_time_pt/2i_time_pt_export_all.csv')
ser = read.csv('Serum_time_pt/Serum_time_pt_export_all.csv')


get.val = function(data, te) {
	sub = data[data$all_data2 %in% te,]
	sub$group = factor(as.character(sub$all_data2), levels=te)
	
	t0 = sub[grep('Ctrl', sub$all_data1),]
	t24 = sub[grep('24hrs', sub$all_data1),]
	t48 = sub[grep('48hrs', sub$all_data1),]
	t1wk = sub[grep('1week', sub$all_data1),]
	
	mat = cbind(
		tapply(t0$all_data3, t0$group, median),
		tapply(t24$all_data3, t24$group, median),
		tapply(t48$all_data3, t48$group, median),
		tapply(t1wk$all_data3, t1wk$group, median))
	colnames(mat) = c('Ctrl', '24h','48h','1wk')
	
	return(mat)
}


te = c('L1 A', 'L1 TF', 'L1 GF', 'ORF1a', 'ORF1b', 'ORF2', 'ORF2b',
	'IAP LTR1', 'IAP LTR2', 'IAP LTR3', 'IAP GAG',
	'MuSD GAG', 'MuLV', 'MuLV GAG', 'MERV', 'MERV GAG')
mat.2i = get.val(x2i, te)
mat.ser = get.val(ser, te)


##Figure 1B

heatmap.2(cbind(mat.2i, mat.ser),
	scale='none', Colv=NA, Rowv=NA, dendrogram='none',
	col=hcl.colors(100),
	trace='none', density.info='none',
	colsep=4)
