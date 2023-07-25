
##plot function

paper.plot = function(data, te.list, g.list, mark, ylab='ChIP signal relative to untreated',
				col=c('orange','darkred')[1:length(g.list)], y.lim=NULL,
				width=4, height=4) {

	sel = data$target %in% te.list & data$cells %in% g.list & data$mark==mark
	val = data$norm2[sel]
	te = factor(data$target[sel], levels=te.list)
	gr = factor(data$cells[sel], levels=g.list)
	
	vlist = list()
	for (i in 1:nlevels(te)) {
		for (j in 1:nlevels(gr)) {
			vlist[[i*nlevels(gr)-nlevels(gr)+j]] = val[te==levels(te)[i] & gr==levels(gr)[j]]
		}
	}
	
	if (nlevels(gr)>1) {
		avg = matrix(unlist(lapply(vlist,mean)), nrow=nlevels(gr))
	} else {
		avg = unlist(lapply(vlist,mean))
	}
	sd = unlist(lapply(vlist,sd))
	if (length(y.lim)==0) y.lim = c(0,max(val)*1.05)
	
	quartz(w=width, h=height)
	par(mar=c(7,4,2,2))
	h = barplot(avg,
		beside=TRUE,
		ylim=y.lim,
		las=2, ylab=ylab,
		col=col,
		names.arg=te.list)
	points(rep(h,unlist(lapply(vlist,length))) + rnorm(length(val), sd=0.05),
		unlist(vlist), pch=19,cex=0.4)
	for (i in 1:length(h)) {
		lines(c(h[i],h[i]), c(avg[i]-sd[i],avg[i]+sd[i]))
	}
	abline(h=1, lty=2)
	
	return(data.frame(gr,te,val))
}


##KDM4 DKD (Figure 3E)

dkd = read.delim('DKD_ChIP.txt')

me3 = paper.plot(dkd, te.list=c('L1Tf','L1A','L1Gf'),
				g.list=c('SCR vitC','DKD vitC'), mark='H3K9me3')
me2 = paper.plot(dkd, te.list=c('L1Tf','L1A','L1Gf'),
				g.list=c('SCR vitC','DKD vitC'), mark='H3K9me2')


summary(aov(val~gr, data=me3))
summary(aov(val~gr, data=me2))


