

##plot function

paper.plot = function(data, te.list, g.list, ylab='Expression relative to Ctrl',
				col=c('orange','grey','darkred')[1:length(g.list)], y.lim=NULL,
				width=5, height=4) {

	group = gsub('[[:digit:]]+$','',data$all_data1)
	sel = data$all_data2 %in% te.list & group %in% g.list
	val = data$all_data3[sel]
	te = factor(data$all_data2[sel], levels=te.list)
	gr = factor(group[sel], levels=g.list)
	
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



##24h TEs (Figure 1A)

x24 = read.csv('2i_24hrs_VC/2i_24hrs_VC_export_all.csv')

x24.df = paper.plot(data=x24,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a', 'ORF1b', 'ORF2', 'ORF2b', 
				'IAP LTR1', 'IAP LTR2', 'IAP LTR3','IAP GAG', 
				'MuLV', 'MuLV GAG', 'EtnI', 'MuSD','MuSD GAG', 
				'MERV', 'MERV GAG'),
		g.list='2i-24hrs-VC',
		width=5.5, height=4.5,
		ylab='Fold change',
		col=rep(c('orange','grey'),c(7,11)))

p = tapply(x24.df$val, x24.df$te, function(x) t.test(x, mu=1)$p.value)
p.adj = p.adjust(p, method='BH')


##aKG (Supplementary Figure 1A)

akg = read.csv('Serum_aKG/Serum_aKG_export_all.csv', as.is=T)
akg = akg[akg$all_data3>0.1,]

akg.df = paper.plot(data=akg,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a'),
		g.list=c('ser-24hrs-VC','ser-Ctrl-aKG','ser-VC-aKG'),
		y.lim=c(0,7.5))

model = aov(akg.df$val ~ akg.df$te * akg.df$gr)
summary(model)
TukeyHSD(model)


##succinate (Supplementary Figure 1B)

succ = read.csv('2i_Succ/2i_Succ_export_all.csv', as.is=T)

succ.df = paper.plot(data=succ,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a'),
		g.list=c('2i-24hrs-VC','2i-Ctrl-Succ','2i-VC-Succ'),
		y.lim=c(0,7.5))

model = aov(succ.df$val ~ succ.df$te * succ.df$gr)
summary(model)


##Tet TKO L1 (Figure 2B)

tet = read.csv('Tet_TKO/Tet_TKO_export_all.csv', as.is=T)

tet.l1 = paper.plot(data=tet,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a', 'ORF1b', 'ORF2', 'ORF2b'),
		g.list=c('2iTTWT-VC','2iTTKO-Ctrl','2iTTKO-VC'),
		ylab='Expression relative to WT')

model = aov(tet.l1$val ~ tet.l1$te * tet.l1$gr)
summary(model)
TukeyHSD(model)


##Tet TKO ERVs (Figure 2C)

tet.erv = paper.plot(data=tet,	
		te.list=c('ETnERV3', 'RLTR13B2', 'RLTR10D2', 'MuLV', 'MuLV GAG'),
		g.list=c('2iTTWT-VC','2iTTKO-Ctrl','2iTTKO-VC'),
		width=4,
		ylab='Expression relative to WT')

model = aov(tet.erv$val ~ tet.erv$te * tet.erv$gr)
summary(model)
TukeyHSD(model)


##DMOG (Supplementary Figure 2A)

dmog = read.csv('2i_DMOG/2i_DMOG_export_all.csv', as.is=T)

dmog.df = paper.plot(data= dmog,	
		te.list=c('ORF1a', 'ORF1b'),
		g.list=c('2i-VC-DMSO','2i-Ctrl-DMOG','2i-VC-DMOG'),
		width=3,
		ylab='Expression relative to DMSO',
		y.lim=c(0,4))

model = aov(dmog.df$val ~ dmog.df$te * dmog.df$gr)
summary(model)


##24h TETs (Supplementary Figure 2B)

x24.tet = paper.plot(data=x24,	
		te.list=c('mTET1', 'mTET2'),
		g.list='2i-24hrs-VC',
		width=2,
		ylab='Expression relative to control')

p = tapply(x24.tet$val, x24.tet$te, function(x) t.test(x, mu=1)$p.value)
p.adj = p.adjust(p, method='BH')


##Tet TKO TETs (Supplementary Figure 2C)

tet.tet = paper.plot(data=tet,	
		te.list=c('mTET1','mTET2'),
		g.list=c('2iTTWT-VC','2iTTKO-Ctrl','2iTTKO-VC'),
		width=3,
		ylab='Expression relative to WT')

model = aov(tet.tet$val ~ tet.tet$te * tet.tet$gr)
summary(model)
TukeyHSD(model)


##Tet TKO IAPs (Supplementary Figure 2D)

tet.iap = paper.plot(data=tet,	
		te.list=c('IAP LTR1', 'IAP LTR2', 'IAP LTR3', 'IAP GAG'),
		g.list=c('2iTTWT-VC','2iTTKO-Ctrl','2iTTKO-VC'),
		y.lim=c(0,4),
		width=4,
		ylab='Expression relative to WT')

summary(aov(tet.iap$val ~ tet.iap$te * tet.iap$gr))


##Dnmt TKO (Supplementary Figure 2E)

dnmt = read.csv('Dnmt_TKO/Dnmt_TKO_export_all.csv', as.is=T)

dnmt.df = paper.plot(data=dnmt,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a', 'ORF2a','ETnERV3'),
		g.list=c('2i-DTKO-Ctrl','2i-DTKO-VC'),
		width=4,
		ylab='Fold change',
		col=c('grey','darkred'))

model = aov(dnmt.df$val ~ dnmt.df$te * dnmt.df$gr)
summary(model)
TukeyHSD(model)


##Kdm3 DKD (Figure 3A)

kdm3 = read.csv('Kdm3_DKD/Kdm3_DKD_export_all.csv', as.is=T)

kdm3.df = paper.plot(data=kdm3,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a'),
		g.list=c('shScr-VC','shKdm3a/b-Ctrl','shKdm3a/b-VC'),
		ylab='Expression relative to shScr')

model = aov(kdm3.df$val ~ kdm3.df$te * kdm3.df$gr)
summary(model)
TukeyHSD(model)


##Kdm3 DKD genes (Supplementary Figure 3A)

kdm3.gene = paper.plot(data=kdm3,	
		te.list=c('Kdm3a', 'Kdm3b', 'Oct4', 'Sox2', 'Nanog', 'Dazl'),
		g.list=c('shScr-VC','shKdm3a/b-Ctrl','shKdm3a/b-VC'),
		ylab='Expression relative to shScr')

p = tapply(kdm3.gene$val, factor(paste(kdm3.gene$te,kdm3.gene$gr)), function(x) t.test(x, mu=1)$p.value)
p.adj = p.adjust(p, method='BH')


##Kdm4 DKD (Figure 3D)

kdm4 = read.csv('Kdm4_DKD/Kdm4_DKD_export_all.csv', as.is=T)

kdm4.df = paper.plot(data=kdm4,	
		te.list=c('L1 TF', 'L1 A', 'L1 GF', 'ORF1a'),
		g.list=c('shScr-VC','shKdm4a/c-Ctrl','shKdm4a/c-VC'),
		ylab='Expression relative to shScr')

model = aov(kdm4.df$val ~ kdm4.df$te * kdm4.df$gr)
summary(model)
TukeyHSD(model)


##Kdm4 DKD genes (Supplementary Figure 3C)

kdm4.gene = paper.plot(data=kdm4,	
		te.list=c('Kdm4a', 'Kdm4c', 'Oct4', 'Sox2', 'Nanog', 'Dazl'),
		g.list=c('shScr-VC','shKdm4a/c-Ctrl','shKdm4a/c-VC'),
		ylab='Expression relative to shScr')

kdm4.gene = kdm4.gene[kdm4.gene$te!='Sox2',] #Sox2 has only one replicate
p = tapply(kdm4.gene$val, factor(paste(kdm4.gene$te,kdm4.gene$gr)), function(x) t.test(x, mu=1)$p.value)
p.adj = p.adjust(p, method='BH')
		

##Kdm4 TKO (Supplementary Figure 3C)

kdm4.tko = read.csv('Kdm4_TKO/Kdm4_TKO_export_all.csv', as.is=T)

kdm4.tko1 = paper.plot(data=kdm4.tko,	
		te.list=c('KDM4A-GT','KDM4B-GT','KDM4C','Nanog'),
		g.list=c('2i-Jmjd-WT-VC','2i-Jmjd-KO-Ctrl','2i-Jmjd-KO-VC'),
		ylab='Expression relative to WT',
		width=3.5, y.lim=c(0,2))
	
p = tapply(kdm4.tko1$val, factor(paste(kdm4.tko1$te,kdm4.tko1$gr)), function(x) t.test(x, mu=1)$p.value)
p.adj = p.adjust(p, method='BH')


#Figure 4A

kdm4.tko2 = paper.plot(data=kdm4.tko,	
		te.list=c('Dazl', 'L1 TF', 'L1 A', 'L1 GF', 'ORF1a', 'ORF1b', 'ORF2', 'ORF2b'),
		g.list=c('2i-Jmjd-WT-VC','2i-Jmjd-KO-Ctrl','2i-Jmjd-KO-VC'),
		ylab='Expression relative to WT')

model = aov(kdm4.tko2$val ~ kdm4.tko2$te * kdm4.tko2$gr)
summary(model)
TukeyHSD(model)


##hESC primed (Figure 6B)

primed = read.csv('hESC_prime/hESC_prime_export_all.csv', as.is=T)

primed.df = paper.plot(data=primed,	
		te.list=c('L1-5UTR', 'L1-HS-1', 'L1-HS-2', 'hORF1', 'hORF2'),
		g.list=c('primed-VC-24hrs', 'primed-VC-48hrs'),
		ylab='Fold change', y.lim=c(0,3),
		width=4,
		col=c('orange','darkred'))


##Setdb1 KD (Figure 4D)

setdb1.kd = read.csv('Setdb1_KD/Setdb1_KD.csv', as.is=T)

setdb1.df = paper.plot(data=setdb1.kd,	
		te.list=c( 'L1 TF', 'L1 A', 'L1 GF', 'ORF1a', 'ORF2a'),
		g.list=c('shScr-VC','shSetdb1-Ctrl','shSetdb1-VC'),
		ylab='Expression relative to shScr',
		width=3.5, y.lim=c(0,35))
	
model = aov(setdb1.df$val ~ setdb1.df$te * setdb1.df$gr)
summary(model)
TukeyHSD(model)

