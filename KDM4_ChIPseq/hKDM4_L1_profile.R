
x = read.delim('hKDM4_L1_profile.txt')


##KDM4A (Figure 6A)


quartz(w=5.5,h=3.5)
plot(x$Position,x$L1Hs.ENCFF000AZD.bam,
	type='l',lwd=2,col='blue',
	bty='n',las=1,
	xlab='Relative position',ylab='Normalised read count')
lines(x$Position,x$L1Hs.ENCFF303YPV.bam,lwd=2,lty=2,col='blue')

lines(x$Position,x$L1PA2.ENCFF000AZD.bam,lwd=2,col='red')
lines(x$Position,x$L1PA2.ENCFF303YPV.bam,lwd=2,lty=2,col='red')

lines(x$Position,x$L1PA3.ENCFF000AZD.bam,lwd=2,col='grey')
lines(x$Position,x$L1PA3.ENCFF303YPV.bam,lwd=2,lty=2,col='grey')


