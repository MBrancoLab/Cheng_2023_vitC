
##Supplementary Figure 4A

p = c(1.012, 0.956)
n = c(0.841, 0.917)
h = c(0.289, 0.312, 0.065, 0.121)



quartz(w=4, h=3)
par(mar=c(5,4,2,2))

plot(c(1,1), p,
	xlim=c(0.5,3.5), ylim=c(0,1.5),
	xaxt='n', xlab='', ylab='Fold change',
	pch=19, bty='n', las=1,
	col=c('blue','red'))
axis(1,1:3,labels=c('Primed','Naive','HEF'))

points(c(2,2),n,pch=19,col=c('blue','red'))
points(rep(3,4),h,pch=19,col=c('blue','red'))

lines(c(0.8,1.2),c(mean(p),mean(p)), lwd=2)
lines(c(1.8,2.2),c(mean(n),mean(n)), lwd=2)
lines(c(2.8,3.2),c(mean(h),mean(h)), lwd=2)
abline(h=1, lty=2)

