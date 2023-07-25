##Supplementary Figure 3B


l1 = c('Lx','L1_Mus4','L1_Mus3','L1_Mus2','L1_Mus1','L1Md_F','L1Md_F2','L1Md_F3','L1Md_Gf','L1Md_A','L1Md_T')

all = c(24769,11854,48280,21044,27127,3988,64716,16093,1077,16802,23639)
all.inc = c(20,8,71,27,34,116,142,114,319,1239,3292)/all*100
all.uni = c(8,3,19,10,5,37,57,69,35,66,373)/all*100

fl = c(170,130,449,349,705,157,3941,841,110,3295,5046)
fl.inc = c(4,2,1,7,11,12,69,83,71,1070,2984)/fl*100
fl.uni = c(2,1,0,3,4,4,36,54,8,47,335)/fl*100

par(mfrow=c(1,2))
plot(1:length(fl),fl.uni,pch=15,col='orange',xaxt='n',ylab='% bound',main='Unique',xlab='',las=1)
lines(1:length(fl),fl.uni,lwd=2,col='orange')
points(1:length(all),all.uni,pch=15,col='black')
lines(1:length(all),all.uni,lwd=2,col='black')
axis(1,1:length(all),labels=l1,las=2)

plot(1:length(fl),fl.inc,pch=15,col='orange',xaxt='n',ylab='% bound',main='Ambiguous',xlab='',las=1)
lines(1:length(fl),fl.inc,lwd=2,col='orange')
points(1:length(all),all.inc,pch=15,col='black')
lines(1:length(all),all.inc,lwd=2,col='black')
axis(1,1:length(all),labels=l1,las=2)

