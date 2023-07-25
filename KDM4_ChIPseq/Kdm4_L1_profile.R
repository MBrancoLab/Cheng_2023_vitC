##Figure 1C


x = read.delim('Kdm4_L1_profile.txt')


##KDM4A

wt = x[,grep('2a_Ctrl',colnames(x))]
ko = x[,grep('2a_OHT',colnames(x))]

quartz(w=5.5,h=3.5)
plot(x$Position,wt$L1Gf.Jmjd2a_2a_Ctrl,
	type='l',lwd=2,col='blue',
	bty='n',las=1,
	xlab='Relative position',ylab='Normalised read count')
lines(x$Position,ko$L1Gf.Jmjd2a_2a_OHT,lwd=2,lty=2,col='blue')

lines(x$Position,wt$L1A.Jmjd2a_2a_Ctrl,lwd=2,col='grey')
lines(x$Position,ko$L1A.Jmjd2a_2a_OHT,lwd=2,lty=2,col='grey')

lines(x$Position,wt$L1Tf.Jmjd2a_2a_Ctrl,lwd=2,col='red')
lines(x$Position,ko$L1Tf.Jmjd2a_2a_OHT,lwd=2,lty=2,col='red')




##KDM4C


ip = x[,grep('2i_ChIP',colnames(x))]
input = x[,grep('2i_input',colnames(x))]

quartz(w=5.5,h=3.5)
plot(x$Position,ip$L1Gf.Jmjd2c_2i_ChIP,
	type='l',lwd=2,col='blue',
	bty='n',las=1,
	xlab='Relative position',ylab='Normalised read count')
lines(x$Position,input$L1Gf.Jmjd2c_2i_input,lwd=2,lty=2,col='blue')

lines(x$Position,ip$L1A.Jmjd2c_2i_ChIP,lwd=2,col='grey')
lines(x$Position,input$L1A.Jmjd2c_2i_input,lwd=2,lty=2,col='grey')

lines(x$Position,ip$L1Tf.Jmjd2c_2i_ChIP,lwd=2,col='red')
lines(x$Position,input$L1Tf.Jmjd2c_2i_input,lwd=2,lty=2,col='red')

