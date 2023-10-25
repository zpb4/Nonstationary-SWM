#plot raw errors vs corrected errors

png('h:/oroville_non-stationary/paper/figs_rev1/si/sep-distribution-example_si.png',width=768,height=312)
par(mfrow=c(1,2),mar=c(2,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

plot(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=1,xi=1),type='l',col='black',lwd=3,
     xlab='',ylab='Density',ylim=c(0,1),xlim=c(-3,3))
lines(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=2,xi=1),col='seagreen4',lwd=3)
lines(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=0.5,xi=1),col='gray50',lwd=2,lty=2)
text(0,0.1,bquote(~xi==.(1)),cex=2.5,col='black',adj=0.5)
text(2.2,0.2,bquote(~beta==.(0)),cex=2.5,col='seagreen4')
text(1,0.5,bquote(~beta==.(1)),cex=2.5,col='black')
text(0.8,0.8,bquote(~beta==.(3)),cex=2.5,col='gray50')

plot(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=1,xi=1),type='l',col='black',lwd=3,
     xlab='',ylab='Density',ylim=c(0,1),xlim=c(-3,3))
lines(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=1,xi=0.1),col='gray50',lwd=2,lty=2)
lines(seq(-5,5,0.1),dsged(seq(-5,5,0.1),mean=0,sd=1,nu=1,xi=10),col='sienna4',lwd=2,lty=2)
text(0,0.1,bquote(~beta==.(1)),cex=2.5,col='black',adj=0.5)
text(0,0.95,bquote(log[10]~xi==.(0)),cex=1.75,col='black',adj=0.5)
text(0,0.85,bquote(~xi==.(1)),cex=1.75,col='black',adj=0.5)
text(2.1,0.5,bquote(log[10]~xi==.(-1)),cex=1.75,col='gray50')
text(-2.1,0.5,bquote(log[10]~xi==.(1)),cex=1.75,col='sienna4')
text(2.1,0.4,bquote(~xi==.(0.1)),cex=1.75,col='gray50')
text(-2.1,0.4,bquote(~xi==.(10)),cex=1.75,col='sienna4')

dev.off()


################################END###################################
