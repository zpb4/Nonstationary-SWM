#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(Kendall)
sma_sites<-c('ORO','SHA','NHG')


ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

yrs<-89:117

tavg_vec<-array(NA,c(3,length(yrs)))
precip_vec<-array(NA,c(3,length(yrs)))
mk_tavg<-c()
mk_precip<-c()

for(k in 1:3){
sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_sites[k],'.rds',sep=''))

#avg temp
for(i in 1:length(yrs)){
  yr_idx<-which(ix2$year==yrs[i])
  tavg_vec[k,i]<-mean(sma_predmat_hist[yr_idx,'tavg 0'])
  precip_vec[k,i]<-sum(sma_predmat_hist[yr_idx,'precip 0'])
}
mk_tavg[k]<-MannKendall(tavg_vec[k,])$sl
mk_precip[k]<-MannKendall(precip_vec[k,])$sl
}


yr2<-1989:2017
yupr_lim<-max(tavg_vec)+3
ylwr_lim<-min(tavg_vec)-3

png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/obs-warm_annual_S9.png',sep=''),width=768,height=1024)
par(mfrow=c(3,1),mar=c(3,4,3,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75,cex.main=3)
for(k in 1:3){
plot(yr2,tavg_vec[k,],main=paste('Mean Annual Temperature - ',sma_sites[k],sep=''),xlab='',ylab='Temp (C)',
     ylim=c(ylwr_lim,yupr_lim))
lfit<-lm(tavg_vec[k,]~yr2)
abline(a=lfit$coefficients[1],b=lfit$coefficients[2],col='red')
text(1998,yupr_lim-0.5,paste('Slope: ',round(lfit$coefficients[2],digits=5),'C/year',sep=''),col='red',cex=2,adj=0)
text(1998,yupr_lim-1.5,paste('Total 1988-2018: ',round(lfit$coefficients[2]*length(yrs),digits=5),'C',sep=''),col='red',cex=2,adj=0)
text(1998,yupr_lim-2.5,paste('Mann-Kendall p: ',round(mk_tavg[k],digits=5),sep=''),col='red',cex=2,adj=0)
}

dev.off()



yupr_lim<-max(precip_vec)+3
ylwr_lim<-min(precip_vec)-3

png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/obs-precip_annual_S9.png',sep=''),width=768,height=1024)
par(mfrow=c(3,1),mar=c(3,4,3,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75,cex.main=3)
for(k in 1:3){
plot(yr2,precip_vec[k,],main=paste('Annual Precipitation - ',sma_sites[k],sep=''),xlab='',ylab='Precip (mm)',
     ylim=c(ylwr_lim,yupr_lim))
lfit<-lm(precip_vec[k,]~yr2)
abline(a=lfit$coefficients[1],b=lfit$coefficients[2],col='red')
text(1998,yupr_lim-50,paste('Slope: ',round(lfit$coefficients[2],digits=5),'mm/year',sep=''),col='red',cex=2,adj=0)
text(1998,yupr_lim-200,paste('Total 1988-2018: ',round(lfit$coefficients[2]*length(yrs),digits=5),' mm',sep=''),col='red',cex=2,adj=0)
text(1998,yupr_lim-350,paste('Mann-Kendall p: ',round(mk_precip[k],digits=5),sep=''),col='red',cex=2,adj=0)
}

dev.off()


########################################END###################################################
