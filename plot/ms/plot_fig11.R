setwd('z:/oro_nonstat/')
library(raster)
library(RColorBrewer)
source('mm-cfs_conversion.R')

hym_site<-'ORO'
sma_site<-'ORO'
vers<-'err13'
gen_period1<-'hist-all'
gen_period2<-'4c-all'
samp_type<-'split'
noise_reg=T
n<-1000

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')
ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))
sim_hymod_tst<-hym_predmat_hist[idx_tst,'sim 0']*hym_kcfs_conv
sim_hymod_4c<-hym_predmat_4c[idx_tst,'sim 0']*hym_kcfs_conv

sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))
sma_predmat_4c<-readRDS(paste('data_rev1/sma_predmat_4c_',sma_site,'.rds',sep=''))
sim_sma_tst<-sma_predmat_hist[idx_tst,'sim 0']*hym_kcfs_conv
sim_sma_4c<-sma_predmat_4c[idx_tst,'sim 0']*hym_kcfs_conv

swm_out_hist<-readRDS(paste('out_rev1/hymod_syn-flow_',hym_site,'_',gen_period1,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))*hym_kcfs_conv
swm_out_4c<-readRDS(paste('out_rev1/hymod_syn-flow_',hym_site,'_',gen_period2,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))*hym_kcfs_conv

swm_out_tst<-swm_out_hist[idx_tst,]
swm_out_4c<-swm_out_4c[idx_tst,]

pcntile<-function(x){
  p<-sort(x)
  l025<-p[round(0.025*length(x))]
  l25<-p[round(0.25*length(x))]
  u75<-p[round(0.75*length(x))]
  u975<-p[round(0.975*length(x))]
  return(c(l025,l25,u75,u975))
}

swm_pcnt_tst<-apply(swm_out_tst,1,pcntile)
swm_pcnt_4c<-apply(swm_out_4c,1,pcntile)

swm_med_tst<-apply(swm_out_tst,1,median)
swm_med_4c<-apply(swm_out_4c,1,median)

#plots - pick year WY2004-2013
yr<-2011  #march

sc<-1.05
#sset_mth<-3

wy<-which(ixx_tst$year==(yr-1900)&ixx_tst$mo>0&ixx_tst$mo<7)
#wy_sset<-which(ixx_tst$year==(yr-1900)&ixx_tst$mo==(sset_mth-1))
ymax<-max(sim_sma_tst[wy],swm_med_tst[wy],sim_sma_4c[wy],swm_med_4c[wy])

wy_lab<-c('Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep')

png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig11/swm-ts_95_fig11_feb-jul_',yr,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),height=900,width=768)

#layout_mat<-matrix(c(1,1,2,2,3,4),ncol=2,byrow=T)
#layout(layout_mat)
par(mfrow=c(2,1),mar=c(2,4,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.5)

#tst timeseries plot
plot(1:length(wy),sim_sma_tst[wy],type='l',lwd=1,xlim=c(7,length(wy)-6),ylim=c(0,sc*ymax),axes=F,xlab='',ylab='Flow (kcfs)')
axis(1,at=seq(15,length(wy),30),labels=F)
axis(2,at=seq(0,sc*ymax,5),labels=seq(0,sc*ymax,5))
polygon(c(1:length(wy),rev(1:length(wy))),c(swm_pcnt_tst[1,wy],rev(swm_pcnt_tst[4,wy])),col='gray90',border=NA)
lines(1:length(wy),sim_hymod_tst[wy],lwd=1,col='skyblue4')
lines(1:length(wy),sim_sma_tst[wy],lwd=2)
lines(1:length(wy),swm_med_tst[wy],lwd=1,col='turquoise3')
legend('topright',c('truth','process','hybrid SWM median','95%'),lwd=c(2,1,1,6),col=c('black','skyblue4','turquoise3','gray90'),cex=2,bty='n')
box(which='plot')
#mtext('Test',side=3,line=1,outer=F,font=2,cex=2)

#4c plot at same dates
plot(1:length(wy),sim_sma_4c[wy],type='l',lwd=1,xlim=c(7,length(wy)-6),ylim=c(0,sc*ymax),axes=F,xlab='',ylab='Flow (kcfs)')
axis(1,at=seq(15,length(wy),30),labels=wy_lab[5:10])
axis(2,at=seq(0,sc*ymax,5),labels=seq(0,sc*ymax,5))
polygon(c(1:length(wy),rev(1:length(wy))),c(swm_pcnt_4c[1,wy],rev(swm_pcnt_4c[4,wy])),col='gray90',border=NA)
lines(1:length(wy),sim_hymod_4c[wy],lwd=1,col='darkorange3')
lines(1:length(wy),sim_sma_4c[wy],lwd=2)
lines(1:length(wy),swm_med_4c[wy],lwd=1,col='coral1')
legend('topright',c('truth','process','hybrid SWM median','95%'),lwd=c(2,1,1,6),col=c('black','darkorange3','coral1','gray90'),cex=2,bty='n')
box(which='plot')
#mtext('Test',side=3,line=1,outer=F,font=2,cex=2)

#plot(1:length(wy_sset),sim_sma_4c[wy_sset],type='l',lwd=1,xlim=c(7,length(wy)-6),ylim=c(0,sc*ymax),axes=F,xlab='',ylab='Flow (kcfs)')
#axis(1,at=seq(15,length(wy_sset),30),labels=wy_lab[3])
#axis(2,at=seq(0,sc*ymax,5),labels=seq(0,sc*ymax,5))
#polygon(c(1:length(wy),rev(1:length(wy))),c(swm_pcnt_4c[1,wy],rev(swm_pcnt_4c[4,wy])),col='gray90',border=NA)
#lines(1:length(wy),sim_hymod_4c[wy],lwd=1,col='darkorange3')
#lines(1:length(wy),sim_sma_4c[wy],lwd=2)
#lines(1:length(wy),swm_med_4c[wy],lwd=1,col='coral1')
#legend('topright',c('truth','process','hybrid SWM median','95%'),lwd=c(2,1,1,6),col=c('black','darkorange3','coral1','gray90'),cex=2,bty='n')
#box(which='plot')

dev.off()

################################END##########################