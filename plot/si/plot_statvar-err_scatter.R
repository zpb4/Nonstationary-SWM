setwd('z:/oro_nonstat/')
library(raster)
library(RColorBrewer)
library(ranger)
library(stringr)
source('mm-cfs_conversion.R')

hym_site<-'ORO'
vers<-'err13'
seed<-1

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
ix_trn<-seq(as.Date('1988-10-01'),as.Date('2004-09-30'),'day')
ixx_trn<-as.POSIXlt(ix_trn)

idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
ix_val<-seq(as.Date('1998-10-01'),as.Date('2004-09-30'),'day')
ixx_val<-as.POSIXlt(ix_val)

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))
rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))

rf_pred<-hym_predmat_hist[,rf_idx]
rf_err<-hym_predmat_hist[,'err 0']

rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#3) Fit multivariate, dynamic GL-SEP model for residuals
#debias calibration errors with RF error correction model
db_trn<-predict(rf_err_corr,data=rf_pred[idx_val,])$predictions
err_db_trn<-rf_err[idx_val]-db_trn


png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/si/error-statvar_scatter.png',sep=''),height=768,width=768)

par(mfrow=c(3,3),mar=c(4,4.5,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.lab=3,cex.axis=2)
for(i in 1:9){
  plot(hym_predmat_hist[idx_val,'err 0'],rf_pred[idx_val,i],xlab='error',ylab=str_remove(colnames(rf_pred)[i],' 0'))
}

dev.off()

png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/si/resid-statvar_scatter.png',sep=''),height=768,width=768)

par(mfrow=c(3,3),mar=c(4,4.5,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.lab=3,cex.axis=2)
for(i in 1:9){
  plot(err_db_trn,rf_pred[idx_val,i],xlab='residual',ylab=str_remove(colnames(rf_pred)[i],' 0'))
}

dev.off()

################################END##########################