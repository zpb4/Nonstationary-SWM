#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(ranger)
library(stringr)

#specifications
hym_site<-'ORO'
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
noise_reg=F

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_cal<-seq(as.Date('1988-10-01'),as.Date('1998-09-30'),'day')
ixx_cal<-as.POSIXlt(ix_cal)
idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')

#fit RF model for error correction
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
err<-hym_predmat_hist[idx_cal,'err 0']

rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}
#--------------------------------------------------------------
#variable importance
var_vec<-colnames(hym_predmat_hist[,rf_idx])
var_vec<-str_remove(var_vec,' 0')
var_imp<-rf_err_corr$variable.importance / sum(rf_err_corr$variable.importance)
srt_var_imp<-sort(var_imp,index.return=T,decreasing=T)

png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig6/fig6_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=640)
par(las=2,mar=c(6.5,5.5,1,0),mgp=c(4,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)
barplot(srt_var_imp$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp$x[1],digits=2),adj=0.5,cex=1.5)

dev.off()

##########################################END######################################
