#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')
library(ranger)


#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
hym_site<-'ORO'
rerun_search=F

noise_reg=T

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

#load predictor arrays
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}


ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

#1) Prepare data
rf_err<-hym_predmat_hist[,'err 0']
rf_pred<-hym_predmat_hist[,rf_idx]

if(rerun_search==T){
  set.seed(seed)
  ntree_grid<-rep(seq(10,500,10),each=12)
  mtry_grid<-rep(1:12,length(seq(10,500,10)))
  pred_err<-c()

  for(i in 1:length(ntree_grid)){
    rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],num.trees = ntree_grid[i],mtry=mtry_grid[i],importance = 'impurity')
    pred_err[i]<-rf_err_corr$prediction.error
  }
  idx<-which.min(pred_err)
  ntree<-ntree_grid[idx]
  mtr<-mtry_grid[idx]

  rf_err_corr_hyp_tune<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],num.trees = ntree,mtry=mtr,importance = 'impurity')
  saveRDS(rf_err_corr_hyp_tune,paste('fit_rev1/rf-err-corr_hyp-tune_v-',vers,'.rds',sep=''))
}

if(rerun_search==F){
  rf_err_corr_hyp_tune<-readRDS(paste('fit_rev1/rf-err-corr_hyp-tune_v-',vers,'.rds',sep=''))
}

rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#default vs hyper parameter tuned combined
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/hyp-tuned-comp_',vers,'_nreg=',noise_reg,'_figsS2.png',sep=''),width=768,height=900)
par(mfrow=c(2,1),mar=c(1,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

db_err_def<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
db_err_hyp<-predict(rf_err_corr_hyp_tune,data=hym_predmat_hist[idx_comp,rf_idx])$predictions

err_db_def<-hym_predmat_hist[idx_comp,'err 0']-db_err_def
err_db_hyp<-hym_predmat_hist[idx_comp,'err 0']-db_err_hyp

err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*3-2)]]<-hym_predmat_hist[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*3-1)]]<-err_db_def[seas]*hym_kcfs_conv
  err_lst[[(i*3)]]<-err_db_hyp[seas]*hym_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=3)+rep(c(0,.3,.5),12),
        names = 1:36,
        las = 2,
        boxwex = rep(c(.3,.15,.15),12),
        notch = F,
        col = rep(c('skyblue4','skyblue','lightblue1'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Raw','Default RF-corrected','Tuned RF-corrected'),lwd=c(6,6,6),col=c('skyblue4','skyblue','lightblue1'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,'a)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

par(mar=c(2,4,0,1))
#4c
db_err_def<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,rf_idx])$predictions
db_err_hyp<-predict(rf_err_corr_hyp_tune,data=hym_predmat_4c[idx_comp,rf_idx])$predictions

err_db_def<-hym_predmat_4c[idx_comp,'err 0']-db_err_def
err_db_hyp<-hym_predmat_4c[idx_comp,'err 0']-db_err_hyp

err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*3-2)]]<-hym_predmat_4c[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*3-1)]]<-err_db_def[seas]*hym_kcfs_conv
  err_lst[[(i*3)]]<-err_db_hyp[seas]*hym_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=3)+rep(c(0,.3,0.5),12),
        names = 1:36,
        las = 2,
        boxwex = rep(c(.3,.15,.15),12),
        notch = F,
        col = rep(c('darkorange3','goldenrod1','gold1'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Raw','Default RF-corrected','Tuned RF-corrected'),lwd=c(6,6,6),col=c('darkorange3','goldenrod1','gold1'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,'b)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

dev.off()

#figs2b
#variable importance
var_vec<-colnames(hym_predmat_hist[,rf_idx])
var_vec<-str_remove(var_vec,' 0')
var_imp<-rf_err_corr_hyp_tune$variable.importance / sum(rf_err_corr_hyp_tune$variable.importance)
srt_var_imp<-sort(var_imp,index.return=T,decreasing=T)

var_imp_def<-rf_err_corr$variable.importance / sum(rf_err_corr$variable.importance)
srt_var_imp_def<-sort(var_imp_def,index.return=T,decreasing=T)

png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/rf-hyp-tune-comp_var-imp_v-',vers,'_nreg=',noise_reg,'_figS2.png',sep=''),width=768,height=384)
par(mfrow=c(1,2),las=2,mar=c(6.5,5.5,1,0),mgp=c(4,0.5,0),tcl=-0.2,cex.lab=1.75,cex.axis=1.75)
barplot(srt_var_imp$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp$x[1],digits=2),adj=0.5,cex=1.5)
text(length(var_vec),0.19,'a)',cex=3,font=2)

barplot(srt_var_imp_def$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp_def$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp_def$x[1],digits=2),adj=0.5,cex=1.5)
text(length(var_vec),0.19,'b)',cex=3,font=2)

dev.off()

#---------------------------------------------------------------------------------
#S2.3
gen_period1<-'hist-all'
gen_period2<-'4c-all'
samp_type<-'split'
n<-10

syn_tst_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period1,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_hyp-tune.rds',sep=''))
syn_4c_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period2,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_hyp-tune.rds',sep=''))

syn_tst_err_hybrid<-syn_tst_err_hyb[idx_comp,]*hym_kcfs_conv
syn_4c_err_hybrid<-syn_4c_err_hyb[idx_comp,]*hym_kcfs_conv

#plot raw errors vs corrected errors
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/rf-hyp-tune-comp_swm-err_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'_figS2.png',sep=''),width=768,height=1024)
par(mfrow=c(2,1),oma=c(1,2,2,0),mar=c(0.5,1.5,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

#a) hybrid test
err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-hym_predmat_hist[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-syn_tst_err_hybrid[seas,]
}

mths<-c('J','F','M','A','M','J','J','A','S','O','N','D')
#mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(-0.05,.305),12),
        names = 1:24,
        las = 2,
        boxwex = 0.25,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red')
text(0.75,4.75,'a)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)

#b) 4c hybrid

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-hym_predmat_4c[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-syn_4c_err_hybrid[seas,]
}

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(-0.05,.305),12),
        names = 1:24,
        las = 2,
        boxwex = 0.25,
        notch = F,
        col = rep(c('darkorange3','coral1'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('darkorange3','coral1'),cex=2)
abline(h=0,lty=2,col='red')
text(0.75,4.75,'b)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)

dev.off()



##########################################END####################################################