#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')
source('GL_maineqs_rev1.R')
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

noise_reg=F

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

seed<-4
#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'

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

res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))

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

rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

norm_vec<-readRDS(paste('fit_rev1/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))
db_err<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-hym_predmat_hist[idx_comp,'err 0']-db_err

#val-fit
if(noise_reg==T){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
}
if(noise_reg==F){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_pred-mat-zero_val_',hym_site,'.rds',sep=''))
}

#test
dyn_res_preds<-hym_predmat_hist[,res_idx]
pred_mat_scale<-(dyn_res_preds[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale)[2],byrow=T)

set.seed(seed)
syn_res<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=err_db)*hym_kcfs_conv
syn_res[is.na(syn_res)==T]<-0

#test4c
dyn_res_preds_4c<-hym_predmat_4c[,res_idx]
pred_mat_scale_4c<-(dyn_res_preds_4c[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T)
pred_mat_zero_min_4c<-pred_mat_scale_4c-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale_4c)[2],byrow=T)

set.seed(seed)
syn_res_4c<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min_4c,beta_var=pred_mat_scale_4c,xi_var=pred_mat_scale_4c,phi_var=pred_mat_zero_min_4c,et=err_db)*hym_kcfs_conv
syn_res_4c[is.na(syn_res_4c)==T]<-0


#cal-fit
if(noise_reg==T){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'_cal.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'_cal.rds',sep=''))
}
if(noise_reg==F){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'_cal.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_pred-mat-zero_val_',hym_site,'.rds',sep=''))
}

#test
dyn_res_preds<-hym_predmat_hist[,res_idx]
pred_mat_scale<-(dyn_res_preds[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale)[2],byrow=T)

set.seed(seed)
syn_res_cal<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=err_db)*hym_kcfs_conv
syn_res_cal[is.na(syn_res_cal)==T]<-0

#test4c
dyn_res_preds_4c<-hym_predmat_4c[,res_idx]
pred_mat_scale_4c<-(dyn_res_preds_4c[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T)
pred_mat_zero_min_4c<-pred_mat_scale_4c-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale_4c)[2],byrow=T)

set.seed(seed)
syn_res_4c_cal<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min_4c,beta_var=pred_mat_scale_4c,xi_var=pred_mat_scale_4c,phi_var=pred_mat_zero_min_4c,et=err_db)*hym_kcfs_conv
syn_res_4c_cal[is.na(syn_res_4c_cal)==T]<-0

#default vs hyper parameter tuned combined
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/cal-val-resid-comp_',vers,'_nreg=',noise_reg,'_figsS3.png',sep=''),width=768,height=900)
par(mfrow=c(2,1),mar=c(1,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

db_err<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-hym_predmat_hist[idx_comp,'err 0']-db_err


err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*3-2)]]<-err_db[seas]*hym_kcfs_conv
  err_lst[[(i*3-1)]]<-syn_res[seas]*hym_kcfs_conv
  err_lst[[(i*3)]]<-syn_res_cal[seas]*hym_kcfs_conv
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
        ylab='Residual (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Sim_val-fit','Sim_cal-fit'),lwd=c(6,6,6),col=c('skyblue4','skyblue','lightblue1'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,'a)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

par(mar=c(2,4,0,1))

#4c
db_err<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,rf_idx])$predictions
err_db<-hym_predmat_4c[idx_comp,'err 0']-db_err

err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*3-2)]]<-err_db[seas]*hym_kcfs_conv
  err_lst[[(i*3-1)]]<-syn_res_4c[seas]*hym_kcfs_conv
  err_lst[[(i*3)]]<-syn_res_4c_cal[seas]*hym_kcfs_conv
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
        ylab='Residual (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Sim_val-fit','Sim_cal-fit'),lwd=c(6,6,6),col=c('darkorange3','goldenrod1','gold1'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,'b)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

dev.off()


##########################################END####################################################