
#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
source('GL_maineqs_rev1.R')

hym_site<-'ORO'
sma_site<-'ORO'
#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
sig_ylms<-c(0,0.8)
bet_ylms<-c(0,2.75)
xi_ylms<-c(-.3,.3)
phi_ylms<-c(0,1)

noise_reg=T

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)}

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))
rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
norm_vec<-readRDS(paste('fit_rev1/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

if(noise_reg==T){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
}
if(noise_reg==F){
  dyn_res_coef<-readRDS(paste('fit_rev1/hymod_dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/hymod_pred-mat-zero_val_',hym_site,'.rds',sep=''))
}

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}

#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))
#--------------------------------------------------------------
#plot raw errors vs corrected errors
db_err<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-hym_predmat_hist[idx_comp,'err 0']-db_err

db_err_4c<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,rf_idx])$predictions
err_db_4c<-hym_predmat_4c[idx_comp,'err 0']-db_err_4c
#----------------------------------------------------------------------
pred_mat_rf<-hym_predmat_hist[idx_comp,rf_idx]
dyn_res_preds<-hym_predmat_hist[idx_comp,res_idx]

pred_mat_scale<-(dyn_res_preds-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale)[2],byrow=T)

param_out_tst<-GL_fun_mv_ar1_lin_params(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=err_db)


#test4c
pred_mat_4c<-hym_predmat_4c[idx_comp,rf_idx]
dyn_res_preds_4c<-hym_predmat_4c[idx_comp,res_idx]

pred_mat_scale_4c<-(dyn_res_preds_4c-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T)
pred_mat_zero_min_4c<-pred_mat_scale_4c-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale_4c)[2],byrow=T)

param_out_4c<-GL_fun_mv_ar1_lin_params(dyn_res_coef,sig_var=pred_mat_zero_min_4c,beta_var=pred_mat_scale_4c,xi_var=pred_mat_scale_4c,phi_var=pred_mat_zero_min_4c,et=err_db_4c)

sig_lst_tst<-vector('list',12)
beta_lst_tst<-vector('list',12)
xi_lst_tst<-vector('list',12)
phi_lst_tst<-vector('list',12)

sig_lst_4c<-vector('list',12)
beta_lst_4c<-vector('list',12)
xi_lst_4c<-vector('list',12)
phi_lst_4c<-vector('list',12)

for(i in 1:12){
  seas_tst<-which(ixx_comp$mon==(i-1))
  sig_lst_tst[[i]]<-param_out_tst[[1]][seas_tst]
  beta_lst_tst[[i]]<-param_out_tst[[2]][seas_tst]
  xi_lst_tst[[i]]<-param_out_tst[[3]][seas_tst]
  phi_lst_tst[[i]]<-param_out_tst[[4]][seas_tst]
  
  seas_4c<-which(ixx_comp$mon==(i-1))
  sig_lst_4c[[i]]<-param_out_4c[[1]][seas_4c]
  beta_lst_4c[[i]]<-param_out_4c[[2]][seas_4c]
  xi_lst_4c[[i]]<-param_out_4c[[3]][seas_4c]
  phi_lst_4c[[i]]<-param_out_4c[[4]][seas_4c]
}

upr<-function(x){
  srt<-sort(x)
  out<-srt[round(0.95*length(x))]
  return(out)}

lwr<-function(x){
  srt<-sort(x)
  out<-srt[round(0.05*length(x))]
  return(out)}

sig_mn_tst<-lapply(sig_lst_tst,mean)
sig_max_tst<-lapply(sig_lst_tst,upr)
sig_min_tst<-lapply(sig_lst_tst,lwr)

sig_mn_4c<-lapply(sig_lst_4c,mean)
sig_max_4c<-lapply(sig_lst_4c,upr)
sig_min_4c<-lapply(sig_lst_4c,lwr)

beta_mn_tst<-lapply(beta_lst_tst,mean)
beta_max_tst<-lapply(beta_lst_tst,upr)
beta_min_tst<-lapply(beta_lst_tst,lwr)

beta_mn_4c<-lapply(beta_lst_4c,mean)
beta_max_4c<-lapply(beta_lst_4c,upr)
beta_min_4c<-lapply(beta_lst_4c,lwr)

xi_mn_tst<-lapply(xi_lst_tst,mean)
xi_max_tst<-lapply(xi_lst_tst,upr)
xi_min_tst<-lapply(xi_lst_tst,lwr)

xi_mn_4c<-lapply(xi_lst_4c,mean)
xi_max_4c<-lapply(xi_lst_4c,upr)
xi_min_4c<-lapply(xi_lst_4c,lwr)

phi_mn_tst<-lapply(phi_lst_tst,mean)
phi_max_tst<-lapply(phi_lst_tst,upr)
phi_min_tst<-lapply(phi_lst_tst,lwr)

phi_mn_4c<-lapply(phi_lst_4c,mean)
phi_max_4c<-lapply(phi_lst_4c,upr)
phi_min_4c<-lapply(phi_lst_4c,lwr)

mos<-c('J','F','M','A','M','J','J','A','S','O','N','D')

png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig9/fig9_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=768)
par(mfrow=c(2,2),mar=c(2,4.5,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=3,cex.axis=1.75,cex.main=3,font.main=4)

plot(1:12,sig_mn_tst,type='l',xlab='',ylab=bquote(~sigma[t]),col='skyblue4',lwd=3,xaxt='n',ylim=sig_ylms)
axis(1,at=1:12,labels = mos)
polygon(c(1:12,12:1),c(sig_max_tst,rev(sig_min_tst)),col='light cyan',border = NA)
lines(1:12,sig_mn_tst,col='skyblue4',lwd=3)
lines(1:12,sig_mn_4c,col='darkorange3',lwd=3)
lines(1:12,sig_max_4c,col='darkorange3',lty=2)
lines(1:12,sig_min_4c,col='darkorange3',lty=2)
box(which='plot')
legend('topright',c('Test','Test+4C'),lwd=c(3,3),col=c('skyblue4','darkorange3'),bty='n',cex=2.5)

plot(1:12,beta_mn_tst,type='l',xlab='',ylab=bquote(~beta[t]),col='skyblue4',lwd=3,xaxt='n',ylim=bet_ylms)
axis(1,at=1:12,labels = mos)
polygon(c(1:12,12:1),c(beta_max_tst,rev(beta_min_tst)),col='light cyan',border = NA)
lines(1:12,beta_mn_tst,col='skyblue4',lwd=3)
lines(1:12,beta_mn_4c,col='darkorange3',lwd=3)
lines(1:12,beta_max_4c,col='darkorange3',lty=2)
lines(1:12,beta_min_4c,col='darkorange3',lty=2)
box(which='plot')

plot(1:12,xi_mn_tst,type='l',xlab='',ylab=bquote(~log[10]~xi[t]),col='skyblue4',lwd=3,xaxt='n',ylim=xi_ylms)
axis(1,at=1:12,labels = mos)
polygon(c(1:12,12:1),c(xi_max_tst,rev(xi_min_tst)),col='light cyan',border = NA)
lines(1:12,xi_mn_tst,col='skyblue4',lwd=3)
lines(1:12,xi_mn_4c,col='darkorange3',lwd=3)
lines(1:12,xi_max_4c,col='darkorange3',lty=2)
lines(1:12,xi_min_4c,col='darkorange3',lty=2)
box(which='plot')

plot(1:12,phi_mn_tst,type='l',xlab='',ylab=bquote(~phi[t]),col='skyblue4',lwd=3,xaxt='n',ylim=phi_ylms)
axis(1,at=1:12,labels = mos)
polygon(c(1:12,12:1),c(phi_max_tst,rev(phi_min_tst)),col='light cyan',border = NA)
lines(1:12,phi_mn_tst,col='skyblue4',lwd=3)
lines(1:12,phi_mn_4c,col='darkorange3',lwd=3)
lines(1:12,phi_max_4c,col='darkorange3',lty=2)
lines(1:12,phi_min_4c,col='darkorange3',lty=2)
box(which='plot')

dev.off()

#####################################END#############################