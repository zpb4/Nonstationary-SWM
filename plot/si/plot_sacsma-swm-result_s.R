setwd('z:/oro_nonstat/')
library(ranger)
source('mm-cfs_conversion.R')


sma_site<-'SHA'
#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
noise_reg=T

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/sacsma_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/sacsma_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'

sma_kcfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs/1000

sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))
sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(sma_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','store 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(sma_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','store 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(sma_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','store 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(sma_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','store 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}

#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(sma_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','store 0')))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

rf_err_corr<-readRDS(paste('fit_rev1/sacsma_rf-err-corr_cal_',sma_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
norm_vec<-readRDS(paste('fit_rev1/sacsma_norm-vec_hist_',sma_site,'.rds',sep=''))

if(noise_reg==T){
  dyn_res_coef<-readRDS(paste('fit_rev1/sacsma_noise-dyn-res-coef_',sma_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/sacsma_noise-pred-mat-zero_val_',sma_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
}
if(noise_reg==F){
  dyn_res_coef<-readRDS(paste('fit_rev1/sacsma_dyn-res-coef_',sma_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('fit_rev1/sacsma_pred-mat-zero_val_',sma_site,'.rds',sep=''))
}
#--------------------------------------------------------------
#plot raw errors vs corrected errors
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/sacsma-swm-result_',sma_site,'_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'_figS12.png',sep=''),width=768,height=900)
par(mar=c(4,4,3,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75,cex.main=2)

layout_mat<-matrix(c(rep(1,8),2:5),ncol=4,byrow=T)
layout(layout_mat)

ylm=c(-5,5)
yaxs=seq(-5,5,2.5)
if(sma_site=='NHG'){ylm=c(-.2,.2);yaxs=seq(-0.2,0.2,0.1)}

err_hist<-sma_predmat_hist[idx_comp,'err 0']

db_err<-predict(rf_err_corr,data=sma_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-err_hist-db_err

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_hist[seas]*sma_kcfs_conv
  err_lst[[(i*2)]]<-err_db[seas]*sma_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.3),12),
        names = 1:24,
        las = 2,
        boxwex = rep(c(.3,.2),12),
        notch = F,
        col = rep(c('skyblue4','skyblue'),12),
        ylim = ylm,
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=yaxs,labels = yaxs )
box(which='plot')
legend('topright',c('Raw','RF-corrected'),lwd=c(6,6),col=c('skyblue4','skyblue'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,sma_site,cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)
if(sma_site=='NHG'){
  text(0.75,0.18,sma_site,cex=3,font=2,adj=0)
  text(11,-0.18,wy_tag,cex=3,font=2,adj=0)
}


source('GL_maineqs_rev1.R')

#test
dyn_res_preds<-sma_predmat_hist[,res_idx]
pred_mat_scale<-(dyn_res_preds[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale)[2],byrow=T)

set.seed(seed)
syn_res<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=err_db)*sma_kcfs_conv
syn_res[is.na(syn_res)==T]<-0

brks<-c(-100,seq(-5,5,0.4),100)
xlm<-c(-4,4)
x_axis<-seq(-4,4,2)
yupr<-1
pad<-0.015
if(sma_site=='NHG'){
  brks<-c(-100,seq(-5,5,0.1),100)
  xlm<-c(-1,1)
  x_axis<-seq(-2,2,1)
  yupr<-5
  pad<-0.01
}

mths<-c(1,4,7,10)

mth<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

for(i in 1:4){
  seas<-which(ixx_comp$mon==(mths[i]-1))
  kd_syn_err<-density(syn_res[seas],n=512)
  if(i==1){dens_lab='Density'}else{dens_lab=''}
  hist(err_db[seas]*sma_kcfs_conv,breaks=brks,xlim=xlm,ylim=c(0,yupr),freq=F,main=mth[mths[i]],col='skyblue',xlab='Residual (kcfs)',ylab=dens_lab,axes=F)
  axis(1,at=x_axis,labels=x_axis)
  if(i==1){axis(2,at=seq(0,5,.5),labels=seq(0,5,.5))}
  if(i>1){axis(2,at=seq(0,5,.5),labels=F)}
  kd_y<-kd_syn_err$y
  kd_y[kd_y<=pad[i]]<-pad[i]
  lines(kd_syn_err$x,kd_y,col='red',lwd=3)
  if(i==1){mtext('Test',side=2,cex=2.5,line=5)
    legend('topright',c('emp','sim'),col=c('skyblue','red'),lwd=c(6,3),cex=1.75,bty='n')}
}


dev.off()


###################################END#########################################