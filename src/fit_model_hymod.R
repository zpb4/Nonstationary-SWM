#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)
library(mco)

#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
seed_length<-10
hym_site<-'ORO'
sma_site<-'ORO'

ll_vec<-c()
noise_ll_vec<-c()

#load predictor arrays
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))

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

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

#set.seed(seed)
#cal_yrs<-sort(sample(1988:2003,10))
#val_yrs<-1988:2003;val_yrs<-val_yrs[!(val_yrs%in%cal_yrs)==T]

#idx_cal<-c()
#for(i in 1:length(cal_yrs)){
  #idx_cal<-c(idx_cal,which(ix==paste(cal_yrs[i],'-10-01',sep='')):which(ix==paste(cal_yrs[i]+1,'-09-30',sep='')))
#}

#idx_val<-c()
#for(i in 1:length(val_yrs)){
  #idx_val<-c(idx_val,which(ix==paste(val_yrs[i],'-10-01',sep='')):which(ix==paste(val_yrs[i]+1,'-09-30',sep='')))
#}

#obs_cal<-sma_hist_vars[idx_cal,'sim']

#wts<-obs_cal/max(obs_cal)

#1) Prepare data
rf_err<-hym_predmat_hist[,'err 0']

rf_pred<-hym_predmat_hist[,rf_idx]

for(s in 1:seed_length){
seed<-s
set.seed(seed)
#2) Fit Random Forest error correction model


#ntree_grid<-rep(seq(10,500,10),each=12)
#mtry_grid<-rep(1:12,length(seq(10,500,10)))
#pred_mn<-c()
#for(i in 1:length(ntree_grid)){
  #rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],num.trees = ntree_grid[i],mtry=mtry_grid[i],importance = 'impurity')
  #db_val<-predict(rf_err_corr,data=rf_pred[idx_val,])$predictions
  #err_db_val<-rf_err[idx_val]-db_val
  #pred_mn[i]<-mean(err_db_val)
#}
#s<-sort(abs(pred_mn),index.return=T)
#head(ntree_grid[s$ix],10)
#head(mtry_grid[s$ix],10)
#head(s$x,20)

#rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],num.trees = 30,mtry=4,importance = 'impurity')
#rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],importance = 'impurity',case.weights = wts)

rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],importance = 'impurity')
saveRDS(rf_err_corr,paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#3) Fit multivariate, dynamic GL-SEP model for residuals
#debias calibration errors with RF error correction model
db_cal<-predict(rf_err_corr,data=rf_pred[idx_cal,])$predictions
err_db_cal<-rf_err[idx_cal]-db_cal

#debias calibration errors with RF error correction model
db_val<-predict(rf_err_corr,data=rf_pred[idx_val,])$predictions
err_db_val<-rf_err[idx_val]-db_val

#see GL_maineqs and GL_subeqs code for more details
source('GL_maineqs_rev1.R')

#define lower bound for sigma intercept
srt_err<-sort(abs(err_db_val))
min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
max_sig<-max(abs(err_db_val))

#define bounds for sigma [par1], beta [par2], xi [par3], and phi [par4]
st_arr<-array(NA,c(4,(length(res_idx)+1)))
st_arr[1,]<-c(0.5,rep(0.5,length(res_idx)))
st_arr[2,]<-c(0,rep(0,length(res_idx)))
st_arr[3,]<-c(0,rep(0,length(res_idx)))
st_arr[4,]<-c(0.5,rep(0,length(res_idx)))

lb_arr<-array(NA,c(4,(length(res_idx)+1)))
lb_arr[1,]<-c(min_sig,rep(0,length(res_idx)))
lb_arr[2,]<-c(-0.99,rep(-5,length(res_idx)))
lb_arr[3,]<-c(-1,rep(-1,length(res_idx)))
lb_arr[4,]<-c(0,rep(-1,length(res_idx)))

ub_arr<-array(NA,c(4,(length(res_idx)+1)))
ub_arr[1,]<-c(max_sig,rep(max_sig,length(res_idx)))
ub_arr[2,]<-c(5,rep(5,length(res_idx)))
ub_arr[3,]<-c(1,rep(1,length(res_idx)))
ub_arr[4,]<-c(1,rep(1,length(res_idx)))

start<-as.vector(st_arr)
lb<-as.vector(lb_arr)
ub<-as.vector(ub_arr)

#normalize predictors
pred_mat<-hym_predmat_hist[,res_idx]
pred_mat_scale<-scale(pred_mat[idx_trn,])
ctr_vec<-attributes(pred_mat_scale)['scaled:center']$'scaled:center'
scale_vec<-attributes(pred_mat_scale)['scaled:scale']$'scaled:scale'

#save normalization vector to reuse on Test and Test+4C data
norm_vec<-rbind(ctr_vec,scale_vec)
saveRDS(norm_vec,paste('fit_rev1/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))

#check with manual calculation
pred_mat_comp<-(pred_mat[idx_trn,]-matrix(rep(ctr_vec,length(idx_trn)),ncol=dim(pred_mat)[2],byrow=T))/matrix(rep(scale_vec,length(idx_trn)),ncol=dim(pred_mat)[2],byrow=T)
max(abs(pred_mat_comp-pred_mat_scale)) #should return 0
noise<-rep(0,length(err_db_val))

pred_mat_zero<-apply(pred_mat_scale,2,min)
saveRDS(pred_mat_zero,paste('fit_rev1/hymod_pred-mat-zero_val_',hym_site,'.rds',sep=''))

pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_trn)),ncol=dim(pred_mat_scale)[2],byrow=T)

dyn_res_preds<-pred_mat_scale[idx_val,]
dyn_res_preds_zero<-pred_mat_zero_min[idx_val,]

#Constrained optimization for MLE
dyn_res_mle<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=noise,neg=F,
    lower = lb,upper = ub,method = 'Rvmmin',
    control = list(maximize=T,all.methods=F))


#opts<-c('lbfgsb','nlminb','spg','Rcgmin','Rvmmin','bobyqa','nmkb','hjkb')

#retrieve parameters
dyn_res_coef<-c()
for(i in 1:length(start)){
  dyn_res_coef[i]<-dyn_res_mle[[i]][1]
}
dyn_res_coef

coef_mat<-matrix(dyn_res_coef,nrow=4,byrow=F)
colnames(coef_mat)<-c('intcpt',res_idx)
rownames(coef_mat)<-c('sigma','beta','xi','phi')

ll<-GL_fun_mv_ar1_lin(dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=noise,neg=F)
param_out<-GL_fun_mv_ar1_lin_params(dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val)

print(ll)
print(coef_mat)
ll_vec[s]<-ll

saveRDS(dyn_res_coef,paste('fit_rev1/hymod_dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(coef_mat,paste('fit_rev1/hymod_dyn-res-coef-mat_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(param_out,paste('fit_rev1/hymod_dyn-res-params_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#----------------------------------------------------------------------------------
#as above but with noise regularization

#inputs/predictand -> z* = z + z*norm(0,sd) [Klotz et al. 2022, p16]
#recommended setting [Rothsfuss et al. 2019, p15] sd_x=0.2, sd_y=0.1 

#Constrained optimization for MLE
noise_pred_mat<-pred_mat_scale + pred_mat_scale * matrix(rnorm(length(idx_trn)*length(res_idx),sd=0.2),nrow=length(idx_trn))
noise_err_db_val<-rnorm(length(err_db_val),sd=0.1)

noise_pred_mat_zero<-apply(noise_pred_mat,2,min)
saveRDS(noise_pred_mat_zero,paste('fit_rev1/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

noise_pred_mat_zero_min<-noise_pred_mat-matrix(rep(noise_pred_mat_zero,length(idx_trn)),ncol=dim(noise_pred_mat)[2],byrow=T)

dyn_res_preds<-noise_pred_mat[idx_val,]
dyn_res_preds_zero<-noise_pred_mat_zero_min[idx_val,]

noise_dyn_res_mle<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,
                          xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=noise_err_db_val,neg=F,
                          lower = lb,upper = ub,method = 'Rvmmin',
                          control = list(maximize=T,all.methods=F))


#opts<-c('lbfgsb','nlminb','spg','Rcgmin','Rvmmin','bobyqa','nmkb','hjkb')

#retrieve parameters
noise_dyn_res_coef<-c()
for(i in 1:length(start)){
  noise_dyn_res_coef[i]<-noise_dyn_res_mle[[i]][1]
}
noise_dyn_res_coef

noise_coef_mat<-matrix(noise_dyn_res_coef,nrow=4,byrow=F)
colnames(noise_coef_mat)<-c('intcpt',res_idx)
rownames(noise_coef_mat)<-c('sigma','beta','xi','phi')

noise_ll<-GL_fun_mv_ar1_lin(noise_dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=noise_err_db_val,neg=F)
noise_param_out<-GL_fun_mv_ar1_lin_params(noise_dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=noise_err_db_val)

print(noise_ll)
print(noise_coef_mat)
noise_ll_vec[s]<-noise_ll

saveRDS(noise_dyn_res_coef,paste('fit_rev1/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(noise_coef_mat,paste('fit_rev1/hymod_noise-dyn-res-coef-mat_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(noise_param_out,paste('fit_rev1/hymod_noise-dyn-res-params_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
}

saveRDS(ll_vec,paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
saveRDS(noise_ll_vec,paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))

rm(list=ls());gc()

##########################################END#################################################
