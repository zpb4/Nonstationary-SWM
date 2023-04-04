#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

seq_cal<-seq(1987,2013,3)
ix_cal<-c()
ixx_cal<-c()
idx_cal<-c()
for(i in 1:length(seq_cal)){
  ix_cal<-as.Date(c(ix_cal,seq(as.Date(paste(seq_cal[i],'-10-01',sep='')),as.Date(paste(seq_cal[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
  ixx_cal<-c(ixx_cal,as.POSIXlt(ix_cal))
  idx_cal<-c(idx_cal,which(ix==paste(seq_cal[i],'-10-01',sep='')):which(ix==paste(seq_cal[i]+1,'-09-30',sep='')))
}

seq_val<-seq(1988,2013,3)
ix_val<-c()
ixx_val<-c()
idx_val<-c()
for(i in 1:length(seq_val)){
  ix_val<-as.Date(c(ix_val,seq(as.Date(paste(seq_val[i],'-10-01',sep='')),as.Date(paste(seq_val[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
  ixx_val<-c(ixx_val,as.POSIXlt(ix_val))
  idx_val<-c(idx_val,which(ix==paste(seq_val[i],'-10-01',sep='')):which(ix==paste(seq_val[i]+1,'-09-30',sep='')))
}

seq_tst<-seq(1989,2012,3)
ix_tst<-c()
ixx_tst<-c()
idx_tst<-c()
for(i in 1:length(seq_tst)){
  ix_tst<-as.Date(c(ix_tst,seq(as.Date(paste(seq_tst[i],'-10-01',sep='')),as.Date(paste(seq_tst[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
  ixx_tst<-c(ixx_tst,as.POSIXlt(ix_tst))
  idx_tst<-c(idx_tst,which(ix==paste(seq_tst[i],'-10-01',sep='')):which(ix==paste(seq_tst[i]+1,'-09-30',sep='')))
}

#fit RF model for error correction
predictor_mat_hist<-readRDS('fit/predictor_mat_hist_sacsma.rds')
err_hist<-readRDS('fit/err_hist_sacsma.rds')

#smoothing as needed
sm_err_vec<-ksmooth(1:length(idx_cal),err_hist[idx_cal],bandwidth = 1)$y

rg_db_def<-ranger(x=predictor_mat_hist[idx_cal,],y=sm_err_vec,classification = F,
                  replace=T,importance = 'impurity',splitrule = 'variance',quantreg = F)

saveRDS(rg_db_def,'fit/rf_corr_sacsma_cal_skip-samp_def.rds')

#fit AR model to corrected errors
db_val<-predict(rg_db_def,data=predictor_mat_hist[idx_val,])$predictions
err_db_val<-err_hist[idx_val]-db_val

pred_scale<-scale(predictor_mat_hist[idx_val,1:11],center=F)
sc_est<-predictor_mat_hist[idx_val,1:11]/pred_scale

sc_factor<-apply(sc_est,2,function(x){mean(x,na.rm=T)})

saveRDS(sc_factor,'fit/sacsma_sc_factor_val_skip-samp.rds')

pred_scale_val<-pred_scale

for(k in 1:11){
  pred_scale_val[,k]<-(pred_scale[,k]-min(0,pred_scale[,k]))
}

#see GL_maineqs and GL_subeqs code for more details
source('GL_maineqs_mv.R')

srt_err<-sort(abs(err_db_val))
srt_err<-srt_err[srt_err>0]
min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])

st_arr<-array(NA,c(4,12))
st_arr[1,]<-c(0.5,rep(0.5,11))
st_arr[2,]<-c(0,rep(0,11))
st_arr[3,]<-c(0,rep(0,11))
st_arr[4,]<-c(0.5,rep(0,11))

lb_arr<-array(NA,c(4,12))
lb_arr[1,]<-c(min_sig,rep(0,11))
#lb_arr[1,]<-c(min_sig,rep(-1,11))
lb_arr[2,]<-c(-0.99,rep(-5,11))
lb_arr[3,]<-c(-1,rep(-1,11))
lb_arr[4,]<-c(0,rep(-1,11))

ub_arr<-array(NA,c(4,12))
ub_arr[1,]<-c(5,rep(1,11))
ub_arr[2,]<-c(5,rep(5,11))
ub_arr[3,]<-c(1,rep(1,11))
ub_arr[4,]<-c(1,rep(1,11))

start<-as.vector(st_arr)
lb<-as.vector(lb_arr)
ub<-as.vector(ub_arr)

gl_mle_val<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val,
                   lower = lb,upper = ub,method = 'Rvmmin',
                   control = list(maximize=T,all.methods=F))

#opts<-c('lbfgsb','nlminb','spg','Rcgmin','Rvmmin','bobyqa','nmkb','hjkb')
gl_mle_val$value

#vectorize validation parameters
par_val<-c()

for(i in 1:48){
  par_val[i]<-gl_mle_val[[i]][1]
}
par_val

GL_fun_mv_ar1_lin(par_val,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val)

saveRDS(par_val,'fit/sacsma-resid-fit_mv-lin-ar_val_skip-samp_def.rds')

rm(list=ls());gc()
##########################################END#################################################

