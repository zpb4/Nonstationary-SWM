#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

ix_cal<-seq(as.Date('1987-10-01'),as.Date('1997-09-30'),'day')
ixx_cal<-as.POSIXlt(ix_cal)

ix_val<-seq(as.Date('1997-10-01'),as.Date('2003-09-30'),'day')
ixx_val<-as.POSIXlt(ix_val)

ix_4c<-seq(as.Date('2003-10-01'),as.Date('2013-09-30'),'day')
ixx_4c<-as.POSIXlt(ix_4c)

idx_cal<-which(ix=='1987-10-01'):which(ix=='1997-09-30')
idx_val<-which(ix=='1997-10-01'):which(ix=='2003-09-30')
idx_4c<-which(ix=='2003-10-01'):which(ix=='2013-09-30')

#fit RF model for error correction
predictor_mat_hist<-readRDS('fit/predictor_mat_hist_sacsma.rds')
err_hist<-readRDS('fit/err_hist_sacsma.rds')

#smoothing as needed
sm_err_vec<-ksmooth(1:length(idx_cal),err_hist[idx_cal],bandwidth = 1)$y

rg_db_def<-ranger(x=predictor_mat_hist[idx_cal,],y=sm_err_vec,classification = F,
                  replace=T,importance = 'impurity',splitrule = 'variance',quantreg = F)

saveRDS(rg_db_def,'fit/rf_corr_sacsma_cal_split-samp_def.rds')

#fit AR model to corrected errors
db_val<-predict(rg_db_def,data=predictor_mat_hist[idx_val,])$predictions
err_db_val<-err_hist[idx_val]-db_val

pred_scale<-scale(predictor_mat_hist[idx_val,1:11],center=F)
sc_est<-predictor_mat_hist[idx_val,1:11]/pred_scale

sc_factor<-apply(sc_est,2,function(x){mean(x,na.rm=T)})

saveRDS(sc_factor,'fit/sacsma_sc_factor_val_split-samp.rds')

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

saveRDS(par_val,'fit/sacsma-resid-fit_mv-lin-ar_val_split-samp_def.rds')

rm(list=ls());gc()
##########################################END#################################################

