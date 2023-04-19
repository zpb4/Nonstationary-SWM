#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)

var_mat_hist<-readRDS('data/var_mat_hist.rds')
var_mat_4c<-readRDS('data/var_mat_4c.rds')

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
#fit RF model for error correction
predictor_mat_hist<-readRDS('fit/predictor_mat_hist.rds')
predictor_mat_4c<-readRDS('fit/predictor_mat_4c.rds')

err_hist<-readRDS('fit/err_hist.rds')
err_4c<-readRDS('fit/err_4c.rds')

rg_db_def<-ranger(x=predictor_mat_hist[idx_cal,],y=sm_err_vec,classification = F,
                  replace=T,importance = 'impurity',splitrule = 'variance',quantreg = F)

saveRDS(rg_db_def,'fit/rf_corr_hymod_cal_87_97_default.rds')

#GL-SEP model for residuals
db_cal<-predict(rg_db,data=predictor_mat_hist[idx_cal,])$predictions
err_db_cal<-var_mat_hist[idx_cal,1]-db_cal

db_val<-predict(rg_db,data=predictor_mat_hist[idx_val,])$predictions
err_db_val<-var_mat_hist[idx_val,1]-db_val

#scaling factor validation
pred_scale_v<-scale(predictor_mat_hist[idx_val,1:10],center=F)
sc_est<-predictor_mat_hist[idx_val,1:10]/pred_scale_v
sc_factor<-apply(sc_est,2,function(x){mean(x,na.rm=T)})
saveRDS(sc_factor,'fit/sc_factor_val.rds')

pred_scale_val<-pred_scale_v
for(k in 1:10){
  pred_scale_val[,k]<-(pred_scale_v[,k]-min(0,pred_scale_v[,k]))
}

#scaling factor calibration
pred_scale_c<-scale(predictor_mat_hist[idx_cal,1:10],center=F)
sc_est<-predictor_mat_hist[idx_val,1:10]/pred_scale_c
sc_factor<-apply(sc_est,2,function(x){mean(x,na.rm=T)})
saveRDS(sc_factor,'fit/sc_factor_cal.rds')

pred_scale_cal<-pred_scale_c
for(k in 1:10){
  pred_scale_cal[,k]<-(pred_scale_c[,k]-min(0,pred_scale_c[,k]))
}

#see GL_maineqs and GL_subeqs code for more details
source('GL_maineqs_mv.R')

srt_err<-sort(abs(err_db_val))
srt_err<-srt_err[srt_err>0]
min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])

st_arr<-array(NA,c(4,11))
st_arr[1,]<-c(0.5,rep(0.5,10))
st_arr[2,]<-c(0,rep(0,10))
st_arr[3,]<-c(0,rep(0,10))
st_arr[4,]<-c(0.5,rep(0,10))

lb_arr<-array(NA,c(4,11))
lb_arr[1,]<-c(min_sig,rep(0,10))
lb_arr[2,]<-c(-0.99,rep(-5,10))
lb_arr[3,]<-c(-1,rep(-1,10))
lb_arr[4,]<-c(0,rep(-1,10))

ub_arr<-array(NA,c(4,11))
ub_arr[1,]<-c(5,rep(1,10))
ub_arr[2,]<-c(5,rep(5,10))
ub_arr[3,]<-c(1,rep(1,10))
ub_arr[4,]<-c(1,rep(1,10))

start<-as.vector(st_arr)
lb<-as.vector(lb_arr)
ub<-as.vector(ub_arr)

gl_mle_val<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val,
    lower = lb,upper = ub,method = 'Rvmmin',
    control = list(maximize=T,all.methods=F,maxit=10000))

gl_mle_cal<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_cal,beta_var=pred_scale_cal,xi_var=pred_scale_cal,phi_var=pred_scale_cal,et=err_db_cal,
    lower = lb,upper = ub,method = 'Rvmmin',
    control = list(maximize=T,all.methods=F))

gl_mle_cal$value
gl_mle_val$value

#vectorize validation parameters
par_val<-c()

for(i in 1:44){
  par_val[i]<-gl_mle_val[[i]][5]
}
par_val

GL_fun_mv_ar1_lin(par_val,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val)

#vectorize calibration parameters
par_cal<-c()

for(i in 1:44){
  par_cal[i]<-gl_mle_cal[[i]][1]
}

par_cal

GL_fun_mv_ar1_lin(par_cal,sig_var=pred_scale_cal,beta_var=pred_scale_cal,xi_var=pred_scale_cal,phi_var=pred_scale_cal,et=err_db_cal)

saveRDS(par_val,'fit/hymod-resid-fit_mv-lin-ar_val_default.rds')
saveRDS(par_cal,'fit/hymod-resid-fit_mv-lin-ar_cal_default.rds')

rm(list=ls());gc()

##########################################END#################################################

