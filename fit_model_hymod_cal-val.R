#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)
seed<-1

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
predictor_mat_hist<-cbind(var_mat_hist[,3:12],c(0,var_mat_hist[-c(length(var_mat_hist[,1])),1]),c(0,0,var_mat_hist[1:(length(var_mat_hist[,1])-2),1]),c(0,0,0,var_mat_hist[1:(length(var_mat_hist[,1])-3),1]))
err_hist<-var_mat_hist[,1]
labs_hist<-c('sim','precip','tavg','et','sm','swe','runoff','baseflow','upr_sm','lwr_sm','lag1','lag2','lag3')
colnames(predictor_mat_hist)<-labs_hist

predictor_mat_4c<-cbind(var_mat_4c[,3:12],c(0,var_mat_4c[-c(length(var_mat_4c[,1])),1]),c(0,0,var_mat_4c[1:(length(var_mat_4c[,1])-2),1]),c(0,0,0,var_mat_4c[1:(length(var_mat_4c[,1])-3),1]))
err_4c<-var_mat_4c[,1]
labs_4c<-c('sim','precip','tavg','et','sm','swe','runoff','baseflow','upr_sm','lwr_sm','lag1','lag2','lag3')
colnames(predictor_mat_4c)<-labs_4c

saveRDS(predictor_mat_hist,'fit/predictor_mat_hist.rds')
saveRDS(predictor_mat_4c,'fit/predictor_mat_4c.rds')

saveRDS(err_hist,'fit/err_hist.rds')
saveRDS(err_4c,'fit/err_4c.rds')

#smoothing as needed
sm_err_vec<-ksmooth(1:length(idx_cal),err_hist[idx_cal],bandwidth = 1)$y

set.seed(seed)
rg_db<-ranger(x=predictor_mat_hist[idx_cal,],y=sm_err_vec,classification = F,
              replace=T,importance = 'impurity',mtry = 5,num.trees = 30,splitrule = 'variance')

saveRDS(rg_db,'fit/rf_corr_hymod_cal_87_97.rds')

#GL-SEP model for residuals
db_cal<-predict(rg_db,data=predictor_mat_hist[idx_cal,])$predictions
err_db_cal<-var_mat_hist[idx_cal,1]-db_cal

db_val<-predict(rg_db,data=predictor_mat_hist[idx_val,])$predictions
err_db_val<-var_mat_hist[idx_val,1]-db_val

pred_scale<-scale(predictor_mat_hist[idx_val,1:10],center=F)
sc_est<-predictor_mat_hist[idx_val,1:10]/pred_scale

sc_factor<-apply(sc_est,2,function(x){mean(x,na.rm=T)})

saveRDS(sc_factor,'fit/sc_factor_val.rds')

pred_scale_val<-pred_scale

for(k in 1:10){
  pred_scale_val[,k]<-(pred_scale[,k]-min(0,pred_scale[,k]))
}

#see GL_maineqs and GL_subeqs code for more details
source('GL_maineqs_mv.R')

srt_err<-sort(abs(err_db_val))
min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])

st_arr<-array(NA,c(4,11))
st_arr[1,]<-c(0.5,rep(0.5,10))
st_arr[2,]<-c(0,rep(0,10))
#st_arr[3,]<-c(1,rep(0,10))
st_arr[3,]<-c(0,rep(0,10))
st_arr[4,]<-c(0.5,rep(0,10))

lb_arr<-array(NA,c(4,11))
lb_arr[1,]<-c(min_sig,rep(0,10))
lb_arr[2,]<-c(-0.99,rep(-5,10))
#lb_arr[3,]<-c(0.1,rep(-2,10))
lb_arr[3,]<-c(-1,rep(-1,10))
lb_arr[4,]<-c(0,rep(-1,10))

ub_arr<-array(NA,c(4,11))
ub_arr[1,]<-c(5,rep(1,10))
ub_arr[2,]<-c(5,rep(5,10))
#ub_arr[3,]<-c(10,rep(2,10))
ub_arr[3,]<-c(1,rep(1,10))
ub_arr[4,]<-c(1,rep(1,10))

start<-as.vector(st_arr)
lb<-as.vector(lb_arr)
ub<-as.vector(ub_arr)

gl_mle_opx<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val,
    lower = lb,upper = ub,method = 'Rvmmin',
    control = list(maximize=T,all.methods=F))
#opts<-c('lbfgsb','nlminb','spg','Rcgmin','Rvmmin','bobyqa','nmkb','hjkb')
#gl_mle_opx$value

par_opx<-c()

for(i in 1:44){
  par_opx[i]<-gl_mle_opx[[i]][1]
}
par_opx

GL_fun_mv_ar1_lin(par_opx,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val)
param_out<-GL_fun_mv_ar1_lin_params(par_opx,sig_var=pred_scale_val,beta_var=pred_scale_val,xi_var=pred_scale_val,phi_var=pred_scale_val,et=err_db_val)
#squrt for beta and xi
gl_mle_opx_sqrt<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_val,beta_var=sqrt(pred_scale_val),xi_var=sqrt(pred_scale_val),phi_var=pred_scale_val,et=err_db_val,
                   lower = lb,upper = ub,#method = 'Rcgmin',
                   control = list(maximize=T,all.methods=T))

par_opx_sqrt<-c()

for(i in 1:44){
  par_opx_sqrt[i]<-gl_mle_opx_sqrt[[i]][4]
}
par_opx_sqrt

GL_fun_mv_ar1_lin(par_opx_sqrt,sig_var=pred_scale_val,beta_var=sqrt(pred_scale_val),xi_var=sqrt(pred_scale_val),phi_var=pred_scale_val,et=err_db_val)

#squrt for beta,xi, phi
gl_mle_opx_sqrt2<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=pred_scale_val,beta_var=sqrt(pred_scale_val),xi_var=sqrt(pred_scale_val),phi_var=sqrt(pred_scale_val),et=err_db_val,
                   lower = lb,upper = ub,#method = 'Rcgmin',
                   control = list(maximize=T,all.methods=T))

par_opx_sqrt2<-c()

for(i in 1:44){
  par_opx_sqrt2[i]<-gl_mle_opx_sqrt2[[i]][5]
}
par_opx_sqrt2

GL_fun_mv_ar1_lin(par_opx_sqrt2,sig_var=pred_scale_val,beta_var=sqrt(pred_scale_val),xi_var=sqrt(pred_scale_val),phi_var=sqrt(pred_scale_val),et=err_db_val)

saveRDS(par_opx,'fit/gl_mle_par_mv_lin-ar.rds')
#saveRDS(par_opx_sqrt,'fit/gl_mle_par_sqrt_mv_lin-ar.rds')
#saveRDS(par_opx_sqrt2,'fit/gl_mle_par_sqrt2_mv_lin-ar.rds')

#rm(list=ls());gc()

##########################################END#################################################

