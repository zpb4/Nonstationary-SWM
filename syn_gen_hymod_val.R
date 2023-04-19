#setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(doParallel)

print(paste('start',Sys.time()))

parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

n<-100

err_hist<-readRDS('fit/err_hist.rds')
err_4c<-readRDS('fit/err_4c.rds')

predictor_mat_hist<-readRDS('fit/predictor_mat_hist.rds')
predictor_mat_4c<-readRDS('fit/predictor_mat_4c.rds')

rg_db<-readRDS('fit/rf_corr_hymod_cal_87_97.rds')

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

ix_cal<-seq(as.Date('1987-10-01'),as.Date('1997-09-30'),'day')
ixx_cal<-as.POSIXlt(ix_cal)

ix_val<-seq(as.Date('1997-10-01'),as.Date('2003-09-30'),'day')
ixx_val<-as.POSIXlt(ix_val)

ix_tst<-seq(as.Date('2003-10-01'),as.Date('2013-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)

ix_4c<-seq(as.Date('2003-10-01'),as.Date('2013-09-30'),'day')
ixx_4c<-as.POSIXlt(ix_4c)

idx_cal<-which(ix=='1987-10-01'):which(ix=='1997-09-30')
idx_val<-which(ix=='1997-10-01'):which(ix=='2003-09-30')
idx_tst<-which(ix=='2003-10-01'):which(ix=='2013-09-30')
idx_4c<-which(ix=='2003-10-01'):which(ix=='2013-09-30')

#-------------------------------------------------------------------------------
#val
idx<-idx_val

db_err<-predict(rg_db,data=predictor_mat_hist[idx,])$predictions
err_db<-err_hist[idx]-db_err

source('GL_maineqs_mv.R')
gl_fit<-readRDS('fit/gl_mle_par_sqrt2_mv_lin-ar.rds')

sc_factor<-readRDS('fit/sc_factor_val.rds')
sc_factor<-matrix(rep(sc_factor,length(idx)),ncol=11,byrow=T)

pred_in<-predictor_mat_hist[idx,1:11]
pred_scale_in<-pred_in/sc_factor
pred_scale<-pred_scale_in

for(k in 1:11){
  pred_scale[,k]<-(pred_scale_in[,k]-min(0,pred_scale_in[,k]))
}


syn_out<-foreach(m = 1:n,.combine='cbind',.packages=c('fGarch','ranger'),.inorder=F) %dopar% {
  
  syn_res<-syn_gen_mv_ar1_lin(gl_fit,sig_var=pred_scale,beta_var=sqrt(pred_scale),xi_var=sqrt(pred_scale),phi_var=sqrt(pred_scale),et=err_db)
  
  syn_gen<-rep(0,(length(idx)+3));syn_gen[1:3]<-sample(syn_res,3)

  vmat<-predictor_mat_hist[idx,]

  for(i in 4:length(syn_gen)){
    dat<-predictor_mat_hist[1:2,]
    dat[1,]<-c(vmat[(i-3),1:11],syn_gen[i-1],syn_gen[i-2],syn_gen[i-3])
    err<-predict(rg_db,data=t(dat[1,]))$predictions
    syn_gen[i]<-err+syn_res[(i-3)]
  }
  
  return(syn_gen[4:length(syn_gen)])
}

saveRDS(syn_out,'out/hymod_val_gl-sep-mv-lin-ar_syn-error.rds')

sim<-matrix(rep(predictor_mat_hist[idx,1],n),ncol=n,byrow=F)
swm_out<-sim+syn_out
swm_out[swm_out<0]<-0

saveRDS(swm_out,'out/hymod_val_gl-sep-mv-lin-ar_syn-flow.rds')

print(paste('end',Sys.time()))

rm(list=ls());gc()

##############################END##################