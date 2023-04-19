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

n<-1000

err_hist<-readRDS('fit/err_hist_sacsma.rds')

predictor_mat_hist<-readRDS('fit/predictor_mat_hist_sacsma.rds')

rg_db<-readRDS('fit/rf_corr_sacsma_cal_skip-samp_def.rds')

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

#-------------------------------------------------------------------------------
#generate
db_err<-predict(rg_db,data=predictor_mat_hist)$predictions
err_db<-err_hist-db_err

source('GL_maineqs_mv.R')
gl_fit<-readRDS('fit/sacsma-resid-fit_mv-lin-ar_val_skip-samp_def.rds')

sc_factor<-readRDS('fit/sacsma_sc_factor_val_skip-samp.rds')
sc_factor<-matrix(rep(sc_factor,length(ix)),ncol=11,byrow=T)

pred_in<-predictor_mat_hist[,1:11]
pred_scale_in<-pred_in/sc_factor
pred_scale<-pred_scale_in

for(k in 1:11){
  pred_scale[,k]<-(pred_scale_in[,k]-min(0,pred_scale_in[,k]))
}


syn_out<-foreach(m = 1:n,.combine='cbind',.packages=c('fGarch','ranger'),.inorder=F) %dopar% {
  
  syn_res<-syn_gen_mv_ar1_lin(gl_fit,sig_var=pred_scale,beta_var=pred_scale,xi_var=pred_scale,phi_var=pred_scale,et=err_db)
  
  syn_gen<-rep(0,(length(ix)+3));syn_gen[1:3]<-sample(syn_res,3)

  vmat<-predictor_mat_hist

  for(i in 4:length(syn_gen)){
    dat<-predictor_mat_hist[1:2,]
    dat[1,]<-c(vmat[(i-3),1:11],syn_gen[i-1],syn_gen[i-2],syn_gen[i-3])
    err<-predict(rg_db,data=t(dat[1,]))$predictions
    syn_gen[i]<-err+syn_res[(i-3)]
  }
  
  return(syn_gen[4:length(syn_gen)])
}

saveRDS(syn_out,'out/sacsma_skip-samp_gl-sep-mv-lin-ar_syn-error_default.rds')

sim<-matrix(rep(predictor_mat_hist[,1],n),ncol=n,byrow=F)
swm_out<-sim+syn_out
swm_out[swm_out<0]<-0

saveRDS(swm_out,'out/sacsma_skip-samp_gl-sep-mv-lin-ar_syn-flow_default.rds')

print(paste('end',Sys.time()))

rm(list=ls());gc()

##############################END##################