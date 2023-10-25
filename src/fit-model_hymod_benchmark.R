#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(optimx)
library(mco)

#specifications
seed<-1
hym_site<-'ORO'

#load predictor arrays
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_hist<-which(ix=='1988-10-01'):which(ix=='2018-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_trn<-seq(as.Date('1988-10-01'),as.Date('2004-09-30'),'day')
ixx_trn<-as.POSIXlt(ix_trn)

ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)

ix_hist<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ixx_hist<-as.POSIXlt(ix_hist)

#fit benchmark model for error correction
err_hist<-hym_predmat_hist[idx_trn,'err 0']
err_hist_db<-err_hist

err_mns<-c()

for(i in 1:12){
  seas<-which(ixx_trn$mon==(i-1))
  err_mns[i]<-mean(err_hist[seas])
  err_hist_db[seas]<-err_hist[seas]-rep(err_mns[i],length(seas))
}

saveRDS(err_mns,'fit_rev1/hymod_benchmark_err-mns.rds')

#GL-SGED model
sim_trn<-hym_predmat_hist[idx_trn,'sim 0']
sim_hist<-hym_predmat_hist[,'sim 0']
sim_4c<-hym_predmat_4c[,'sim 0']

source('GL_maineqs_uv.R')

lb<-c(0.0001,0,0,0,0,-0.99,0.1) #lower bounds for sig0, sig1, phi1, phi2, phi3, beta, xi; -1 for beta gives an error
ub<-c(20,5,1,1,1,5,10) #upper bounds for sig0, sig1, phi1, phi2, phi3, beta, xi
st<-c(.5,.5,.5,.5,.5,0,1) #starting parameters for sig0, sig1, phi1, phi2, phi3, beta, xi

gl_par<-array(NA,c(12,7))

for(i in 1:12){
  seas<-which(ixx_trn$mon==(i-1))
  et_sort<-sort(abs(err_hist_db[seas]))
  min_sig<-max(0.0001,mean(et_sort[1:round(0.1*length(seas))]))
  lb[1]<-min_sig
  gl_mle<-optim(par=st,fn=GL_fun_noscale_ar3_bfgs,inflow=sim_hist[seas],et=err_hist_db[seas],
              method = 'L-BFGS-B',lower = lb, upper = ub,
              control = list(fnscale=-1,maxit=100000))
  gl_par[i,]<-gl_mle$par
}

saveRDS(gl_par,'fit_rev1/hymod_benchmark_gl-par.rds')

#hist
n<-1000

syn_err_hist<-array(NA,c(length(ix_hist),n))
syn_flow_hist<-array(NA,c(length(ix_hist),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_hist$mon==(i-1))
    et<-et_syn(sim_hist[seas],sim_hist[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err_hist[seas,m]<-et+rep(err_mns[i],length(seas))
    syn_flow_hist[seas,m]<-sim_hist[seas]+syn_err_hist[seas,m]
  }
}

saveRDS(syn_err_hist,'out_rev1/hymod_benchmark_syn-err_hist.rds')
saveRDS(syn_flow_hist,'out_rev1/hymod_benchmark_syn-flow_hist.rds')

#4c
n<-1000

syn_err_4c<-array(NA,c(length(ix_hist),n))
syn_flow_4c<-array(NA,c(length(ix_hist),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_hist$mon==(i-1))
    et<-et_syn(sim_4c[seas],sim_4c[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err_4c[seas,m]<-et+rep(err_mns[i],length(seas))
    syn_flow_4c[seas,m]<-sim_4c[seas]+syn_err_4c[seas,m]
  }
}

saveRDS(syn_err_4c,'out_rev1/hymod_benchmark_syn-err_4c.rds')
saveRDS(syn_flow_4c,'out_rev1/hymod_benchmark_syn-flow_4c.rds')

rm(list=ls());gc()

##########################################END#################################################

