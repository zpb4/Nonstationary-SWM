#fit model between hymod and sacsma
#setwd('z:/oro_nonstat/')
library(fGarch)

print(paste('start',Sys.time()))

var_mat_hist<-readRDS('data/var_mat_hist.rds')
var_mat_4c<-readRDS('data/var_mat_4c.rds')

predictor_mat_hist<-readRDS('fit/predictor_mat_hist.rds')
predictor_mat_4c<-readRDS('fit/predictor_mat_4c.rds')

err_hist<-readRDS('fit/err_hist.rds')
err_4c<-readRDS('fit/err_4c.rds')

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

ix_cal<-seq(as.Date('1987-10-01'),as.Date('1997-09-30'),'day')
ixx_cal<-as.POSIXlt(ix_cal)

ix_val<-seq(as.Date('1997-10-01'),as.Date('2003-09-30'),'day')
ixx_val<-as.POSIXlt(ix_val)

ix_hst<-seq(as.Date('1987-10-01'),as.Date('2003-09-30'),'day')
ixx_hst<-as.POSIXlt(ix_hst)

ix_tst<-seq(as.Date('2003-10-01'),as.Date('2013-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)

ix_4c<-seq(as.Date('2003-10-01'),as.Date('2013-09-30'),'day')
ixx_4c<-as.POSIXlt(ix_4c)

idx_cal<-which(ix=='1987-10-01'):which(ix=='1997-09-30')
idx_val<-which(ix=='1997-10-01'):which(ix=='2003-09-30')
idx_hst<-which(ix=='1987-10-01'):which(ix=='2003-09-30')
idx_tst<-which(ix=='2003-10-01'):which(ix=='2013-09-30')
idx_4c<-which(ix=='2003-10-01'):which(ix=='2013-09-30')

#fit benchmark model for error correction
err_hist<-var_mat_hist[idx_hst,1]
err_hist_db<-err_hist

err_mns<-c()

for(i in 1:12){
  seas<-which(ixx_hst$mon==(i-1))
  err_mns[i]<-mean(err_hist[seas])
  err_hist_db[seas]<-err_hist[seas]-rep(err_mns[i],length(seas))
}

saveRDS(err_mns,'fit/hymod_benchmark_err-mns.rds')

#GL-SGED model
sim_hist<-predictor_mat_hist[idx_hst,1]
sim_tst<-predictor_mat_hist[idx_tst,1]
sim_4c<-predictor_mat_4c[,1]

source('GL_maineqs_uv.R')

lb<-c(0.0001,0,0,0,0,-0.99,0.1) #lower bounds for sig0, sig1, phi1, phi2, phi3, beta, xi; -1 for beta gives an error
ub<-c(20,5,1,1,1,1,10) #upper bounds for sig0, sig1, phi1, phi2, phi3, beta, xi
st<-c(.5,.5,.5,.5,.5,0,1) #starting parameters for sig0, sig1, phi1, phi2, phi3, beta, xi

gl_par<-array(NA,c(12,7))

for(i in 1:12){
  seas<-which(ixx_hst$mon==(i-1))
  gl_mle<-optim(par=st,fn=GL_fun_noscale_ar3_bfgs,inflow=sim_hist[seas],et=err_hist_db[seas],
              method = 'L-BFGS-B',lower = lb, upper = ub,
              control = list(fnscale=-1,maxit=100000))
  gl_par[i,]<-gl_mle$par
}

saveRDS(gl_par,'fit/hymod_benchmark_gl-par.rds')


n<-1000

syn_err_hist<-array(NA,c(length(ix_hst),n))
syn_flow_hist<-array(NA,c(length(ix_hst),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_hst$mon==(i-1))
    et<-et_syn(sim_hist[seas],sim_hist[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err_hist[seas,m]<-et+rep(err_mns[i],length(seas))
    syn_flow_hist[seas,m]<-sim_hist[seas]-syn_err_hist[seas,m]
  }
}

saveRDS(syn_err_hist,'out/hymod_benchmark_syn-err_hist.rds')
saveRDS(syn_flow_hist,'out/hymod_benchmark_syn-flow_hist.rds')


syn_err_tst<-array(NA,c(length(ix_tst),n))
syn_flow_tst<-array(NA,c(length(ix_tst),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_tst$mon==(i-1))
    et<-et_syn(sim_tst[seas],sim_tst[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err_tst[seas,m]<-et+rep(err_mns[i],length(seas))
    syn_flow_tst[seas,m]<-sim_tst[seas]-syn_err_tst[seas,m]
  }
}

saveRDS(syn_err_tst,'out/hymod_benchmark_syn-err_tst.rds')
saveRDS(syn_flow_tst,'out/hymod_benchmark_syn-flow_tst.rds')


syn_err_4c<-array(NA,c(length(ix_4c),n))
syn_flow_4c<-array(NA,c(length(ix_4c),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_4c$mon==(i-1))
    et<-et_syn(sim_4c[seas],sim_4c[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err_4c[seas,m]<-et+rep(err_mns[i],length(seas))
    syn_flow_4c[seas,m]<-sim_4c[seas]-syn_err_4c[seas,m]
  }
}

saveRDS(syn_err_4c,'out/hymod_benchmark_syn-err_4c.rds')
saveRDS(syn_flow_4c,'out/hymod_benchmark_syn-flow_4c.rds')

print(paste('end',Sys.time()))

#--------------------------------------------------------------

#ar3_fits<-vector('list',12)
#err_hist_db_uc<-err_hist_db

#for(i in 1:12){
  #seas<-which(ixx_hst$mon==(i-1))
  #ar3_fits[[i]]<-arima(err_hist_db[seas],order=c(3,0,0),include.mean = F)
  #err_hist_db_uc[seas]<-ar3_fits[[i]]$residuals
#}

#saveRDS(ar3_fits,'fit/hymod_benchmark_ar3-fits.rds')

#sgd_fits<-array(NA,c(12,4))

#for(i in 1:12){
  #seas<-which(ixx_hst$mon==(i-1))
  #sgd_fits[i,]<-sgedFit(err_hist_db_uc[seas])$par
#}

#saveRDS(sgd_fits,'fit/hymod_benchmark_sged-fits.rds')

#seas_4c<-which(ixx_4c$mon==2)

#syn_bench_err<-array(NA,c(length(seas_4c),n))
#syn_bench_flow<-array(NA,c(length(seas_4c),n))

#for(m in 1:n){
  #innov_new<-rsged(length(seas_4c),mean=sgd_fits[3,1],sd=sgd_fits[3,2],nu=sgd_fits[3,3],xi=sgd_fits[3,4])
  #sim<-arima.sim(ar3_fits,n=length(seas_4c),innov = innov_new)
  #syn_bench_err[,i]<-sim+err_mns[3]
  #syn_bench_flow[,i]<-var_mat_4c[seas_4c,3]-sim
#}

#saveRDS(syn_bench_err,'out/hymod_benchmark_syn-err.rds')
#saveRDS(sgd_fits,'out/hymod_benchmark_syn-flow.rds')

rm(list=ls());gc()

##########################################END#################################################

