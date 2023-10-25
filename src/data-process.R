#Process and arrange data into labeled matrices

setwd('z:/oro_nonstat/')

sma_site<-'ORO'
hym_site<-'ORO'

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

#simulation data
sma_hist<-read.table(paste('data_rev1/sacsma_',sma_site,'.txt',sep=''))
hym_hist<-read.table(paste('data_rev1/hymod_',hym_site,'.txt',sep=''))
sma_4c<-read.table(paste('data_rev1/sacsma_',sma_site,'_4C.txt',sep=''))
hym_4c<-read.table(paste('data_rev1/hymod_',hym_site,'_4C.txt',sep=''))

if(sma_site=='ORO'&hym_site=='ORO'){
  sma_hist<-read.table(paste('data_rev1/sacsma_',sma_site,'_NEW.txt',sep=''))
  hym_hist<-read.table(paste('data_rev1/hymod_',hym_site,'_NEW.txt',sep=''))
  sma_4c<-read.table(paste('data_rev1/sacsma_',sma_site,'_4C_NEW.txt',sep=''))
  hym_4c<-read.table(paste('data_rev1/hymod_',hym_site,'_4C_NEW.txt',sep=''))
}

#hymod output columns
cnames_hym<-c('year','mo','day','precip','tavg','obs*','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm')
cnames_sma<-c('year','mo','day','precip','tavg','obs','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store')

#obs*: HYMOD obs are actually the SACSMA flows (truth model)

#store raw data - HYMOD
hym_hist_vars<-hym_hist
hym_4c_vars<-hym_4c

colnames(hym_hist_vars)<-cnames_hym
colnames(hym_4c_vars)<-cnames_hym
rownames(hym_hist_vars)<-as.character(ix)
rownames(hym_4c_vars)<-as.character(ix)

saveRDS(hym_hist_vars,paste('data_rev1/hym_hist_vars_',hym_site,'.rds',sep=''))
saveRDS(hym_4c_vars,paste('data_rev1/hym_4c_vars_',hym_site,'.rds',sep=''))

#store raw data - SACSMA
sma_hist_vars<-sma_hist
sma_4c_vars<-sma_4c

colnames(sma_hist_vars)<-cnames_sma
colnames(sma_4c_vars)<-cnames_sma
rownames(sma_hist_vars)<-as.character(ix)
rownames(sma_4c_vars)<-as.character(ix)

obs_correct_hist<-sma_hist_vars[,'obs']
obs_correct_hist[obs_correct_hist<0]<-0
sma_hist_vars[,'obs']<-obs_correct_hist

obs_correct_4c<-sma_4c_vars[,'obs']
obs_correct_4c[obs_correct_4c<0]<-0
sma_4c_vars[,'obs']<-obs_correct_4c

saveRDS(sma_hist_vars,paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))
saveRDS(sma_4c_vars,paste('data_rev1/sma_4c_vars_',sma_site,'.rds',sep=''))

#----------------------------------------------------
#calculate lagged SWE and precip averages
avg_upto<-function(x,wdow){
  out<-c()
  for(i in 1:length(x)){
    samp<-max(1,(i-wdow)):min(length(x),i)
    out[i]<-mean(x[samp])}
  return(out)}

trend<-function(x,wdow){
  out<-c()
  for(i in 1:length(x)){
    samp<-max(1,(i-wdow)):min(length(x),i)
    out[i]<-lm(x[samp]~c(1:length(samp)))$coefficients[2]
    if(is.na(out[i])==T){out[i]<-x[samp][length(samp)]-x[samp][1]}
  }
  return(out)}

#--------------------------------------------------
##HYMOD
#long lag average from t-window to t.
window<-360

hym_hist_swe_avg<-avg_upto(hym_hist_vars[,'swe'],window)
hym_hist_et_avg<-avg_upto(hym_hist_vars[,'et'],window)
hym_hist_tavg_avg<-avg_upto(hym_hist_vars[,'tavg'],window)
hym_hist_sim_avg<-avg_upto(hym_hist_vars[,'sim'],window)
hym_hist_bflow_avg<-avg_upto(hym_hist_vars[,'baseflow'],window)
hym_hist_roff_avg<-avg_upto(hym_hist_vars[,'runoff'],window)
hym_hist_lsm_avg<-avg_upto(hym_hist_vars[,'lwr_sm'],window)
hym_hist_usm_avg<-avg_upto(hym_hist_vars[,'upr_sm'],window)

hym_hist_avg_llag<-cbind(hym_hist_tavg_avg,hym_hist_sim_avg,hym_hist_roff_avg,hym_hist_bflow_avg,hym_hist_et_avg,hym_hist_swe_avg,hym_hist_usm_avg,hym_hist_lsm_avg)

hym_4c_swe_avg<-avg_upto(hym_4c_vars[,'swe'],window)
hym_4c_et_avg<-avg_upto(hym_4c_vars[,'et'],window)
hym_4c_tavg_avg<-avg_upto(hym_4c_vars[,'tavg'],window)
hym_4c_sim_avg<-avg_upto(hym_4c_vars[,'sim'],window)
hym_4c_bflow_avg<-avg_upto(hym_4c_vars[,'baseflow'],window)
hym_4c_roff_avg<-avg_upto(hym_4c_vars[,'runoff'],window)
hym_4c_lsm_avg<-avg_upto(hym_4c_vars[,'lwr_sm'],window)
hym_4c_usm_avg<-avg_upto(hym_4c_vars[,'upr_sm'],window)

hym_4c_avg_llag<-cbind(hym_4c_tavg_avg,hym_4c_sim_avg,hym_4c_roff_avg,hym_4c_bflow_avg,hym_4c_et_avg,hym_4c_swe_avg,hym_4c_usm_avg,hym_4c_lsm_avg)

#t-window to t linear model trend in variable
window<-30

hym_hist_swe_avg<-trend(hym_hist_vars[,'swe'],window)
hym_hist_et_avg<-trend(hym_hist_vars[,'et'],window)
hym_hist_tavg_avg<-trend(hym_hist_vars[,'tavg'],window)
hym_hist_sim_avg<-trend(hym_hist_vars[,'sim'],window)
hym_hist_bflow_avg<-trend(hym_hist_vars[,'baseflow'],window)
hym_hist_roff_avg<-trend(hym_hist_vars[,'runoff'],window)
hym_hist_lsm_avg<-trend(hym_hist_vars[,'lwr_sm'],window)
hym_hist_usm_avg<-trend(hym_hist_vars[,'upr_sm'],window)

hym_hist_trend_llag<-cbind(hym_hist_tavg_avg,hym_hist_sim_avg,hym_hist_roff_avg,hym_hist_bflow_avg,hym_hist_et_avg,hym_hist_swe_avg,hym_hist_usm_avg,hym_hist_lsm_avg)

hym_4c_swe_avg<-trend(hym_4c_vars[,'swe'],window)
hym_4c_et_avg<-trend(hym_4c_vars[,'et'],window)
hym_4c_tavg_avg<-trend(hym_4c_vars[,'tavg'],window)
hym_4c_sim_avg<-trend(hym_4c_vars[,'sim'],window)
hym_4c_bflow_avg<-trend(hym_4c_vars[,'baseflow'],window)
hym_4c_roff_avg<-trend(hym_4c_vars[,'runoff'],window)
hym_4c_lsm_avg<-trend(hym_4c_vars[,'lwr_sm'],window)
hym_4c_usm_avg<-trend(hym_4c_vars[,'upr_sm'],window)

hym_4c_trend_llag<-cbind(hym_4c_tavg_avg,hym_4c_sim_avg,hym_4c_roff_avg,hym_4c_bflow_avg,hym_4c_et_avg,hym_4c_swe_avg,hym_4c_usm_avg,hym_4c_lsm_avg)

#-----------------------------------------
##SACSMA
#long lag average from t-window to t.
window<-360

sma_hist_swe_avg<-avg_upto(sma_hist_vars[,'swe'],window)
sma_hist_et_avg<-avg_upto(sma_hist_vars[,'et'],window)
sma_hist_tavg_avg<-avg_upto(sma_hist_vars[,'tavg'],window)
sma_hist_sim_avg<-avg_upto(sma_hist_vars[,'sim'],window)
sma_hist_bflow_avg<-avg_upto(sma_hist_vars[,'baseflow'],window)
sma_hist_roff_avg<-avg_upto(sma_hist_vars[,'runoff'],window)
sma_hist_lsm_avg<-avg_upto(sma_hist_vars[,'lwr_sm'],window)
sma_hist_usm_avg<-avg_upto(sma_hist_vars[,'upr_sm'],window)
sma_hist_stor_avg<-avg_upto(sma_hist_vars[,'store'],window)

sma_hist_avg_llag<-cbind(sma_hist_tavg_avg,sma_hist_sim_avg,sma_hist_roff_avg,sma_hist_bflow_avg,sma_hist_et_avg,sma_hist_swe_avg,sma_hist_usm_avg,sma_hist_lsm_avg,sma_hist_stor_avg)

sma_4c_swe_avg<-avg_upto(sma_4c_vars[,'swe'],window)
sma_4c_et_avg<-avg_upto(sma_4c_vars[,'et'],window)
sma_4c_tavg_avg<-avg_upto(sma_4c_vars[,'tavg'],window)
sma_4c_sim_avg<-avg_upto(sma_4c_vars[,'sim'],window)
sma_4c_bflow_avg<-avg_upto(sma_4c_vars[,'baseflow'],window)
sma_4c_roff_avg<-avg_upto(sma_4c_vars[,'runoff'],window)
sma_4c_lsm_avg<-avg_upto(sma_4c_vars[,'lwr_sm'],window)
sma_4c_usm_avg<-avg_upto(sma_4c_vars[,'upr_sm'],window)
sma_4c_stor_avg<-avg_upto(sma_4c_vars[,'store'],window)

sma_4c_avg_llag<-cbind(sma_4c_tavg_avg,sma_4c_sim_avg,sma_4c_roff_avg,sma_4c_bflow_avg,sma_4c_et_avg,sma_4c_swe_avg,sma_4c_usm_avg,sma_4c_lsm_avg,sma_4c_stor_avg)

#t-window to t linear model trend in variable
window<-30

sma_hist_swe_avg<-trend(sma_hist_vars[,'swe'],window)
sma_hist_et_avg<-trend(sma_hist_vars[,'et'],window)
sma_hist_tavg_avg<-trend(sma_hist_vars[,'tavg'],window)
sma_hist_sim_avg<-trend(sma_hist_vars[,'sim'],window)
sma_hist_bflow_avg<-trend(sma_hist_vars[,'baseflow'],window)
sma_hist_roff_avg<-trend(sma_hist_vars[,'runoff'],window)
sma_hist_lsm_avg<-trend(sma_hist_vars[,'lwr_sm'],window)
sma_hist_usm_avg<-trend(sma_hist_vars[,'upr_sm'],window)
sma_hist_stor_avg<-trend(sma_hist_vars[,'store'],window)

sma_hist_trend_llag<-cbind(sma_hist_tavg_avg,sma_hist_sim_avg,sma_hist_roff_avg,sma_hist_bflow_avg,sma_hist_et_avg,sma_hist_swe_avg,sma_hist_usm_avg,sma_hist_lsm_avg,sma_hist_stor_avg)

sma_4c_swe_avg<-trend(sma_4c_vars[,'swe'],window)
sma_4c_et_avg<-trend(sma_4c_vars[,'et'],window)
sma_4c_tavg_avg<-trend(sma_4c_vars[,'tavg'],window)
sma_4c_sim_avg<-trend(sma_4c_vars[,'sim'],window)
sma_4c_bflow_avg<-trend(sma_4c_vars[,'baseflow'],window)
sma_4c_roff_avg<-trend(sma_4c_vars[,'runoff'],window)
sma_4c_lsm_avg<-trend(sma_4c_vars[,'lwr_sm'],window)
sma_4c_usm_avg<-trend(sma_4c_vars[,'upr_sm'],window)
sma_4c_stor_avg<-trend(sma_4c_vars[,'store'],window)

sma_4c_trend_llag<-cbind(sma_4c_tavg_avg,sma_4c_sim_avg,sma_4c_roff_avg,sma_4c_bflow_avg,sma_4c_et_avg,sma_4c_swe_avg,sma_4c_usm_avg,sma_4c_lsm_avg,sma_4c_stor_avg)

#-----------------------------------------------------------------------
#2) Define arrays for predictions

#define predictive errors (truth - process)
#hym_err_hist<-hym_hist_vars[,'obs*']-hym_hist_vars[,'sim']
#hym_err_4c<-hym_4c_vars[,'obs*']-hym_4c_vars[,'sim']

hym_err_hist<-sma_hist_vars[,'sim']-hym_hist_vars[,'sim']
hym_err_4c<-sma_4c_vars[,'sim']-hym_4c_vars[,'sim']

sma_err_hist<-sma_hist_vars[,'obs']-sma_hist_vars[,'sim']
sma_err_4c<-sma_4c_vars[,'obs']-sma_4c_vars[,'sim']

#2a) HYMOD historical
#remove date columns
hym_predmat_hist<-hym_hist_vars[,4:dim(hym_hist_vars)[2]]
#remove obs (truth model) values (not used for prediction)
hym_predmat_hist<-hym_predmat_hist[-c(which(colnames(hym_predmat_hist)=='obs*'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
hym_predmat_hist<-cbind(hym_err_hist,hym_predmat_hist)
labs_hym<-paste(rep(c('err','precip','tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),3),rep(0:(-3),each=dim(hym_predmat_hist)[2]))
lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')
labs_hym<-c(labs_hym,lab_ll_avg,lab_ll_trend)
hym_predmat_hist<-cbind(hym_predmat_hist,
                      rbind(matrix(0,nrow=1,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-1),])),
                      rbind(matrix(0,nrow=2,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-2),])),
                      rbind(matrix(0,nrow=3,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-3),])),
                      hym_hist_avg_llag,
                      hym_hist_trend_llag
                      ) 
#relabel and save
colnames(hym_predmat_hist)<-labs_hym


saveRDS(hym_predmat_hist,paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))

#2b) HYMOD 4C
#remove date columns
hym_predmat_4c<-hym_4c_vars[,4:dim(hym_4c_vars)[2]]
#remove obs (truth model) values (not used for prediction)
hym_predmat_4c<-hym_predmat_4c[-c(which(colnames(hym_predmat_4c)=='obs*'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
hym_predmat_4c<-cbind(hym_err_4c,hym_predmat_4c)
hym_predmat_4c<-cbind(hym_predmat_4c,
                      rbind(matrix(0,nrow=1,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-1),])),
                      rbind(matrix(0,nrow=2,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-2),])),
                      rbind(matrix(0,nrow=3,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-3),])),
                      hym_4c_avg_llag,
                      hym_4c_trend_llag
                      ) 
                      
#relabel and save
colnames(hym_predmat_4c)<-labs_hym

saveRDS(hym_predmat_4c,paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

#2c) SACSMA historical
#remove date columns
sma_predmat_hist<-sma_hist_vars[,4:dim(sma_hist_vars)[2]]
#remove obs (truth model) values (not used for prediction)
sma_predmat_hist<-sma_predmat_hist[-c(which(colnames(sma_predmat_hist)=='obs'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
sma_predmat_hist<-cbind(sma_err_hist,sma_predmat_hist)
labs_sma<-paste(rep(c('err','precip','tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),3),rep(0:(-3),each=dim(sma_predmat_hist)[2]))
lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),'llag-trend')
labs_sma<-c(labs_sma,lab_ll_avg,lab_ll_trend)
sma_predmat_hist<-cbind(sma_predmat_hist,
                        rbind(matrix(0,nrow=1,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-1),])),
                        rbind(matrix(0,nrow=2,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-2),])),
                        rbind(matrix(0,nrow=3,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-3),])),
                        sma_hist_avg_llag,
                        sma_hist_trend_llag
) 

#relabel and save
colnames(sma_predmat_hist)<-labs_sma

saveRDS(sma_predmat_hist,paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))

#2b) SACSMA 4C
#remove date columns
sma_predmat_4c<-sma_4c_vars[,4:dim(sma_4c_vars)[2]]
#remove obs (truth model) values (not used for prediction)
sma_predmat_4c<-sma_predmat_4c[-c(which(colnames(sma_predmat_4c)=='obs'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
sma_predmat_4c<-cbind(sma_err_4c,sma_predmat_4c)
sma_predmat_4c<-cbind(sma_predmat_4c,
                        rbind(matrix(0,nrow=1,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-1),])),
                        rbind(matrix(0,nrow=2,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-2),])),
                        rbind(matrix(0,nrow=3,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-3),])),
                        sma_4c_avg_llag,
                        sma_4c_trend_llag
)
#relabel and save
colnames(sma_predmat_4c)<-labs_sma

saveRDS(sma_predmat_4c,paste('data_rev1/sma_predmat_4c_',hym_site,'.rds',sep=''))

rm(list=ls());gc()

####################################################END#####################################################