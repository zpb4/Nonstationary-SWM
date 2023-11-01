#Script to process raw simulated and observed streamflow data and state-variable predictors for 3 sites

#setup path to current directory
current_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)

#specify sites
sma_site<-'ORO'
hym_site<-'ORO'

#date indices
ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

#retrieve raw simulation data
sma_hist<-read.table(paste('../data/sacsma_',sma_site,'.txt',sep=''))
hym_hist<-read.table(paste('../data/hymod_',hym_site,'.txt',sep=''))
sma_4c<-read.table(paste('../data/sacsma_',sma_site,'_4C.txt',sep=''))
hym_4c<-read.table(paste('../data/hymod_',hym_site,'_4C.txt',sep=''))

if(sma_site=='ORO'&hym_site=='ORO'){
  sma_hist<-read.table(paste('../data/sacsma_',sma_site,'.txt',sep=''))
  hym_hist<-read.table(paste('../data/hymod_',hym_site,'.txt',sep=''))
  sma_4c<-read.table(paste('../data/sacsma_',sma_site,'_4C.txt',sep=''))
  hym_4c<-read.table(paste('../data/hymod_',hym_site,'_4C.txt',sep=''))
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

saveRDS(hym_hist_vars,paste('../data/hym_hist_vars_',hym_site,'.rds',sep=''))
saveRDS(hym_4c_vars,paste('../data/hym_4c_vars_',hym_site,'.rds',sep=''))

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

saveRDS(sma_hist_vars,paste('../data/sma_hist_vars_',sma_site,'.rds',sep=''))
saveRDS(sma_4c_vars,paste('../data/sma_4c_vars_',sma_site,'.rds',sep=''))

#-----------------------------------------------------------------------
#2) Define arrays for predictions

#define predictive errors (truth - process)
hym_err_hist<-sma_hist_vars[,'sim']-hym_hist_vars[,'sim']
hym_err_4c<-sma_4c_vars[,'sim']-hym_4c_vars[,'sim']

sma_err_hist<-sma_hist_vars[,'obs']-sma_hist_vars[,'sim']
sma_err_4c<-sma_4c_vars[,'obs']-sma_4c_vars[,'sim']

#2a) HYMOD historical
#remove date columns
hym_predmat_hist<-hym_hist_vars[,4:dim(hym_hist_vars)[2]]
#remove obs (truth model) values (not used for prediction)
hym_predmat_hist<-hym_predmat_hist[-c(which(colnames(hym_predmat_hist)=='obs*'))]
#add predictive errors (predictand) and state variable predictors
hym_predmat_hist<-cbind(hym_err_hist,hym_predmat_hist)

#create labels for final array
labs_hym<-paste(rep(c('err','precip','tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),3),rep(0:(-3),each=dim(hym_predmat_hist)[2]))

#create super-array of state variables out to lag3 and including long-lag predictors
hym_predmat_hist<-cbind(hym_predmat_hist,
                      rbind(matrix(0,nrow=1,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-1),])),
                      rbind(matrix(0,nrow=2,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-2),])),
                      rbind(matrix(0,nrow=3,ncol=dim(hym_predmat_hist)[2]),as.matrix(hym_predmat_hist[1:(length(ix)-3),]))
                      ) 
#relabel and save
colnames(hym_predmat_hist)<-labs_hym
saveRDS(hym_predmat_hist,paste('../data/hym_predmat_hist_',hym_site,'.rds',sep=''))

#2b) HYMOD 4C
#remove date columns
hym_predmat_4c<-hym_4c_vars[,4:dim(hym_4c_vars)[2]]
#remove obs (truth model) values (not used for prediction)
hym_predmat_4c<-hym_predmat_4c[-c(which(colnames(hym_predmat_4c)=='obs*'))]
#add predictive errors (predictand) and state variable predictors
hym_predmat_4c<-cbind(hym_err_4c,hym_predmat_4c)

#create super-array for 4C predictors
hym_predmat_4c<-cbind(hym_predmat_4c,
                      rbind(matrix(0,nrow=1,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-1),])),
                      rbind(matrix(0,nrow=2,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-2),])),
                      rbind(matrix(0,nrow=3,ncol=dim(hym_predmat_4c)[2]),as.matrix(hym_predmat_4c[1:(length(ix)-3),]))
                      ) 
                      
#relabel and save
colnames(hym_predmat_4c)<-labs_hym
saveRDS(hym_predmat_4c,paste('../data/hym_predmat_4c_',hym_site,'.rds',sep=''))

#2c) SACSMA historical
#remove date columns
sma_predmat_hist<-sma_hist_vars[,4:dim(sma_hist_vars)[2]]
#remove obs (truth model) values (not used for prediction)
sma_predmat_hist<-sma_predmat_hist[-c(which(colnames(sma_predmat_hist)=='obs'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
sma_predmat_hist<-cbind(sma_err_hist,sma_predmat_hist)

#create labels for super-array of SAC-SMA predictors
labs_sma<-paste(rep(c('err','precip','tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm','store'),3),rep(0:(-3),each=dim(sma_predmat_hist)[2]))

#create 
sma_predmat_hist<-cbind(sma_predmat_hist,
                        rbind(matrix(0,nrow=1,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-1),])),
                        rbind(matrix(0,nrow=2,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-2),])),
                        rbind(matrix(0,nrow=3,ncol=dim(sma_predmat_hist)[2]),as.matrix(sma_predmat_hist[1:(length(ix)-3),]))
) 

#relabel and save
colnames(sma_predmat_hist)<-labs_sma
saveRDS(sma_predmat_hist,paste('../data/sma_predmat_hist_',sma_site,'.rds',sep=''))

#2b) SACSMA 4C
#remove date columns
sma_predmat_4c<-sma_4c_vars[,4:dim(sma_4c_vars)[2]]
#remove obs (truth model) values (not used for prediction)
sma_predmat_4c<-sma_predmat_4c[-c(which(colnames(sma_predmat_4c)=='obs'))]
#add predictive errors (predictand) and lag1:3 errors as predictors
sma_predmat_4c<-cbind(sma_err_4c,sma_predmat_4c)

#create super-array of predictors
sma_predmat_4c<-cbind(sma_predmat_4c,
                        rbind(matrix(0,nrow=1,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-1),])),
                        rbind(matrix(0,nrow=2,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-2),])),
                        rbind(matrix(0,nrow=3,ncol=dim(sma_predmat_4c)[2]),as.matrix(sma_predmat_4c[1:(length(ix)-3),]))
)

#relabel and save
colnames(sma_predmat_4c)<-labs_sma
saveRDS(sma_predmat_4c,paste('../data/sma_predmat_4c_',hym_site,'.rds',sep=''))

rm(list=ls());gc()

####################################################END#####################################################