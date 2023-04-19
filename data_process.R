#Process and arrange data into labeled matrices

setwd('d:/Nonstationary-SWM/')

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1987-10-01'):which(ix=='1997-09-30')
idx_val<-which(ix=='1997-10-01'):which(ix=='2003-09-30')

#simulation data
sma_hist<-read.table('data/simflow_sacsma_ORO.txt')
hym_hist<-read.table('data/simflow_hymod_ORO.txt')
sma_4c<-read.table('data/simflow_sacsma_ORO_4C.txt')
hym_4c<-read.table('data/simflow_hymod_ORO_4C.txt')

#define errors
err_hist<-sma_hist[,4]-hym_hist[,4]
err_4c<-sma_4c[,4]-hym_4c[,4]

lr_hist<-log(hym_hist[,4]/sma_hist[,4])
lr_4c<-log(hym_4c[,4]/sma_4c[,4])

#historical state variables
meteo_hist<-read.table('data/meteo_livneh_CDEC_watershed_ORO_fix.txt')
ix3<-seq(as.Date('1950-01-01'),as.Date('2013-12-31'),'day')
ix3_sub<-which(ix3=='1987-10-01'):which(ix3=='2013-09-30')
precip_hist<-meteo_hist[ix3_sub,4]
tavg_hist<- (meteo_hist[ix3_sub,5] + meteo_hist[ix3_sub,6]) / 2


et_hist<-read.table('data/simaet_hymod_ORO.txt')[,4]
#swe_hist<-read.table('data/simswe_hymod_ORO.txt')[,4]
sm_hist<-read.table('data/simsm_hymod_ORO.txt')[,4]

#hist output variables
addvars_hist<-read.table('data/simvars_hymod_ORO.txt')
roff_hist<-addvars_hist[,5]
bflow_hist<-addvars_hist[,6]
swe_hist<-addvars_hist[,7]
usoil_hist<-addvars_hist[,8]
lsoil_hist<-addvars_hist[,9]


#4c state variables
idx4<-which(ix=='2003-10-01'):which(ix=='2013-09-30')
meteo_4c<-read.table('data/hymod_pr_tas_sm_et_swe_ORO.txt')

precip_4c<-meteo_4c[idx4,4]
tavg_4c<- meteo_4c[idx4,5]
et_4c<-meteo_4c[idx4,7]
#swe_4c<-meteo_4c[idx4,8]
sm_4c<-meteo_4c[idx4,6]

#4c output variables
addvars_4c<-read.table('data/simvars_hymod_ORO_4C.txt')
roff_4c<-addvars_4c[idx4,5]
bflow_4c<-addvars_4c[idx4,6]
swe_4c<-addvars_4c[idx4,7]
usoil_4c<-addvars_4c[idx4,8]
lsoil_4c<-addvars_4c[idx4,9]

#simulations
sim_hist<-hym_hist[,4]
sim_4c<-hym_4c[idx4,4]

labs<-c('err','lr_err','sim','precip','tavg','et','sm','swe','runoff','baseflow','upr_sm','lwr_sm')

var_mat_hist<-cbind(err_hist,lr_hist,sim_hist,precip_hist,tavg_hist,et_hist,sm_hist,swe_hist,roff_hist,bflow_hist,usoil_hist,lsoil_hist)
var_mat_4c<-cbind(err_4c[idx4],lr_4c[idx4],sim_4c,precip_4c,tavg_4c,et_4c,sm_4c,swe_4c,roff_4c,bflow_4c,usoil_4c,lsoil_4c)

colnames(var_mat_hist)<-labs
colnames(var_mat_4c)<-labs

rownames(var_mat_hist)<-as.character(ix)
rownames(var_mat_4c)<-as.character(ix[idx4])

saveRDS(var_mat_hist,'data/var_mat_hist.rds')
saveRDS(var_mat_4c,'data/var_mat_4c.rds')

rm(list=ls());gc()

####################################################END#####################################################