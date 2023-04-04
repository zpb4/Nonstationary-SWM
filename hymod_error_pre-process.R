#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')

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

rm(list=ls());gc()

################################END###############################