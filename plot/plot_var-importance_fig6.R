#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
library(ranger)

ix<-seq(as.Date('1987-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

ix_cal<-seq(as.Date('1987-10-01'),as.Date('1997-09-30'),'day')
ixx_cal<-as.POSIXlt(ix_cal)
idx_cal<-which(ix=='1987-10-01'):which(ix=='1997-09-30')

#fit RF model for error correction
predictor_mat_hist<-readRDS('fit/predictor_mat_hist.rds')
err_hist<-readRDS('fit/err_hist.rds')

sm_err_vec<-ksmooth(1:length(idx_cal),err_hist[idx_cal],bandwidth = 1)$y

hyp_param_err<-readRDS('fit/rf_hyp-param-err.rds')

srt_idx_err<-sort(hyp_param_err[3,],index.return=T)

ntree<-hyp_param_err[1,srt_idx_err$ix[1]]
mtry<-hyp_param_err[2,srt_idx_err$ix[1]]
seed<-hyp_param_err[4,srt_idx_err$ix[1]]

set.seed(seed)
rg_db<-ranger(x=predictor_mat_hist[idx_cal,],y=sm_err_vec,classification = F,
                  num.trees = ntree,mtry = mtry,
                  replace=T,importance = 'permutation',splitrule = 'variance',quantreg = F)

rg_db<-readRDS('fit/rf_corr_hymod_cal_87_97_default.rds')
rg_db_tune<-readRDS('fit/rf_corr_hymod_cal_87_97_hyp-tune.rds')


#--------------------------------------------------------------
#variable importance
var_vec<-c('sim','precip','tavg','et','sm','swe','runoff','baseflow','upr_sm','lwr_sm','lag1','lag2','lag3')
var_imp<-rg_db$variable.importance / sum(rg_db$variable.importance)
srt_var_imp<-sort(var_imp,index.return=T,decreasing=T)

png('h:/oroville_non-stationary/paper/figs/fig6_def.png',width=768,height=640)
par(las=2,mar=c(6.5,5.5,1,0),mgp=c(4,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)
barplot(srt_var_imp$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp$x[1],digits=2),adj=0.5,cex=1.5)

dev.off()

#plot SI figure S2
var_imp_tune<-rg_db_tune$variable.importance / sum(rg_db_tune$variable.importance)
srt_var_imp_tune<-sort(var_imp_tune,index.return=T,decreasing=T)

png('h:/oroville_non-stationary/paper/figs/default-tune-compare_figS2.png',width=768,height=448)
par(mfrow=c(1,2),las=2,mar=c(6.5,5.5,1,0),mgp=c(4,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)
barplot(srt_var_imp_tune$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp_tune$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp_tune$x[1],digits=2),adj=0.5,cex=1.5)

barplot(srt_var_imp$x,ylim=c(0,0.2),names.arg = var_vec[srt_var_imp$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')
text(0.75,0.18,round(srt_var_imp$x[1],digits=2),adj=0.5,cex=1.5)

dev.off()

##########################################END######################################
