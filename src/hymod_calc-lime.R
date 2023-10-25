#setwd('z:/oro_nonstat/')
library(ranger)
library(lime)

print(paste('start',Sys.time()))

#specifications
hym_site<-'ORO'
vers<-'err13'
noise_reg<-F  #use noise regularized coefficients?

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)}

#load predictor arrays
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))
rf_err_corr<-readRDS(paste('fit_rev1/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

idx_cal<-which(ix=='1988-10-01'):which(ix=='1989-09-30')

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}

pred_df<-as.data.frame(hym_predmat_hist[idx_tst,rf_idx])
pred_df_4c<-as.data.frame(hym_predmat_4c[idx_tst,rf_idx])

#lime setup
explainer <- lime(pred_df,rf_err_corr)
explainer_4c <- lime(pred_df_4c,rf_err_corr)

explanation <- explain(pred_df, explainer,n_labels=1,n_features = length(rf_idx))
explanation_4c <- explain(pred_df_4c, explainer_4c,n_labels=1,n_features = length(rf_idx))

saveRDS(explanation,paste('out_rev1/explanation_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(explanation_4c,paste('out_rev1/explanation_4c_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))

feat_wt<-matrix(explanation$feature_weight,ncol=length(rf_idx),byrow=T)
feat_wt_4c<-matrix(explanation_4c$feature_weight,ncol=length(rf_idx),byrow=T)

lime_feat_wt<-feat_wt
colnames(lime_feat_wt)<-colnames(hym_predmat_hist[,rf_idx])
rownames(lime_feat_wt)<-rownames(hym_predmat_hist[idx_tst,])

lime_feat_wt_4c<-feat_wt_4c
colnames(lime_feat_wt_4c)<-colnames(hym_predmat_4c[,rf_idx])
rownames(lime_feat_wt_4c)<-rownames(hym_predmat_4c[idx_tst,])

saveRDS(lime_feat_wt,paste('out_rev1/lime_feat_wt_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(lime_feat_wt_4c,paste('out_rev1/lime_feat_wt_4c_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))

print(paste('end',Sys.time()))

rm(list=ls());gc()


##########################################END################################