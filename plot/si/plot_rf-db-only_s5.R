setwd('z:/oro_nonstat/')
library(ranger)
source('mm-cfs_conversion.R')

hym_site<-'ORO'
sma_site<-'ORO'
#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
noise_reg=T

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000
sma_kcfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs/1000

sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0',
                                                                          lab_ll_avg,lab_ll_trend)))}

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

obs_tst<-sma_hist_vars[idx_tst,'obs']
tst_wy<-2005:2018
obs_sum<-c()

for(i in 1:length(tst_wy)){
  obs_sum[i]<-sum(obs_tst[which(ix_tst==as.Date(paste(tst_wy[i]-1,'-10-01',sep=''))):which(ix_tst==as.Date(paste(tst_wy[i],'-09-30',sep='')))])
}

srt_obs<-sort(obs_sum,index.return=T)
tst_wy[tail(srt_obs$ix,5)] #5 wettest
#[1] 2016 2005 2011 2006 2017
tst_wy[head(srt_obs$ix,5)] #5 driest
#[1] 2014 2015 2008 2007 2012

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

rf_err<-hym_predmat_hist[,'err 0']
rf_pred<-hym_predmat_hist[,rf_idx]

rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],importance = 'impurity')

#--------------------------------------------------------------
#plot raw errors vs corrected errors
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/rf-db-only_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'_figS5.png',sep=''),width=768,height=900)
par(mfrow=c(2,1),mar=c(1,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

err_hist<-hym_predmat_hist[idx_comp,'err 0']

db_err<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-err_hist-db_err

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_hist[seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-err_db[seas]*hym_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.3),12),
        names = 1:24,
        las = 2,
        boxwex = rep(c(.3,.2),12),
        notch = F,
        col = rep(c('skyblue4','skyblue'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Raw','RF-corrected'),lwd=c(6,6),col=c('skyblue4','skyblue'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(0.75,4.75,'a)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

par(mar=c(2,4,0,1))

#4c
err_4c<-hym_predmat_4c[idx_comp,'err 0']

db_err_4c<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,rf_idx])$predictions
err_db_4c<-err_4c-db_err_4c
#db_err_4c<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,2:10])$predictions
#err_db_4c<-err_4c-db_err_4c

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_4c[seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-err_db_4c[seas]*hym_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.3),12),
        names = 1:24,
        las = 2,
        boxwex = rep(c(.3,.2),12),
        notch = F,
        col = rep(c('darkorange3','goldenrod1'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Raw','RF-corrected'),lwd=c(6,6),col=c('darkorange3','goldenrod1'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(.75,4.75,'b)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)

dev.off()


###################################END#########################################