setwd('z:/oro_nonstat/')
library(lime)
library(Hmisc)

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
plot_mth<-3

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'

if(noise_reg==F){
  ll_vec<-readRDS(paste('fit_rev1/hymod_llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(ll_vec)
  print(ll_vec[seed])}
if(noise_reg==T){
  noise_ll_vec<-readRDS(paste('fit_rev1/hymod_noise-llvec_v-',vers,'.rds',sep=''))
  seed<-which.max(noise_ll_vec)
  print(noise_ll_vec[seed])}

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000
sma_kcfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

lime_feat_wt_tst<-readRDS(paste('out_rev1/lime_feat_wt_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
lime_feat_wt_4c<-readRDS(paste('out_rev1/lime_feat_wt_4c_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))

lab_ll_avg<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-avg')
lab_ll_trend<-paste(c('tavg','sim','runoff','baseflow','et','swe','upr_sm','lwr_sm'),'llag-trend')

if(vers=='all'){rf_idx<-which(colnames(hym_predmat_hist)!='err 0')}
if(vers=='err13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))}
if(vers=='sim13'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3')))}
if(vers=='err_sim13'){{rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','sim -1','sim -2','sim -3','err -1','err -2','err -3')))}}
if(vers=='err13_llag'){rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3',
                                                                          lab_ll_avg,lab_ll_trend)))}

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

ix2<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
idx_comp2<-c()
for(i in 1:length(wy_comp)){
  idx_comp2<-c(idx_comp2,which(ix2==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix2==paste(wy_comp[i],'-09-30',sep='')))
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

#calculate subsets based on positive or negative biased errors
err_lst<-c()
seas<-which(ixx_comp$mon==(plot_mth-1))
  
err_hst<-hym_predmat_hist[idx_comp,'err 0'][seas]
pos_idx_hst<-which(err_hst>=0)
neg_idx_hst<-which(err_hst<0)
err_hst_pos<-err_hst[pos_idx_hst]
err_hst_neg<-err_hst[neg_idx_hst]
  
err_wm<-hym_predmat_4c[idx_comp,'err 0'][seas]
pos_idx_4c<-which(err_wm>=0)
neg_idx_4c<-which(err_wm<0)
err_4c_pos<-err_wm[pos_idx_4c]
err_4c_neg<-err_wm[neg_idx_4c]
  
err_lst[1]<-median(err_hst_neg)*hym_kcfs_conv
err_lst[2]<-median(err_4c_pos)*hym_kcfs_conv
  
for(k in 2:length(rf_idx)){
  lime_comp<-lime_feat_wt_tst[idx_comp2,]
  lime_hst<-lime_comp[seas,(k-1)]
  lime_comp_4c<-lime_feat_wt_4c[idx_comp2,]
  lime_4c<-lime_comp_4c[seas,(k-1)]
  err_lst[(k*2-1)]<-median(lime_hst[neg_idx_hst])*4.5
  err_lst[(k*2)]<-median(lime_4c[pos_idx_4c])*4.5
}


png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig7/fig7_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),height=512,width=768)
my.col<-alpha('peachpuff',alpha = 0.3)

cnames<-str_remove(colnames(lime_feat_wt_tst)[1:9],' 0')

par(cex.lab=2,cex.axis=1.5,mar=c(6,5,1.5,4.5))
barplot(err_lst[1:20],space=c(0,0.1,1.5,rep(c(0.1,0.75),8),0.1),col=rep(c('skyblue4','darkorange3'),12),
        ylab='Median error (kcfs)',ylim=c(-1.5,1.5),border=c(rep('black',2),rep('gray50',22)))
box(which='plot')
axis(1,at=c(1,seq(4.5,30,2.85)),las=2,labels = c('error',cnames))
axis(4,at=seq(-1,1,0.25)*4.5,labels=seq(-1,1,0.25),col='gray30',col.axis='gray30')
mtext('Median Feature Weight',side=4,cex=2,line=3,col='gray30')
polygon(c(-1.3,-1.3,3,3),c(-1.5,1.5,1.5,-1.5),col=my.col)
text(25,1.25,'Mar',font=4,cex=3)
abline(h=0,lty=2)
legend('bottomright',c('Test','Test+4C'),lwd=c(6,6),col=c('skyblue4','darkorange3'),cex=2,bty='n')

dev.off()

#------------------------------------------------------------------------------------------------------
#SI figure with AR terms

for(k in 2:length(rf_idx)){
  lime_hst<-lime_feat_wt_tst[seas,(k-1)]
  lime_4c<-lime_feat_wt_4c[seas,(k-1)]
  err_lst[(k*2-1)]<-median(lime_hst[neg_idx_hst])*1.75
  err_lst[(k*2)]<-median(lime_4c[pos_idx_4c])*1.75
}

png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig7/fig7_ar-terms_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),height=512,width=768)
my.col<-alpha('peachpuff',alpha = 0.3)

cnames<-str_remove(colnames(lime_feat_wt_tst)[1:12],' 0')

par(cex.lab=2,cex.axis=1.5,mar=c(6,5,1,4.5))
barplot(err_lst,space=c(0,0.1,1.5,rep(c(0.1,0.75),11),0.1),col=rep(c('skyblue4','darkorange3'),15),
        ylab='Median error (kcfs)',ylim=c(-2.75,2.75),border=c(rep('black',2),rep('gray50',28)))
box(which='plot')
axis(1,at=c(1,seq(4.5,36,2.85)),las=2,labels = c('error',cnames))
axis(4,at=seq(-2,2,0.5)*1.75,labels=seq(-2,2,0.5),col='gray30',col.axis='gray30')
mtext('Median Feature Weight',side=4,cex=2,line=3,col='gray30')
polygon(c(-1.6,-1.6,3,3),c(-2.75,2.75,2.75,-2.75),col=my.col)
text(35,2.25,'Mar',font=4,cex=3)
abline(h=0,lty=2)
legend(5,-1.5,c('Test','Test+4C'),lwd=c(6,6),col=c('skyblue4','darkorange3'),cex=2,bty='n')

dev.off()

#############################END##############################