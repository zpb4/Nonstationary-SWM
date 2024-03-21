setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')

hym_site<-'ORO'
sma_site<-'ORO'
gen_period1<-'hist-all'
gen_period2<-'4c-all'
samp_type<-'split'
n<-1000
#specifications
vers<-'err13' 
# 'err13' - SV + lag1:3 errors
# 'sim13' - SV + lag1:3 sim
# 'err_sim13' - SV + lag1:3 errors and sim
# 'err13_llag' - SV + lag1:3 errors and long lagged variables
# 'all' - all SV (including lag 1:3 terms) and long lagged variables
noise_reg=T
pcnt<-99

kcfs_to_mafy = 2.29568411*10**-5 * 86400 / 1000


#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))
sma_4c_vars<-readRDS(paste('data_rev1/sma_4c_vars_',sma_site,'.rds',sep=''))

syn_tst_flow_static<-readRDS('out_rev1/hymod_benchmark_syn-flow_hist.rds')
syn_4c_flow_static<-readRDS('out_rev1/hymod_benchmark_syn-flow_4c.rds')

syn_tst_flow_hyb<-readRDS(paste('out_rev1/hymod_syn-flow_',hym_site,'_',gen_period1,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
syn_4c_flow_hyb<-readRDS(paste('out_rev1/hymod_syn-flow_',hym_site,'_',gen_period2,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

hym_sim_hist<-hym_predmat_hist[idx_comp,'sim 0']*hym_kcfs_conv
hym_sim_4c<-hym_predmat_4c[idx_comp,'sim 0']*hym_kcfs_conv

sma_sim_hist<-sma_hist_vars[idx_comp,'sim']*hym_kcfs_conv
sma_sim_4c<-sma_4c_vars[idx_comp,'sim']*hym_kcfs_conv

syn_tst_flow_static<-syn_tst_flow_static[idx_comp,]*hym_kcfs_conv
syn_4c_flow_static<-syn_4c_flow_static[idx_comp,]*hym_kcfs_conv

syn_tst_flow_hybrid<-syn_tst_flow_hyb[idx_comp,]*hym_kcfs_conv
syn_4c_flow_hybrid<-syn_4c_flow_hyb[idx_comp,]*hym_kcfs_conv

pct<-pcnt/100
idx<-round(pct*length(idx_comp))
  
hyb_tst_dist<-c()
hyb_4c_dist<-c()
stat_tst_dist<-c()
stat_4c_dist<-c()
for(i in 1:n){
  hyb_tst_dist[i]<-sort(syn_tst_flow_hybrid[,i])[idx]
  hyb_4c_dist[i]<-sort(syn_4c_flow_hybrid[,i])[idx]
  stat_tst_dist[i]<-sort(syn_tst_flow_static[,i])[idx]
  stat_4c_dist[i]<-sort(syn_4c_flow_static[,i])[idx]
}

yrs<-2005:2017

sma_tst_annloflow<-c()
sma_4c_annloflow<-c()
hym_tst_annloflow<-c()
hym_4c_annloflow<-c()

for(i in 1:length(yrs)){
  sset = which(ixx_comp$year==(yrs[i]-1900))
  sma_tst_annloflow[i]<-sum(sma_sim_hist[sset])
  sma_4c_annloflow[i]<-sum(sma_sim_4c[sset])
  hym_tst_annloflow[i]<-sum(hym_sim_hist[sset])
  hym_4c_annloflow[i]<-sum(hym_sim_4c[sset])
}
  
sma_tst_annlow<-min(sma_tst_annloflow)
sma_4c_annlow<-min(sma_tst_annloflow)
hym_tst_annlow<-min(hym_tst_annloflow)
hym_4c_annlow<-min(hym_4c_annloflow)



hyb_tst_annloflow<-c()
hyb_4c_annloflow<-c()
stat_tst_annloflow<-c()
stat_4c_annloflow<-c()

for(s in 1:n){
  hyb_tst<-c()
  hyb_4c<-c()
  stat_tst<-c()
  stat_4c<-c()
  for(i in 1:length(yrs)){
    sset = which(ixx_comp$year==(yrs[i]-1900))
    hyb_tst[i]<-sum(syn_tst_flow_hybrid[sset,s])
    hyb_4c[i]<-sum(syn_4c_flow_hybrid[sset,s])
    stat_tst[i]<-sum(syn_tst_flow_static[sset,s])
    stat_4c[i]<-sum(syn_4c_flow_static[sset,s])
  }
  hyb_tst_annloflow[s]<-min(hyb_tst)
  hyb_4c_annloflow[s]<-min(hyb_4c)
  stat_tst_annloflow[s]<-min(stat_tst)
  stat_4c_annloflow[s]<-min(stat_4c)
}

#--------------------------------------------------------------
#plot raw errors vs corrected errors
#png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/fig10_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=640)
par(mfrow=c(2,2),oma=c(1,2,2,0),mar=c(0.5,1.5,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)
breaks<-seq(0,200,1)

par(mfrow=c(2,2))
hist(hyb_tst_dist,breaks=breaks,xlim=c(30,60))
abline(v=sort(sma_sim_hist)[idx])
abline(v=sort(hym_sim_hist)[idx],col='green')

hist(stat_tst_dist,breaks=breaks,xlim=c(30,60))
abline(v=sort(sma_sim_hist)[idx])
abline(v=sort(hym_sim_hist)[idx],col='green')

hist(hyb_4c_dist,breaks=breaks,xlim=c(30,60))
abline(v=sort(sma_sim_4c)[idx])
abline(v=sort(hym_sim_4c)[idx],col='green')

hist(stat_4c_dist,breaks=breaks,xlim=c(30,60))
abline(v=sort(sma_sim_4c)[idx])
abline(v=sort(hym_sim_4c)[idx],col='green')


#dev.off()

png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/si/si_extremes-comp_',pcnt,'pct_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=384)
par(mfrow=c(1,2),mar=c(3,4,3,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75,cex.main=2)
breaks<-seq(0,200,1)

hist(hyb_tst_dist,breaks=breaks,xlim=c(30,45),main='Hybrid SWM - Test',xlab='Flow (kcfs)')
abline(v=sort(sma_sim_hist)[idx],col='black',lwd=2)
abline(v=sort(hym_sim_hist)[idx],col='skyblue4',lwd=2)
legend('topleft',c('truth','process','hybrid SWM'),col=c('black','skyblue4','gray'),lwd=c(2,2,6),bty='n')

hist(hyb_4c_dist,breaks=breaks,xlim=c(30,45),main='Hybrid SWM - Test+4C',xlab='Flow (kcfs)')
abline(v=sort(sma_sim_4c)[idx],col='black',lwd=2)
abline(v=sort(hym_sim_4c)[idx],col='darkorange3',lwd=2)
legend('topleft',c('truth','process','hybrid SWM'),col=c('black','darkorange3','gray'),lwd=c(2,2,6),bty='n')

dev.off()

#-----------------------------------------------------------------------------------

#kde plots 2 x 2 with cumulative annual flow stat

png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/fig12/fig12_',pcnt,'pct_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=768)
par(mfrow=c(2,2),mar=c(3,4,3,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75,cex.main=2)

kd_hyb_tst<-density(hyb_tst_dist,n=512)
hyb_y<-kd_hyb_tst$y
hyb_x<-kd_hyb_tst$x
kd_stat_tst<-density(stat_tst_dist,n=512)
stat_y<-kd_stat_tst$y
stat_x<-kd_stat_tst$x
plot(hyb_x,hyb_y,col='palevioletred3',lwd=3,xlim=c(32,52),type='l',xlab='Flow (kcfs)', ylab='Density',ylim=c(0,0.4))
lines(stat_x,stat_y,col='gray70',lwd=3)
abline(v=sort(sma_sim_hist)[idx],col='black',lwd=1,lty=2)
abline(v=sort(hym_sim_hist)[idx],col='dodgerblue',lwd=1,lty=2)
legend('topright',c('truth','process','hyb-SWM','stat-SWM'),col=c('black','dodgerblue','palevioletred3','gray70'),lwd=c(1,1,3,3),lty=c(2,2,1,1),bty='n',cex=1.75)
text(33,0.38,'a)',font=2,cex=3)

kd_hyb_4c<-density(hyb_4c_dist,n=512)
hyb_y<-kd_hyb_4c$y
hyb_x<-kd_hyb_4c$x
kd_stat_4c<-density(stat_4c_dist,n=512)
stat_y<-kd_stat_4c$y
stat_x<-kd_stat_4c$x
plot(hyb_x,hyb_y,col='darkorange3',lwd=3,xlim=c(32,52),type='l',xlab='Flow (kcfs)', ylab='',ylim=c(0,0.4))
lines(stat_x,stat_y,col='gray70',lwd=3)
abline(v=sort(sma_sim_4c)[idx],col='black',lwd=1,lty=2)
abline(v=sort(hym_sim_4c)[idx],col='seagreen3',lwd=1,lty=2)
legend('topright',c('truth','process','hyb-SWM','stat-SWM'),col=c('black','seagreen3','darkorange3','gray70'),lwd=c(1,1,3,3),lty=c(2,2,1,1),bty='n',cex=1.75)
text(33,0.38,'b)',font=2,cex=3)

kd_hyb_tst<-density(hyb_tst_annloflow * kcfs_to_mafy,n=512)
hyb_y<-kd_hyb_tst$y
hyb_x<-kd_hyb_tst$x
kd_stat_tst<-density(stat_tst_annloflow * kcfs_to_mafy,n=512)
stat_y<-kd_stat_tst$y
stat_x<-kd_stat_tst$x
plot(hyb_x,hyb_y,col='palevioletred3',lwd=3,xlim=c(0.75,1.75),type='l',xlab='Cumulative flow (MAF)', ylab='Density',ylim=c(0,8))
lines(stat_x,stat_y,col='gray70',lwd=3)
abline(v=sma_tst_annlow * kcfs_to_mafy,col='black',lwd=1,lty=2)
abline(v=hym_tst_annlow * kcfs_to_mafy,col='dodgerblue',lwd=1,lty=2)
#legend('topleft',c('truth','process','hybrid SWM','static SWM'),col=c('black','skyblue4','palevioletred3','gray70'),lwd=c(1,1,3,3),bty='n')
text(0.8,7.5,'c)',font=2,cex=3)

kd_hyb_4c<-density(hyb_4c_annloflow * kcfs_to_mafy,n=512)
hyb_y<-kd_hyb_4c$y
hyb_x<-kd_hyb_4c$x
kd_stat_4c<-density(stat_4c_annloflow * kcfs_to_mafy,n=512)
stat_y<-kd_stat_4c$y
stat_x<-kd_stat_4c$x
plot(hyb_x,hyb_y,col='darkorange3',lwd=3,xlim=c(0.75,1.75),type='l',xlab='Cumulative flow (MAF)', ylab='',ylim=c(0,8))
lines(stat_x,stat_y,col='gray70',lwd=3)
abline(v=sma_4c_annlow * kcfs_to_mafy,col='black',lwd=1,lty=2)
abline(v=hym_4c_annlow * kcfs_to_mafy,col='seagreen3',lwd=1,lty=2)
#legend('topleft',c('truth','process','hybrid SWM','static SWM'),col=c('black','darkorange3','seagreen4','gray70'),lwd=c(1,1,3,3),bty='n')
text(0.8,7.5,'d)',font=2,cex=3)

dev.off()
###################################END#########################################

