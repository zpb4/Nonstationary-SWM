setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')

hym_site<-'ORO'
sma_site<-'ORO'
#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

wy_comp1<-c(2007,2008,2012,2014,2015)
wy_tag1<-'5dry'

hym_cfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs
sma_cfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs

sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))
hym_hist_vars<-readRDS(paste('data_rev1/hym_hist_vars_',sma_site,'.rds',sep=''))
hym_4c_vars<-readRDS(paste('data_rev1/hym_4c_vars_',sma_site,'.rds',sep=''))
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))
sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))
sma_predmat_4c<-readRDS(paste('data_rev1/sma_predmat_4c_',sma_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

#comparison indices
#tag0
ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

#tag1
ix_comp1<-c()
idx_comp1<-c()
for(i in 1:length(wy_comp)){
  ix_comp1<-c(ix_comp1,seq(as.Date(paste(wy_comp1[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp1[i],'-09-30',sep='')),'day'))
  idx_comp1<-c(idx_comp1,which(ix==paste(wy_comp1[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp1[i],'-09-30',sep='')))
}
ix_comp1<-as.Date(ix_comp1,origin = '1970-01-01')
ixx_comp1<-as.POSIXlt(ix_comp1)

#data process for plot panel a
err_tst<-hym_predmat_hist[idx_comp,'err']*hym_cfs_conv/1000 #convert to cfs
err_4c<-hym_predmat_4c[idx_comp,'err']*hym_cfs_conv/1000    #convert to kcfs
err_tst1<-hym_predmat_hist[idx_comp1,'err']*hym_cfs_conv/1000 #convert to cfs
err_4c1<-hym_predmat_4c[idx_comp1,'err']*hym_cfs_conv/1000    #convert to kcfs

#data process for plot panel b
sim_tst<-sma_predmat_hist[idx_comp,'sim']
sim_4c<-sma_predmat_4c[idx_comp,'sim']
obs_tst<-sma_hist_vars[idx_comp,'obs']
sim_tst1<-sma_predmat_hist[idx_comp1,'sim']
sim_4c1<-sma_predmat_4c[idx_comp1,'sim']
obs_tst1<-sma_hist_vars[idx_comp1,'obs']

ix_leap<-seq(as.Date('2008-01-01'),as.Date('2008-12-31'),'day')
ixx_leap<-as.POSIXlt(ix_leap)
dy_idx<-ixx_leap$mday
mth_idx<-ixx_leap$mo
mov_avg_window<-15 #size of smoothing window for climatology

#define moving average function
mavg<-function(x,wdow){
  out<-c()
  for(i in 1:length(x)){
    samp<-max(1,(i-wdow)):min(length(x),(i+wdow))
    out[i]<-mean(x[samp])}
  return(out)}

#calculate daily means across 1989:2019 for each day w/Feb 29 included
dly_vals_sma_tst<-array(NA,c(3,366))
dly_vals_sma_4c<-array(NA,c(3,366))
dly_vals_obs_tst<-array(NA,c(3,366))

for(i in 1:length(dy_idx)){
  dly_vals_obs_tst[2,i]<-mean(obs_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_obs_tst[1,i]<-min(obs_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_obs_tst[3,i]<-max(obs_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_tst[2,i]<-mean(sim_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_tst[1,i]<-min(sim_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_tst[3,i]<-max(sim_tst[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_4c[2,i]<-mean(sim_4c[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_4c[1,i]<-min(sim_4c[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
  dly_vals_sma_4c[3,i]<-max(sim_4c[which(ixx_comp$mo==mth_idx[i]&ixx_comp$mday==dy_idx[i])])
}

#define moving average climatology for leap and non-leap years
mavg_leap_sma_tst<-mavg(dly_vals_sma_tst[2,],mov_avg_window)*hym_cfs_conv/1000 #convert to kcfs
mavg_leap_sma_4c<-mavg(dly_vals_sma_4c[2,],mov_avg_window)*hym_cfs_conv/1000   #convert to kcfs
mavg_leap_obs_tst<-mavg(dly_vals_obs_tst[2,],mov_avg_window)*hym_cfs_conv/1000 #convert to kcfs

#tag1
dly_vals_sma_tst1<-array(NA,c(3,366))
dly_vals_sma_4c1<-array(NA,c(3,366))
dly_vals_obs_tst1<-array(NA,c(3,366))

for(i in 1:length(dy_idx)){
  dly_vals_obs_tst1[2,i]<-mean(obs_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_obs_tst1[1,i]<-min(obs_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_obs_tst1[3,i]<-max(obs_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_tst1[2,i]<-mean(sim_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_tst1[1,i]<-min(sim_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_tst1[3,i]<-max(sim_tst1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_4c1[2,i]<-mean(sim_4c1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_4c1[1,i]<-min(sim_4c1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
  dly_vals_sma_4c1[3,i]<-max(sim_4c1[which(ixx_comp1$mo==mth_idx[i]&ixx_comp1$mday==dy_idx[i])])
}

#define moving average climatology for leap and non-leap years
mavg_leap_sma_tst1<-mavg(dly_vals_sma_tst1[2,],mov_avg_window)*hym_cfs_conv/1000 #convert to kcfs
mavg_leap_sma_4c1<-mavg(dly_vals_sma_4c1[2,],mov_avg_window)*hym_cfs_conv/1000   #convert to kcfs
mavg_leap_obs_tst1<-mavg(dly_vals_obs_tst1[2,],mov_avg_window)*hym_cfs_conv/1000 #convert to kcfs


#--------------------------------------------------------------
png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig4_',wy_tag,'_',wy_tag1,'.png',sep=''),width=768,height=1024)

pos_mat<-matrix(1:3,ncol=1,byrow=T)
layout(pos_mat,heights = c(1.75,1.75,1))
par(mar=c(3,6,1,1),mgp=c(3.5,1,0),tcl=-0.3,cex.lab=3,cex.axis=2.75)

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_tst[seas]
  err_lst[[(i*2)]]<-err_4c[seas]
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','darkorange3'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = NA)#mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Test','Test+4C'),lwd=c(6,6),col=c('skyblue4','darkorange3'),cex=3,bty='n')
abline(h=0,lty=2,col='red')
text(1,4.75,'a)',cex=4,font=2)
text(11.5,-4.75,wy_tag,font=2,cex=4)

#wy_tag2
err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp1$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_tst1[seas]
  err_lst[[(i*2)]]<-err_4c1[seas]
}

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','darkorange3'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = NA)#mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
#legend('topright',c('Test','Test+4C'),lwd=c(6,6),col=c('skyblue4','darkorange3'),cex=3,bty='n')
abline(h=0,lty=2,col='red')
text(1,4.75,'b)',cex=4,font=2)
text(11.5,-4.75,wy_tag1,font=2,cex=4)

#smoothed SMA means plot
plot(1:366,mavg_leap_obs_tst,col='black',type='l',ylim=c(0,25),xlim=c(12,354),lwd=2,ylab='Flow (kcfs)',xlab='',
     axes=F)
lines(1:366,mavg_leap_sma_tst,col='skyblue4',lwd=4)
lines(1:366,mavg_leap_sma_4c,col='darkorange3',lwd=4)

lines(1:366,mavg_leap_obs_tst1,col='black',lwd=2,lty=3)
lines(1:366,mavg_leap_sma_tst1,col='skyblue4',lwd=4,lty=3)
lines(1:366,mavg_leap_sma_4c1,col='darkorange3',lwd=4,lty=3)
axis(1,at=seq(15,366,30.5),labels=mths)
axis(2,at=seq(0,25,10),labels=seq(0,25,10))
text(10,23,'c)',cex=4,font=2)
legend(175,25,c('obs','SAC-SMA, Test','SAC-SMA, Test+4C'),col=c('black','skyblue4','darkorange3'),lwd=c(2,4,4),cex=3,bty='n')
box(which='plot')
text(90,23,wy_tag,cex=3,font=2)
text(90,9,wy_tag1,cex=3,font=2)

dev.off()

################################END###########################################
