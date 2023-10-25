#fit model between hymod and sacsma
setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')
library(ranger)

hym_site<-'ORO'
#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

#load predictor arrays
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

#---------------------------------------------------------------------------------
#S7

syn_tst_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_hist.rds')
syn_4c_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_4c.rds')

syn_tst_err_static<-syn_tst_err_static[idx_comp,]*hym_kcfs_conv
syn_4c_err_static<-syn_4c_err_static[idx_comp,]*hym_kcfs_conv

#plot raw errors vs corrected errors
png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/static-swm_',wy_tag,'_figS7.png',sep=''),width=768,height=1024)
par(mfrow=c(2,1),oma=c(1,2,2,0),mar=c(0.5,1.5,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

#a) hybrid test
err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-hym_predmat_hist[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-syn_tst_err_static[seas,]
}

mths<-c('J','F','M','A','M','J','J','A','S','O','N','D')
#mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(-0.05,.305),12),
        names = 1:24,
        las = 2,
        boxwex = 0.25,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red')
text(0.75,4.75,'a)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)

#b) 4c hybrid

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-hym_predmat_4c[idx_comp,'err 0'][seas]*hym_kcfs_conv
  err_lst[[(i*2)]]<-syn_4c_err_static[seas,]
}

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(-0.05,.305),12),
        names = 1:24,
        las = 2,
        boxwex = 0.25,
        notch = F,
        col = rep(c('darkorange3','coral1'),12),
        ylim = c(-5,5),
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(-5,5,2.5),labels = seq(-5,5,2.5) )
box(which='plot')
legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('darkorange3','coral1'),cex=2)
abline(h=0,lty=2,col='red')
text(0.75,4.75,'b)',cex=3,font=2,adj=0)
text(11,-4.75,wy_tag,cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)

dev.off()



##########################################END####################################################