setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')
library(scales)

hym_site<-'ORO'
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


#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

syn_tst_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_hist.rds')
syn_4c_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_4c.rds')

syn_tst_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period1,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
syn_4c_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period2,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

err_hist<-hym_predmat_hist[idx_comp,'err 0']*hym_kcfs_conv
err_4c<-hym_predmat_4c[idx_comp,'err 0']*hym_kcfs_conv

syn_tst_err_static<-syn_tst_err_static[idx_comp,]*hym_kcfs_conv
syn_4c_err_static<-syn_4c_err_static[idx_comp,]*hym_kcfs_conv

syn_tst_err_hybrid<-syn_tst_err_hyb[idx_comp,]*hym_kcfs_conv
syn_4c_err_hybrid<-syn_4c_err_hyb[idx_comp,]*hym_kcfs_conv
#--------------------------------------------------------------
#plot raw errors vs corrected errors
png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/fig10/fig10_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=640)
par(mfrow=c(2,2),oma=c(1,2,2,0),mar=c(0.5,1.5,0.5,0.5),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)
my.gray = alpha('gray50',alpha=0.5)

#a) hybrid test
err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_hist[seas]
  err_lst[[(i*2)]]<-syn_tst_err_hybrid[seas,]
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
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col=my.gray)
text(11,-4.5,'a)',cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)
mtext('Hybrid SWM',side=3, line=.5,cex=2,font=2)

#c) Test static
err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_hist[seas]
  err_lst[[(i*2)]]<-syn_tst_err_static[seas,]
}

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
axis(2,at=seq(-5,5,2.5),labels = F)#seq(-1,1,0.5) )
box(which='plot')
legend(6.5,5.5,c('empirical','simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2,bty='n')
abline(h=0,lty=2,col=my.gray)
text(11,-4.5,'c)',cex=3,font=2,adj=0)
mtext('Static SWM',side=3,line=.5,cex=2,font=2)

#b) 4c hybrid

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_4c[seas]
  err_lst[[(i*2)]]<-syn_4c_err_hybrid[seas,]
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
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('darkorange3','coral1'),cex=2)
abline(h=0,lty=2,col=my.gray)
text(11,-4.5,'b)',cex=3,font=2,adj=0)
mtext('Error (kcfs)',side=2,line=2,cex=1.75)

#d) 4c Static

err_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_4c[seas]
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
axis(2,at=seq(-5,5,2.5),labels = F)#seq(-1,1,0.5) )
box(which='plot')
legend(6.5,5.5,c('empirical','simulated'),lwd=c(6,6),col=c('darkorange3','coral1'),cex=2,bty='n')
abline(h=0,lty=2,col=my.gray)
text(11,-4.5,'d)',cex=3,font=2,adj=0)

dev.off()


###################################END#########################################