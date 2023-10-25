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
noise_reg=F


#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-2005:2018
wy_tag<-'tst-all'
hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000

hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))
sma_predmat_4c<-readRDS(paste('data_rev1/sma_predmat_4c_',sma_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

cp_tst<-readRDS(paste('out_rev1/cp_tst_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
cp_tst_bench<-readRDS('out_rev1/cp_tst_bench.rds')
cp_4c<-readRDS(paste('out_rev1/cp_4c_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
cp_4c_bench<-readRDS('out_rev1/cp_4c_bench.rds')

syn_tst_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_hist.rds')
syn_4c_err_static<-readRDS('out_rev1/hymod_benchmark_syn-err_4c.rds')

syn_tst_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period1,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
syn_4c_err_hyb<-readRDS(paste('out_rev1/hymod_syn-err_',hym_site,'_',gen_period2,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))

err_hist<-hym_predmat_hist[idx_comp,'err 0']*hym_kcfs_conv
err_4c<-hym_predmat_4c[idx_comp,'err 0']*hym_kcfs_conv

syn_tst_err_static<-syn_tst_err_static[idx_comp,]*hym_kcfs_conv
syn_4c_err_static<-syn_4c_err_static[idx_comp,]*hym_kcfs_conv

syn_tst_err_hybrid<-syn_tst_err_hyb[idx_comp,]*hym_kcfs_conv
syn_4c_err_hybrid<-syn_4c_err_hyb[idx_comp,]*hym_kcfs_conv

#barplot
mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

mth<-4

err_lst_tst<-c()
err_lst_wm<-c()

seas<-which(ixx_comp$mon==(mth-1))

err_tst<-err_hist[seas]
err_wm<-err_4c[seas]

syn_err_tst<-median(as.vector(syn_tst_err_hybrid[seas,]))
bench_syn_err_tst<-median(as.vector(syn_tst_err_static[seas,]))

syn_err_4c<-median(as.vector(syn_4c_err_hybrid[seas,]))
bench_syn_err_4c<-median(as.vector(syn_4c_err_static[seas,]))

scale<-0.5

err_lst_tst[1]<-median(err_tst)*scale
err_lst_tst[2]<-bench_syn_err_tst*scale
err_lst_tst[3]<-syn_err_tst*scale

err_lst_tst[4]<-cp_tst_bench[3,mth]/100
err_lst_tst[5]<-cp_tst[3,mth]/100

err_lst_tst[6]<-cp_tst_bench[2,mth]/100
err_lst_tst[7]<-cp_tst[2,mth]/100

err_lst_tst[8]<-cp_tst_bench[1,mth]/100
err_lst_tst[9]<-cp_tst[1,mth]/100

err_lst_wm[1]<-median(err_wm)*scale
err_lst_wm[2]<-bench_syn_err_4c*scale
err_lst_wm[3]<-syn_err_4c*scale

err_lst_wm[4]<-cp_4c_bench[3,mth]/100
err_lst_wm[5]<-cp_4c[3,mth]/100

err_lst_wm[6]<-cp_4c_bench[2,mth]/100
err_lst_wm[7]<-cp_4c[2,mth]/100

err_lst_wm[8]<-cp_4c_bench[1,mth]/100
err_lst_wm[9]<-cp_4c[1,mth]/100

png(paste('h:/oroville_non-stationary/paper/figs_rev1/si/covg-prob_',wy_tag,'_',mths[mth],'_v-',vers,'_nreg=',noise_reg,'.png',sep=''),height=768,width=768)
my.col<-alpha('papayawhip',alpha = 0.3)

par(mfrow=c(2,1),cex.lab=2,cex.axis=1.5,mar=c(3,4,1,4),mgp=c(2,0.75,0))
barplot(err_lst_tst,space=c(0,0.1,0.1,1,0.1,1,0.1,1,0.1),col=c('skyblue4','aquamarine1','turquoise3','aquamarine1','turquoise3','aquamarine1','turquoise3','aquamarine1','turquoise3'),
        ylab='Median Error (kcfs)',ylim=c(-.5,1),border=c(rep('black',3),rep('gray50',6)),axes=F)
box(which='plot')
axis(1,at=c(1.6,5.25,8.25,11.5),labels = c('error','0.90','0.50','0.10'))
axis(2,at=seq(-5,5,1)*scale,labels = seq(-5,5,1))
axis(4,at=seq(0,1,0.2),las=3,labels=seq(0,1,.2),col='gray50',col.axis='gray50')
mtext('Coverage Probability',side=4,line=2.5,cex=2,col='gray50',adj=1)
polygon(c(-1.3,-1.3,3.75,3.75),c(-0.5,1,1,-0.5),col=my.col)
text(12,0.8,mths[mth],font=4,cex=3,adj=1)
text(8.25,-0.4,'Test',font=2,cex=3,adj=0.5)
text(6.75,0.85,'0.90',col='gray70',cex=1.5)
text(6.75,0.45,'0.50',col='gray70',cex=1.5)
text(6.75,0.05,'0.10',col='gray70',cex=1.5)
abline(h=0,lty=2)
segments(3.75,0.9,13.25,0.9,col='gray70',lty=2,lwd=2)
segments(3.75,0.5,13.25,0.5,col='gray70',lty=2,lwd=2)
segments(3.75,0.1,13.25,0.1,col='gray70',lty=2,lwd=2)
legend('topleft',c('Raw','Static','Hybrid'),lwd=c(6,6,6),col=c('skyblue4','aquamarine1','turquoise3'),cex=2)

barplot(err_lst_wm,space=c(0,0.1,0.1,1,0.1,1,0.1,1,0.1),col=c('darkorange3','palevioletred1','coral1','palevioletred1','coral1','palevioletred1','coral1','palevioletred1','coral1'),
        ylab='Median Error (mm)',ylim=c(-0.5,1),border=c(rep('black',3),rep('gray50',6)),axes=F)
box(which='plot')
axis(1,at=c(1.6,5.25,8.25,11.5),labels = c('error','0.90','0.50','0.10'))
axis(2,at=seq(-5,5,1)*scale,labels = seq(-5,5,1))
axis(4,at=seq(0,1,0.2),las=3,labels=seq(0,1,.2),col='gray50',col.axis='gray50')
mtext('Coverage Probability',side=4,line=2.5,cex=2,col='gray50',adj=1)
polygon(c(-1.3,-1.3,3.75,3.75),c(-0.5,1,1,-0.5),col=my.col)
text(12,0.8,mths[mth],font=4,cex=3,adj=1)
text(8.25,-0.4,'Test+4C',font=2,cex=3,adj=0.5)
text(6.75,0.85,'0.90',col='gray70',cex=1.5)
text(6.75,0.45,'0.50',col='gray70',cex=1.5)
text(6.75,0.05,'0.10',col='gray70',cex=1.5)
abline(h=0,lty=2)
segments(3.75,0.9,13.25,0.9,col='gray70',lty=2,lwd=2)
segments(3.75,0.5,13.25,0.5,col='gray70',lty=2,lwd=2)
segments(3.75,0.1,13.25,0.1,col='gray70',lty=2,lwd=2)
legend('topleft',c('Raw','Static','Hybrid'),lwd=c(6,6,6),col=c('darkorange3','palevioletred1','coral1'),cex=2)

dev.off()

################################END##############################