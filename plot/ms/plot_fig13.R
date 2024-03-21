setwd('z:/oro_nonstat/')
source('mm-cfs_conversion.R')

sma_site1<-'ORO'
sma_site2<-'SHA'
sma_site3<-'NHG'
gen_period<-'hist-all'
samp_type<-'split'
vers<-'err13'
noise_reg<-T
n<-1000

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp1<-c(2005,2006,2011,2016,2017)
wy_tag1<-'5wet'

wy_comp2<-c(2007,2008,2012,2014,2015)
wy_tag2<-'5dry'

ylm_s1<-c(-5,5) 
ylm_s2<-c(-5,5)
ylm_s3<-c(-0.25,0.25)

y_axis_s1<-seq(-5,5,2.5) 
y_axis_s2<-seq(-5,5,2.5)
y_axis_s3<-seq(-0.25,0.25,0.1)

s1_kcfs_conv<-as.numeric(area_calc[sma_site1,2])*mm_to_cfs/1000
s1_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site1,'.rds',sep=''))

s2_kcfs_conv<-as.numeric(area_calc[sma_site2,2])*mm_to_cfs/1000
s2_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site2,'.rds',sep=''))

s3_kcfs_conv<-as.numeric(area_calc[sma_site3,2])*mm_to_cfs/1000
s3_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site3,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_comp1<-c()
idx_comp1<-c()
for(i in 1:length(wy_comp1)){
  ix_comp1<-c(ix_comp1,seq(as.Date(paste(wy_comp1[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp1[i],'-09-30',sep='')),'day'))
  idx_comp1<-c(idx_comp1,which(ix==paste(wy_comp1[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp1[i],'-09-30',sep='')))
}
ix_comp1<-as.Date(ix_comp1,origin = '1970-01-01')
ixx_comp1<-as.POSIXlt(ix_comp1)

ix_comp2<-c()
idx_comp2<-c()
for(i in 1:length(wy_comp2)){
  ix_comp2<-c(ix_comp2,seq(as.Date(paste(wy_comp2[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp2[i],'-09-30',sep='')),'day'))
  idx_comp2<-c(idx_comp2,which(ix==paste(wy_comp2[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp2[i],'-09-30',sep='')))
}
ix_comp2<-as.Date(ix_comp2,origin = '1970-01-01')
ixx_comp2<-as.POSIXlt(ix_comp2)

syn_err_s1<-readRDS(paste('out_rev1/sacsma_syn-err_',sma_site1,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
syn_err_s2<-readRDS(paste('out_rev1/sacsma_syn-err_',sma_site2,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
syn_err_s3<-readRDS(paste('out_rev1/sacsma_syn-err_',sma_site3,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))

err_hist_s1_col1<-s1_predmat_hist[idx_comp1,'err 0']
err_hist_s1_col2<-s1_predmat_hist[idx_comp2,'err 0']
err_hist_s2_col1<-s2_predmat_hist[idx_comp1,'err 0']
err_hist_s2_col2<-s2_predmat_hist[idx_comp2,'err 0']
err_hist_s3_col1<-s3_predmat_hist[idx_comp1,'err 0']
err_hist_s3_col2<-s3_predmat_hist[idx_comp2,'err 0']


syn_err_s1_col1<-syn_err_s1[idx_comp1,]
syn_err_s1_col2<-syn_err_s1[idx_comp2,]
syn_err_s2_col1<-syn_err_s2[idx_comp1,]
syn_err_s2_col2<-syn_err_s2[idx_comp2,]
syn_err_s3_col1<-syn_err_s3[idx_comp1,]
syn_err_s3_col2<-syn_err_s3[idx_comp2,]

#--------------------------------------------------------------------------------------------------------------
#plot raw v synthetic residuals

png(paste('h:/oroville_non-stationary/paper/figs_rev1/fig13/fig13_v-',vers,'_nreg=',noise_reg,'.png',sep=''),width=768,height=1024)

par(mfrow=c(3,2),oma=c(1,4,4,0),mar=c(2,5,1,1),mgp=c(3,1,0),tcl=-0.2,cex.lab=2,cex.axis=2)

  
mths<-c('J','F','M','A','M','J','J','A','S','O','N','D')

#row 1, col1
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp1$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s1_col1[seas]*s1_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s1_col1[seas,]*s1_kcfs_conv
}

boxplot(res_lst,
      main = '',
      at = rep(1:12,each=2)+rep(c(0,.25),12),
      names = 1:24,
      las = 2,
      boxwex = .2,
      notch = F,
      col = rep(c('skyblue4','turquoise3'),12),
      ylim = ylm_s1,
      xlim = c(1,12.2),
      xlab='',
      ylab='Error (kcfs)',
      outline = F,
      axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=y_axis_s1,labels = y_axis_s1 )
box(which='plot')
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 
mtext(sma_site1,side=2,line=5,cex=3,font=2)
mtext('5wet',side=3,line=1.5,cex=3,font=2)

#row 1, col2
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp2$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s1_col2[seas]*s1_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s1_col2[seas,]*s1_kcfs_conv
}

boxplot(res_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = ylm_s1,
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=y_axis_s1,labels = F )
box(which='plot')
legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 
mtext('5dry',side=3,line=1.5,cex=3,font=2)
  
#row 2, col1
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp1$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s2_col1[seas]*s2_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s2_col1[seas,]*s2_kcfs_conv
}

boxplot(res_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = ylm_s2,
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=y_axis_s2,labels = y_axis_s2 )
box(which='plot')
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 
mtext(sma_site2,side=2,line=5,cex=3,font=2)

#row 2, col2
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp2$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s2_col2[seas]*s2_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s2_col2[seas,]*s2_kcfs_conv
}

boxplot(res_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = ylm_s2,
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = F)
axis(2,at=y_axis_s2,labels = F )
box(which='plot')
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 

#row 3, col1
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp1$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s3_col1[seas]*s3_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s3_col1[seas,]*s3_kcfs_conv
}

boxplot(res_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = ylm_s3,
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=y_axis_s3,labels = y_axis_s3 )
box(which='plot')
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 
mtext(sma_site3,side=2,line=5,cex=3,font=2)

#row 3, col2
res_lst<-vector('list',24)

for(i in 1:12){
  seas<-which(ixx_comp2$mon==(i-1))
  res_lst[[(i*2-1)]]<-err_hist_s3_col2[seas]*s3_kcfs_conv
  res_lst[[(i*2)]]<-syn_err_s3_col2[seas,]*s3_kcfs_conv
}

boxplot(res_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = .2,
        notch = F,
        col = rep(c('skyblue4','turquoise3'),12),
        ylim = ylm_s3,
        xlim = c(1,12.2),
        xlab='',
        ylab='',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=y_axis_s3,labels = F )
box(which='plot')
#legend('topright',c('Empirical','Simulated'),lwd=c(6,6),col=c('skyblue4','turquoise3'),cex=2)
abline(h=0,lty=2,col='red') 

dev.off()

################################END###########################################