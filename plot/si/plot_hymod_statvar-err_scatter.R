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

hym_hist_vars<-readRDS(paste('data_rev1/hym_hist_vars_',sma_site,'.rds',sep=''))
hym_4c_vars<-readRDS(paste('data_rev1/hym_4c_vars_',sma_site,'.rds',sep=''))
hym_predmat_hist<-readRDS(paste('data_rev1/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('data_rev1/hym_predmat_4c_',hym_site,'.rds',sep=''))

sma_hist_vars<-readRDS(paste('data_rev1/sma_hist_vars_',sma_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_trn<-seq(as.Date('1988-10-01'),as.Date('2004-09-30'),'day')
ixx_trn<-as.POSIXlt(ix_trn)
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')

obs_trn<-sma_hist_vars[idx_trn,'obs']
trn_wy<-1989:2004
obs_sum<-c()

for(i in 1:length(trn_wy)){
  obs_sum[i]<-sum(obs_trn[which(ix_trn==as.Date(paste(trn_wy[i]-1,'-10-01',sep=''))):which(ix_trn==as.Date(paste(trn_wy[i],'-09-30',sep='')))])
}

srt_obs<-sort(obs_sum,index.return=T)
trn_wy[tail(srt_obs$ix,5)] #5 wettest
#[1] 1996 1993 1997 1998 1995
trn_wy[head(srt_obs$ix,5)] #5 driest
#[1] 1994 1992 1991 1990 2001

wy_comp_trn<-c(1993,1995,1996,1997,1998)
wy_comp1_trn<-c(1990,1991,1992,1994,2001)


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

#tag0_trn
ix_comp_trn<-c()
idx_comp_trn<-c()
for(i in 1:length(wy_comp_trn)){
  ix_comp_trn<-c(ix_comp_trn,seq(as.Date(paste(wy_comp_trn[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp_trn[i],'-09-30',sep='')),'day'))
  idx_comp_trn<-c(idx_comp_trn,which(ix==paste(wy_comp_trn[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp_trn[i],'-09-30',sep='')))
}
ix_comp_trn<-as.Date(ix_comp_trn,origin = '1970-01-01')
ixx_comp_trn<-as.POSIXlt(ix_comp_trn)

#tag1_trn
ix_comp1_trn<-c()
idx_comp1_trn<-c()
for(i in 1:length(wy_comp1_trn)){
  ix_comp1_trn<-c(ix_comp1_trn,seq(as.Date(paste(wy_comp1_trn[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp1_trn[i],'-09-30',sep='')),'day'))
  idx_comp1_trn<-c(idx_comp1_trn,which(ix==paste(wy_comp1_trn[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp1_trn[i],'-09-30',sep='')))
}
ix_comp1_trn<-as.Date(ix_comp1_trn,origin = '1970-01-01')
ixx_comp1_trn<-as.POSIXlt(ix_comp1_trn)

#data process for plot panel a
err_trn<-hym_predmat_hist[idx_comp_trn,'err 0']*hym_cfs_conv/1000 #convert to cfs
err_tst<-hym_predmat_hist[idx_comp,'err 0']*hym_cfs_conv/1000 #convert to cfs
err_4c<-hym_predmat_4c[idx_comp,'err 0']*hym_cfs_conv/1000    #convert to kcfs
err_trn1<-hym_predmat_hist[idx_comp1_trn,'err 0']*hym_cfs_conv/1000 #convert to cfs
err_tst1<-hym_predmat_hist[idx_comp1,'err 0']*hym_cfs_conv/1000 #convert to cfs
err_4c1<-hym_predmat_4c[idx_comp1,'err 0']*hym_cfs_conv/1000    #convert to kcfs



#--------------------------------------------------------------
png(paste('h:/Projects/Nonstationary SWM/oroville_non-stationary/paper/figs_rev1/si/err-comp_trn-tst-tst4c_',wy_tag,'_',wy_tag1,'.png',sep=''),width=768,height=1024)

pos_mat<-matrix(1:2,ncol=1,byrow=T)
layout(pos_mat,heights = c(1,1))
par(mar=c(3,6,1,1),mgp=c(3.5,1,0),tcl=-0.3,cex.lab=3,cex.axis=2.75)

err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  seas_trn<-which(ixx_comp_trn$mon==(i-1))
  err_lst[[(i*3-2)]]<-err_trn[seas_trn]
  err_lst[[(i*3-1)]]<-err_tst[seas]
  err_lst[[(i*3)]]<-err_4c[seas]
}

#mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
mths<-c('J','F','M','A','M','J','J','A','S','O','N','D')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=3)+rep(c(-0.2,0,.2),12),
        names = 1:36,
        las = 2,
        boxwex = .15,
        notch = F,
        col = rep(c('skyblue2','skyblue4','darkorange3'),12),
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
legend('topright',c('Train','Test','Test+4C'),lwd=c(6,6),col=c('skyblue2','skyblue4','darkorange3'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(1,4.75,'a)',cex=3,font=2)
text(11.5,-4.75,wy_tag,font=2,cex=3)

#wy_tag2
err_lst<-vector('list',36)

for(i in 1:12){
  seas<-which(ixx_comp1$mon==(i-1))
  seas_trn<-which(ixx_comp1_trn$mon==(i-1))
  err_lst[[(i*3-2)]]<-err_trn1[seas_trn]
  err_lst[[(i*3-1)]]<-err_tst1[seas]
  err_lst[[(i*3)]]<-err_4c1[seas]
}

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=3)+rep(c(-0.2,0,.2),12),
        names = 1:36,
        las = 2,
        boxwex = .15,
        notch = F,
        col = rep(c('skyblue2','skyblue4','darkorange3'),12),
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
#legend('topright',c('Test','Test+4C'),lwd=c(6,6),col=c('skyblue4','darkorange3'),cex=3,bty='n')
abline(h=0,lty=2,col='red')
text(1,4.75,'b)',cex=3,font=2)
text(11.5,-4.75,wy_tag1,font=2,cex=3)


dev.off()

################################END###########################################
