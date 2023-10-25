setwd('z:/oro_nonstat/')
#setwd('z:/oro_nonstat/')

sma_site<-'NHG'
samp_type<-'split' # 'split' 'skip'

sma_predmat_hist<-readRDS(paste('data_rev1/sma_predmat_hist_',sma_site,'.rds',sep=''))

source('mm-cfs_conversion.R')
sma_kcfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs/1000

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

if(samp_type=='split'){
  idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
  idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
  idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
  idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')
  
  ix_cal<-seq(as.Date('1988-10-01'),as.Date('1998-09-30'),'day')
  ixx_cal<-as.POSIXlt(ix_cal)
  ix_val<-seq(as.Date('1998-10-01'),as.Date('2004-09-30'),'day')
  ixx_val<-as.POSIXlt(ix_val)
  ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
  ixx_tst<-as.POSIXlt(ix_tst)
  ix_trn<-seq(as.Date('1988-10-01'),as.Date('2004-09-30'),'day')
  ixx_trn<-as.POSIXlt(ix_trn)
}

if(samp_type=='skip'){
  seq_cal<-seq(1988,2016,3)
  ix_cal<-c()
  ixx_cal<-c()
  idx_cal<-c()
  for(i in 1:length(seq_cal)){
    ix_cal<-as.Date(c(ix_cal,seq(as.Date(paste(seq_cal[i],'-10-01',sep='')),as.Date(paste(seq_cal[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
    ixx_cal<-c(ixx_cal,as.POSIXlt(ix_cal))
    idx_cal<-c(idx_cal,which(ix==paste(seq_cal[i],'-10-01',sep='')):which(ix==paste(seq_cal[i]+1,'-09-30',sep='')))
  }
  
  seq_val<-seq(1989,2018,3)
  ix_val<-c()
  ixx_val<-c()
  idx_val<-c()
  for(i in 1:length(seq_val)){
    ix_val<-as.Date(c(ix_val,seq(as.Date(paste(seq_val[i],'-10-01',sep='')),as.Date(paste(seq_val[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
    ixx_val<-c(ixx_val,as.POSIXlt(ix_val))
    idx_val<-c(idx_val,which(ix==paste(seq_val[i],'-10-01',sep='')):which(ix==paste(seq_val[i]+1,'-09-30',sep='')))
  }
  
  seq_tst<-seq(1990,2018,3)
  ix_tst<-c()
  ixx_tst<-c()
  idx_tst<-c()
  for(i in 1:length(seq_tst)){
    ix_tst<-as.Date(c(ix_tst,seq(as.Date(paste(seq_tst[i],'-10-01',sep='')),as.Date(paste(seq_tst[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
    ixx_tst<-c(ixx_tst,as.POSIXlt(ix_tst))
    idx_tst<-c(idx_tst,which(ix==paste(seq_tst[i],'-10-01',sep='')):which(ix==paste(seq_tst[i]+1,'-09-30',sep='')))
  }
  
  seq_trn<-c(sort(c(seq_cal,seq_val)))
  ix_trn<-c()
  ixx_trn<-c()
  idx_trn<-c()
  for(i in 1:length(seq_trn)){
    ix_trn<-as.Date(c(ix_trn,seq(as.Date(paste(seq_trn[i],'-10-01',sep='')),as.Date(paste(seq_trn[i]+1,'-09-30',sep='')),'day')),origin='1970-01-01')
    ixx_trn<-c(ixx_trn,as.POSIXlt(ix_trn))
    idx_trn<-c(idx_trn,which(ix==paste(seq_trn[i],'-10-01',sep='')):which(ix==paste(seq_trn[i]+1,'-09-30',sep='')))
  }
}

err_hist<-sma_predmat_hist[,'err 0']
#--------------------------------------------------------------
#cal-val-test
png(paste('h:/oroville_non-stationary/paper/figs_rev1/cal-val-test_',sma_site,'_',samp_type,'.png',sep=''),width=768,height=512)
par(mfrow=c(1,1),mar=c(2,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

#cal
err_lst<-vector('list',36)

for(i in 1:12){
  seas_cal<-which(ixx_cal$mon==(i-1))
  seas_val<-which(ixx_val$mon==(i-1))
  seas_tst<-which(ixx_tst$mon==(i-1))
  err_lst[[(i*3-2)]]<-err_hist[idx_cal][seas_cal]*sma_kcfs_conv
  err_lst[[(i*3-1)]]<-err_hist[idx_val][seas_val]*sma_kcfs_conv
  err_lst[[(i*3)]]<-err_hist[idx_tst][seas_tst]*sma_kcfs_conv
}

ylm=c(-7.5,7.5)
if(sma_site=='NHG'){ylm<-c(-.3,.3)}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=3)+rep(c(-.1,.1,.3),12),
        names = 1:24,
        las = 2,
        boxwex = 0.15,
        notch = F,
        col = rep(c('turquoise1','turquoise3','skyblue4'),12),
        ylim = ylm,
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(ylm[1],ylm[2],ylm[2]/3),labels = seq(ylm[1],ylm[2],ylm[2]/3) )
box(which='plot')
legend(8.5,ylm[1]/2,c('cal','val','test'),lwd=c(6,6,6),col=c('turquoise1','turquoise3','skyblue4'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(11,ylm[2],sma_site,cex=3,font=2,adj=0)

dev.off()

#train-test
png(paste('h:/oroville_non-stationary/paper/figs_rev1/train-test_',sma_site,'_',samp_type,'.png',sep=''),width=768,height=512)
par(mfrow=c(1,1),mar=c(2,4,1,1),mgp=c(2,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=1.75)

err_lst<-vector('list',24)

for(i in 1:12){
  seas_trn<-which(ixx_trn$mon==(i-1))
  seas_tst<-which(ixx_tst$mon==(i-1))
  err_lst[[(i*2-1)]]<-err_hist[idx_trn][seas_trn]*sma_kcfs_conv
  err_lst[[(i*2)]]<-err_hist[idx_tst][seas_tst]*sma_kcfs_conv
}

mths<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

boxplot(err_lst,
        main = '',
        at = rep(1:12,each=2)+rep(c(0,.25),12),
        names = 1:24,
        las = 2,
        boxwex = 0.2,
        notch = F,
        col = rep(c('turquoise3','skyblue4'),12),
        ylim = ylm,
        xlim = c(1,12.2),
        xlab='',
        ylab='Error (kcfs)',
        outline = F,
        axes=F
)
axis(1,at=seq(1.12,12.12,1),labels = mths)
axis(2,at=seq(ylm[1],ylm[2],ylm[2]/3),labels = seq(ylm[1],ylm[2],ylm[2]/3) )
box(which='plot')
legend(8.5,ylm[1]/2,c('train','test'),lwd=c(6,6),col=c('turquoise3','skyblue4'),cex=2,bty='n')
abline(h=0,lty=2,col='red')
text(11,ylm[2],sma_site,cex=3,font=2,adj=0)

dev.off()


###################################END#########################################