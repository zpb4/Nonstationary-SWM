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

sim_hymod_tst<-hym_predmat_hist[idx_comp,'sim 0']*hym_kcfs_conv
sim_hymod_4c<-hym_predmat_4c[idx_comp,'sim 0']*hym_kcfs_conv

sim_sma_tst<-sma_predmat_hist[idx_comp,'sim 0']*hym_kcfs_conv
sim_sma_4c<-sma_predmat_4c[idx_comp,'sim 0']*hym_kcfs_conv

syn_tst_flow_static<-syn_tst_flow_static[idx_comp,]*hym_kcfs_conv
syn_4c_flow_static<-syn_4c_flow_static[idx_comp,]*hym_kcfs_conv

syn_tst_flow_hybrid<-syn_tst_flow_hyb[idx_comp,]*hym_kcfs_conv
syn_4c_flow_hybrid<-syn_4c_flow_hyb[idx_comp,]*hym_kcfs_conv

pcntile<-function(x){
  p<-sort(x)
  l05<-p[round(0.05*length(x))]
  l25<-p[round(0.25*length(x))]
  l45<-p[round(0.45*length(x))]
  u55<-p[round(0.55*length(x))]
  u75<-p[round(0.75*length(x))]
  u95<-p[round(0.95*length(x))]
  return(c(l05,l25,l45,u55,u75,u95))
}

swm_pcnt_tst<-apply(syn_tst_flow_hybrid,1,pcntile)
swm_pcnt_4c<-apply(syn_4c_flow_hybrid,1,pcntile)

bench_swm_pcnt_tst<-apply(syn_tst_flow_static,1,pcntile)
bench_swm_pcnt_4c<-apply(syn_4c_flow_static,1,pcntile)

#--------------------------------------------------------
#coverage probabilities
#test
cp_idx_tst<-array(0,c(3,length(idx_comp)))

for(i in 1:length(idx_comp)){
  if(sim_sma_tst[i]>swm_pcnt_tst[4,i]|sim_sma_tst[i]<swm_pcnt_tst[3,i]){cp_idx_tst[1,i]<-1}
  if(sim_sma_tst[i]>swm_pcnt_tst[5,i]|sim_sma_tst[i]<swm_pcnt_tst[2,i]){cp_idx_tst[2,i]<-1}
  if(sim_sma_tst[i]>swm_pcnt_tst[6,i]|sim_sma_tst[i]<swm_pcnt_tst[1,i]){cp_idx_tst[3,i]<-1}
}

cp_10_tst<-(1-sum(cp_idx_tst[1,])/length(idx_comp))*100
cp_50_tst<-(1-sum(cp_idx_tst[2,])/length(idx_comp))*100
cp_90_tst<-(1-sum(cp_idx_tst[3,])/length(idx_comp))*100

cp_10_tst_mth<-c()
cp_50_tst_mth<-c()
cp_90_tst_mth<-c()

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  cp_10_tst_mth[i]<-(1-sum(cp_idx_tst[1,seas])/length(seas))*100
  cp_50_tst_mth[i]<-(1-sum(cp_idx_tst[2,seas])/length(seas))*100
  cp_90_tst_mth[i]<-(1-sum(cp_idx_tst[3,seas])/length(seas))*100
  print(paste('Month',i,'cp50',cp_50_tst_mth[i],'cp90',cp_90_tst_mth[i]))
}

#test-bench
cp_idx_tst_bench<-array(0,c(3,length(idx_comp)))

for(i in 1:length(idx_comp)){
  if(sim_sma_tst[i]>bench_swm_pcnt_tst[4,i]|sim_sma_tst[i]<bench_swm_pcnt_tst[3,i]){cp_idx_tst_bench[1,i]<-1}
  if(sim_sma_tst[i]>bench_swm_pcnt_tst[5,i]|sim_sma_tst[i]<bench_swm_pcnt_tst[2,i]){cp_idx_tst_bench[2,i]<-1}
  if(sim_sma_tst[i]>bench_swm_pcnt_tst[6,i]|sim_sma_tst[i]<bench_swm_pcnt_tst[1,i]){cp_idx_tst_bench[3,i]<-1}
}

cp_10_tst_bench<-(1-sum(cp_idx_tst_bench[1,])/length(idx_comp))*100
cp_50_tst_bench<-(1-sum(cp_idx_tst_bench[2,])/length(idx_comp))*100
cp_90_tst_bench<-(1-sum(cp_idx_tst_bench[3,])/length(idx_comp))*100

cp_10_tst_bench_mth<-c()
cp_50_tst_bench_mth<-c()
cp_90_tst_bench_mth<-c()

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  cp_10_tst_bench_mth[i]<-(1-sum(cp_idx_tst_bench[1,seas])/length(seas))*100
  cp_50_tst_bench_mth[i]<-(1-sum(cp_idx_tst_bench[2,seas])/length(seas))*100
  cp_90_tst_bench_mth[i]<-(1-sum(cp_idx_tst_bench[3,seas])/length(seas))*100
  print(paste('Month',i,'cp50',cp_50_tst_bench_mth[i],'cp90',cp_90_tst_bench_mth[i]))
}

#4c
cp_idx_4c<-array(0,c(3,length(idx_comp)))

for(i in 1:length(idx_comp)){
  if(sim_sma_4c[i]>swm_pcnt_4c[4,i]|sim_sma_4c[i]<swm_pcnt_4c[3,i]){cp_idx_4c[1,i]<-1}
  if(sim_sma_4c[i]>swm_pcnt_4c[5,i]|sim_sma_4c[i]<swm_pcnt_4c[2,i]){cp_idx_4c[2,i]<-1}
  if(sim_sma_4c[i]>swm_pcnt_4c[6,i]|sim_sma_4c[i]<swm_pcnt_4c[1,i]){cp_idx_4c[3,i]<-1}
}

cp_10_4c<-(1-sum(cp_idx_4c[1,])/length(idx_comp))*100
cp_50_4c<-(1-sum(cp_idx_4c[2,])/length(idx_comp))*100
cp_90_4c<-(1-sum(cp_idx_4c[3,])/length(idx_comp))*100

cp_10_4c_mth<-c()
cp_50_4c_mth<-c()
cp_90_4c_mth<-c()

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  cp_10_4c_mth[i]<-(1-sum(cp_idx_4c[1,seas])/length(seas))*100
  cp_50_4c_mth[i]<-(1-sum(cp_idx_4c[2,seas])/length(seas))*100
  cp_90_4c_mth[i]<-(1-sum(cp_idx_4c[3,seas])/length(seas))*100
  print(paste('Month',i,'cp50',cp_50_4c_mth[i],'cp90',cp_90_4c_mth[i]))
}


#4c - bench
cp_idx_4c_bench<-array(0,c(3,length(idx_comp)))

for(i in 1:length(idx_comp)){
  if(sim_sma_4c[i]>bench_swm_pcnt_4c[4,i]|sim_sma_4c[i]<bench_swm_pcnt_4c[3,i]){cp_idx_4c_bench[1,i]<-1}
  if(sim_sma_4c[i]>bench_swm_pcnt_4c[5,i]|sim_sma_4c[i]<bench_swm_pcnt_4c[2,i]){cp_idx_4c_bench[2,i]<-1}
  if(sim_sma_4c[i]>bench_swm_pcnt_4c[6,i]|sim_sma_4c[i]<bench_swm_pcnt_4c[1,i]){cp_idx_4c_bench[3,i]<-1}
}

cp_10_4c_bench<-(1-sum(cp_idx_4c_bench[1,])/length(idx_comp))*100
cp_50_4c_bench<-(1-sum(cp_idx_4c_bench[2,])/length(idx_comp))*100
cp_90_4c_bench<-(1-sum(cp_idx_4c_bench[3,])/length(idx_comp))*100

cp_10_4c_bench_mth<-c()
cp_50_4c_bench_mth<-c()
cp_90_4c_bench_mth<-c()

for(i in 1:12){
  seas<-which(ixx_comp$mon==(i-1))
  cp_10_4c_bench_mth[i]<-(1-sum(cp_idx_4c_bench[1,seas])/length(seas))*100
  cp_50_4c_bench_mth[i]<-(1-sum(cp_idx_4c_bench[2,seas])/length(seas))*100
  cp_90_4c_bench_mth[i]<-(1-sum(cp_idx_4c_bench[3,seas])/length(seas))*100
  print(paste('Month',i,'cp50',cp_50_4c_bench_mth[i],'cp90',cp_90_4c_bench_mth[i]))
}

cp_tst<-rbind(cp_10_tst_mth,cp_50_tst_mth,cp_90_tst_mth)
cp_tst_bench<-rbind(cp_10_tst_bench_mth,cp_50_tst_bench_mth,cp_90_tst_bench_mth)

cp_4c<-rbind(cp_10_4c_mth,cp_50_4c_mth,cp_90_4c_mth)
cp_4c_bench<-rbind(cp_10_4c_bench_mth,cp_50_4c_bench_mth,cp_90_4c_bench_mth)

saveRDS(cp_tst,paste('out_rev1/cp_tst_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(cp_tst_bench,'out_rev1/cp_tst_bench.rds')
saveRDS(cp_4c,paste('out_rev1/cp_4c_',wy_tag,'_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(cp_4c_bench,'out_rev1/cp_4c_bench.rds')

rm(list=ls());gc()

###################################END#####################################