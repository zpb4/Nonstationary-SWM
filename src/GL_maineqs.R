#Main equations for Generalized Likelihood Function, Schoups & Vrugt (2010)

source('GL_subeqs.R') #bring in subequations in terms of parameters below
library(fGarch)

#1) GL function with linear models for each parameter
#pars in order: 
#par[1]: sigma_0, intercept for linearly scaled SD (heteroscedastic)
#par[2]: sigma_1, linear coefficient for scale SD
#par[3]: beta_0, intercept - kurtosis parameter (-1,1) for normalized SEP
#par[4]: beta_1, slope - kurtosis parameter (-1,1) for normalized SEP
#par[5]: xi_0, intercept -  skewness parameter (0.1,10) for normalized SEP 
#par[6]: xi_1, slope -  skewness parameter (0.1,10) for normalized SEP 
#par[7]: phi_0, intercept - ar1 coefficient (0,1)
#par[8]: phi_1, slope - ar1 coefficient (0,1)

#with ar1 model
GL_fun_mv_ar1_lin<-function(pars,sig_var,beta_var,xi_var,phi_var,et,noise,neg){
  n<-length(et)
  par<-matrix(pars,nrow=4)
  svar_len<-length(pars)/4
  
  srt_err<-sort(abs(et))
  srt_err<-srt_err[srt_err>0]
  min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  #sig_t<-sig_t[1,]
  sig_t_viol<-F
  if(any(sig_t<=0)==T){sig_t_viol<-T}
  sig_t[sig_t<=0]<-1
  
  beta_t<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  #beta[beta<(-0.99)]<-(-0.99)#;beta[beta>1]<-1
  #beta<-beta[1,]
  
  xi_in<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  #xi[xi<0.1]<-0.1;xi[xi>10]<-10
  #xi_int<-xi[1,] #approximately linearizes skewness
  xi_t<-10^xi_in
  #xi<-xi[1,]
  
  phi_t<-par[4,1]+par[4,2:svar_len]%*%t(phi_var)
  #phi[phi<0]<-0;phi[phi>1]<-1
  #phi<-phi[1,]
  
  sig_xi<-sigma_xi(M1(beta_t),M2,xi_t) #Eq A8 from subequations
  om_b<-omega_beta(beta_t) #Eq A2 from subequations
  cb<-c_beta(beta_t) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  at<-a_t(et,phi_t,sig_t)
  at<-at+at*noise #apply noise regularization
  a_xt<-c()
  for(k in 1:n){
    a_xt[k]<-a_xi_t(xi_t[k],mu_xi(M1(beta_t[k]),xi_t[k]),sig_xi[k],at[k])
  }
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-c()
  for(i in 1:n){
    gl_ll[i]<-log((2*sig_xi[i]*om_b[i])/(xi_t[i] + xi_t[i]^-1)) - sum(log(sig_t[i])) - cb[i] * sum(abs(a_xt[i])^(2/(1+beta_t[i])))
  }
 
  ll<-sum(gl_ll)
  
  if(is.na(ll)==T){ll<-(-10e9)} #optimizer can't handle NAs
  if(ll == Inf|ll == -Inf){ll<-(-10e9)} #optimizer can't handle Inf or -Inf values
  if(any(sig_t<min_sig)==T){ll<-(-10e9)}  #penalty for phi values outside acceptable range (0-1)
  if(sig_t_viol==T){ll<-(-10e9)}
  if(any(xi_t<0.1|xi_t>10)==T){ll<-(-10e9)} #penalty for xi values outside acceptable range (0.1-10)
  if(any(beta_t<(-0.99))==T){ll<-(-10e9)} #penalty for beta values outside acceptable range (<(-0.99))
  if(any(phi_t<0|phi_t>1)==T){ll<-(-10e9)}  #penalty for phi values outside acceptable range (0-1)

  if(neg==F){
    return(ll) 
  }
  
  if(neg==T){
    return(-ll) 
  }
}

#with ar1 model, output param values
GL_fun_mv_ar1_lin_params<-function(pars,sig_var,beta_var,xi_var,phi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=4)
  svar_len<-length(pars)/4
  
  srt_err<-sort(abs(et))
  srt_err<-srt_err[srt_err>0]
  min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t[sig_t<min_sig]<-min_sig
  sig_t<-sig_t[1,]
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99)#;beta[beta>1]<-1
  beta<-beta[1,]
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  #xi_int<-xi[1,] #approximately linearizes skewness
  #xi<-10^xi_in
  #xi[xi<0.1]<-0.1;xi[xi>10]<-10
  xi[xi<-1]<-(-1);xi[xi>1]<-1
  xi<-xi[1,]
  
  phi<-par[4,1]+par[4,2:svar_len]%*%t(phi_var)
  phi[phi<0]<-0;phi[phi>1]<-1
  phi<-phi[1,]
  
  return(list(sig_t,beta,xi,phi))
}

#no ar1 model
GL_fun_mv_noar1<-function(pars,sig_var,beta_var,xi_var,et,noise,neg){
  n<-length(et)
  par<-matrix(pars,nrow=3)
  svar_len<-length(pars)/3
  
  srt_err<-sort(abs(et))
  srt_err<-srt_err[srt_err>0]
  min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  #sig_t<-sig_t[1,]
  sig_t_viol<-F
  if(any(sig_t<=0)==T){sig_t_viol<-T}
  sig_t[sig_t<=0]<-1
  
  beta_t<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  #beta[beta<(-0.99)]<-(-0.99)#;beta[beta>1]<-1
  #beta<-beta[1,]
  
  xi_in<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  #xi[xi<0.1]<-0.1;xi[xi>10]<-10
  #xi_int<-xi[1,] #approximately linearizes skewness
  xi_t<-10^xi_in
  #xi<-xi[1,]
  
  #phi_t<-par[4,1]+par[4,2:svar_len]%*%t(phi_var)
  #phi[phi<0]<-0;phi[phi>1]<-1
  #phi<-phi[1,]
  
  sig_xi<-sigma_xi(M1(beta_t),M2,xi_t) #Eq A8 from subequations
  om_b<-omega_beta(beta_t) #Eq A2 from subequations
  cb<-c_beta(beta_t) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  at<-et/
  at<-at+at*noise #apply noise regularization
  a_xt<-c()
  for(k in 1:n){
    a_xt[k]<-a_xi_t(xi_t[k],mu_xi(M1(beta_t[k]),xi_t[k]),sig_xi[k],at[k])
  }
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-c()
  for(i in 1:n){
    gl_ll[i]<-log((2*sig_xi[i]*om_b[i])/(xi_t[i] + xi_t[i]^-1)) - sum(log(sig_t[i])) - cb[i] * sum(abs(a_xt[i])^(2/(1+beta_t[i])))
  }
  
  ll<-sum(gl_ll)
  
  if(is.na(ll)==T){ll<-(-10e9)} #optimizer can't handle NAs
  if(ll == Inf|ll == -Inf){ll<-(-10e9)} #optimizer can't handle Inf or -Inf values
  if(any(sig_t<min_sig)==T){ll<-(-10e9)}  #penalty for phi values outside acceptable range (0-1)
  if(sig_t_viol==T){ll<-(-10e9)}
  if(any(xi_t<0.1|xi_t>10)==T){ll<-(-10e9)} #penalty for xi values outside acceptable range (0.1-10)
  if(any(beta_t<(-0.99))==T){ll<-(-10e9)} #penalty for beta values outside acceptable range (<(-0.99))
  if(any(phi_t<0|phi_t>1)==T){ll<-(-10e9)}  #penalty for phi values outside acceptable range (0-1)
  
  if(neg==F){
    return(ll) 
  }
  
  if(neg==T){
    return(-ll) 
  }
}
#no ar1 model but db
GL_fun_mv_lin_db<-function(pars,sig_var,beta_var,xi_var,mu_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=4)
  svar_len<-length(pars)/4
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t<-sig_t[1,]
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99);#beta[beta>1]<-1
  beta<-beta[1,]
  
  #beta<-exp(par[2,1]+par[2,2:12]%*%t(beta_var))
  #beta[beta<(-0.99)]<-(-0.99);#beta[beta>1]<-1
  #beta<-beta[1,]
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi[xi<0.1]<-0.1;xi[xi>10]<-10
  xi<-xi[1,]
  #xi<-rep(xi,n)
  
  #xi<-exp(par[3,1]+par[3,2:12]%*%t(xi_var))
  #xi[xi<0.1]<-0.1;xi[xi>10]<-10
  #xi<-xi[1,]
  
  mu<-par[4,1]+par[4,2:svar_len]%*%t(mu_var)
  mu<-mu[1,]
  
  et_db<-et-mu
  
  sig_xi<-sigma_xi(M1(beta),M2,xi) #Eq A8 from subequations
  om_b<-omega_beta(beta) #Eq A2 from subequations
  cb<-c_beta(beta) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  at<-et_db/sig_t
  a_xt<-c()
  for(k in 1:n){
    a_xt[k]<-a_xi_t(xi[k],mu_xi(M1(beta[k]),xi[k]),sig_xi[k],at[k])
  }
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-c()
  for(i in 1:n){
    gl_ll[i]<-log((2*sig_xi[i]*om_b[i])/(xi[i] + xi[i]^-1)) - log(sig_t[i]) - cb[i] * abs(a_xt[i])^(2/(1+beta[i]))
  }
  
  ll<-sum(gl_ll)
  
  if (is.na(ll)==T){ll<-(-10e9)}
  if (ll == Inf|ll == -Inf){ll<-(-10e9)} #optimizer can't handle Inf or -Inf values
  
  return(ll) #mult by -1 to enable maximization
}

#no ar1 model
GL_fun_mv_lin<-function(pars,sig_var,beta_var,xi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=3)
  svar_len<-length(pars)/3
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t<-sig_t[1,]
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99);#beta[beta>1]<-1
  beta<-beta[1,]
  
  #beta<-exp(par[2,1]+par[2,2:12]%*%t(beta_var))
  #beta[beta<(-0.99)]<-(-0.99);#beta[beta>1]<-1
  #beta<-beta[1,]
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi[xi<0.1]<-0.1;xi[xi>10]<-10
  xi<-xi[1,]
  #xi<-rep(xi,n)
  
  #xi<-exp(par[3,1]+par[3,2:12]%*%t(xi_var))
  #xi[xi<0.1]<-0.1;xi[xi>10]<-10
  #xi<-xi[1,]
  
  sig_xi<-sigma_xi(M1(beta),M2,xi) #Eq A8 from subequations
  om_b<-omega_beta(beta) #Eq A2 from subequations
  cb<-c_beta(beta) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  at<-et/sig_t
  a_xt<-c()
  for(k in 1:n){
    a_xt[k]<-a_xi_t(xi[k],mu_xi(M1(beta[k]),xi[k]),sig_xi[k],at[k])
  }
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-c()
  for(i in 1:n){
    gl_ll[i]<-log((2*sig_xi[i]*om_b[i])/(xi[i] + xi[i]^-1)) - log(sig_t[i]) - cb[i] * abs(a_xt[i])^(2/(1+beta[i]))
  }

  ll<-sum(gl_ll)
  
  if (is.na(ll)==T){ll<-(-10e9)}
  if (ll == Inf|ll == -Inf){ll<-(-10e9)} #optimizer can't handle Inf or -Inf values
  
  return(ll) #mult by -1 to enable maximization
}


#2) Function to output decorrelated, normalized a_t

#with ar1 model
at_mv_ar1_lin<-function(pars,sig_var,beta_var,xi_var,phi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=4)
  svar_len<-length(pars)/4
  
  srt_err<-sort(abs(et))
  srt_err<-srt_err[srt_err>0]
  min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t[sig_t<min_sig]<-min_sig
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99)
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi<-10**xi
  xi[xi<0.1]<-0.1;xi[xi>10]<-10
  
  phi<-par[4,1]+par[4,2:svar_len]%*%t(phi_var)
  phi[phi<0]<-0;phi[phi>1]<-1
  
  ar_resid<-et-phi*c(0,et[-c(length(et))])
  norm_resid<-t(ar_resid/sig_t)
  
  return(norm_resid)
}

#no ar1 model
at_mv_lin<-function(pars,sig_var,beta_var,xi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=3)
  svar_len<-length(pars)/3
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99)#;beta[beta>1]<-1
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi[xi<0.1]<-0.1;xi[xi>10]<-10

  norm_resid<-et/t(sig_t)
  
  return(norm_resid)
}

#3) Generate from fitted AR1 - linear modeled parameters
#pars in order: 
#par[1]: sigma_0, intercept for linearly scaled SD (heteroscedastic)
#par[2]: sigma_1, linear coefficient for scale SD
#par[3]: beta_0, intercept - kurtosis parameter (-1,1) for normalized SEP
#par[4]: beta_1, slope - kurtosis parameter (-1,1) for normalized SEP
#par[5]: xi_0, intercept -  skewness parameter (0.1,10) for normalized SEP 
#par[6]: xi_1, slope -  skewness parameter (0.1,10) for normalized SEP 
#par[7]: phi_0, intercept - ar1 coefficient (0,1)
#par[8]: phi_1, slope - ar1 coefficient (0,1)

#with ar1 model
syn_gen_mv_ar1_lin<-function(pars,sig_var,beta_var,xi_var,phi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=4)
  svar_len<-length(pars)/4
  
  srt_err<-sort(abs(et))
  srt_err<-srt_err[srt_err>0]
  min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t[sig_t<min_sig]<-min_sig
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99)#;beta[beta>1]<-1
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi<-10^xi
  xi[xi<0.1]<-0.1;xi[xi>10]<-10
  
  phi<-par[4,1]+par[4,2:svar_len]%*%t(phi_var)
  phi[phi<0]<-0;phi[phi>1]<-1
  
  res<-c()
  
  for(k in 1:n){
    res[k]<-rsged(1,mean=0,sd=1,nu=(2 / (1 + beta[k])),xi=xi[k])
  }

  res_scale<-res*sig_t
  
  res_scale<-c(0,res_scale)
  
  syn_et<-rep(0,(n+1))
  syn_et[1]<-sample(res_scale,1)
  
  for(i in 1:n){
    syn_et[i+1]<-syn_et[i]*phi[i]+res_scale[i+1]
  }
  
  return(syn_et[2:(n+1)])
}

#no ar1 model
syn_gen_mv_lin<-function(pars,sig_var,beta_var,xi_var,et){
  n<-length(et)
  par<-matrix(pars,nrow=3)
  svar_len<-length(pars)/3
  
  sig_t<-par[1,1]+par[1,2:svar_len]%*%t(sig_var) #Eq 5 from subequations
  sig_t<-sig_t[1,]
  
  beta<-par[2,1]+par[2,2:svar_len]%*%t(beta_var)
  beta[beta<(-0.99)]<-(-0.99)
  beta<-beta[1,]
  
  xi<-par[3,1]+par[3,2:svar_len]%*%t(xi_var)
  xi<-10^xi
  xi[xi<0.1]<-0.1;xi[xi>10]<-10
  xi<-xi[1,]
  
  res<-c()
  
  for(k in 1:n){
    res[k]<-rsged(1,mean=0,sd=1,nu=(2 / (1 + beta[k])),xi=xi[k])
  }
  
  res_scale<-res*sig_t
  
  return(res_scale)
}

#4) define SEP density as function of calculated xi and beta and random variable axt (a_xi_t)
SEP_dens<-function(xi,beta,axt){
  SEP_dens<-(2*sigma_xi(M1(beta),M2,xi))/(xi + xi^(-1))*omega_beta(beta)*
    exp(-c_beta(beta)*abs(a_xi_t(xi,mu_xi(M1(beta),xi),sigma_xi(M1(beta),M2,xi),axt))^(2/(1+beta)))
}

