#Main equations for Generalized Likelihood Function

source('z:/oro_nonstat/GL_subeqs_uv.R') #bring in subequations in terms of parameters below

#1) Define Generalized Likelihood function per Schoups and Vrugt (2010)
#pars in order: 
#par[1]: sigma_0, intercept for linearly scaled SD (heteroscedastic)
#par[2]: sigma_1, linear coefficient for scale SD
#par[3]: phi_1, AR(1) coefficient
#par[4]: phi_2, AR(2) coefficient
#par[5]: phi_2, AR(3) coefficient
#par[6]: beta, kurtosis parameter (-1,1) for normalized SEP
#par[7]: xi, skewness parameter (0.1,10) for normalized SEP 
#par[8]: mu_h, exponential scaling parameter based on simulated flow, set to zero for no-scaling
#for this version of function, but retained for more general version

#1a) GL function with up to 3 AR orders 
GL_fun_noscale_ar3_bfgs<-function(pars,inflow,et){
  n<-length(et)
  e_t<-et
  sig_xi<-sigma_xi(M1(pars[6]),M2,pars[7]) #Eq A8 from subequations
  om_b<-omega_beta(pars[6]) #Eq A2 from subequations
  sig_t<-sigma_t(pars[1],pars[2],inflow) #Eq 5 from subequations
  cb<-c_beta(pars[6]) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  a_xt<-a_xi_t(pars[7],mu_xi(M1(pars[6]),pars[7]),sig_xi,a_tx(e_t,pars[3],pars[4],pars[5],sig_t))
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-n*log((2*sig_xi*om_b)/(pars[7] + pars[7]^-1)) - sum(log(sig_t)) - cb * sum(abs(a_xt)^(2/(1+pars[6])))
  
  if (gl_ll == Inf|gl_ll == -Inf) {gl_ll<-(-1e10)} #optimizer can't handle Inf or -Inf values
  if (anyNA(gl_ll)==T) {gl_ll<-(-1e10)} #optimizer can't handle Inf or -Inf values

  return(gl_ll) 
}


#2) Synthetically generate residuals
et_syn<-function(sim_inflow,scale_flow,sigma_0,sigma_1,phi_1,phi_2,phi_3,beta,xi){
  at<-rsged(length(sim_inflow),mean=0,sd=1,nu=(2/(1+beta)),xi=xi)
  e_t<-c()
  sigma_t1<-sigma_t(sigma_0,sigma_1,scale_flow[1:3])
  e_t[1:3]<-sigma_t1*at[1:3] #assume initial error of zero
  
  for(i in 4:length(sim_inflow)){
    e_t[i]<-e_t[i-1]*phi_1 + e_t[i-2]*phi_2 + e_t[i-3]*phi_3 + sigma_t(sigma_0,sigma_1,scale_flow[i])*at[i]
  }
  Y_t<-sim_inflow + e_t
  Y_t[Y_t<0]<-0
  return(list(e_t,Y_t)) #returns list with el [[1]] = residuals, el [[2]] = synthetic flow
}



#4) define SEP density as function of calculated xi and beta and random variable axt (a_xi_t)
SEP_dens<-function(xi,beta,axt){
  SEP_dens<-(2*sigma_xi(M1(beta),M2,xi))/(xi + xi^(-1))*omega_beta(beta)*
    exp(-c_beta(beta)*abs(a_xi_t(xi,mu_xi(M1(beta),xi),sigma_xi(M1(beta),M2,xi),axt))^(2/(1+beta)))
}

