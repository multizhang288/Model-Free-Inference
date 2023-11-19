#Load Packages
#install.packages(c("pracma","dplyr","splines","stringr","MASS","ncvreg","hdi","glmnet"))
library(pracma)
library(dplyr)
library(splines)
library(stringr)
library(MASS)
library(ncvreg)
#library(hdi)
library(mvtnorm)
library(sn)
#library(LaplacesDemon)
################
gen_error<-function(N,p,rho,seed){
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
}

data_seq = function(x,h)
{
  value = seq(0,1,by=1/(h))
  then = quantile(x,probs=value)
  return(as.numeric(then))
}
#################################
Highdim.Ad = function(n,p,h,hfix,s,df,Trans,penalty,delta,Inactive,seed)
  {
  ################
  #n: sample size
  #p: dimension of covariates
  #h: degrees of freedom for B-spline Basis
  #hfix: if hfix = -1, h is forced to be 7. Otherwise, set your own h
  #s: sparsity level
  #d: the minimum structural dimension
  #df: degree of freedom of Chi-squared distribution
  #Trans: transformation of Y. Bspline, SIR, SIR2 and Poly
  #active: random: the indices of active variables are evenly set across 1:p. Otherwise, the first 1:s variables are active 
  #penalty: "lasso" or "SCAD"
  #CovMatirx: covariance matrix of X
  #model: data-generation model
  #delta: strength of active variables
  #seed: seed
  ################
  set.seed(seed)
  sigma = 1
  mu = matrix(rep(0,n*p),n,p)#rep(0,p)
  X = mu
  mu[,1:min(floor(p/2),100)] = -df 
  X_1 = gen_error(n,p,0,seed = seed) + mu
  ##########
  mu = matrix(rep(0,n*p),n,p)#rep(0,p)
  mu[,1:min(floor(p/2),100)] = df 
  X_2 = gen_error(n,p,0,seed = seed) + mu
  ###############
  U = runif(n,0,1)
  U_1 = which(U<=0.5)
  U_2 = which(U>0.5)
  X[U_1,] = X_1[U_1,]
  X[U_2,] = X_2[U_2,]
  interest = c(Inactive,1)
 ###########################  
  index.ini = 1:p
    error <- rnorm(n = n,0,sigma) 
    beta.true = rep(0,p)
    beta.true[1] = delta
    Y =  abs(X%*%beta.true + 3) + 1*error
    d=1
    if(hfix==-1)
    {
      h=7
    }
  ##############################Create the basis function with h = max(d+3,5)
  if(Trans == "Bspline")
  {
    if(h>=2)
    {
      temp = bs(Y,df=h,degree = 2)
    }
    if(h==1)
    {
      percens = data_seq(Y,h)
      temp = NULL
      for(k in 1:h)
      {
        fyk = ifelse(Y>=percens[k]&Y<=percens[k+1],Y,0)
        temp = cbind(temp,fyk)
      }
    }
  }
  if(Trans == "SIR")
  {
    percens = data_seq(Y,h)
    temp = NULL
    for(k in 1:(h-1))
    {
      fyk = ifelse(Y>=percens[k]&Y<percens[k+1],1,0)
      temp = cbind(temp,fyk)
    }
    fyk = ifelse(Y>=percens[h]&Y<=percens[h+1],1,0)
    temp = cbind(temp,fyk)
    
  }
  if(Trans == "SIR2")
  {
    percens = data_seq(Y,h)
    temp = NULL
    for(k in 1:(h-1))
    {
      fyk = ifelse(Y>=percens[k]&Y<percens[k+1],Y,0)
      temp = cbind(temp,fyk)
    }
    fyk = ifelse(Y>=percens[h]&Y<=percens[h+1],Y,0)
    temp = cbind(temp,fyk)
  }
  if(Trans == "Poly")
  {
    temp = NULL
    T = max(abs(Y))
    for(k in 1:h)
    {
      fyk1 = (Y/T)^k
      temp = cbind(temp,fyk1)
    }
  }
###############################Estimate gamma_{kj} using algorithm (3.2) in Remark 1
   model.save = rep(list(1),h)
    FYK =rep(list(1),h)
    for(k in 1:h)
    {
      fyk = temp[,k]
      model1 = cv.ncvreg(X,fyk,nlambda = 50,seed = seed,penalty = penalty)
      fit = model1$fit
      beta = fit$beta[,model1$min]
      FYK[[k]] = fyk
      model.save[[k]] = as.numeric(beta)
    }
  #####################Decorrelated score part
    model2 = cv.ncvreg(X,Y,penalty = penalty,seed = seed,nlambda = 50)
    fit = model2$fit
    beta = fit$beta[,model2$min]
    hatdecor = as.numeric(beta)
  ##############################
  record.stat = record.pvalue = record.pvalue1= rep(-1,p)
  ###############################################
  for(i in c(interest))
  {
    time0 = Sys.time()
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    cvfit = cv.ncvreg(Z,Xi,penalty=penalty,seed = seed,nlambda = 50)
    fit = cvfit$fit
    beta <- fit$beta[,cvfit$min]
    theta = as.numeric(beta)
    Z1 = cbind(1,Z)
    res = Xi - Z1%*%theta
    fisher = mean(res^2)
    eta = matrix(rep(0,n*h),n,h) 
    for(k in 1:h)
    {
      fyk = temp[,k]
      beta = model.save[[k]]
      hatgammak = as.numeric(beta)[-c(i+1)]
      eta.temp1 = fyk -Z1%*%hatgammak
      eta[,k] = eta.temp1*res
    }
    omega = t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    #####Calculate Statistic
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean #sometimes omega is not invertible, we use ginv here.
    record.stat[i] = Wni
    record.pvalue[i] = 1-pchisq(Wni,h)#Calculate the P-value
#####################Construction of TNL
    resdecor = Y-Z1%*%hatdecor[-c(i+1)]
    sigmahat = t(resdecor)%*%resdecor/n
    etadecor = resdecor*res
    Un = sum(etadecor)/sqrt(sigmahat*n*fisher)
    record.pvalue1[i] = 2*(1-pnorm(abs(Un)))
  }
  list(pvalue = record.pvalue[interest], pvalue.decor = record.pvalue1[interest],stat = record.stat[interest],h = h,d =d)
}

