#Load Packages
#install.packages(c("pracma","dplyr","splines","stringr","MASS","ncvreg","hdi","glmnet"))
library(pracma)
library(dplyr)
library(splines)
library(stringr)
library(MASS)
library(ncvreg)
library(hdi)
library(glmnet)
library(mvtnorm)
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
################
#n: sample size
#p: dimension of covariates
#h: degrees of freedom for B-spline Basis
#penalty: Penalty function, you can choose one between lasso and SCAD
#CovMatirx: covariance matrix of X
#model: data-generation model
#interest: the variable index of interest
################
Highdim.Fast = function(n,p,rho,h,penalty,model,q,interest)
  {
  index.ini = 1:p
  if(model==1)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    beta.true = c(1,1,rep(0,p-2))
    error <- rnorm(n = n,0,sigma) 
    Y = X%*%beta.true + 1*error
    Y = as.vector(Y)
  }
  if(model==3)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    error <- rnorm(n = n,0,sigma) 
    Y = 3*sin(X[,1])+3*sin(X[,p])+exp(-2*X[,3])*error
  }
  if(model == 2)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    error <- rnorm(n = n,0,sigma) 
    Y = (X[,1]+X[,2])/(0.5+(1.5+X[,3]+X[,4])^2)+0.1*error  
  }
  #####################Decorrelated score part
  model2 = cv.ncvreg(X,Y,penalty = penalty,seed = q,nlambda = 50)
  fit = model2$fit
  beta = fit$beta[,model2$min]
  hatdecor = as.numeric(beta)
  ###############################
  record_pvalue = record_pvalue1 = rep(-1,p)
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
        fyk = ifelse(Y>=percens[k]&Y<percens[k+1],Y,0)
        temp = cbind(temp,fyk)
      }
    }
  ##########################Estimate gamma_{kj} using algorithm (3.2) in Remark 1
    model.save = rep(list(1),h)
    FYK =rep(list(1),h)
    for(k in 1:h)
    {
      fyk = temp[,k]
      model1 = cv.ncvreg(X,fyk,nlambda = 50,seed = q,penalty = penalty)
      fit = model1$fit
      beta = fit$beta[,model1$min]
      FYK[[k]] = fyk
      model.save[[k]] = as.numeric(beta)
    }
  ##########################
  for(i in interest)
  {
    time0 = Sys.time()
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    Z1 = cbind(1,Z)
    cvfit = cv.ncvreg(Z,Xi,penalty=penalty,seed = q,nlambda = 50)
    fit = cvfit$fit
    beta <- fit$beta[,cvfit$min]
    hattheta = as.numeric(beta)
    res = Xi - Z1%*%hattheta
    fisher =mean(res^2)
    eta = matrix(rep(0,n*h),n,h)
    for(k in 1:h)
    {
      fyk = temp[,k]
      beta <- model.save[[k]]
      hatgammak = as.numeric(beta)[-c(i+1)]
      eta.temp = fyk-Z1%*%hatgammak
      eta[,k] = eta.temp*res
    }
    omega =t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    omega = t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    #####Calculate Statistic
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean #sometimes omega is not invertible, we use ginv here.
    record_pvalue[i] = 1-pchisq(Wni,h)#Calculate the P-value
    #####################Construction of TNL
    resdecor = Y-Z1%*%hatdecor[-c(i+1)]
    sigmahat = t(resdecor)%*%resdecor/n
    etadecor = resdecor*res
    Un = sum(etadecor)/sqrt(sigmahat*n*fisher)
    record_pvalue1[i] = 2*(1-pnorm(abs(Un)))
    #############################
  }
  list(pvalue = record_pvalue[interest], pvalue1 = record_pvalue1[interest])
}

################
Highdim.Slow = function(n,p,rho,h,penalty,model,q,interest)
{
  index.ini = 1:p
  if(model==1)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    beta.true = c(1,1,rep(0,p-2))
    error <- rnorm(n = n,0,sigma) 
    Y = X%*%beta.true + 1*error
    Y = as.vector(Y)
  }
  if(model==3)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    error <- rnorm(n = n,0,sigma) 
    Y = 3*sin(X[,1])+3*sin(X[,p])+exp(-2*X[,3])*error
  }
  if(model == 2)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    error <- rnorm(n = n,0,sigma) 
    Y = (X[,1]+X[,2])/(0.5+(1.5+X[,3]+X[,4])^2)+0.1*error  
  }
  #####################Decorrelated score part
  model2 = cv.ncvreg(X,Y,penalty = penalty,seed = q,nlambda = 50)
  fit = model2$fit
  beta = fit$beta[,model2$min]
  hatdecor = as.numeric(beta)
  ###############################
  record_pvalue = record_pvalue1 = rep(-1,p)
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
      fyk = ifelse(Y>=percens[k]&Y<percens[k+1],Y,0)
      temp = cbind(temp,fyk)
    }
  }
  ##########################
  for(i in interest)
  {
    time0 = Sys.time()
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    Z1 = cbind(1,Z)
    cvfit = cv.ncvreg(Z,Xi,penalty=penalty,seed = q,nlambda = 50)
    fit = cvfit$fit
    beta <- fit$beta[,cvfit$min]
    hattheta = as.numeric(beta)
    res = Xi - Z1%*%hattheta
    fisher =mean(res^2)
    eta = matrix(rep(0,n*h),n,h)
    for(k in 1:h)
    {
      fyk = temp[,k]
      cvfit = cv.ncvreg(Z,fyk,penalty=penalty,seed = q,nlambda = 50)
      fit = cvfit$fit
      beta <- fit$beta[,cvfit$min]
      hatgammak = as.numeric(beta)
      eta.temp = fyk-Z1%*%hatgammak
      eta[,k] = eta.temp*res
    }
    omega =t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    omega = t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    #####Calculate Statistic
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean #sometimes omega is not invertible, we use ginv here.
    record_pvalue[i] = 1-pchisq(Wni,h)#Calculate the P-value
    #####################Construction of TNL
    resdecor = Y-Z1%*%hatdecor[-c(i+1)]
    sigmahat = t(resdecor)%*%resdecor/n
    etadecor = resdecor*res
    Un = sum(etadecor)/sqrt(sigmahat*n*fisher)
    record_pvalue1[i] = 2*(1-pnorm(abs(Un)))
    #############################
  }
  list(pvalue = record_pvalue[interest], pvalue1 = record_pvalue1[interest])
}
