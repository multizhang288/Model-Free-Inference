#Load Packages
#install.packages(c("pracma","dplyr","splines","stringr","MASS","ncvreg","hdi","glmnet"))
library(pracma)
library(dplyr)
library(splines)
library(stringr)
library(MASS)
library(ncvreg)
library(hdi)
library(mvtnorm)
library(sn)
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
Highdim.Study2 = function(n,p,h,hfix,s,d,Dist,Trans,active,penalty,CovMatrix,model,delta,Inactive,seed)
  {
  ################
  #n: sample size
  #p: dimension of covariates
  #h: degrees of freedom for B-spline Basis
  #hfix: if hfix = -1, h is forced to be 7. Otherwise, set your own h
  #s: sparsity level
  #d: the minimum structural dimension
  #Dist: distribution of X. MVN, Mixture, MSN, Tdist and Chi
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
  if(Dist == "MVN")
  {
    X = gen_error(n,p,0.5,seed = seed)
  }
  if(Dist == "Tdist")
  {
    X = rmvt(n,CovMatrix,df = 5)
  }
  if(Dist == "Chi")
  {
    X = NULL
    for(kk in 1:n)
    {
      X0 = rchisq(p,df=5) - 5
      X = rbind(X,X0)
    }
  }
  if(Dist == "Mixture")
  {
    chose = sample(c(1,2,3),1,prob=c(0.4,0.2,0.4))
    if(chose==1)
    {
      CovMatrix2 = matrix(rep(0,p^2),p,p)
      for(i in 1:p)
      {
        for(j in 1:p)
        {
          CovMatrix2[i,j] = 0.2^(abs(i-j))  
        }
      }
      mu = rep(0,p)
      mu[1:6] = -1 
      X = rmvnorm(n,mu,CovMatrix2)
    }
    if(chose==2)
    {
      CovMatrix2 = matrix(rep(0,p^2),p,p)
      for(i in 1:p)
      {
        for(j in 1:p)
        {
          CovMatrix2[i,j] = 0.5^(abs(i-j))  
        }
      }
      mu = rep(0,p)
      X = rmvnorm(n,mu,CovMatrix2)
    }
    if(chose==3)
    {
      CovMatrix2 = matrix(rep(0,p^2),p,p)
      for(i in 1:p)
      {
        for(j in 1:p)
        {
          CovMatrix2[i,j] = 0.8^(abs(i-j))  
        }
      }
      mu = rep(0,p)
      mu[1:6] = 1 
      X = rmvnorm(n,mu,CovMatrix2)
    }
  }
  if(Dist == "MSN")
  {
    CovMatrix2 = diag(1,p)
    for(i in 1:p)
    {
      for(j in 1:p)
      {
        CovMatrix2[i,j] = 0.5^(abs(i-j))  
      }
    }
    mu = rep(0,p)
    mu[c(1,2,5,6)]=1
    X = rmsn(n,xi= rep(0,p),Omega=CovMatrix2,alpha = mu)
  }
###########################Determine the location of actice variables
  intv = floor(p/s)
  if(active == "random")
  {
    location = intv*(0:(s-1))+1
  }
  if(active == "consecutive")
  {
    location = 1:s
    s = length(location)
  }
  interest = c(Inactive,location)
  ###########################Determine the index we used (No overlap)
  Index = NULL
  d=10
  intd = ceiling(s/d)
  for(l in 1:d)
  {
    sindex = ((l-1)*intd+1):(l*intd)
    sindex[sindex>s] = sindex[sindex>s]-s
    sindex1 = location[sindex]
    thetal = rep(delta,intd)
    if(intd==1)
    {
      Index = cbind(Index,X[,sindex1]*thetal)
    }
    if(intd>1)
    {
      Index = cbind(Index,X[,sindex1]%*%thetal)
    }
  }
  ##########################Create the mean function .i.e, require model==-1
  if(model==6)
  {
    mean_fun = NULL
    mean_fun = (Index[,1]) + sinh(Index[,2]) + exp(Index[,3]) + abs(Index[,4]+3) + 4*sin(Index[,5]) +sinh(Index[,6]) + exp(Index[,7]) + abs(Index[,8]+3) + 4*sin(Index[,9]) 
    d = 10
    if(hfix==-1)
    {
      h= 15
    }
  }
 ###########################  
  index.ini = 1:p
  error = rnorm(n,0,sigma)
  if(model==6)
  {
    Y = mean_fun + exp(-2*Index[,10])*error 
  }
  if(model==1)
  {
    error <- rnorm(n = n,0,sigma) 
    beta.true = rep(0,p)
    beta.true[location] = delta
    Y =(X%*%beta.true)+ 1*error
    Y = as.vector(Y)
    
    d=1
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model==2)
  {
    error <- rnorm(n = n,0,sigma) 
    beta1 = beta2 = rep(0,p)
    ll = length(location)
    beta1[location[1:(floor(ll*0.8))]] = delta
    beta2[location[-c(1:(floor(ll*0.8)))]] = delta
    Y = (X%*%beta1)/(0.5+(1.5+X%*%beta2)^2)+0.1*error  
    
    d=2
    if(hfix==-1)
    {
      h=7
    }
    }
  if(model == 3)
  {
    error <- rnorm(n = n,0,sigma) 
    beta1 = beta2 = beta3 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/3)]] = delta
    beta2[location[c((floor(ll/3)+1):(2*floor(ll/3)))]] = delta
    beta3[location[-c(1:(2*floor(ll/3)))]] = delta
    Y = sinh(X%*%beta1) + sinh(X%*%beta2) + sinh(X%*%beta3) +0.5*error

    d=3
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model == 4)
  {
    error <- rnorm(n = n,0,sigma) 
    beta1 = beta2 = beta3 = beta4 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/4)]] = delta
    beta2[location[c((floor(ll/4)+1):(2*floor(ll/4)))]] = delta
    beta3[location[c((2*floor(ll/4)+1):(3*floor(ll/4)))]] = delta
    beta4[location[-c(1:(3*floor(ll/4)))]] = delta
    Y = X%*%beta1 + sinh(X%*%beta2) + 0.3*exp(X%*%beta3) + abs(X%*%beta4+3) +1*error
    
    d=4
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model == 5)
  {
    error <- rnorm(n = n,0,sigma) 
    beta1 = beta2 = beta3 = beta4 = beta5 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/5)]] = delta
    beta2[location[c((floor(ll/5)+1):(2*floor(ll/5)))]] = delta
    beta3[location[c((2*floor(ll/5)+1):(3*floor(ll/5)))]] = delta
    beta4[location[c((3*floor(ll/5)+1):(4*floor(ll/5)))]] = delta
    beta5[location[-c(1:(4*floor(ll/5)))]] = delta
    Y = X%*%beta1 + sinh(X%*%beta2) + exp(X%*%beta3) + abs(X%*%beta4+3) + 4*sin(X%*%beta5) + 1*error
    
    d=5
    if(hfix==-1)
    {
      h=7
    }
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

